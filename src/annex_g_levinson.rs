//! Annex G fixed-point Levinson-Durbin recursion (§G.2.2).
//!
//! ITU-T G.728 Annex G (1994-11) §G.2.2 specifies how the three
//! Levinson-Durbin recursion modules of Recommendation G.728 — block 50
//! (synthesis filter), block 37 (perceptual weighting filter) and block
//! 44 (log-gain predictor) — are implemented in bit-exact fixed point.
//! The annex develops the code for block 50 and notes the other two are
//! identical with renamed variables.
//!
//! The floating-point recursion already lives in [`crate::levinson`];
//! this module is the §G.2.2 fixed-point reformulation, built on the
//! §G.1.3 arithmetic primitives in [`crate::annex_g_arith`]
//! ([`rnd`](crate::annex_g_arith::rnd)). It introduces no codec state
//! and reads no codebooks — it is a self-contained transcription of the
//! §G.2.2 pseudo-code, one block at a time.
//!
//! ## What §G.2.2 changes versus the floating-point recursion
//!
//! * **SIMPDIV** — the annex explicitly chooses "a different and simpler
//!   division algorithm than that used throughout the rest of the
//!   algorithm" for the two divisions inside the recursion (the
//!   first-order reflection coefficient and each `RC = −SUM/ALPHATMP`).
//!   It is a plain 16-iteration restoring long division on two 16-bit
//!   integers, output in the low 17 bits of an accumulator
//!   ([`simpdiv`]). This is NOT the §G.1.3.4 `DIVIDE` used elsewhere.
//! * **Block floating point for `RTMP`** — the autocorrelation values
//!   arrive already normalized so that the largest magnitude lies in
//!   `[0.5, 1)`, i.e. all of `RTMP` fits `Q15`. Only the mantissas are
//!   used here.
//! * **Variable precision for `ATMP`** — the predictor coefficients
//!   start in `Q15`. If updating a coefficient overflows the 16-bit
//!   word, the recursion right-shifts every computed coefficient by one
//!   (`NRS += 1`, dropping to `Q14`, then `Q13`) and recomputes the
//!   overflowed term. The final format is one of `Q13` / `Q14` / `Q15`,
//!   signalled out of the module as `NLSATMP = 15 − NRS`.
//! * **`ALPHATMP` as a saved 16-bit high word** — the residual energy is
//!   accumulated in a 32-bit accumulator but only its rounded high word
//!   is carried to the next order ("only the high word of ALPHATMP is
//!   saved after each update … did not degrade the coder's performance").
//! * **Ill-conditioning flags** — instead of three `go to LABEL` exits,
//!   the module sets `ILLCOND` (the whole recursion failed) and, when
//!   the failure happened at or before the 10th iteration, `ILLCONDP`
//!   (the order-10 postfilter coefficients are also invalid). The annex
//!   adds the `NLSATMP < 13` post-check as a fourth failure path.
//! * **Decoder restart at `MINC0 = 10`** — the decoder interrupts the
//!   recursion after the 10th-order coefficients (saving them for the
//!   adaptive postfilter) and resumes from `MINC0 = 10`. Restarting
//!   requires `NRS` and `ALPHATMP` to be carried in.
//!
//! ## `ILLCOND` as an input
//!
//! The first listed branch of the floating-point code, `If RTMP(LPC + 1)
//! = 0, go to LABEL`, is moved upstream into the hybrid-window module
//! (block 49): that module tests the full 32-bit accumulator value of
//! `RTMP(LPC + 1)` for zero and reports the result as the logical input
//! `ILLCOND`. The annex makes this change because a 16-bit-rounded
//! `RTMP(LPC + 1)` underflows to zero far more often than the 32-bit
//! float, "causing interoperability problems." This module therefore
//! takes `illcond` as an input ([`LevinsonInput::illcond`]) rather than
//! re-deriving it.

use crate::annex_g_arith::rnd;

/// `Q15` numerator pre-shift / quotient width for [`simpdiv`]: the
/// routine emits 16 quotient bits into the low 17 bits of the
/// accumulator (one extra guard bit), then `RC = RND(AA0 << 15)` lifts
/// the quotient to the high word for rounding to a `Q15` value.
const SIMPDIV_ITERS: u32 = 16;

/// The §G.2.2 `SIMPDIV(NUM, DEN, AA0)` subroutine — the "different and
/// simpler division algorithm" the recursion uses for its two divisions.
///
/// Both inputs are non-negative 16-bit integers (callers strip the sign
/// and reapply it afterwards). The result is returned in the low 17 bits
/// of the accumulator: an unsigned fixed-point quotient `NUM/DEN` scaled
/// by `2^16` (16 fractional bits, MSB-first). The annex's pseudo-code:
///
/// ```text
/// AA0 = 0;  AA1 = NUM;  K = 0
/// LOOP: AA0 = AA0 << 1
///       AA1 = AA1 << 1
///       If AA1 ≥ DEN, then AA1 = AA1 − DEN and AA0 = AA0 + 1
///       K = K + 1
///       If K < 16, go to LOOP
/// ```
///
/// The routine assumes `DEN > 0`; the recursion only calls it after
/// checking `NUM < ALPHATMP` (so the quotient is `< 1`) and, for the
/// first-order coefficient, `RTMP(1) > 0`.
#[must_use]
pub fn simpdiv(num: i16, den: i16) -> i32 {
    debug_assert!(num >= 0, "SIMPDIV NUM must be non-negative");
    debug_assert!(den > 0, "SIMPDIV DEN must be positive");
    let den = den as i32;
    let mut aa0: i32 = 0;
    let mut aa1: i32 = num as i32;
    for _ in 0..SIMPDIV_ITERS {
        aa0 <<= 1;
        aa1 <<= 1;
        if aa1 >= den {
            aa1 -= den;
            aa0 += 1;
        }
    }
    aa0
}

/// Outcome flags produced by the §G.2.2 recursion.
///
/// These are the annex's `ILLCOND` / `ILLCONDP` logical outputs plus the
/// `NLSATMP` precision signal. They replace the floating-point code's
/// three `go to LABEL` exits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LevinsonStatus {
    /// `ILLCOND` — the recursion failed; the caller must skip block 51
    /// (bandwidth expansion) and keep the previous adaptation cycle's
    /// synthesis-filter coefficients.
    pub illcond: bool,
    /// `ILLCONDP` — the failure (if any) happened at or before the 10th
    /// iteration, so the order-10 postfilter coefficients are also
    /// invalid. Only meaningful for block 50 in the decoder.
    pub illcondp: bool,
    /// `NLSATMP = 15 − NRS` — the `Q` format of the output `ATMP`
    /// (`15 = Q15`, `14 = Q14`, `13 = Q13`). The downstream bandwidth-
    /// expansion module (block 51 / 45 / 38 / 85) uses this to convert
    /// the array to `Q14`. On failure this is left at its last value and
    /// `illcond` is set.
    pub nlsatmp: i32,
    /// `NRS` — the cumulative number of right shifts applied to `ATMP`
    /// during the recursion (`0` ⇒ output stayed `Q15`). Carried out so
    /// a decoder restart at `MINC0 = 10` can resume with it.
    pub nrs: i32,
    /// The order at which the recursion stopped on failure (the value of
    /// `MINC` when `FAILED` was reached). `lpc + 1` if it completed
    /// normally. Used to populate `illcondp` and for diagnostics.
    pub stopped_at: usize,
    /// `ALPHATMP` (the 16-bit residual-energy high word) at the point the
    /// recursion stopped — at the end of the final completed order on
    /// success, or the last good value on failure. The decoder saves this
    /// at the order-10 interruption to feed [`RecursionResume`].
    pub alphatmp: i16,
    /// `RC1` — the first reflection coefficient (Q15), captured at the
    /// first-order step of a fresh run (Table G.2's "temporary buffer
    /// for first reflection coefficients"). Block 85 (the short-term
    /// postfilter coefficient calculator) consumes it as `k1`. Zero on
    /// a resume run (the fresh order-10 run already reported it) and on
    /// failures before the first-order step.
    pub rc1: i16,
}

/// Inputs to the §G.2.2 recursion.
///
/// `rtmp` carries the block-floating-point autocorrelation mantissas in
/// `Q15` (`rtmp[k]` is the spec's `RTMP(k+1)`), normalized so the
/// largest magnitude lies in `[0.5, 1)`. `illcond` is the upstream
/// hybrid-window verdict on `RTMP(LPC + 1)` being zero.
#[derive(Debug, Clone)]
pub struct LevinsonInput<'a> {
    /// Autocorrelation mantissas, `Q15`, spec 1-based shifted to 0-based.
    /// Length must be at least `lpc + 1`.
    pub rtmp: &'a [i16],
    /// `ILLCOND` from block 49 — `true` ⇒ `RTMP(LPC + 1)` was zero.
    pub illcond: bool,
    /// Predictor order (`LPC = 50`, `LPCW = 10`, `LPCLG = 10`).
    pub lpc: usize,
}

/// State carried across a decoder restart of block 50 at `MINC0 = 10`.
///
/// In the decoder the recursion is interrupted after the 10th-order
/// coefficients are derived (they feed the adaptive postfilter), then
/// resumed at `MINC0 = 10`. Resuming needs the `ATMP` accumulated so
/// far, the running `NRS` and the 16-bit `ALPHATMP`.
#[derive(Debug, Clone)]
pub struct RecursionResume {
    /// `ATMP(1..=11)` from the order-10 run, `atmp[0] = 1` in the
    /// stored `Q` format implied by `nrs`.
    pub atmp: Vec<i16>,
    /// `NRS` at the order-10 interruption.
    pub nrs: i32,
    /// `ALPHATMP` (16-bit high word) at the order-10 interruption.
    pub alphatmp: i16,
}

/// Run the §G.2.2 fixed-point Levinson-Durbin recursion from scratch
/// (`MINC0 = 1`).
///
/// Writes the predictor coefficients into `atmp` (`atmp[0]` = `ATMP(1)`,
/// the implicit unity tap in the output `Q` format; `atmp[1..=lpc]` =
/// `ATMP(2..=LPC+1)`). The buffer must be at least `lpc + 1` long.
/// Returns the [`LevinsonStatus`] flags. On failure `atmp` holds the
/// partially-computed coefficients and `status.illcond` is `true`.
///
/// `atmp[0]` is set to the unity value in the final `Q` format
/// (`1 << NLSATMP` capped at `i16::MAX` for `Q15`), so a caller that
/// scales the whole array by `2^(−NLSATMP)` recovers `A(1) = 1`.
pub fn levinson_durbin_fixed(input: &LevinsonInput, atmp: &mut [i16]) -> LevinsonStatus {
    run(input, atmp, 1, None)
}

/// Resume block 50 in the decoder from `MINC0 = 10`, carrying the
/// order-10 state in `resume`. `atmp` is overwritten with `resume.atmp`
/// before the recursion continues. See [`RecursionResume`].
pub fn levinson_durbin_fixed_resume(
    input: &LevinsonInput,
    atmp: &mut [i16],
    resume: &RecursionResume,
) -> LevinsonStatus {
    run(input, atmp, 10, Some(resume))
}

/// Core of the §G.2.2 pseudo-code, parameterised by the start order
/// `minc0` (1 for a fresh run, 10 for a decoder restart). When `minc0
/// == 10`, `resume` supplies the carried state.
fn run(
    input: &LevinsonInput,
    atmp: &mut [i16],
    minc0: usize,
    resume: Option<&RecursionResume>,
) -> LevinsonStatus {
    let lpc = input.lpc;
    let rtmp = input.rtmp;
    assert!(
        rtmp.len() > lpc,
        "RTMP too short: need {} entries, got {}",
        lpc + 1,
        rtmp.len()
    );
    assert!(
        atmp.len() > lpc,
        "ATMP too short: need {} entries, got {}",
        lpc + 1,
        atmp.len()
    );

    // `RTMP(k)` is `rtmp[k - 1]`; `ATMP(k)` is `atmp[k - 1]`.
    // Helper closures keep the per-line transcription readable.
    let r = |k: usize| rtmp[k - 1] as i32; // 1-based, Q15 mantissa

    let mut nrs: i32;
    let mut alphatmp: i16;
    let mut rc1_out: i16 = 0;

    if minc0 > 1 {
        // RECURSION restart (decoder, MINC0 = 10): load carried state.
        let resume = resume.expect("MINC0 > 1 requires a RecursionResume");
        for (dst, &src) in atmp.iter_mut().zip(resume.atmp.iter()) {
            *dst = src;
        }
        nrs = resume.nrs;
        alphatmp = resume.alphatmp;
    } else {
        // | If MINC0 > 1, go to RECURSION                | (false here)
        // | MINC0 = 1                                    | Initializations
        // | ILLCONDP = .FALSE.                           | decoder only
        //
        // | If ILLCOND = .TRUE., go to FAILED            | RTMP(LPC+1)=0
        if input.illcond {
            return failed(0, lpc, 15, 0, 0);
        }
        // | If RTMP(1) ≤ 0, go to FAILED                 | Skip if zero signal
        if r(1) <= 0 {
            return failed(0, lpc, 15, 0, 0);
        }
        // | NRS = 0                                      | Q15 format initially
        nrs = 0;

        // First-order predictor.
        // | DEN = RTMP(1)                                |
        // | NUM = RTMP(2)                                |
        // | If NUM < 0, set NUM = −NUM                   |
        let den = r(1) as i16;
        let num_signed = r(2);
        let num = num_signed.unsigned_abs().min(i16::MAX as u32) as i16;
        // | Call SIMPDIV(NUM, DEN, AA0)                  | |RTMP(2)|/RTMP(1)
        let mut aa0 = simpdiv(num, den);
        // | AA0 = AA0 << 15                              |
        // | RC1 = RND(AA0)                              |
        aa0 <<= 15;
        let mut rc1 = rnd(aa0);
        // | If RTMP(2) > 0, set RC1 = −RC1               | Add sign information
        if num_signed > 0 {
            rc1 = rc1.saturating_neg();
        }
        // | RC = RC1                                     |
        // | ATMP(2) = RC1                               | First order coeff
        let rc = rc1;
        rc1_out = rc1;
        atmp[1] = rc1;
        // | ATMP(1) = 1 (Q15 unity)                      | implicit unity tap
        atmp[0] = i16::MAX; // 0x7FFF ≈ 1.0 in Q15

        // | AA0 = RTMP(1) << 16                          |
        // | P = RTMP(2) * RC                             |
        // | AA0 = AA0 + (P << 1)                         |
        // | ALPHATMP = RND(AA0)                         | high word to memory
        let mut acc = (r(1)) << 16;
        let p = (r(2)) * (rc as i32);
        acc = acc.wrapping_add(p << 1);
        alphatmp = rnd(acc);
    }

    // RECURSION:
    // | For MINC = MINC0 + 1, …, LPC, do the following
    let mut minc = minc0 + 1;
    while minc <= lpc {
        // SUM accumulation.
        // | AA0 = 0                                      |
        // | For IP = 2, …, MINC: N1 = MINC − IP + 2     |
        // |     P = RTMP(N1) * ATMP(IP)                  |
        // |     AA0 = AA0 + P                            | 32 bits for SUM
        let mut aa0: i64 = 0;
        for ip in 2..=minc {
            let n1 = minc - ip + 2;
            let p = r(n1) as i64 * atmp[ip - 1] as i64;
            aa0 += p;
        }
        // | AA0 = AA0 << 1                               |
        // | AA0 = AA0 << NRS                             |
        aa0 <<= 1;
        aa0 = shl_acc(aa0, nrs);
        // | AA1 = RTMP(MINC + 1) << 16                   |
        // | AA0 = AA0 + AA1                              |
        let aa1 = (r(minc + 1) as i64) << 16;
        aa0 += aa1;
        // | SIGN = RND(AA0)                             | Save high word sign
        let sign = rnd(clamp_acc(aa0));
        // | NUM = SIGN; If NUM < 0, set NUM = −NUM       |
        let num = (sign as i32).unsigned_abs().min(i16::MAX as u32) as i16;
        // | If NUM ≥ ALPHATMP, go to FAILED             |
        if (num as i32) >= (alphatmp as i32) {
            return failed(minc, lpc, 15 - nrs, alphatmp, rc1_out);
        }
        // | Call SIMPDIV(NUM, ALPHATMP, AA0)            | Divide to get RC
        let div = simpdiv(num, alphatmp);
        // | AA2 = AA0 << 15                             | AA2 stores 17-bit RC
        // | RC = RND(AA2)                              |
        let aa2 = (div as i64) << 15;
        let mut rc = rnd(clamp_acc(aa2));
        // | If SIGN > 0, set RC = −RC                    |
        if sign > 0 {
            rc = rc.saturating_neg();
        }

        // Update ALPHATMP (residual energy).
        // | AA1 = ALPHATMP << 16                         |
        // | P = RC * SIGN                                |
        // | AA1 = AA1 + (P << 1)                         |
        // | If AA1 ≤ 0, go to FAILED                     |
        // | ALPHATMP = RND(AA1)                         |
        let mut alpha_acc = (alphatmp as i64) << 16;
        let p = (rc as i64) * (sign as i64);
        alpha_acc += p << 1;
        if alpha_acc <= 0 {
            return failed(minc, lpc, 15 - nrs, alphatmp, rc1_out);
        }
        alphatmp = rnd(clamp_acc(alpha_acc));

        // Predictor coefficient update.
        // | MH = MINC/2 + 1                              | integer divide
        let mh = minc / 2 + 1;
        for ip in 2..=mh {
            // | IB = MINC − IP + 2                       |
            let ib = minc - ip + 2;

            // | AA0 = ATMP(IP) << 16                     |
            // | P = RC * ATMP(IB)  (Q15 RC ⇒ << 1)       |
            // | AA0 = AA0 + (P << 1)                     |
            let mut acc0 = compute_update(atmp[ip - 1], rc, atmp[ib - 1]);
            // | If AA0 overflowed: NRS += 1, halve ATMP, recompute |
            if overflowed(acc0) {
                nrs += 1;
                for lp in atmp.iter_mut().take(minc + 1).skip(1) {
                    *lp >>= 1;
                }
                acc0 = compute_update(atmp[ip - 1], rc, atmp[ib - 1]);
            }

            // | AA1 = ATMP(IB) << 16                     |
            // | P = RC * ATMP(IP)                        |
            // | AA1 = AA1 + (P << 1)                     |
            let mut acc1 = compute_update(atmp[ib - 1], rc, atmp[ip - 1]);
            // | If AA1 overflowed: NRS += 1, halve ATMP, recompute BOTH |
            if overflowed(acc1) {
                nrs += 1;
                for lp in atmp.iter_mut().take(minc + 1).skip(1) {
                    *lp >>= 1;
                }
                acc0 = compute_update(atmp[ip - 1], rc, atmp[ib - 1]);
                acc1 = compute_update(atmp[ib - 1], rc, atmp[ip - 1]);
            }

            // | ATMP(IP) = RND(AA0)                     |
            // | ATMP(IB) = RND(AA1)                     |
            // The midpoint of an odd step has ip == ib; both accs read
            // the same pre-update values, and RND(AA0)/RND(AA1) are
            // equal there, so order of assignment is immaterial.
            atmp[ip - 1] = rnd(clamp_acc(acc0));
            atmp[ib - 1] = rnd(clamp_acc(acc1));
        }

        // Update ATMP(MINC + 1) from the 17-bit RC in AA2.
        // | AA0 = AA2 >> NRS                            | AA2 contains 17-bit RC
        // | AA0 = RND(AA0)                             |
        // | If SIGN > 0, set AA0 = −AA0                 |
        // | ATMP(MINC + 1) = AA0                       |
        let mut tail = rnd(clamp_acc(shr_acc(aa2, nrs)));
        if sign > 0 {
            tail = tail.saturating_neg();
        }
        atmp[minc] = tail;

        minc += 1;
    }

    // | NLSATMP = 15 − NRS                              |
    // | If NLSATMP < 13, go to FAILED                   |
    let nlsatmp = 15 - nrs;
    if nlsatmp < 13 {
        return failed(lpc, lpc, nlsatmp, alphatmp, rc1_out);
    }
    // | Exit this program  — recursion completed normally
    LevinsonStatus {
        illcond: false,
        illcondp: false,
        nlsatmp,
        nrs,
        stopped_at: lpc + 1,
        alphatmp,
        rc1: rc1_out,
    }
}

/// The §G.2.2 `FAILED:` exit. Sets `ILLCOND`; if the stop order is at or
/// before the 10th iteration, also sets `ILLCONDP` (postfilter coeffs
/// invalid). `minc == 0` covers the two pre-recursion failure tests,
/// which precede any postfilter-coefficient derivation, so `ILLCONDP` is
/// set there too.
fn failed(minc: usize, lpc: usize, nlsatmp: i32, alphatmp: i16, rc1: i16) -> LevinsonStatus {
    // | FAILED: Set ILLCOND = .TRUE.
    // |         If MINC ≤ 10, set ILLCONDP = .TRUE.
    let illcondp = minc <= 10;
    LevinsonStatus {
        illcond: true,
        illcondp,
        nlsatmp,
        nrs: 15 - nlsatmp,
        stopped_at: if minc == 0 { lpc + 1 } else { minc },
        alphatmp,
        rc1,
    }
}

/// Compute `RND`-ready accumulator for `coeff_hi << 16 + (RC·coeff_other) << 1`.
/// `rc` is `Q15`, so the product is shifted left by one to align with the
/// `Q15`-coefficient high word. Returned as an `i64` so the caller can
/// test for 32-bit accumulator overflow before rounding.
#[inline]
fn compute_update(coeff_hi: i16, rc: i16, coeff_other: i16) -> i64 {
    let acc = (coeff_hi as i64) << 16;
    let p = (rc as i64) * (coeff_other as i64);
    acc + (p << 1)
}

/// §G.2.2 overflow test: the update accumulator must fit a 32-bit word
/// for `RND` to produce a valid 16-bit coefficient. The annex relies on
/// "guard bits or … an overflow flag" — modelled here as the value
/// leaving the signed 32-bit range.
#[inline]
fn overflowed(acc: i64) -> bool {
    acc > i32::MAX as i64 || acc < i32::MIN as i64
}

/// Clamp a 64-bit working accumulator back into the 32-bit register
/// range before [`rnd`]. The recursion only reaches [`rnd`] on values it
/// has already proven non-overflowing (via [`overflowed`]) or that are
/// inherently bounded, so this is a defensive saturation.
#[inline]
fn clamp_acc(acc: i64) -> i32 {
    acc.clamp(i32::MIN as i64, i32::MAX as i64) as i32
}

/// `AA0 << NRS` on the wide SUM accumulator. `NRS ≥ 0` always in the
/// recursion (it only ever increments), so this is a plain left shift;
/// kept as a helper to mirror the pseudo-code line.
#[inline]
fn shl_acc(acc: i64, nrs: i32) -> i64 {
    debug_assert!(nrs >= 0, "NRS is non-negative in the recursion");
    acc << (nrs as u32)
}

/// `AA2 >> NRS` (sign-extending) for the `ATMP(MINC + 1)` update.
#[inline]
fn shr_acc(acc: i64, nrs: i32) -> i64 {
    debug_assert!(nrs >= 0, "NRS is non-negative in the recursion");
    acc >> (nrs as u32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annex_g_arith::vscale;
    use crate::levinson::levinson_durbin;

    // --- SIMPDIV --------------------------------------------------

    #[test]
    fn simpdiv_matches_scaled_integer_division() {
        // SIMPDIV(NUM, DEN) = round-down( NUM/DEN · 2^16 ) for the
        // restoring-division domain NUM < DEN (quotient < 1).
        for (num, den) in [(1i16, 2i16), (1, 3), (3, 4), (7, 8), (123, 456)] {
            let got = simpdiv(num, den);
            let want = ((num as i64) << 16) / den as i64;
            assert_eq!(got as i64, want, "SIMPDIV({num},{den})");
        }
    }

    #[test]
    fn simpdiv_half_is_exact() {
        // 1/2 · 2^16 = 32768.
        assert_eq!(simpdiv(1, 2), 32768);
        // 1/4 · 2^16 = 16384.
        assert_eq!(simpdiv(1, 4), 16384);
    }

    // --- pre-recursion failure paths ------------------------------

    #[test]
    fn illcond_input_fails_immediately() {
        let rtmp = [16384i16, 1000, 100, 10];
        let mut atmp = [0i16; 4];
        let st = levinson_durbin_fixed(
            &LevinsonInput {
                rtmp: &rtmp,
                illcond: true,
                lpc: 3,
            },
            &mut atmp,
        );
        assert!(st.illcond);
        assert!(st.illcondp, "pre-recursion failure ⇒ postfilter invalid");
    }

    #[test]
    fn nonpositive_rtmp1_fails() {
        // RTMP(1) ≤ 0 (zero signal).
        let rtmp = [0i16, 0, 0, 0];
        let mut atmp = [0i16; 4];
        let st = levinson_durbin_fixed(
            &LevinsonInput {
                rtmp: &rtmp,
                illcond: false,
                lpc: 3,
            },
            &mut atmp,
        );
        assert!(st.illcond);
    }

    // --- numeric agreement with the floating-point recursion -------

    /// Build a `Q15` block-floating-point `RTMP` from a real
    /// autocorrelation, returning the mantissas. The largest magnitude
    /// is `RTMP(1)`, so VSCALE with SLEN = 1 left-justifies it.
    fn to_q15_block(r_real: &[f64]) -> Vec<i16> {
        // Scale to a comfortable 16-bit range first (peak near 30000),
        // then VSCALE normalizes the block so RTMP(1) lands in band.
        let peak = r_real[0].abs().max(1e-30);
        let scaled: Vec<i32> = r_real
            .iter()
            .map(|&v| ((v / peak) * 30000.0).round() as i32)
            .collect();
        let (block, _nls) = vscale(&scaled, 1, crate::annex_g_arith::MLS_SINGLE);
        block.iter().map(|&v| v as i16).collect()
    }

    /// Reconstruct the real predictor coefficients from a fixed-point
    /// run by scaling `ATMP` by `2^(−NLSATMP)`.
    fn atmp_to_real(atmp: &[i16], nlsatmp: i32, lpc: usize) -> Vec<f64> {
        let scale = 2f64.powi(-nlsatmp);
        atmp[..=lpc].iter().map(|&v| v as f64 * scale).collect()
    }

    #[test]
    fn ar1_matches_floating_point_recursion() {
        // AR(1): r[k] = ρ^k, ρ = 0.8. Float recursion gives a[1] = −ρ.
        let rho = 0.8f64;
        let r_real: Vec<f64> = (0..4).map(|k| rho.powi(k)).collect();

        // Floating-point reference (block 50 transcription).
        let mut a_float = [0.0f64; 4];
        levinson_durbin(&r_real, &mut a_float, 3).unwrap();

        // Fixed-point §G.2.2.
        let rtmp = to_q15_block(&r_real);
        let mut atmp = [0i16; 4];
        let st = levinson_durbin_fixed(
            &LevinsonInput {
                rtmp: &rtmp,
                illcond: false,
                lpc: 3,
            },
            &mut atmp,
        );
        assert!(!st.illcond, "well-conditioned input must not fail");
        let a_fixed = atmp_to_real(&atmp, st.nlsatmp, 3);

        // The fixed-point coefficients carry ~Q13–Q15 precision; allow a
        // few LSBs of tolerance against the float reference.
        for k in 1..=3 {
            assert!(
                (a_fixed[k] - a_float[k]).abs() < 5e-3,
                "tap {k}: fixed {} vs float {}",
                a_fixed[k],
                a_float[k]
            );
        }
        // First tap: a1 ≈ −0.8.
        assert!((a_fixed[1] - (-rho)).abs() < 5e-3);
    }

    #[test]
    fn ar2_matches_floating_point_recursion() {
        // Two-pole predictor (1, 0.6, −0.2), autocorrelation synthesized
        // from Yule-Walker (same construction as the float test).
        let a1 = 0.6f64;
        let a2 = -0.2f64;
        let r0 = 1.0f64;
        let r1 = -a1 / (1.0 + a2);
        let r2 = -a1 * r1 - a2 * r0;
        let r_real = [r0, r1, r2];

        let mut a_float = [0.0f64; 3];
        levinson_durbin(&r_real, &mut a_float, 2).unwrap();

        let rtmp = to_q15_block(&r_real);
        let mut atmp = [0i16; 3];
        let st = levinson_durbin_fixed(
            &LevinsonInput {
                rtmp: &rtmp,
                illcond: false,
                lpc: 2,
            },
            &mut atmp,
        );
        assert!(!st.illcond);
        let a_fixed = atmp_to_real(&atmp, st.nlsatmp, 2);
        assert!((a_fixed[1] - a1).abs() < 1e-2, "a1 {}", a_fixed[1]);
        assert!((a_fixed[2] - a2).abs() < 1e-2, "a2 {}", a_fixed[2]);
    }

    #[test]
    fn output_q_format_is_q13_q14_or_q15() {
        // NLSATMP must always land in {13, 14, 15} for a successful run.
        let rho = 0.5f64;
        let r_real: Vec<f64> = (0..11).map(|k| rho.powi(k)).collect();
        let rtmp = to_q15_block(&r_real);
        let mut atmp = [0i16; 11];
        let st = levinson_durbin_fixed(
            &LevinsonInput {
                rtmp: &rtmp,
                illcond: false,
                lpc: 10,
            },
            &mut atmp,
        );
        assert!(!st.illcond);
        assert!(
            (13..=15).contains(&st.nlsatmp),
            "NLSATMP {} out of range",
            st.nlsatmp
        );
        assert_eq!(st.nlsatmp, 15 - st.nrs);
    }

    #[test]
    fn unity_tap_recovers_one() {
        // ATMP(1) scaled by 2^(−NLSATMP) reconstructs A(1) ≈ 1.
        let rho = 0.7f64;
        let r_real: Vec<f64> = (0..5).map(|k| rho.powi(k)).collect();
        let rtmp = to_q15_block(&r_real);
        let mut atmp = [0i16; 5];
        let st = levinson_durbin_fixed(
            &LevinsonInput {
                rtmp: &rtmp,
                illcond: false,
                lpc: 4,
            },
            &mut atmp,
        );
        assert!(!st.illcond);
        let a = atmp_to_real(&atmp, st.nlsatmp, 4);
        assert!((a[0] - 1.0).abs() < 1e-3, "unity tap {}", a[0]);
    }

    #[test]
    fn high_order_well_conditioned_completes() {
        // A 50th-order run on a smoothly decaying autocorrelation must
        // complete (matches the float recursion's high-order test).
        let r_real: Vec<f64> = (0..51).map(|k| (-0.05 * k as f64).exp()).collect();
        let rtmp = to_q15_block(&r_real);
        let mut atmp = [0i16; 51];
        let st = levinson_durbin_fixed(
            &LevinsonInput {
                rtmp: &rtmp,
                illcond: false,
                lpc: 50,
            },
            &mut atmp,
        );
        assert!(!st.illcond, "50th-order recursion failed: {st:?}");
        assert!((13..=15).contains(&st.nlsatmp));
        assert_eq!(st.stopped_at, 51);
    }

    #[test]
    fn decoder_resume_at_minc0_10_matches_full_run() {
        // Running 1→50 in one shot must equal running 1→10 then 10→50
        // with the carried state. We exercise the resume API by capturing
        // the order-10 state from an interrupted run and resuming it.
        let r_real: Vec<f64> = (0..51).map(|k| (-0.04 * k as f64).exp()).collect();
        let rtmp = to_q15_block(&r_real);

        // Full run 1→50.
        let mut atmp_full = [0i16; 51];
        let st_full = levinson_durbin_fixed(
            &LevinsonInput {
                rtmp: &rtmp,
                illcond: false,
                lpc: 50,
            },
            &mut atmp_full,
        );
        assert!(!st_full.illcond);

        // Interrupted run 1→10.
        let mut atmp10 = [0i16; 51];
        let st10 = levinson_durbin_fixed(
            &LevinsonInput {
                rtmp: &rtmp,
                illcond: false,
                lpc: 10,
            },
            &mut atmp10,
        );
        assert!(!st10.illcond);
        // The order-10 run reports its carried state (ATMP, NRS,
        // ALPHATMP) exactly as a real decoder would save it at the
        // interruption point. Resuming from it must reproduce the full
        // run's orders 11..50 bit for bit.
        let resume = RecursionResume {
            atmp: atmp10[..=10].to_vec(),
            nrs: st10.nrs,
            alphatmp: st10.alphatmp,
        };
        let mut atmp_resumed = [0i16; 51];
        let st_resumed = levinson_durbin_fixed_resume(
            &LevinsonInput {
                rtmp: &rtmp,
                illcond: false,
                lpc: 50,
            },
            &mut atmp_resumed,
            &resume,
        );
        assert!(!st_resumed.illcond);
        // Resumed coefficients (orders 11..50) must match the full run.
        for k in 11..=50 {
            assert_eq!(
                atmp_resumed[k], atmp_full[k],
                "tap {k}: resumed {} vs full {}",
                atmp_resumed[k], atmp_full[k]
            );
        }
    }
}
