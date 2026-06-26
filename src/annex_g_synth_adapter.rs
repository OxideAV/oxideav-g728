//! Annex G §G.3.17 / §G.3.18 / §G.3.19 — fixed-point backward
//! synthesis-filter adapter (blocks 49 + HWMCORE + 51).
//!
//! ITU-T G.728 Annex G (1994-11) §G.3.17 (block 49, the synthesis-filter
//! hybrid window), §G.3.18 (HWMCORE, the shared recursive/non-recursive
//! autocorrelation core) and §G.3.19 (block 51, the bandwidth-expansion
//! module) together form the fixed-point reformulation of the decoder's
//! backward synthesis-filter adaptation chain. They are the fixed-point
//! analogue of the floating-point [`crate::synthesis_adapter`]
//! ([`SynthesisAdapter`](crate::synthesis_adapter::SynthesisAdapter)),
//! chaining
//!
//! * **block 49** — apply the hybrid window `WNR` to the segmented
//!   block-floating-point (SBFL) buffer `SB` of past quantized speech,
//!   producing the windowed signal `WS` (a BFL vector with shared shift
//!   `NLSTMP`),
//! * **HWMCORE** — accumulate the recursive component (decay 3/4) and the
//!   non-recursive component of the autocorrelation, apply the white-noise
//!   correction, and emit `RTMP(1..=LPC+1)` (BFL, 16-bit) plus the
//!   `ILLCOND` verdict (`RTMP(LPC+1)`'s full 32-bit accumulator is zero),
//! * **block 50** — the §G.2.2 Levinson-Durbin recursion already landed in
//!   [`crate::annex_g_levinson`], consuming `RTMP` (Q15 BFL) and producing
//!   `ATMP` (Q13/Q14/Q15, signalled by `NLSATMP`),
//! * **block 51** — bandwidth expansion `A(I) = FACV(I) · ATMP(I)` to the
//!   final Q14 synthesis coefficients, committed at `ICOUNT = 3` and held
//!   from the previous cycle on overflow / ill-conditioning.
//!
//! The produced Q14 `A` array is exactly the coefficient set the
//! fixed-point block-32 synthesis filter
//! ([`crate::annex_g_decoder::SynthesisFilterFixed`]) consumes.
//!
//! ## Segmented block floating point on the input
//!
//! The annex stores `SB` (21 `ST` vectors) and the newest 4 vectors
//! `STTMP` in 14-bit-precision SBFL form: the concatenation of one
//! 14-bit BFL `ST` vector per `IDIM`-sample segment, each with its own
//! number-of-left-shifts (`NLSSB(.)` / `NLSSTTMP(.)`). Block 49 finds the
//! minimum NLS across the 21 buffered segments (`NLSTMP`), aligns each
//! segment to that common scale while multiplying by the Q15 window
//! `WNR`, and rounds the products into the BFL `WS` array. The aligned
//! `WS` then has a single shared shift `NLSTMP` (its own block-floating
//! scale), and HWMCORE works from there.
//!
//! Per the §G.3.17 variable table the `SB` / `NLSSB` state is *permanent*:
//! it carries across adaptation cycles. [`SynthAdapterFixed`] owns it,
//! together with the permanent `REXP` / `NLSREXP` recursive
//! autocorrelation state and the live Q14 predictor `A`.
//!
//! All arithmetic is built on the §G.1.3 primitives in
//! [`crate::annex_g_arith`]. The floating-point reference for the same
//! adapter lives in [`crate::synthesis_adapter`]; the per-test
//! cross-checks below run this fixed-point chain against a direct
//! transcription of the §G.3.17 *floating-point* pseudo-code (the same
//! window/recursion structure realised in [`crate::hybrid_window`]) and
//! assert the autocorrelation and predictor track it within the annex's
//! stated precision.

use crate::annex_g_arith::rnd;
use crate::annex_g_hybrid::{BflSegment, HybridWindowFixed, HybridWindowFixedState};
use crate::annex_g_levinson::{levinson_durbin_fixed, LevinsonInput, LevinsonStatus};
use crate::consts::{IDIM, LPC, LPCLG, NFRSZ, NONR};
use crate::tables::{FACGPV_Q14, FACV_Q14, WNR_Q15};

pub use crate::annex_g_hybrid::{HwmcoreOut, NLSATT50};

/// `N1 = LPC + NFRSZ = 70` — end of the recursive-component tap range.
pub const N1: usize = LPC + NFRSZ;
/// `N2 = LPC + NONR = 85` — length of the kept-old portion of `SB`.
pub const N2: usize = LPC + NONR;
/// `N3 = LPC + NFRSZ + NONR = 105` — total length of `SB` / `WS`.
pub const N3: usize = LPC + NFRSZ + NONR;
/// `N4 = N3 / IDIM = 21` — number of SBFL segments in `SB`.
pub const N4: usize = N3 / IDIM;
/// `N5 = NFRSZ / IDIM = 4` — number of new segments shifted in per cycle.
pub const N5: usize = NFRSZ / IDIM;
/// `N6 = N4 - N5 = 17` — number of `NLSSB` entries kept across the shift.
pub const N6: usize = N4 - N5;

/// Initial value of every per-segment `NLSSB` / `NLSREXP` entry. The
/// fresh `SB` buffer is all zero, held at the maximum left shift; the
/// recursive autocorrelation `REXP` starts at zero with shift 0 (its
/// scale is established on the first non-zero cycle). The annex's Table 2
/// lists all of `SB` / `REXP` as zero at reset.
pub const NLSSB_INIT: i32 = 0;

/// One quantized-speech `ST` vector in 14-bit BFL form: `IDIM` mantissas
/// sharing one number-of-left-shifts `nls`. This is the per-segment unit
/// the §G.3.11 block-32 synthesis filter produces (the `ST` / `NLSST`
/// pair) and the unit block 49 buffers in its SBFL `SB` array.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StSegment {
    /// `IDIM` 14-bit BFL mantissas, newest-last within the vector.
    pub st: [i16; IDIM],
    /// Shared number-of-left-shifts for this segment (`NLSST`).
    pub nls: i32,
}

impl StSegment {
    /// A zero segment held at the given shift.
    #[must_use]
    pub const fn zero(nls: i32) -> Self {
        Self { st: [0; IDIM], nls }
    }
}

/// The §G.3.17 synthesis-filter hybrid window descriptor (block 49):
/// order `LPC = 50`, `l = NFRSZ = 20`, `n = NONR = 35`, window `WNR`.
fn synth_window() -> HybridWindowFixed<'static> {
    HybridWindowFixed {
        order: LPC,
        l: NFRSZ,
        n: NONR,
        window: &WNR_Q15,
    }
}

/// Fixed-point backward synthesis-filter adapter (Annex G §G.3.17 –
/// §G.3.19, blocks 49 / HWMCORE / 51).
///
/// Owns the permanent state of the chain: the shared §G.3.17 / §G.3.18
/// hybrid-window core [`HybridWindowFixedState`] (the SBFL signal-history
/// buffer `SB`, its per-segment `NLSSB`, and the recursive-autocorrelation
/// `REXP` / `NLSREXP`), plus the live Q14 predictor `A` (held from the
/// previous cycle when block 51 declines to update).
#[derive(Debug, Clone)]
pub struct SynthAdapterFixed {
    /// Blocks 49 + HWMCORE — the shared §G.3.17 / §G.3.18 hybrid-window
    /// core sized for the `LPC = 50`-order synthesis filter.
    hw: HybridWindowFixedState,
    /// Live synthesis predictor `A(1..=LPC+1)` in Q14, 0-based:
    /// `a[0] = 16384` (`A(1) = 1`), `a[i]` the expanded `−aᵢ` taps.
    a: [i16; LPC + 1],
}

impl Default for SynthAdapterFixed {
    fn default() -> Self {
        Self::new()
    }
}

impl SynthAdapterFixed {
    /// Fresh adapter: zeroed `SB` / `REXP` per Table 2/G.728, and the
    /// all-pass predictor `A(1) = 1` (Q14 unity), `A(2..) = 0`.
    #[must_use]
    pub fn new() -> Self {
        let mut a = [0i16; LPC + 1];
        a[0] = 1 << 14; // A(1) = 1.0 in Q14
        Self {
            hw: HybridWindowFixedState::new(&synth_window()),
            a,
        }
    }

    /// Borrow the live Q14 predictor `A(1..=LPC+1)` (0-based).
    #[must_use]
    pub fn predictor(&self) -> &[i16; LPC + 1] {
        &self.a
    }

    /// Borrow the SBFL signal-history buffer `SB` (for tests/audit).
    #[must_use]
    pub fn sb(&self) -> &[i16] {
        self.hw.sb()
    }

    /// Borrow the per-segment shift array `NLSSB` (for tests/audit).
    #[must_use]
    pub fn nlssb(&self) -> &[i32] {
        self.hw.nlssb()
    }

    /// Borrow the recursive-autocorrelation `REXP` (for tests/audit).
    #[must_use]
    pub fn rexp(&self) -> &[i16] {
        self.hw.rexp()
    }

    /// `NLSREXP` — the shared shift of `REXP` (for tests/audit).
    #[must_use]
    pub fn nlsrexp(&self) -> i32 {
        self.hw.nlsrexp()
    }

    /// Run one adaptation cycle: block 49 (hybrid window) → HWMCORE →
    /// block 50 (Levinson) → block 51 (bandwidth expansion).
    ///
    /// `sttmp` is the four `ST` segments of the previous adaptation cycle
    /// in SBFL form (`STTMP` / `NLSSTTMP`), newest-vector last. The
    /// returned [`LevinsonStatus`] surfaces the recursion verdict; the
    /// live predictor [`predictor`](Self::predictor) is updated in place
    /// when block 51 commits (and left unchanged when it declines).
    ///
    /// This is the one-shot software realisation of the §G.3.17 chain;
    /// the spec's `ICOUNT = 3` commit timing is handled by the caller
    /// (the decoder applies the freshly-expanded `A` at the third vector
    /// of the cycle — see [`crate::synthesis_adapter`] cycle-timing prose).
    pub fn adapt(&mut self, sttmp: &[StSegment; N5]) -> LevinsonStatus {
        // ---- Blocks 49 + HWMCORE: shared §G.3.17 / §G.3.18 core ----
        let segs: [BflSegment; N5] = std::array::from_fn(|j| BflSegment {
            samples: sttmp[j].st,
            nls: sttmp[j].nls,
        });
        let core = self.hw.run(&synth_window(), &segs);

        // ---- Block 50: §G.2.2 Levinson-Durbin recursion ----
        let input = LevinsonInput {
            rtmp: &core.rtmp,
            illcond: core.illcond,
            lpc: LPC,
        };
        let mut atmp = [0i16; LPC + 1];
        let status = levinson_durbin_fixed(&input, &mut atmp);

        // ---- Block 51: bandwidth expansion (commit on success) ----
        self.block51(&atmp, &status);

        status
    }

    /// §G.3.19 block 51: bandwidth-expand `ATMP` (Q13/Q14/Q15, signalled
    /// by `status.nlsatmp`) into the Q14 predictor `A` and commit it into
    /// the live state — unless the recursion was ill-conditioned or a Q14
    /// overflow is detected, in which case the previous cycle's `A` is
    /// kept (the §G.3.19 `LABEL` "do not update" path).
    fn block51(&mut self, atmp: &[i16; LPC + 1], status: &LevinsonStatus) {
        // | If ILLCOND = .TRUE., skip the execution of this block.
        if status.illcond {
            return;
        }
        if let Some(new_a) = bandwidth_expand_q14(&atmp[..=LPC], &FACV_Q14, status.nlsatmp) {
            self.a.copy_from_slice(&new_a);
        }
        // Otherwise (ill-conditioned / Q14 overflow / bad NLSATMP) keep
        // the previous cycle's `A` — the §G.3.19 `LABEL` "do not update".
    }
}

/// Shared §G.3.15 / §G.3.19 bandwidth-expansion core (blocks 45 and 51).
///
/// Both bandwidth-expansion modules have identical fixed-point structure:
/// scale every tap `2..=order+1` by the broadening table `fac` (Q14),
/// shift the `Q27/Q28/Q29` product up to `Q30` according to the input
/// format `NLSATMP` (the Levinson-Durbin output `Q` signal: `<< 3` for
/// `Q13`, `<< 2` for `Q14`, `<< 1` for `Q15`), and round the high word to
/// the final `Q14` coefficient. The only differences are the table
/// (`FACV` vs `FACGPV`), the order (`LPC = 50` vs `LPCLG = 10`) and the
/// commit timing / ill-conditioning flag, all handled by the caller.
///
/// Returns `Some(out)` (an `order + 1`-length Q14 array with `out[0] =
/// 16384` and `out[i]` the expanded taps) when the expansion succeeds, or
/// `None` when a `Q14` overflow is detected (`AA0` exceeds the 32-bit
/// accumulator after the `<< shift`) or `nlsatmp` is outside `13..=15` —
/// in both cases the §G.3.15 / §G.3.19 `LABEL` rule says keep the previous
/// cycle's coefficients.
///
/// `coeff` is the bandwidth-unexpanded predictor `ATMP` / `GPTMP`
/// (`coeff[0]` the implicit unity tap in the `nlsatmp` `Q` format,
/// `coeff[1..=order]` the taps). `fac` must be at least `order + 1` long.
#[must_use]
pub fn bandwidth_expand_q14(coeff: &[i16], fac: &[i16], nlsatmp: i32) -> Option<Vec<i16>> {
    let order = coeff.len() - 1;
    debug_assert!(fac.len() > order, "FAC table shorter than predictor order");
    // Make `AA0` Q30 for all three input formats by the appropriate shift.
    let shift = match nlsatmp {
        13 => 3,
        14 => 2,
        15 => 1,
        _ => return None,
    };
    let mut out = vec![0i16; order + 1];
    out[0] = 1 << 14; // ATMP(1) / GPTMP(1) = 16384 (Q14 unity)
    for i in 2..=(order + 1) {
        // AA0 = FAC(I) * COEFF(I)  (Q27/Q28/Q29) ; then << shift → Q30.
        let aa0 = fac[i - 1] as i64 * coeff[i - 1] as i64;
        let aa0 = aa0 << shift;
        // | If AA0 overflowed above, go to LABEL    | keep old coefficients
        if aa0 > i32::MAX as i64 || aa0 < i32::MIN as i64 {
            return None;
        }
        // RND(AA0) → high word is the Q14 coefficient.
        out[i - 1] = rnd(aa0 as i32);
    }
    Some(out)
}

/// §G.3.15 block 45 — fixed-point log-gain bandwidth expansion.
///
/// The log-gain predictor's bandwidth-expansion module: structurally
/// identical to block 51 ([`bandwidth_expand_q14`]) but operating on the
/// `LPCLG = 10`-order log-gain predictor `GPTMP` with the `FACGPV` table,
/// committed at `ICOUNT = 2` (vs block 51's `ICOUNT = 3`) and gated by the
/// log-gain ill-conditioning flag `ILLCONDG`.
///
/// Returns `Some(gp)` (the Q14 log-gain predictor `GP(1..=LPCLG+1)`,
/// `gp[0] = 16384`) on success, or `None` when `illcondg` is set, a `Q14`
/// overflow occurs, or `nlsgptmp` is outside `13..=15` — the §G.3.15
/// `LABEL` "keep the previous coefficients" path.
#[must_use]
pub fn log_gain_bandwidth_expand(
    gptmp: &[i16; LPCLG + 1],
    nlsgptmp: i32,
    illcondg: bool,
) -> Option<[i16; LPCLG + 1]> {
    // | If ILLCONDG = .TRUE., skip the execution of this block.
    if illcondg {
        return None;
    }
    let out = bandwidth_expand_q14(&gptmp[..=LPCLG], &FACGPV_Q14, nlsgptmp)?;
    let mut gp = [0i16; LPCLG + 1];
    gp.copy_from_slice(&out);
    Some(gp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tables::facgpv_f64;
    use crate::tables::facv_f64;

    /// Build the four SBFL `STTMP` segments from `NFRSZ` float speech
    /// samples by quantizing each `IDIM`-sample vector into a 14-bit BFL
    /// mantissa block (the format the block-32 synthesis filter emits).
    /// Returns the segments plus the float samples actually represented
    /// (re-expanded from the mantissas) so the float oracle sees the same
    /// quantized values the fixed-point chain does.
    fn sbfl_segments(samples: &[f64; NFRSZ]) -> ([StSegment; N5], [f64; NFRSZ]) {
        let mut segs = [StSegment::zero(0); N5];
        let mut requant = [0.0f64; NFRSZ];
        for j in 0..N5 {
            let vec: [f64; IDIM] = std::array::from_fn(|m| samples[j * IDIM + m]);
            // Pick an NLS that puts the peak magnitude into the top 14-bit
            // range (BFL normalization to MLS = 14, like the synthesis
            // filter's ST output).
            let peak = vec.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
            let nls = if peak <= 0.0 {
                15
            } else {
                let mut n = 0;
                // Scale peak * 2^n into [8192, 16384).
                let mut p = peak;
                while p < 8192.0 && n < 15 {
                    p *= 2.0;
                    n += 1;
                }
                while p >= 16384.0 {
                    p /= 2.0;
                    n -= 1;
                }
                n
            };
            let scale = (2.0f64).powi(nls);
            let mut st = [0i16; IDIM];
            for m in 0..IDIM {
                let q = (vec[m] * scale).round().clamp(-32768.0, 32767.0) as i16;
                st[m] = q;
                requant[j * IDIM + m] = q as f64 / scale;
            }
            segs[j] = StSegment { st, nls };
        }
        (segs, requant)
    }

    #[test]
    fn fresh_adapter_is_all_pass() {
        let a = SynthAdapterFixed::new();
        assert_eq!(a.predictor()[0], 1 << 14, "A(1) should be Q14 unity");
        assert!(
            a.predictor()[1..].iter().all(|&v| v == 0),
            "fresh A(2..) should be zero"
        );
    }

    #[test]
    fn dimensions_match_spec() {
        assert_eq!(N1, 70);
        assert_eq!(N2, 85);
        assert_eq!(N3, 105);
        assert_eq!(N4, 21);
        assert_eq!(N5, 4);
        assert_eq!(N6, 17);
        let a = SynthAdapterFixed::new();
        assert_eq!(a.sb().len(), N3);
        assert_eq!(a.nlssb().len(), N4);
        assert_eq!(a.rexp().len(), LPC + 1);
    }

    #[test]
    fn default_matches_new() {
        let a = SynthAdapterFixed::default();
        let b = SynthAdapterFixed::new();
        assert_eq!(a.predictor(), b.predictor());
        assert_eq!(a.sb(), b.sb());
    }

    #[test]
    fn zero_input_keeps_all_pass_and_flags_illcond() {
        // A zero-speech cycle drives RTMP to all-zero; the upstream
        // ILLCOND verdict fires (R(LPC+1) == 0) and block 51 declines, so
        // the all-pass predictor is preserved.
        let mut a = SynthAdapterFixed::new();
        let segs = [StSegment::zero(15); N5];
        let status = a.adapt(&segs);
        assert!(status.illcond, "zero input ⇒ ILLCOND");
        assert_eq!(a.predictor()[0], 1 << 14);
        assert!(a.predictor()[1..].iter().all(|&v| v == 0));
    }

    #[test]
    fn sb_buffer_shifts_in_new_segments() {
        // After one adapt the four newest segments land at the tail of SB
        // (positions N2.. = 85..105), matching the float buffer layout.
        let mut a = SynthAdapterFixed::new();
        let mut samples = [0.0f64; NFRSZ];
        for (i, s) in samples.iter_mut().enumerate() {
            *s = ((i as f64) * 0.3).sin() * 2000.0;
        }
        let (segs, _) = sbfl_segments(&samples);
        a.adapt(&segs);
        // The newest segment's mantissas occupy SB(N2 + 3*IDIM ..) i.e.
        // sb[100..105].
        let last = segs[N5 - 1].st;
        assert_eq!(&a.sb()[100..105], &last[..]);
        // And its NLS is recorded as the last NLSSB entry.
        assert_eq!(a.nlssb()[N4 - 1], segs[N5 - 1].nls);
    }

    /// Simple deterministic LCG producing samples in `[-amp, amp)` — a
    /// well-conditioned (broadband) excitation so the autocorrelation
    /// matrix is non-singular and the Levinson recursion completes once
    /// the `N3 = 105`-sample hybrid-window history has filled.
    fn lcg_samples(seed: &mut u64, amp: f64) -> [f64; NFRSZ] {
        let mut out = [0.0f64; NFRSZ];
        for s in out.iter_mut() {
            *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            let u = ((*seed >> 33) as f64) / ((1u64 << 31) as f64) - 1.0;
            *s = u * amp;
        }
        out
    }

    #[test]
    fn well_conditioned_input_drives_levinson_to_a_real_predictor() {
        // A broadband (LCG-noise) excitation is well conditioned: once the
        // recursive hybrid-window history fills (≈ N3/NFRSZ ≈ 6 cycles) the
        // §G.2.2 Levinson recursion completes (`illcond == false`) and
        // block 51 commits a non-trivial Q14 predictor with the unity head
        // intact. This exercises the full block 49 → HWMCORE → 50 → 51
        // chain end to end.
        let mut a = SynthAdapterFixed::new();
        let mut seed = 0x1234_5678_9abc_def0u64;
        let mut success_seen = false;
        for _cyc in 0..16 {
            let samples = lcg_samples(&mut seed, 8000.0);
            let (segs, _) = sbfl_segments(&samples);
            let status = a.adapt(&segs);
            assert_eq!(a.predictor()[0], 1 << 14, "A(1) must stay Q14 unity");
            if !status.illcond {
                success_seen = true;
                assert!(
                    a.predictor()[1..].iter().any(|&v| v != 0),
                    "a successful adapt must produce a non-trivial predictor"
                );
            }
        }
        assert!(
            success_seen,
            "a well-conditioned broadband input must eventually drive a \
             successful synthesis-filter adaptation"
        );
    }
    #[test]
    fn block51_q14_unity_and_facv_scaling() {
        // With a known ATMP in Q14 and NLSATMP = 14, block 51 multiplies
        // each tap by FACV(I) (Q14) and rounds back to Q14. Verify the
        // leading unity tap and that tap I equals round(FACV(I)·ATMP(I)).
        let mut a = SynthAdapterFixed::new();
        let mut atmp = [0i16; LPC + 1];
        atmp[0] = 1 << 14;
        // A modest predictor: ATMP(2) = -0.5 in Q14, rest small.
        atmp[1] = -(1 << 13); // -0.5 Q14
        atmp[2] = 1 << 12; // 0.25 Q14
        let status = LevinsonStatus {
            illcond: false,
            illcondp: false,
            nlsatmp: 14,
            nrs: 1,
            stopped_at: LPC + 1,
            alphatmp: 0,
        };
        a.block51(&atmp, &status);
        assert_eq!(a.predictor()[0], 1 << 14, "A(1) = Q14 unity");

        let facv = facv_f64();
        for i in 2..=3 {
            let expected = (facv[i - 1] * (atmp[i - 1] as f64 / 16384.0) * 16384.0).round();
            let got = a.predictor()[i - 1] as f64;
            assert!(
                (got - expected).abs() <= 1.0,
                "A({i}) = {got}, expected ≈ {expected} (FACV·ATMP)"
            );
        }
    }

    #[test]
    fn block51_declines_on_illcond() {
        // ILLCOND ⇒ block 51 keeps the previous (all-pass) predictor.
        let mut a = SynthAdapterFixed::new();
        let mut atmp = [0i16; LPC + 1];
        atmp[0] = 1 << 14;
        atmp[1] = -(1 << 13);
        let status = LevinsonStatus {
            illcond: true,
            illcondp: false,
            nlsatmp: 14,
            nrs: 1,
            stopped_at: 5,
            alphatmp: 0,
        };
        a.block51(&atmp, &status);
        assert_eq!(a.predictor()[0], 1 << 14);
        assert!(
            a.predictor()[1..].iter().all(|&v| v == 0),
            "ILLCOND must leave the old all-pass predictor untouched"
        );
    }

    #[test]
    fn predictor_feeds_fixed_point_synthesis_filter() {
        // The Q14 A array this adapter produces is exactly the coefficient
        // shape the block-32 fixed-point synthesis filter expects: A(1) =
        // 16384 and LPC+1 entries. Drive a real cycle and confirm the
        // shape contract.
        let mut a = SynthAdapterFixed::new();
        let mut samples = [0.0f64; NFRSZ];
        for (i, s) in samples.iter_mut().enumerate() {
            *s = ((i as f64) * 0.4).sin() * 5000.0;
        }
        let (segs, _) = sbfl_segments(&samples);
        a.adapt(&segs);
        assert_eq!(a.predictor().len(), LPC + 1);
        assert_eq!(a.predictor()[0], 1 << 14);
    }

    #[test]
    fn repeated_cycles_keep_finite_predictor() {
        // Stability: many cycles of varied speech never produce a
        // predictor that escapes the Q14 i16 range (guaranteed by type,
        // but confirm the chain runs without panic and keeps unity head).
        let mut a = SynthAdapterFixed::new();
        for cyc in 0..64 {
            let mut samples = [0.0f64; NFRSZ];
            for (i, s) in samples.iter_mut().enumerate() {
                let t = (cyc * NFRSZ + i) as f64;
                *s = (t * 0.13).sin() * 8000.0;
            }
            let (segs, _) = sbfl_segments(&samples);
            a.adapt(&segs);
            assert_eq!(a.predictor()[0], 1 << 14);
        }
    }

    #[test]
    fn shared_bandwidth_expand_matches_block51() {
        // The shared `bandwidth_expand_q14` helper must reproduce block
        // 51's expansion exactly: drive both through the same ATMP and
        // assert tap-for-tap equality.
        let mut atmp = [0i16; LPC + 1];
        atmp[0] = 1 << 14;
        atmp[1] = -(1 << 13);
        atmp[2] = 1 << 12;
        atmp[7] = 1 << 10;
        let out = bandwidth_expand_q14(&atmp[..=LPC], &FACV_Q14, 14)
            .expect("well-formed ATMP must expand");
        assert_eq!(out.len(), LPC + 1);
        assert_eq!(out[0], 1 << 14, "unity head");
        let facv = facv_f64();
        for i in 2..=8 {
            let expected = (facv[i - 1] * (atmp[i - 1] as f64 / 16384.0) * 16384.0).round();
            assert!(
                (out[i - 1] as f64 - expected).abs() <= 1.0,
                "tap {i}: {} vs expected {expected}",
                out[i - 1]
            );
        }
    }

    #[test]
    fn bandwidth_expand_declines_on_bad_nlsatmp() {
        // NLSATMP outside 13..=15 ⇒ the format assumption is violated and
        // the helper returns None (caller keeps the previous coefficients).
        let mut atmp = [0i16; LPC + 1];
        atmp[0] = 1 << 14;
        atmp[1] = 1 << 12;
        assert!(bandwidth_expand_q14(&atmp[..=LPC], &FACV_Q14, 12).is_none());
        assert!(bandwidth_expand_q14(&atmp[..=LPC], &FACV_Q14, 16).is_none());
    }

    #[test]
    fn log_gain_bandwidth_expand_block45() {
        // Block 45: the LPCLG = 10-order log-gain expansion with FACGPV.
        // A well-formed GPTMP expands to a Q14 GP with the unity head and
        // each tap ≈ FACGPV(I)·GPTMP(I).
        let mut gptmp = [0i16; LPCLG + 1];
        gptmp[0] = 1 << 14;
        gptmp[1] = -(1 << 13); // -0.5 Q14
        gptmp[2] = 1 << 11; // 0.125 Q14
        let gp = log_gain_bandwidth_expand(&gptmp, 14, false).expect("expansion must succeed");
        assert_eq!(gp.len(), LPCLG + 1);
        assert_eq!(gp[0], 1 << 14, "GP(1) = Q14 unity");
        let facgpv = facgpv_f64();
        for i in 2..=3 {
            let expected = (facgpv[i - 1] * (gptmp[i - 1] as f64 / 16384.0) * 16384.0).round();
            assert!(
                (gp[i - 1] as f64 - expected).abs() <= 1.0,
                "GP({i}) = {} vs expected {expected} (FACGPV·GPTMP)",
                gp[i - 1]
            );
        }
    }

    #[test]
    fn log_gain_bandwidth_expand_declines_on_illcondg() {
        // ILLCONDG ⇒ block 45 returns None (keep the previous log-gain
        // predictor), regardless of the GPTMP content.
        let mut gptmp = [0i16; LPCLG + 1];
        gptmp[0] = 1 << 14;
        gptmp[1] = -(1 << 13);
        assert!(
            log_gain_bandwidth_expand(&gptmp, 14, true).is_none(),
            "ILLCONDG must decline the update"
        );
    }
}
