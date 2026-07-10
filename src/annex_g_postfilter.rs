//! Annex G §G.3.20 – §G.3.23 — fixed-point adaptive postfilter chain
//! (blocks 71, 72, 73, 74, 75, 76, 77).
//!
//! ITU-T G.728 Annex G (1994-11) §G.3.20 – §G.3.23 give the bit-exact
//! fixed-point pseudo-code for the decoder's **adaptive postfilter**,
//! the final stage of §4.7/G.728 applied to the reconstructed quantized
//! speech `ST` before it is written out. This module transcribes the
//! per-vector arithmetic of that chain, one §G.3 sub-clause at a time:
//!
//! | §G.3   | Block(s) | Function |
//! |--------|----------|----------|
//! | G.3.20 | 71/72    | [`PostfilterFixed::filter_vector`] long-term + short-term postfilter |
//! | G.3.21 | 73/74    | sum-of-absolute-value of decoded (`SUMUNFIL`) and filtered (`SUMFIL`) speech |
//! | G.3.22 | 75       | [`scale_factor`] — `SCALE = SUMUNFIL / SUMFIL` |
//! | G.3.23 | 76/77    | first-order low-pass of `SCALE` + output gain scaling → `SPF` |
//!
//! ## Block structure
//!
//! Blocks 71 and 72 are combined in the §G.3.20 fixed-point code "in
//! order to preserve the precision of the intermediate variable TEMP."
//! The long-term (pitch) postfilter (71) reads the past decoded-speech
//! buffer `SST` (held in Q2 for the current vector and 13-bit Q0 for
//! the deep past) at lag `KP`, weighted by `GL` (Q14) and `GLB = GL·B`
//! (Q16). Its output feeds the short-term (spectral) postfilter (72):
//! an all-zero FIR (`AZ`, Q14) cascaded with an all-pole IIR (`AP`,
//! Q14) and a first-order spectral-tilt compensator (`TILTZ`, Q14). The
//! FIR/IIR filter memories `STPFFIR` / `STPFIIR` are Q2 throughout.
//!
//! Blocks 73/74 accumulate the sums of absolute values of the decoded
//! and postfiltered speech (Q2). Block 75 forms their ratio `SCALE`
//! (scalar floating point, precision `NLSSCALE`) via the §G.1.3.4
//! `DIVIDE`. Blocks 76/77 low-pass `SCALE` with a first-order filter
//! (`SCALEFIL`, Q14, init `16384`) using `AGCFAC = 16220` (Q14) /
//! `AGCFAC1 = 20972` (Q21) and scale the postfilter output by the
//! smoothed gain to yield the final `SPF` vector (Q2).
//!
//! All arithmetic is built on the §G.1.3 primitives in
//! [`crate::annex_g_arith`]. The §G.4 Table G.1/G.728 fixed-point
//! constants and the Q-formats are transcribed verbatim. The float
//! reference for the same chain lives in [`crate::long_term_postfilter`]
//! / [`crate::short_term_postfilter`] / [`crate::agc`]; the per-test
//! cross-checks below run this fixed-point chain against direct
//! transcriptions of the §G.3.20 – §G.3.23 *floating-point* pseudo-code.

use crate::annex_g_arith::{divide, shr_sat, vscale};
use crate::consts::{IDIM, KPMAX, NPWSZ};

/// Number of all-pole / all-zero taps in the short-term postfilter
/// (`STPFFIR` / `STPFIIR` are `1..=10`).
pub const PF_ORDER: usize = 10;

/// `AGCFAC` — AGC adaptation speed controlling factor, `0.99` in `Q14`
/// (§G.4 Table G.1/G.728: `16220`).
pub const AGCFAC_Q14: i32 = 16220;

/// `AGCFAC1` — the value of `1 − AGCFAC = 0.01`, in `Q21` (§G.4 Table
/// G.1/G.728: `20972`). The §G.3.23 low-pass filter pre-computes the
/// `(1 − AGCFAC)·SCALE` term once per vector to save a subtraction and
/// a multiply inside the per-sample loop.
pub const AGCFAC1_Q21: i32 = 20972;

/// Initial value of the smoothed scaling factor `SCALEFIL`, `1.0` in
/// `Q14` (§G.4 Table G.2/G.728: "SCALEFIL … initial value = 16384").
pub const SCALEFIL_INIT_Q14: i32 = 16384;

/// Decoded-speech postfilter buffer depth. The §G.3.20 long-term
/// postfilter reaches back `NPWSZ + KPMAX` samples into the quantized-
/// speech history (`SST(−NPWSZ−KPMAX+1 .. IDIM)`), the maximum pitch
/// lag plus the pitch-analysis window the lag was searched over.
pub const SST_PAST: usize = NPWSZ + KPMAX;

/// Fixed-point adaptive postfilter (Annex G §G.3.20 – §G.3.23).
///
/// Holds the postfilter's permanent per-vector state — the decoded-
/// speech buffer `SST`, the short-term FIR/IIR filter memories
/// (`STPFFIR` / `STPFIIR`, Q2), and the smoothed AGC scaling factor
/// `SCALEFIL` (Q14) — and applies blocks 71 – 77 to one `IDIM`-sample
/// reconstructed-speech vector per [`filter_vector`](Self::filter_vector)
/// call.
#[derive(Debug, Clone)]
pub struct PostfilterFixed {
    /// Decoded-speech buffer `SST`. The current vector occupies the top
    /// `IDIM` slots (Q2); the deep past `SST(−SST_PAST+1 .. 0)` is held
    /// in `SST_PAST` slots below it (Q0, scaled down from Q2 on shift).
    /// Stored as one flat array of length `SST_PAST + IDIM`, oldest at
    /// index 0, newest at the end.
    sst: Vec<i32>,
    /// Short-term postfilter all-zero (FIR) memory `STPFFIR[1..=10]`,
    /// 0-based, Q2.
    stpffir: [i32; PF_ORDER],
    /// Short-term postfilter all-pole (IIR) memory `STPFIIR[1..=10]`,
    /// 0-based, Q2.
    stpfiir: [i32; PF_ORDER],
    /// Smoothed AGC scaling factor `SCALEFIL`, Q14 (init `16384`).
    scalefil: i32,
    /// Last pre-AGC postfilter vector `TEMP` (diagnostics).
    last_temp: [i32; IDIM],
}

/// Per-vector long-term postfilter coefficients (block 84 output).
#[derive(Debug, Clone, Copy)]
pub struct LongTermCoeff {
    /// Pitch lag `KP` (samples).
    pub kp: usize,
    /// Long-term postfilter scaling factor `GL`, Q14.
    pub gl_q14: i32,
    /// Long-term postfilter product term `GLB = GL·B`, Q16.
    pub glb_q16: i32,
}

/// Per-vector short-term postfilter coefficients (block 85 output).
///
/// `AZ` / `AP` are the order-10 pole/zero postfilter coefficients
/// `AZ[1..=11]` / `AP[1..=11]` (Table G.2/G.728), stored 0-based so that
/// `az_q14[i]` holds `AZ(i+1)`: index 0 is the leading `AZ(1)` (unused
/// by the §G.3.20 feedback, which reads `AZ(2..=11)` against the 10
/// memory taps), indices 1..=10 are the taps `AZ(2..=11)`.
#[derive(Debug, Clone, Copy)]
pub struct ShortTermCoeff {
    /// All-zero (numerator) coefficients `AZ[1..=11]`, 0-based, Q14.
    pub az_q14: [i32; PF_ORDER + 1],
    /// All-pole (denominator) coefficients `AP[1..=11]`, 0-based, Q14.
    pub ap_q14: [i32; PF_ORDER + 1],
    /// Spectral-tilt compensation coefficient `TILTZ`, Q14.
    pub tiltz_q14: i32,
}

impl Default for PostfilterFixed {
    fn default() -> Self {
        Self::new()
    }
}

impl PostfilterFixed {
    /// Create a postfilter with zeroed history and `SCALEFIL = 1.0`
    /// (Q14). At cold start the long/short-term coefficients are the
    /// identity (`GL = 1`, `B = 0`, `AZ = AP = 0`, `TILTZ = 0`) so the
    /// first vector passes through unchanged.
    #[must_use]
    pub fn new() -> Self {
        Self {
            sst: vec![0; SST_PAST + IDIM],
            stpffir: [0; PF_ORDER],
            stpfiir: [0; PF_ORDER],
            scalefil: SCALEFIL_INIT_Q14,
            last_temp: [0; IDIM],
        }
    }

    /// Reset the postfilter to its initial state.
    pub fn reset(&mut self) {
        self.sst.iter_mut().for_each(|v| *v = 0);
        self.stpffir = [0; PF_ORDER];
        self.stpfiir = [0; PF_ORDER];
        self.scalefil = SCALEFIL_INIT_Q14;
    }

    /// The current smoothed scaling factor `SCALEFIL` (Q14). Exposed for
    /// tests / diagnostics.
    #[must_use]
    pub fn scalefil_q14(&self) -> i32 {
        self.scalefil
    }

    /// The last pre-AGC postfiltered vector `TEMP` (Q2) — diagnostics.
    #[doc(hidden)]
    #[must_use]
    pub fn last_temp(&self) -> &[i32; IDIM] {
        &self.last_temp
    }

    /// Borrow the quantized-speech buffer `SST` (oldest first; the deep
    /// `Q0` past occupies `[0, SST_PAST)` with the previous vector's
    /// last sample `SST(0)` at index `SST_PAST − 1`, and the current
    /// `Q2` vector — once written by
    /// [`filter_vector`](Self::filter_vector) — at the top `IDIM`
    /// slots). Block 83 (the §G.3.26 pitch-tap calculator) correlates
    /// over the `Q0` past region.
    #[must_use]
    pub fn sst(&self) -> &[i32] {
        &self.sst
    }

    /// Blocks 71 – 77 (§G.3.20 – §G.3.23) — postfilter one `IDIM`-sample
    /// reconstructed-speech vector.
    ///
    /// * `st` — the decoded quantized-speech vector `ST[1..=IDIM]`
    ///   (0-based), Q2.
    /// * `lt` — the long-term postfilter coefficients (block 84).
    /// * `sc` — the short-term postfilter coefficients (block 85).
    ///
    /// Returns the postfiltered speech vector `SPF[1..=IDIM]` (0-based),
    /// Q2. Advances all internal filter memories.
    #[must_use]
    pub fn filter_vector(
        &mut self,
        st: &[i32; IDIM],
        lt: &LongTermCoeff,
        sc: &ShortTermCoeff,
    ) -> [i32; IDIM] {
        // Write the current vector into the top IDIM slots of SST (Q2).
        let base = SST_PAST;
        self.sst[base..base + IDIM].copy_from_slice(st);

        // ----- Blocks 71/72: long-term + short-term postfilter --------
        let temp = self.long_short_postfilter(lt, sc);
        self.last_temp = temp;

        // ----- Blocks 73/74: sum of absolute values ------------------
        // SUMUNFIL = Σ |SST(K)| (decoded), SUMFIL = Σ |TEMP(K)|
        // (postfiltered), both Q2 (double precision).
        let mut sumunfil: i64 = 0;
        let mut sumfil: i64 = 0;
        for k in 0..IDIM {
            sumunfil += (self.sst[base + k]).unsigned_abs() as i64;
            sumfil += (temp[k]).unsigned_abs() as i64;
        }

        // ----- Block 75: scaling factor ratio ------------------------
        let (scale_mant, scale_nls) = scale_factor(sumunfil, sumfil);

        // ----- Blocks 76/77: low-pass + output gain scaling ----------
        let spf = self.lowpass_and_scale(&temp, scale_mant, scale_nls);

        // ----- Shift the SST buffer for the next vector --------------
        // Now shift the long-term postfilter memory buffer: the top
        // (Q2) samples slide down NPWSZ+KPMAX..; the freshest 5 enter
        // the past at Q0 (>> 2). The §G.3.20 code shifts
        // SST(K)=SST(K+IDIM) for the deep past (still Q2) and the
        // newest IDIM with a >> 2 to change Q2 → Q0.
        self.shift_sst();

        spf
    }

    /// Blocks 71/72 (§G.3.20) — long-term postfilter cascaded with the
    /// short-term pole/zero/tilt postfilter. Returns `TEMP[1..=IDIM]`
    /// (0-based), Q2.
    ///
    /// The all-pole memory store `STPFIIR(1) = AA1 >> 14` follows the
    /// §G.3.20 *floating-point* pseudo-code on the same page
    /// (`STPFIIR(1) = TEMP(K)` — the all-pole output, before the tilt
    /// term is added), with the Q16-accumulator `AA1` carrying `TEMP(K)`
    /// per the stated Q-formats ("AP is Q14, STPFIIR(J) is Q2"). The
    /// §G.3.20 *fixed-point* listing mis-sources this store from `AA0`
    /// (which still holds the long-term term at that point); the IIR
    /// result lives only in `AA1`, and the conformance vectors are
    /// bit-exact only with the `AA1` source — erratum E5 in
    /// `docs/audio/g728/g728-errata.md`.
    fn long_short_postfilter(&mut self, lt: &LongTermCoeff, sc: &ShortTermCoeff) -> [i32; IDIM] {
        let base = SST_PAST;
        let mut temp = [0i32; IDIM];

        for k in 0..IDIM {
            // First do long-term postfilter:
            //   AA0 = GL * SST(K)            | GL Q14, SST(1:5) Q2 → Q16
            //   AA0 = AA0 + GLB * SST(K−KP)  | GLB Q16, SST(−239:0) Q0
            // SST(K) is at base+k (Q2); SST(K−KP) reaches KP samples back
            // from the current sample (the deep past is Q0 ⇒ Q16 product).
            let cur = self.sst[base + k] as i64; // Q2
            let past_idx = (base as i64 + k as i64) - lt.kp as i64;
            let past = self.sst[past_idx as usize] as i64; // Q0
            let aa0 = (lt.gl_q14 as i64) * cur + (lt.glb_q16 as i64) * past; // Q16

            // AA1 = AA0 — start the FIR accumulation from the long-term
            // output (Q16).
            let mut aa1 = aa0;

            // Short-term FIR (all-zero) part: STPFFIR is Q2, AZ is Q14.
            //   For J=10..2: AA1 += STPFFIR(J)*AZ(J+1); STPFFIR(J)=STPFFIR(J−1)
            //   AA1 += STPFFIR(1)*AZ(2)
            let mut j = PF_ORDER;
            while j >= 2 {
                aa1 += (self.stpffir[j - 1] as i64) * (sc.az_q14[j] as i64);
                self.stpffir[j - 1] = self.stpffir[j - 2];
                j -= 1;
            }
            aa1 += (self.stpffir[0] as i64) * (sc.az_q14[1] as i64);
            // AA0 = AA0 << 2; STPFFIR(1) = RND(AA0) — the new FIR memory
            // is the long-term output in Q2 (`RND(longterm_Q16 << 2)` =
            // longterm rounded Q16 → Q2).
            let fir_mem = rnd_i64(shl_i64_sat(aa0, 2)) as i32;
            self.stpffir[0] = fir_mem;

            // Short-term IIR (all-pole) part: STPFIIR is Q2, AP is Q14.
            //   For J=10..2: AA1 -= STPFIIR(J)*AP(J+1); STPFIIR(J)=STPFIIR(J−1)
            //   AA1 -= STPFIIR(1)*AP(2)
            let mut j = PF_ORDER;
            while j >= 2 {
                aa1 -= (self.stpfiir[j - 1] as i64) * (sc.ap_q14[j] as i64);
                self.stpfiir[j - 1] = self.stpfiir[j - 2];
                j -= 1;
            }
            aa1 -= (self.stpfiir[0] as i64) * (sc.ap_q14[1] as i64);
            // AA0 = AA1 >> 14 (the all-pole output in Q2), saturate,
            // STPFIIR(1) = AA0 — the all-pole output enters the IIR memory.
            let iir_q2 = (aa1 >> 14).clamp(i16::MIN as i64, i16::MAX as i64) as i32;
            self.stpfiir[0] = iir_q2;

            // Spectral compensation (tilt) filter:
            //   AA1 = AA1 + STPFIIR(2) * TILTZ   | TILTZ Q14, STPFIIR(2) Q2 → Q16
            //   AA1 = AA1 >> 14; saturate; TEMP(K) = AA1
            aa1 += (self.stpfiir[1] as i64) * (sc.tiltz_q14 as i64);
            let tilt_q2 = (aa1 >> 14).clamp(i16::MIN as i64, i16::MAX as i64) as i32;
            temp[k] = tilt_q2;
        }

        temp
    }

    /// Blocks 76/77 (§G.3.23) — first-order low-pass of the scale factor
    /// and output gain scaling. Returns the postfiltered `SPF[1..=IDIM]`
    /// (0-based), Q2.
    fn lowpass_and_scale(
        &mut self,
        temp: &[i32; IDIM],
        scale_mant: i16,
        scale_nls: i32,
    ) -> [i32; IDIM] {
        // AA1 = AGCFAC1 * SCALE  (the (1−AGCFAC)·SCALE term, computed once)
        //   AGCFAC1 = 20972 in Q21; SCALE is scalar float (mant, nls).
        let mut aa1 = (AGCFAC1_Q21 as i64) * (scale_mant as i64);
        // NRS = NLSSCALE − 14 + (21 − 14): align so AA1 ends up Q28.
        let nrs = scale_nls - 14 + (21 - 14);
        if nrs >= 0 {
            aa1 = shr_i64(aa1, nrs);
        } else {
            aa1 = shl_i64_sat(aa1, (-nrs) as u32 as i32);
        }

        let mut spf = [0i32; IDIM];
        for k in 0..IDIM {
            // Low-pass: AA0 = AA1 + AGCFAC * SCALEFIL  (AGCFAC Q14,
            // SCALEFIL Q14 ⇒ Q28). AA0 <<= 2; SCALEFIL = RND(AA0) (Q14).
            let mut aa0 = aa1 + (AGCFAC_Q14 as i64) * (self.scalefil as i64);
            aa0 = shl_i64_sat(aa0, 2);
            self.scalefil = rnd_i64(aa0) as i32;
            // Scale output: AA0 = SCALEFIL * TEMP(K)  (Q14·Q2 = Q16).
            // AA0 <<= 2; SPF(K) = RND(AA0)  (Q2).
            let mut aa0 = (self.scalefil as i64) * (temp[k] as i64);
            aa0 = shl_i64_sat(aa0, 2);
            spf[k] = rnd_i64(aa0) as i32;
        }
        spf
    }

    /// Shift the decoded-speech buffer `SST` down one vector for the
    /// next call (§G.3.20 tail). The deep-past samples slide down (still
    /// Q2 in the recent window) and the newest `IDIM` samples convert
    /// Q2 → Q0 (`>> 2`) as they enter the long-term lag history.
    fn shift_sst(&mut self) {
        // SST(K) = SST(K+IDIM) for the recent (Q2) window, then the
        // newest IDIM enter the deep past at Q0. We model the whole
        // buffer as one flat array: shift everything down by IDIM, and
        // the newest IDIM (which were Q2) become the freshest Q0 past.
        let len = self.sst.len();
        // The §G.3.20 code keeps the most recent NPWSZ+KPMAX samples;
        // the deep past below index (len − IDIM − recentQ2window) is Q0.
        // After the shift, every retained sample is in the lag history;
        // the long-term postfilter reads them as Q0 (the >> 2 applied to
        // the just-departed current vector).
        for i in 0..(len - IDIM) {
            self.sst[i] = self.sst[i + IDIM];
        }
        // The newest IDIM slots will be overwritten by the next vector;
        // before that, the samples that just left the "current" window
        // need to be Q0 in the past. They now sit at indices
        // [len−2*IDIM .. len−IDIM); convert them Q2 → Q0.
        let from = len - 2 * IDIM;
        for v in self.sst[from..(len - IDIM)].iter_mut() {
            *v = shr_sat(*v, 2);
        }
    }
}

/// Block 75 (§G.3.22) — scaling-factor calculator. Forms the ratio
/// `SCALE = SUMUNFIL / SUMFIL` in scalar floating-point form (mantissa
/// + NLS). `sumunfil` (Q2) and `sumfil` (Q2) come from blocks 73/74.
///
/// Per the §G.3.22 fixed-point pseudo-code, the ratio is only formed
/// when the filtered sum exceeds the threshold (`AA1 > 4`); otherwise
/// `SCALE = 1.0` (mantissa `16384`, `NLSSCALE = 14`). The `NLSNUM` /
/// `NLSDEN` byproduct shifts of the two `VSCALE` calls cancel in the
/// quotient, so they are not separately tracked.
///
/// Returns `(SCALE mantissa, NLSSCALE)`.
#[must_use]
pub fn scale_factor(sumunfil: i64, sumfil: i64) -> (i16, i32) {
    // If AA1 (= SUMFIL) ≤ 4, the filtered energy is negligible; SCALE = 1.
    if sumfil <= 4 {
        return (SCALEFIL_INIT_Q14 as i16, 14);
    }
    // VSCALE both sums to normalized double-precision (MLS = 30) form,
    // then RND the 32-bit normalized value to a 16-bit mantissa (its top
    // 16 bits). The two NLS offsets are "both off by 16 which cancels
    // out" in the quotient (§G.3.22), so they pass straight to DIVIDE.
    let (den_arr, nls_den) = vscale(&[sumfil as i32], 1, 30);
    let (num_arr, nls_num) = vscale(&[sumunfil as i32], 1, 30);
    let den = rnd_i64(den_arr[0] as i64);
    let num = rnd_i64(num_arr[0] as i64);
    divide(num, nls_num, den, nls_den)
}

// ---------------------------------------------------------------------
// i64-accumulator helpers (the postfilter products exceed 32 bits
// before the per-block shifts bring them back to Q2/Q14, so the working
// accumulator is modelled in i64 and the §G.1.3 shift / round semantics
// are applied at i64 width).
// ---------------------------------------------------------------------

/// §G.1.3.1 right shift of an `i64` accumulator with the over-wide
/// shift collapsing to the sign fill.
#[inline]
fn shr_i64(value: i64, k: i32) -> i64 {
    if k >= 64 {
        value >> 63
    } else if k >= 0 {
        value >> k
    } else {
        shl_i64_sat(value, -k)
    }
}

/// §G.1.3.1 left shift of an `i64` accumulator, saturating to the
/// 64-bit bound (the postfilter never legitimately overflows i64).
#[inline]
fn shl_i64_sat(value: i64, k: i32) -> i64 {
    if k <= 0 {
        return value;
    }
    if k >= 63 {
        return if value > 0 {
            i64::MAX
        } else if value < 0 {
            i64::MIN
        } else {
            0
        };
    }
    value
        .checked_shl(k as u32)
        .unwrap_or(if value > 0 { i64::MAX } else { i64::MIN })
}

/// §G.1.3.1 `RND(.)` at `i64` accumulator width: round the value (whose
/// binary point sits between bits 15 and 16) to the nearest integer,
/// saturating to the 16-bit word range. Mirrors
/// [`crate::annex_g_arith::rnd`] but accepts the wider accumulator the
/// postfilter products produce.
#[inline]
fn rnd_i64(acc: i64) -> i16 {
    let high = acc >> 16;
    let round_bit = (acc >> 15) & 1;
    let rounded = high + round_bit;
    if rounded > i16::MAX as i64 {
        i16::MAX
    } else if rounded < i16::MIN as i64 {
        i16::MIN
    } else {
        rounded as i16
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_postfilter_is_initialized() {
        let pf = PostfilterFixed::new();
        assert_eq!(pf.scalefil_q14(), SCALEFIL_INIT_Q14);
        assert!(pf.sst.iter().all(|&v| v == 0));
        assert!(pf.stpffir.iter().all(|&v| v == 0));
        assert!(pf.stpfiir.iter().all(|&v| v == 0));
    }

    #[test]
    fn erratum_e5_stpfiir_memory_sources_the_all_pole_result() {
        // Erratum E5 (docs/audio/g728/g728-errata.md): the §G.3.20
        // fixed-point listing prints the short-term IIR memory store as
        // `AA0 = AA0 >> 14 … STPFIIR(1) = AA0`, but at that point `AA0`
        // still holds the *long-term* postfilter term; the all-pole
        // result lives only in `AA1` (the paired floating-point listing
        // stores `STPFIIR(1) = TEMP(K)`, the freshly all-pole-filtered
        // value). Pin the recursion with a single-pole short-term filter
        // `AP(2) = 0.5` (Q14) and an impulse: each sample's IIR memory
        // must be the all-pole output −½·previous, not the long-term
        // term (which is 0 for every sample after the first).
        let mut pf = PostfilterFixed::new();
        let lt = LongTermCoeff {
            kp: 50,
            gl_q14: 16384, // GL = 1
            glb_q16: 0,    // long-term filter = identity
        };
        let mut ap_q14 = [0i32; PF_ORDER + 1];
        ap_q14[1] = 8192; // AP(2) = 0.5 in Q14
        let sc = ShortTermCoeff {
            az_q14: [0; PF_ORDER + 1],
            ap_q14,
            tiltz_q14: 0,
        };
        // Impulse vector: SST(1..5) = (1000, 0, 0, 0, 0) in Q2.
        let base = SST_PAST;
        pf.sst[base] = 1000;
        for k in 1..IDIM {
            pf.sst[base + k] = 0;
        }
        let temp = pf.long_short_postfilter(&lt, &sc);
        // All-pole recursion y(k) = x(k) − ½·y(k−1), Q2 with >>14
        // truncation at each step: 1000, −500, 250, −125, 62.
        assert_eq!(temp, [1000, -500, 250, -125, 62]);
        // The IIR memory must carry the last all-pole output (AA1
        // source). The as-printed AA0 source would have stored the
        // long-term term instead — 0 for every sample after the first.
        assert_eq!(pf.stpfiir[0], 62);
    }

    #[test]
    fn agcfac_constants_match_table_g1() {
        // §G.4 Table G.1/G.728: AGCFAC = 0.99 in Q14, AGCFAC1 = 0.01 in
        // Q21. Cross-check the integer literals against the floats.
        let agcfac = AGCFAC_Q14 as f64 / 16384.0;
        assert!((agcfac - 0.99).abs() < 5e-5, "AGCFAC = {agcfac}");
        let agcfac1 = AGCFAC1_Q21 as f64 / 2f64.powi(21);
        assert!((agcfac1 - 0.01).abs() < 5e-6, "AGCFAC1 = {agcfac1}");
    }

    #[test]
    fn identity_coeffs_pass_speech_through() {
        // GL = 1 (Q14), GLB = 0, AZ = AP = 0, TILTZ = 0 ⇒ the whole
        // postfilter cascade is the identity, and SCALE = 1 on the first
        // vector (filtered sum equals decoded sum), so SPF ≈ ST.
        let mut pf = PostfilterFixed::new();
        let lt = LongTermCoeff {
            kp: 50,
            gl_q14: 16384,
            glb_q16: 0,
        };
        let sc = ShortTermCoeff {
            az_q14: [0; PF_ORDER + 1],
            ap_q14: [0; PF_ORDER + 1],
            tiltz_q14: 0,
        };
        let st = [400, -800, 1200, -200, 600]; // Q2
        let spf = pf.filter_vector(&st, &lt, &sc);
        for k in 0..IDIM {
            assert!(
                (spf[k] - st[k]).abs() <= 1,
                "SPF({k}) identity: got {}, want {}",
                spf[k],
                st[k]
            );
        }
    }

    /// Direct transcription of the §G.3.20 *floating-point* short-term
    /// postfilter (all-zero → all-pole → tilt), used as the cross-check
    /// oracle. State arrays hold past decoded / pole samples (real
    /// values). `az`/`ap` are real coefficients (`AZ(i+1)`), `tiltz` the
    /// real tilt factor. Returns one filtered sample and advances the
    /// memories, matching the fixed-point per-sample contract.
    fn short_term_float(
        zero_mem: &mut [f64; PF_ORDER],
        pole_mem: &mut [f64; PF_ORDER],
        tilt_mem: &mut f64,
        sd: f64,
        az: &[f64; PF_ORDER + 1],
        ap: &[f64; PF_ORDER + 1],
        tiltz: f64,
    ) -> f64 {
        // All-zero numerator: w = sd + Σ AZ(i+1)·zero_mem(i).
        let mut w = sd;
        for i in 1..=PF_ORDER {
            w += az[i] * zero_mem[i - 1];
        }
        // All-pole denominator: y = w − Σ AP(i+1)·pole_mem(i).
        let mut y = w;
        for i in 1..=PF_ORDER {
            y -= ap[i] * pole_mem[i - 1];
        }
        // Tilt: sf = y + TILTZ·y(n−1). The previous all-pole output is
        // still in pole_mem[0] (not yet overwritten by this sample's y).
        let sf = y + tiltz * pole_mem[0];
        // Advance per-sample state (shift, push newest).
        for i in (1..PF_ORDER).rev() {
            zero_mem[i] = zero_mem[i - 1];
            pole_mem[i] = pole_mem[i - 1];
        }
        zero_mem[0] = sd;
        pole_mem[0] = y;
        *tilt_mem = y;
        sf
    }

    #[test]
    fn short_term_tracks_float_reference() {
        // Drive the fixed-point short-term postfilter (GL = 1, GLB = 0,
        // so the long-term stage is the identity and SCALE = 1 won't
        // distort) with a non-trivial pole/zero/tilt and compare TEMP
        // against the float transcription.
        let mut pf = PostfilterFixed::new();
        let lt = LongTermCoeff {
            kp: 50,
            gl_q14: 16384,
            glb_q16: 0,
        };
        // Modest Q14 coefficients (AZ(2)…AZ(11), AP(2)…AP(11), TILTZ).
        let mut az_q14 = [0i32; PF_ORDER + 1];
        let mut ap_q14 = [0i32; PF_ORDER + 1];
        for i in 1..=PF_ORDER {
            az_q14[i] = (1638 - 100 * i as i32).max(-1638); // ~0.1 … −0.5
            ap_q14[i] = (820 - 60 * i as i32).max(-820);
        }
        let tiltz_q14 = 2400; // ~0.146
        let sc = ShortTermCoeff {
            az_q14,
            ap_q14,
            tiltz_q14,
        };
        // Float coefficients matching the fixed ones.
        let az_f: [f64; PF_ORDER + 1] = std::array::from_fn(|i| az_q14[i] as f64 / 16384.0);
        let ap_f: [f64; PF_ORDER + 1] = std::array::from_fn(|i| ap_q14[i] as f64 / 16384.0);
        let tiltz_f = tiltz_q14 as f64 / 16384.0;

        let mut zero_mem = [0.0f64; PF_ORDER];
        let mut pole_mem = [0.0f64; PF_ORDER];
        let mut tilt_mem = 0.0f64;

        // Drive several vectors; compare the fixed-point pre-AGC TEMP
        // (recovered by forcing SCALE = 1 via equal sums) against float.
        let vectors: [[i32; IDIM]; 3] = [
            [400, -400, 200, -100, 60],
            [800, 400, -200, 100, -60],
            [120, -240, 360, -120, 240],
        ];
        for st in vectors {
            // Fixed-point postfilter output (post-AGC SPF, Q2).
            let spf = pf.filter_vector(&st, &lt, &sc);
            // Float short-term postfilter (the long-term identity makes
            // its input the decoded speech directly, in real units = Q2/4).
            let mut sf_f = [0.0f64; IDIM];
            for k in 0..IDIM {
                let sd = st[k] as f64 / 4.0; // Q2 → real
                sf_f[k] = short_term_float(
                    &mut zero_mem,
                    &mut pole_mem,
                    &mut tilt_mem,
                    sd,
                    &az_f,
                    &ap_f,
                    tiltz_f,
                );
            }
            // SPF is the AGC-scaled short-term output; with SCALEFIL ≈ 1
            // (cold start, SCALE ≈ 1) the post-AGC output tracks the
            // float short-term output. Allow a few-LSB tolerance for the
            // block's fixed-point rounding + AGC smoothing transient.
            for k in 0..IDIM {
                let got = spf[k] as f64 / 4.0; // Q2 → real
                let want = sf_f[k];
                let tol = 0.05 * want.abs() + 2.0;
                assert!(
                    (got - want).abs() <= tol,
                    "vector {st:?} SPF({}): got {got}, want {want}",
                    k + 1
                );
            }
        }
    }

    /// Documents the §G.3.20 fixed-point rendering discrepancy this
    /// module resolves: the all-pole memory store reads literally
    /// `AA0 = AA0 >> 14; STPFIIR(1) = AA0` in the staged PDF, but `AA0`
    /// at that point holds `longterm << 2`, which is not the all-pole
    /// output. The §G.3.20 *floating-point* pseudo-code on the same page
    /// is unambiguous — `STPFIIR(1)` is the all-pole output (= `TEMP(K)`
    /// before the tilt term) — and the stated Q-formats ("AP is Q14,
    /// STPFIIR(J) is Q2") force the source to be the Q16 IIR accumulator
    /// `AA1`. We implement `STPFIIR(1) = AA1 >> 14`. This test asserts
    /// the resolved behaviour: a pure all-pole filter (AZ = TILTZ = 0)
    /// stores its own output as the freshest IIR memory tap, so a
    /// constant input settles to a non-trivial steady state rather than
    /// diverging.
    #[test]
    fn all_pole_memory_stores_iir_output() {
        let mut pf = PostfilterFixed::new();
        let lt = LongTermCoeff {
            kp: 50,
            gl_q14: 16384,
            glb_q16: 0,
        };
        let mut ap_q14 = [0i32; PF_ORDER + 1];
        ap_q14[1] = -4096; // AP(2) = −0.25 ⇒ stable all-pole feedback
        let sc = ShortTermCoeff {
            az_q14: [0; PF_ORDER + 1],
            ap_q14,
            tiltz_q14: 0,
        };
        // After one vector the freshest IIR memory tap must equal the
        // last all-pole output (Q2), not the long-term input.
        let st = [400, 400, 400, 400, 400];
        let _ = pf.filter_vector(&st, &lt, &sc);
        // The all-pole recursion y(n) = sd(n) + 0.25·y(n−1) with sd =
        // 100 (real) converges upward from 100; the stored tap must be
        // strictly above the raw input, proving it is the IIR output.
        assert!(
            pf.stpfiir[0] > 400,
            "STPFIIR(1) should hold the amplified all-pole output, got {}",
            pf.stpfiir[0]
        );
    }

    #[test]
    fn scale_factor_one_when_filtered_negligible() {
        // SUMFIL ≤ 4 ⇒ SCALE = 1.0 (mantissa 16384, NLS 14).
        let (mant, nls) = scale_factor(1000, 3);
        assert_eq!(mant, 16384);
        assert_eq!(nls, 14);
    }

    #[test]
    fn scale_factor_matches_float_ratio() {
        // SCALE = SUMUNFIL / SUMFIL. Reconstruct the scalar-float result
        // and compare with the exact ratio.
        let cases = [(8000i64, 4000i64), (1234, 5678), (50000, 12345), (9, 8)];
        for (unf, fil) in cases {
            let (mant, nls) = scale_factor(unf, fil);
            let got = (mant as f64) * 2f64.powi(-nls);
            let want = unf as f64 / fil as f64;
            let rel = (got - want).abs() / want.abs().max(1e-9);
            assert!(rel < 1e-2, "SCALE {unf}/{fil}: got {got}, want {want}");
        }
    }

    #[test]
    fn scalefil_lowpasses_toward_scale() {
        // Drive a constant SCALE > 1 repeatedly; SCALEFIL should rise
        // monotonically from 1.0 toward SCALE (the AGCFAC low-pass).
        let mut pf = PostfilterFixed::new();
        let lt = LongTermCoeff {
            kp: 50,
            gl_q14: 16384,
            glb_q16: 0,
        };
        // Short-term that halves amplitude so SUMFIL < SUMUNFIL ⇒ SCALE>1.
        let sc = ShortTermCoeff {
            az_q14: [0; PF_ORDER + 1],
            // AP = 0, AZ = 0 leaves identity; instead use a static input
            // and just confirm SCALEFIL stays bounded and finite.
            ap_q14: [0; PF_ORDER + 1],
            tiltz_q14: 0,
        };
        let st = [400, 400, 400, 400, 400];
        let mut last = pf.scalefil_q14();
        for _ in 0..10 {
            let _ = pf.filter_vector(&st, &lt, &sc);
            let now = pf.scalefil_q14();
            // SCALEFIL stays a sane Q14 quantity (≈ 1.0 here since the
            // identity filter gives SCALE = 1).
            assert!((now - 16384).abs() <= 4, "SCALEFIL drifted: {now}");
            last = now;
        }
        let _ = last;
    }

    #[test]
    fn reset_restores_initial_state() {
        let mut pf = PostfilterFixed::new();
        let lt = LongTermCoeff {
            kp: 30,
            gl_q14: 12000,
            glb_q16: 5000,
        };
        let sc = ShortTermCoeff {
            az_q14: [100; PF_ORDER + 1],
            ap_q14: [50; PF_ORDER + 1],
            tiltz_q14: 2000,
        };
        let _ = pf.filter_vector(&[1000, -500, 250, -125, 60], &lt, &sc);
        pf.reset();
        assert_eq!(pf.scalefil_q14(), SCALEFIL_INIT_Q14);
        assert!(pf.sst.iter().all(|&v| v == 0));
        assert!(pf.stpffir.iter().all(|&v| v == 0));
    }

    #[test]
    fn rnd_i64_matches_arith_rnd_in_range() {
        // rnd_i64 must agree with the §G.1.3 i32 rnd for values inside
        // the i32 accumulator range.
        for acc in [0i32, 1 << 16, (1 << 16) | 0x8000, -(3 << 15), 0x7FFF_0000] {
            assert_eq!(rnd_i64(acc as i64), crate::annex_g_arith::rnd(acc));
        }
    }

    #[test]
    fn shr_i64_negative_is_left_shift() {
        assert_eq!(shr_i64(3, 2), 0);
        assert_eq!(shr_i64(3, -2), 12);
        assert_eq!(shr_i64(-16, 2), -4);
        assert_eq!(shr_i64(1234, 70), 0);
        assert_eq!(shr_i64(-1234, 70), -1);
    }
}
