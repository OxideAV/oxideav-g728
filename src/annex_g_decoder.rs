//! Annex G §G.3.11 — fixed-point decoder synthesis filter (block 32).
//!
//! ITU-T G.728 Annex G (1994-11) §G.3.11 gives the bit-exact
//! fixed-point pseudo-code for block 32, the 50th-order LPC **decoder
//! synthesis filter**. Where the coder front-end blocks (11 – 21,
//! landed in [`crate::annex_g_codebook`]) turn a target vector into a
//! 10-bit channel index, block 32 is the decoder's inverse step: it
//! drives the gain-scaled excitation vector `ET` through the backward-
//! adapted synthesis filter `1/A(z)` to reconstruct the quantized
//! speech vector `ST`.
//!
//! ## Why this is the hardest §G.3 block
//!
//! The annex states block 32 "follows the same methodology used in
//! block 9 except that there is no memory update for the perceptual
//! weighting filter." Like block 9, it splits the IIR filtering into a
//! **zero-input response** (the filter "ringing" off its stored memory
//! `STATELPC`) and a **zero-state response** (the response to the new
//! input `ET` alone), then sums them and writes the result back to the
//! filter memory.
//!
//! The fixed-point subtlety is precision management. `STATELPC`
//! (`LPC = 50` taps) is held in **segmented block floating point**: the
//! 50 taps are grouped into `IDIM = 5`-tap segments, and the array
//! `NLSSTATE[1..=10]` carries one shared left-shift count per segment
//! (`NLSSTATE[11]` is scratch for the running minimum). The synthesis
//! coefficients `A` are `Q14`; `ET` enters as a 15-bit block-floating
//! value with shift `NLSET`; the output `ST` is a 14-bit block-floating
//! value with shift `NLSST`. The pseudo-code threads the running NLS of
//! the accumulator through every multiply-add so the zero-input and
//! zero-state responses can be aligned to a common scale before being
//! summed (LABEL2), with an overflow-driven down-scaling retry on the
//! zero-state side (LABEL1).
//!
//! ## What this module computes
//!
//! [`SynthesisFilterFixed`] owns the two permanent state arrays the
//! §G.3.11 pseudo-code names — `STATELPC` (the 50-tap filter memory,
//! 14-bit precision) and `NLSSTATE[1..=11]` (its segmented NLS, init
//! 16) — and runs one `IDIM`-sample vector through the filter per call
//! to [`SynthesisFilterFixed::synthesize`]. The synthesis coefficients
//! `A[1..=LPC+1]` (Q14) are supplied per call by the backward adapter
//! (block 49/50/51); this module does not own them.
//!
//! All arithmetic is built on the §G.1.3 primitives in
//! [`crate::annex_g_arith`]. The floating-point reference for the same
//! filter lives in the decoder's synthesis path; the per-test
//! cross-checks below run this fixed-point block against a direct
//! transcription of the §G.3.11 *floating-point* pseudo-code and assert
//! the reconstructed speech tracks it within the annex's stated
//! precision.

use crate::annex_g_arith::{shl_sat, shr_sat, vscale};
use crate::consts::{IDIM, LPC};

/// Number of segments the 50-tap `STATELPC` memory is divided into for
/// segmented block-floating-point representation: `LPC / IDIM = 10`.
/// `NLSSTATE[1..=NSEG]` holds one shared left-shift count per segment;
/// `NLSSTATE[NSEG+1]` (index 11, 1-based) is scratch for the running
/// minimum used to align the zero-input-response accumulation.
pub const NSEG: usize = LPC / IDIM;

/// Initial value of every `NLSSTATE` entry (§G.4 Table G.2/G.728:
/// "NLSSTATE … initial value = 16"). At reset the synthesis filter
/// memory is all zero, and a zero block-floating-point value is held
/// with the maximum left shift.
pub const NLSSTATE_INIT: i32 = 16;

/// `STATELPC` clipping bound. The §G.3.11 pseudo-code clips the summed
/// zero-input + zero-state filter memory to `±4095` *at the segment
/// scale* (`AA1 = 4095` shifted by `NLSSTATE[10]`), mirroring the
/// floating-point range limit on `STATELPC` (the synthesis-filter
/// memory is bounded so it can never overflow a later multiplier
/// input). `4095` is the §G.3.11 constant ("`4095 = STATELPC clipping
/// level`").
pub const STATELPC_CLIP: i32 = 4095;

/// Fixed-point decoder synthesis filter (Annex G §G.3.11, block 32).
///
/// Holds the 50-tap synthesis-filter memory `STATELPC` in segmented
/// block-floating-point form together with its per-segment NLS array
/// `NLSSTATE`. One call to [`synthesize`](Self::synthesize) consumes
/// one `IDIM`-sample gain-scaled excitation vector `ET` and produces
/// the reconstructed `IDIM`-sample quantized-speech vector `ST`.
#[derive(Debug, Clone)]
pub struct SynthesisFilterFixed {
    /// Synthesis filter memory, 1-based `STATELPC[1..=LPC]` stored
    /// 0-based in `[0..LPC]`. 14-bit segmented block-floating-point
    /// mantissas; the shared shift of segment `s` (0-based) is
    /// `nlsstate[s]`.
    statelpc: [i32; LPC],
    /// Segmented NLS array, 1-based `NLSSTATE[1..=NSEG+1]` stored
    /// 0-based in `[0..=NSEG]`. `nlsstate[0..NSEG]` are the per-segment
    /// shifts; `nlsstate[NSEG]` is the scratch running minimum.
    nlsstate: [i32; NSEG + 1],
    /// `NLSST` — the NLS of the most recently produced `ST` vector
    /// (= `NLSSTATE[10]` after the filtering). Carried for the decoder
    /// stages downstream of block 32 (the LPC inverse filter of block
    /// 81 reads it).
    nlsst: i32,
}

impl Default for SynthesisFilterFixed {
    fn default() -> Self {
        Self::new()
    }
}

impl SynthesisFilterFixed {
    /// Create a synthesis filter with zeroed memory and every
    /// `NLSSTATE` entry at the §G.4 initial value 16.
    #[must_use]
    pub fn new() -> Self {
        Self {
            statelpc: [0; LPC],
            nlsstate: [NLSSTATE_INIT; NSEG + 1],
            nlsst: NLSSTATE_INIT,
        }
    }

    /// Reset the filter to its initial (all-zero memory) state.
    pub fn reset(&mut self) {
        self.statelpc = [0; LPC];
        self.nlsstate = [NLSSTATE_INIT; NSEG + 1];
        self.nlsst = NLSSTATE_INIT;
    }

    /// The NLS of the most recently produced quantized-speech vector
    /// (`NLSST`). Valid after the first [`synthesize`](Self::synthesize)
    /// call; the §G.3.24 LPC inverse filter consumes it.
    #[must_use]
    pub fn nlsst(&self) -> i32 {
        self.nlsst
    }

    /// Read the current segmented NLS array (`NLSSTATE[1..=NSEG+1]`),
    /// 0-based. Exposed for cross-stage diagnostics and tests.
    #[must_use]
    pub fn nlsstate(&self) -> &[i32; NSEG + 1] {
        &self.nlsstate
    }

    /// Block 32 (§G.3.11) — drive one `IDIM`-sample gain-scaled
    /// excitation vector through the synthesis filter `1/A(z)`.
    ///
    /// * `et` — the gain-scaled excitation vector `ET[1..=IDIM]`
    ///   (0-based), a 15-bit block-floating mantissa array.
    /// * `nlset` — the shared NLS of `et` (`NLSET`).
    /// * `a` — the synthesis-filter coefficients `A[1..=LPC+1]`
    ///   (0-based, `a[0] = A(1) = 1` in `Q14 = 16384`), in `Q14`.
    ///
    /// Returns the reconstructed quantized-speech vector `ST[1..=IDIM]`
    /// (0-based) as a 14-bit block-floating mantissa array; its shared
    /// NLS is available afterwards via [`nlsst`](Self::nlsst).
    ///
    /// Mutates the internal `STATELPC` / `NLSSTATE` memory in place, as
    /// the §G.3.11 pseudo-code does (the updated memory becomes the
    /// starting point for the next vector).
    #[must_use]
    pub fn synthesize(&mut self, et: &[i32; IDIM], nlset: i32, a: &[i32; LPC + 1]) -> [i32; IDIM] {
        // The pseudo-code is transcribed in four phases:
        //   (A) zero-input response (filter ringing), with STATELPC
        //       memory shift, accumulated per output sample K;
        //   (B) re-normalize STATELPC, refresh NLSSTATE;
        //   (C) zero-state response of ET (LABEL1, with overflow retry);
        //   (D) align (LABEL2), sum ZIR+ZSR, clip, store, VSCALE, emit.

        let mut temp = self.zero_input_response(a);
        // After phase A the per-segment STATELPC has been shifted; phase
        // B re-normalizes it to 15 bits and refreshes NLSSTATE.
        self.renormalize_state();

        // Phase C: zero-state response of ET into TEMP (overwrites it),
        // returning the NLS of the resulting ZSR (NLSET as carried).
        let nlset_zsr = self.zero_state_response(et, nlset, a, &mut temp);

        // Phase D: combine ZIR (STATELPC) + ZSR (TEMP), clip, store,
        // VSCALE, and emit ST.
        self.combine_and_emit(&mut temp, nlset_zsr)
    }

    /// Phase A — zero-input response of the synthesis filter (the
    /// §G.3.11 fixed-point "first" loop). Produces `TEMP[1..=IDIM]`
    /// (0-based) and shifts the `STATELPC` memory in place. Returns the
    /// zero-input response vector `TEMP`.
    fn zero_input_response(&mut self, a: &[i32; LPC + 1]) -> [i32; IDIM] {
        // NLSSTATE(11) = min over the per-segment NLSs.
        let mut nls_min = self.nlsstate[0];
        for &nls in &self.nlsstate[1..NSEG] {
            if nls < nls_min {
                nls_min = nls;
            }
        }
        self.nlsstate[NSEG] = nls_min;

        let mut temp = [0i32; IDIM];

        // For K = 1, 2, …, IDIM.
        for k in 1..=IDIM {
            // I = 1; L = 6 − K; J = LPC; AA0 = 0.
            // First inner loop over LL = 1..L (segment I = 1), with
            // memory shift. Segment 1 is STATELPC[1..=IDIM].
            let mut j = LPC; // 1-based index into STATELPC / A.
            let mut aa0: i64 = 0;
            let l1 = 6 - k; // number of multiply-adds for the first segment
            for _ in 0..l1 {
                // AA0 = AA0 − STATELPC(J) * A(J + 1).
                aa0 -= (self.statelpc[j - 1] as i64) * (a[j] as i64);
                // STATELPC(J) = STATELPC(J − 1).  (J ≥ 2 here since
                // L = 6 − K ≤ 5 and J starts at LPC = 50.)
                self.statelpc[j - 1] = self.statelpc[j - 2];
                j -= 1;
            }
            // NLS = NLSSTATE(I) − NLSSTATE(11); AA1 = AA0 >> NLS.
            // I = 1 (segment index, 1-based) ⇒ nlsstate[0].
            let nls = self.nlsstate[0] - self.nlsstate[NSEG];
            let mut aa1: i64 = shr_i64(aa0, nls);

            // For I = 2, …, NSEG: another IDIM multiply-adds each,
            // shifted to the segment NLS, accumulated into AA1.
            for seg in 2..=NSEG {
                aa0 = 0;
                for _ in 0..IDIM {
                    aa0 -= (self.statelpc[j - 1] as i64) * (a[j] as i64);
                    // STATELPC(0) is garbage if J = 1; it is OK (the
                    // annex notes this — the garbage tap is never read).
                    if j >= 2 {
                        self.statelpc[j - 1] = self.statelpc[j - 2];
                    }
                    j -= 1;
                }
                let nls = self.nlsstate[seg - 1] - self.nlsstate[NSEG];
                aa0 = shr_i64(aa0, nls);
                aa1 += aa0;
            }

            // If K = 1, go to SHIFT2. Otherwise do the final partial
            // segment of L = K − 1 multiply-adds for the newest segment.
            if k != 1 {
                let l_last = k - 1;
                aa0 = 0;
                for _ in 0..l_last {
                    aa0 -= (self.statelpc[j - 1] as i64) * (a[j] as i64);
                    if j >= 2 {
                        self.statelpc[j - 1] = self.statelpc[j - 2];
                    }
                    j -= 1;
                }
                // No shift necessary for this time.
                aa1 += aa0;
            }

            // SHIFT2: AA1 = AA1 >> 14 (A() was Q14, NLS of AA1 is now
            // NLSSTATE(11)). Clip to 16 bits.
            aa1 = shr_i64(aa1, 14);
            let mut aa1 = aa1.clamp(i16::MIN as i64, i16::MAX as i64) as i32;
            // STATELPC(1) = AA1 — save the new lowest-index memory word.
            self.statelpc[0] = aa1;
            // IR = NLSSTATE(11) − 2 — make TEMP Q2 format.
            let ir = self.nlsstate[NSEG] - 2;
            if ir > 0 {
                aa1 = shr_sat(aa1, ir as u32);
            } else if ir < 0 {
                aa1 = shl_sat(aa1, (-ir) as u32);
            }
            temp[k - 1] = aa1;
        }

        temp
    }

    /// Phase B — re-normalize `STATELPC` to 15 bits and refresh the
    /// segmented `NLSSTATE` array (the §G.3.11 `VSCALE(STATELPC, …, 13,
    /// …)` call plus the `NLSSTATE` shift-down).
    fn renormalize_state(&mut self) {
        // Call VSCALE(STATELPC, IDIM, IDIM, 13, STATELPC, NLS).  The
        // annex passes LEN = IDIM and SLEN = IDIM: only the *first*
        // segment (the newly-written 5 words) is re-scaled and searched,
        // its NLS adjusted, and the change propagated. MLS = 13 (15-bit
        // headroom: "Re-normalize new STATELPC to 15 bits").
        let mut first_seg: [i32; IDIM] = [0; IDIM];
        first_seg.copy_from_slice(&self.statelpc[0..IDIM]);
        let (scaled, nls) = vscale(&first_seg, IDIM, 13);
        self.statelpc[0..IDIM].copy_from_slice(&scaled);
        // NLSSTATE(11) = NLSSTATE(11) + NLS.
        self.nlsstate[NSEG] += nls;
        // For L = 1, …, NSEG: NLSSTATE(L) = NLSSTATE(L + 1) — shift the
        // segment NLSs down one position (the freshest segment inherits
        // the running-minimum-plus-rescale shift at index NSEG, and the
        // oldest segment falls off).
        for l in 0..NSEG {
            self.nlsstate[l] = self.nlsstate[l + 1];
        }
    }

    /// Phase C — zero-state response of the excitation `ET` (the
    /// §G.3.11 `LABEL1` loop). Overwrites `temp` with the zero-state
    /// response and returns its NLS (`NLSET`, possibly decreased by the
    /// overflow-retry path).
    fn zero_state_response(
        &self,
        et: &[i32; IDIM],
        nlset: i32,
        a: &[i32; LPC + 1],
        temp: &mut [i32; IDIM],
    ) -> i32 {
        // The zero-state response runs the new input ET through the
        // all-pole synthesis filter with zero initial memory. A(1) = 1
        // in Q14 = 16384, so the K = 1 term is ET(1) directly (after the
        // Q14 alignment). The loop may overflow if a sample needs 16+
        // bits; the annex's retry halves every ET sample, bumps NLSET
        // down, and restarts ("GO TO LABEL1").
        let mut nlset = nlset;
        let mut et = *et;
        loop {
            // LABEL1: TEMP(1) = ET(1).
            temp[0] = et[0];
            let mut overflowed = false;
            // For K = 2, 3, …, IDIM.
            for k in 2..=IDIM {
                // AA0 = ET(K) << 14 (because A(1) = 1 in Q14 = 16384).
                let mut aa0: i64 = (et[k - 1] as i64) << 14;
                // For I = K, K − 1, …, 2: TEMP(I) = TEMP(I − 1);
                // AA0 = AA0 − A(I) * TEMP(I).
                let mut i = k;
                while i >= 2 {
                    temp[i - 1] = temp[i - 2];
                    aa0 -= (a[i - 1] as i64) * (temp[i - 1] as i64);
                    i -= 1;
                }
                // AA1 = AA0 << 3 — overflow check: after AA0 >> 14
                // later, the result must not exceed 15 bits.
                let aa1 = shl_i64_checked(aa0);
                if aa1.overflowed {
                    overflowed = true;
                    break;
                }
                // AA0 = AA0 >> 14 (compensate for A() being Q14); keep
                // lowest 16 bits.
                let v = shr_i64(aa0, 14);
                temp[0] = (v & 0xFFFF) as i16 as i32;
            }
            if overflowed {
                // For I = 1, …, IDIM: ET(I) = ET(I) >> 1; NLSET -= 1;
                // GO TO LABEL1.
                for v in et.iter_mut() {
                    *v >>= 1;
                }
                nlset -= 1;
                continue;
            }
            return nlset;
        }
    }

    /// Phase D — align the zero-input response (`STATELPC`) and
    /// zero-state response (`temp`) to a common scale, sum them, clip,
    /// store the result back into `STATELPC`, re-normalize, and emit the
    /// reconstructed `ST` vector (the §G.3.11 `LABEL2` tail).
    fn combine_and_emit(&mut self, temp: &mut [i32; IDIM], nlset_zsr: i32) -> [i32; IDIM] {
        // NLSET here is the NLS of the zero-state response (TEMP). The
        // zero-input response sits in STATELPC[1..=IDIM] (the freshest
        // segment) at NLS = NLSSTATE[10] (0-based index NSEG-1, the
        // freshest after the phase-B shift).  Align the two to whichever
        // is the *smaller* NLS (more headroom) by right-shifting the
        // higher-NLS operand, exactly as the annex's three cascaded
        // tests do.
        let nlset = nlset_zsr;
        let nlsstate10 = self.nlsstate[NSEG - 1];

        if nlset == nlsstate10 {
            // No changes necessary.
        } else if nlset < nlsstate10 {
            // Lose precision in STATELPC by NLSD bits.
            let nlsd = nlsstate10 - nlset;
            for v in self.statelpc[0..IDIM].iter_mut() {
                *v = shr_sat(*v, nlsd as u32);
            }
            self.nlsstate[NSEG - 1] = nlset;
        } else {
            // nlset > nlsstate10: lose precision in TEMP by NLSD bits.
            let nlsd = nlset - nlsstate10;
            for v in temp.iter_mut() {
                *v = shr_sat(*v, nlsd as u32);
            }
        }

        // LABEL2: AA1 = 4095 (STATELPC clipping level) shifted to align
        // with STATELPC's current segment NLS.
        let nls10 = self.nlsstate[NSEG - 1];
        let clip = if nls10 >= 0 {
            shl_sat(STATELPC_CLIP, nls10 as u32)
        } else {
            shr_sat(STATELPC_CLIP, (-nls10) as u32)
        };

        // For K = 1, …, IDIM: STATELPC(K) = STATELPC(K) + TEMP(K),
        // clipped to ±clip then to 16 bits.
        for k in 0..IDIM {
            let mut aa0 = (self.statelpc[k] as i64) + (temp[k] as i64);
            if aa0 > clip as i64 {
                aa0 = clip as i64;
            } else if aa0 < -(clip as i64) {
                aa0 = -(clip as i64);
            }
            aa0 = aa0.clamp(i16::MIN as i64, i16::MAX as i64);
            self.statelpc[k] = aa0 as i32;
        }

        // Call VSCALE(STATELPC, IDIM, IDIM, 12, STATELPC, NLS) — scale
        // STATELPC to 14 bits to avoid overflow in the zero-input
        // response calculation later. NLSSTATE(10) += NLS.
        let mut first_seg: [i32; IDIM] = [0; IDIM];
        first_seg.copy_from_slice(&self.statelpc[0..IDIM]);
        let (scaled, nls) = vscale(&first_seg, IDIM, 12);
        self.statelpc[0..IDIM].copy_from_slice(&scaled);
        self.nlsstate[NSEG - 1] += nls;

        // Obtain quantized speech by reversing the order of the top
        // IDIM synthesis-filter memory locations: ST(K) = STATELPC(I −
        // K) with I = IDIM + 1.
        let mut st = [0i32; IDIM];
        for k in 1..=IDIM {
            st[k - 1] = self.statelpc[(IDIM + 1 - k) - 1];
        }
        // NLSST = NLSSTATE(10).
        self.nlsst = self.nlsstate[NSEG - 1];
        st
    }
}

/// Right-shift an `i64` accumulator by a possibly-negative count,
/// interpreting a negative count as a left shift, mirroring the
/// §G.1.3.1 negative-shift convention. Saturates an over-wide shift to
/// the sign fill, as the annex requires.
#[inline]
fn shr_i64(value: i64, k: i32) -> i64 {
    if k >= 0 {
        if k >= 64 {
            value >> 63
        } else {
            value >> k
        }
    } else {
        let s = (-k) as u32;
        if s >= 64 {
            // Over-wide left shift: saturate to the 32-bit accumulator
            // bound (block 32 never legitimately produces this).
            if value > 0 {
                i32::MAX as i64
            } else if value < 0 {
                i32::MIN as i64
            } else {
                0
            }
        } else {
            value << s
        }
    }
}

/// Result of the §G.3.11 `AA1 = AA0 << 3` overflow probe: the shifted
/// value (unused when it overflows) and whether the `<< 3` would have
/// pushed the value past the 32-bit accumulator bound, signalling that
/// the zero-state response needs a down-scale retry.
struct ShlProbe {
    overflowed: bool,
}

/// Compute the §G.3.11 `AA1 = AA0 << 3` overflow probe. The annex left
/// shifts by 3 to "make sure after AA0 >> 14 later, the result does not
/// exceed 15 bits"; if the shift overflows the 32-bit accumulator the
/// zero-state response is too large and must be halved.
#[inline]
fn shl_i64_checked(aa0: i64) -> ShlProbe {
    let shifted = aa0 << 3;
    // Overflow if the high bits don't agree with the sign after the
    // shift (i.e. the value left the 32-bit range).
    let overflowed = shifted > i32::MAX as i64 || shifted < i32::MIN as i64;
    ShlProbe { overflowed }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Q14 unit (A(1) = 1).
    const Q14_ONE: i32 = 16384;

    /// Direct transcription of the §G.3.11 *floating-point* block-32
    /// pseudo-code, used as the cross-check oracle for the fixed-point
    /// path. `statelpc` is the 50-tap memory (real values), `et` the
    /// excitation (real values), `a` the Q14-scaled coefficients
    /// expressed as real numbers (`a_real[i] = A(i+1)/16384`). Returns
    /// the reconstructed speech vector `st` (real values) and leaves
    /// `statelpc` updated, matching the fixed-point memory contract.
    fn block32_float(
        statelpc: &mut [f64; LPC],
        et: &[f64; IDIM],
        a_real: &[f64; LPC + 1],
    ) -> [f64; IDIM] {
        // Zero-input response with memory shift.
        let mut temp = [0.0f64; IDIM];
        for k in 1..=IDIM {
            let mut acc = 0.0;
            // For J = LPC, LPC−1, …, 3, 2.
            let mut j = LPC;
            while j >= 2 {
                acc -= statelpc[j - 1] * a_real[j];
                statelpc[j - 1] = statelpc[j - 2];
                j -= 1;
            }
            // Handle last one differently: J = 1.
            acc -= statelpc[0] * a_real[1];
            statelpc[0] = acc;
            temp[k - 1] = acc;
            let _ = k;
        }
        // Zero-state response of ET (zero initial memory).
        let mut tmp2 = [0.0f64; IDIM];
        tmp2[0] = et[0];
        for k in 2..=IDIM {
            let mut a0 = et[k - 1];
            let mut i = k;
            while i >= 2 {
                tmp2[i - 1] = tmp2[i - 2];
                a0 -= a_real[i - 1] * tmp2[i - 1];
                i -= 1;
            }
            tmp2[0] = a0;
        }
        // Update filter memory by adding zero-state to zero-input.
        for k in 0..IDIM {
            statelpc[k] += tmp2[k];
            // Range limit (clip) as in the float spec.
            let bound = STATELPC_CLIP as f64;
            if statelpc[k] > bound {
                statelpc[k] = bound;
            } else if statelpc[k] < -bound {
                statelpc[k] = -bound;
            }
        }
        // Obtain quantized speech by reversing the top IDIM memory
        // locations.
        let mut st = [0.0f64; IDIM];
        for k in 1..=IDIM {
            st[k - 1] = statelpc[(IDIM + 1 - k) - 1];
        }
        st
    }

    /// Reconstruct the real value of a block-floating ST sample.
    fn st_value(mantissa: i32, nlsst: i32) -> f64 {
        (mantissa as f64) * 2f64.powi(-nlsst)
    }

    #[test]
    fn new_filter_starts_zeroed() {
        let f = SynthesisFilterFixed::new();
        assert!(f.statelpc.iter().all(|&v| v == 0));
        assert!(f.nlsstate.iter().all(|&v| v == NLSSTATE_INIT));
        assert_eq!(f.nlsst(), NLSSTATE_INIT);
    }

    #[test]
    fn nseg_is_ten() {
        // LPC = 50, IDIM = 5 ⇒ 10 segments of segmented block floating
        // point, with an 11th scratch NLS slot.
        assert_eq!(NSEG, 10);
        assert_eq!(LPC / IDIM, NSEG);
    }

    #[test]
    fn zero_input_zero_excitation_stays_zero() {
        // A(1) = 1, all other coefficients zero: with zero memory and
        // zero excitation the reconstructed speech must be zero and the
        // memory must stay zero.
        let mut f = SynthesisFilterFixed::new();
        let mut a = [0i32; LPC + 1];
        a[0] = Q14_ONE;
        let st = f.synthesize(&[0; IDIM], 0, &a);
        assert!(st.iter().all(|&v| v == 0));
        assert!(f.statelpc.iter().all(|&v| v == 0));
    }

    #[test]
    fn pure_gain_passes_excitation_through() {
        // With A(z) = 1 (only A(1) = 1, all taps zero), the synthesis
        // filter 1/A(z) is the identity: ST should equal ET (up to the
        // block-floating scale).  The zero-state response of the
        // identity filter just stacks the input, and the §G.3.11 output
        // reversal of the filter memory unstacks it, so ST(K) = ET(K)
        // after one vector from zero memory.
        let mut f = SynthesisFilterFixed::new();
        let mut a = [0i32; LPC + 1];
        a[0] = Q14_ONE;
        // ET in a modest Q2-ish block-floating form (NLSET = 2).
        let et = [40, -80, 120, -20, 60];
        let nlset = 2;
        let st = f.synthesize(&et, nlset, &a);
        let nlsst = f.nlsst();
        // Reconstructed real values: ST(K) ≈ ET(K) / 2^2.
        for k in 1..=IDIM {
            let want = (et[k - 1] as f64) / 4.0;
            let got = st_value(st[k - 1], nlsst);
            assert!(
                (got - want).abs() <= 0.5,
                "ST({k}) identity: got {got}, want {want}"
            );
        }
    }

    #[test]
    fn tracks_float_reference_first_order_filter() {
        // A first-order all-pole filter A(z) = 1 + a1 z^-1, driven by a
        // small excitation, with zero initial memory. Compare the
        // fixed-point reconstruction against the float pseudo-code.
        let a1_q14 = 8192; // 0.5 in Q14
        let mut a_fix = [0i32; LPC + 1];
        a_fix[0] = Q14_ONE;
        a_fix[1] = a1_q14;
        let mut a_real = [0.0f64; LPC + 1];
        a_real[0] = 1.0;
        a_real[1] = a1_q14 as f64 / 16384.0;

        let mut f = SynthesisFilterFixed::new();
        let mut statelpc_f = [0.0f64; LPC];

        // Drive several vectors and check tracking each time. Use a
        // larger excitation scale (NLSET = 0) so the block-floating
        // quantum is a meaningful fraction of the value, keeping the
        // tolerance check honest.
        let vectors: [[i32; IDIM]; 3] = [
            [4000, 0, 0, 0, 0],
            [2000, -2000, 1000, -500, 250],
            [0, 0, 4000, 0, 0],
        ];
        let nlset = 0;
        for et in vectors {
            let st = f.synthesize(&et, nlset, &a_fix);
            let nlsst = f.nlsst();
            let et_f: [f64; IDIM] = std::array::from_fn(|i| (et[i] as f64) / 2f64.powi(nlset));
            let st_f = block32_float(&mut statelpc_f, &et_f, &a_real);
            for k in 0..IDIM {
                let got = st_value(st[k], nlsst);
                let want = st_f[k];
                // Block-floating reconstruction tracks the float
                // pseudo-code to ~0.5 % plus a few LSB at the working
                // scale.
                let tol = 0.01 * want.abs() + 4.0;
                assert!(
                    (got - want).abs() <= tol,
                    "vector {et:?} ST({}): got {got}, want {want}",
                    k + 1
                );
            }
        }
    }

    #[test]
    fn reset_clears_state() {
        let mut f = SynthesisFilterFixed::new();
        let mut a = [0i32; LPC + 1];
        a[0] = Q14_ONE;
        a[1] = 4096;
        let _ = f.synthesize(&[100, -50, 25, -10, 5], 2, &a);
        assert!(f.statelpc.iter().any(|&v| v != 0));
        f.reset();
        assert!(f.statelpc.iter().all(|&v| v == 0));
        assert!(f.nlsstate.iter().all(|&v| v == NLSSTATE_INIT));
    }

    #[test]
    fn shr_i64_negative_count_is_left_shift() {
        assert_eq!(shr_i64(3, 2), 0);
        assert_eq!(shr_i64(3, -2), 12);
        assert_eq!(shr_i64(-8, 1), -4);
        // Over-wide right shift collapses to sign fill.
        assert_eq!(shr_i64(12345, 70), 0);
        assert_eq!(shr_i64(-12345, 70), -1);
    }
}
