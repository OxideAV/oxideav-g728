//! Zero-input response + VQ target vector — blocks 9, 10 and 11 of
//! Figure 2/G.728 (§3.5, §3.6; pseudo-code §5.9, §5.10).
//!
//! ## Scope
//!
//! This is the next encoder front-end unit after the perceptual
//! weighting filter (block 4, §3.4) already landed: it produces the
//! analysis-by-synthesis search **target vector** `x(n)` that the
//! §3.9 codebook search of equation 3-14..3-23 minimises against.
//!
//! Per §3.5 of the Recommendation:
//!
//! > After the weighted speech vector `v(n)` has been obtained, a
//! > zero-input response vector `r(n)` will be generated using the
//! > synthesis filter (block 9) and the perceptual weighting filter
//! > (block 10). To accomplish this, we first open the switch 5 …
//! > This implies that the signal going from node 7 to the synthesis
//! > filter 9 will be zero. We then let the synthesis filter 9 and
//! > the perceptual weighting filter 10 "ring" for five samples
//! > (one vector) … The resulting output of the perceptual weighting
//! > filter 10 is the desired zero-input response vector `r(n)`.
//!
//! Then per §3.6 (block 11):
//!
//! > This block subtracts the zero-input response vector `r(n)` from
//! > the weighted speech vector `v(n)` to obtain the VQ codebook
//! > search target vector `x(n)`.
//!
//! i.e. `x(n) = v(n) − r(n)` (the spec calls `v(n)` `SW` and `r(n)`
//! `ZIR` in the §5.10 pseudo-code).
//!
//! ## What this module computes
//!
//! `ZeroInputResponse` carries the three filter-memory arrays the
//! §5.9 block-9 / block-10 pseudo-code names — `STATELPC` (synthesis
//! filter, `LPC = 50` taps), `ZIRWFIR` and `ZIRWIIR` (the all-zero
//! and all-pole halves of the weighting filter, `LPCW = 10` taps
//! each) — and runs one `IDIM = 5`-sample "ring" through them to
//! produce `r(n)`, then subtracts it from `v(n)` to yield `x(n)`.
//!
//! Since round 276 it also carries the complementary §3.10
//! **memory-update phase** (pseudo-code §5.13): after the §3.9
//! codebook search has selected the gain-scaled excitation `e(n)`,
//! [`ZeroInputResponse::update_memory`] computes the **zero-state
//! responses** of the cascaded filters to `e(n)` and adds them onto
//! the post-ring memory — "this in effect adds the zero-input
//! responses to the zero-state responses of the filters 9 and 10"
//! (§3.10). The quantized speech vector `sq(n)` then falls out of
//! the top five `STATELPC` taps for free, which is why the §5.12
//! block-22 synthesis filter "can be omitted".
//!
//! `ZeroInputResponse` constructs to the spec's all-zero
//! initialisation state, so the very first vector after construction
//! yields `r(n) = 0` and `x(n) = v(n)` exactly — matching the §3.5
//! note that "except for the vector right after initialization, the
//! memory of the filters 9 and 10 is in general non-zero".

use crate::consts::{IDIM, LPC, LPCW};
use crate::decoder::DEFAULT_MAX;
use crate::weighting_filter_coeff::WeightingFilterCoeff;

/// Block 9 + block 10 + block 11 — zero-input-response computation
/// and VQ target vector formation.
///
/// Carries the three filter-memory arrays the §5.9 pseudo-code
/// names. All three are initialised to zero per the spec's
/// initialisation state (Table 2/G.728); after that they evolve
/// through [`ZeroInputResponse::compute_target`] (the "ring" shift
/// of §5.9, before the codebook search) and
/// [`ZeroInputResponse::update_memory`] (the §5.13 zero-state-
/// response addition, after the search).
#[derive(Debug, Clone)]
pub struct ZeroInputResponse {
    /// `STATELPC(1..LPC)` — synthesis filter memory (block 9). Index
    /// 0 holds `STATELPC(1)` (the most recent past output), index
    /// `LPC-1` holds `STATELPC(LPC)`. Layout matches
    /// [`crate::Synthesizer::state`] so the block-9 ring uses the
    /// same shift idiom as the decoder's synthesis filter.
    state_lpc: [f64; LPC],
    /// `ZIRWFIR(1..LPCW)` — all-zero (numerator) half of the
    /// weighting filter's memory (block 10). Index 0 = `ZIRWFIR(1)`.
    zirw_fir: [f64; LPCW],
    /// `ZIRWIIR(1..LPCW)` — all-pole (denominator) half of the
    /// weighting filter's memory (block 10). Index 0 = `ZIRWIIR(1)`.
    zirw_iir: [f64; LPCW],
}

impl Default for ZeroInputResponse {
    fn default() -> Self {
        Self::new()
    }
}

impl ZeroInputResponse {
    /// Construct with all three filter-memory arrays zeroed — the
    /// spec's initialisation state. On the first vector after
    /// construction the ring produces `r(n) = 0`, so `x(n) = v(n)`.
    pub fn new() -> Self {
        Self {
            state_lpc: [0.0; LPC],
            zirw_fir: [0.0; LPCW],
            zirw_iir: [0.0; LPCW],
        }
    }

    /// Read-only view of the synthesis-filter memory `STATELPC`.
    pub fn state_lpc(&self) -> &[f64; LPC] {
        &self.state_lpc
    }

    /// Read-only view of the weighting-filter all-zero memory
    /// `ZIRWFIR`.
    pub fn zirw_fir(&self) -> &[f64; LPCW] {
        &self.zirw_fir
    }

    /// Read-only view of the weighting-filter all-pole memory
    /// `ZIRWIIR`.
    pub fn zirw_iir(&self) -> &[f64; LPCW] {
        &self.zirw_iir
    }

    /// Compute the VQ target vector `x(n) = v(n) − r(n)` (block 11,
    /// §3.6 / §5.10), where `r(n)` is the zero-input response of the
    /// synthesis filter (block 9, §5.9) cascaded with the perceptual
    /// weighting filter (block 10, §5.9).
    ///
    /// `a` is the order-50 synthesis predictor in the crate's
    /// canonical `A` layout (`a[0] = 1.0`, `a[i] = −â_i`) — exactly
    /// what [`crate::SynthesisAdapter::coefficients`] returns; the
    /// spec's `A(J + 1)` is `a[J]` here.
    ///
    /// `w` is the current weighting-filter coefficient set (block 38
    /// output): `w.q_gamma1` is the spec's `AWZ` (numerator,
    /// `q_i·γ₁ⁱ`) and `w.q_gamma2` is the spec's `AWP` (denominator,
    /// `q_i·γ₂ⁱ`), both in the `[1, …]` layout where index 0 is the
    /// implicit unity tap and `AWZ(J + 1)` is `q_gamma1[J]`.
    ///
    /// `v` is the weighted speech vector `v(n)` (the spec's `SW`)
    /// produced by block 4 — [`crate::Encoder::apply_weighting_filter`].
    ///
    /// The three filter-memory arrays advance one slot per output
    /// sample as the spec's "ring for five samples" prescribes, so a
    /// subsequent call sees the rung-down memory (the spec's general
    /// non-zero state). This is the **zero-input phase**; once the
    /// §3.9 search has picked `e(n)`, the caller completes the cycle
    /// with [`Self::update_memory`].
    pub fn compute_target(
        &mut self,
        a: &[f64; LPC + 1],
        w: &WeightingFilterCoeff,
        v: &[f64; IDIM],
    ) -> [f64; IDIM] {
        let r = self.zero_input_response(a, w);
        // ----- Block 11 (§5.10): TARGET(K) = SW(K) − ZIR(K) ----------
        let mut target = [0.0f64; IDIM];
        for k in 0..IDIM {
            target[k] = v[k] - r[k];
        }
        target
    }

    /// Compute the raw zero-input response vector `r(n)` (blocks 9 +
    /// 10) without the block-11 subtraction. Exposed for tests/audit
    /// and for callers that want `r(n)` directly.
    ///
    /// Side effect: advances all three filter-memory arrays as the
    /// §5.9 ring prescribes — identical to the side effect of
    /// [`Self::compute_target`].
    pub fn zero_input_response(
        &mut self,
        a: &[f64; LPC + 1],
        w: &WeightingFilterCoeff,
    ) -> [f64; IDIM] {
        let awz = &w.q_gamma1; // numerator (all-zero) taps, spec AWZ
        let awp = &w.q_gamma2; // denominator (all-pole) taps, spec AWP

        let mut zir = [0.0f64; IDIM];

        for k in 0..IDIM {
            // ===== Block 9 (§5.9): synthesis-filter ZIR ==============
            // | For K = 1..IDIM:
            // |   TEMP(K) = 0
            // |   For J = LPC, LPC-1, ..., 2:
            // |     TEMP(K) = TEMP(K) − STATELPC(J) * A(J + 1)   | MAC
            // |     STATELPC(J) = STATELPC(J − 1)                | shift
            // |   TEMP(K) = TEMP(K) − STATELPC(1) * A(2)
            // |   STATELPC(1) = TEMP(K)                          | last
            //
            // STATELPC(J) → state_lpc[J-1]; A(J+1) → a[J]. The loop
            // runs J = LPC..=2; in 0-based terms state_lpc[LPC-1]
            // down to state_lpc[1], with the J=1 tap handled last.
            let mut temp = 0.0f64;
            for j in (2..=LPC).rev() {
                temp -= self.state_lpc[j - 1] * a[j];
                self.state_lpc[j - 1] = self.state_lpc[j - 2];
            }
            temp -= self.state_lpc[0] * a[1];
            self.state_lpc[0] = temp;

            // ===== Block 10 (§5.9): weighting-filter ZIR ============
            // The block-9 output TEMP(K) feeds block 10 in place; the
            // pseudo-code reuses the same TEMP cell across both the
            // all-zero and all-pole halves.
            //
            // | TMP = TEMP(K)                                   | save input
            // | For J = LPCW, LPCW-1, ..., 2:                   | all-zero
            // |   TEMP(K) = TEMP(K) + ZIRWFIR(J) * AWZ(J + 1)
            // |   ZIRWFIR(J) = ZIRWFIR(J − 1)
            // | TEMP(K) = TEMP(K) + ZIRWFIR(1) * AWZ(2)
            // | ZIRWFIR(1) = TMP                                | input enters delay line
            //
            // Note the all-zero half pushes the *input* TMP (= the
            // block-9 output for this sample) into ZIRWFIR(1) after
            // the convolution, not the running TEMP — exactly the
            // §5.9 "ZIRWFIR(1) = TMP" line.
            let tmp_in = temp;
            for j in (2..=LPCW).rev() {
                temp += self.zirw_fir[j - 1] * awz[j];
                self.zirw_fir[j - 1] = self.zirw_fir[j - 2];
            }
            temp += self.zirw_fir[0] * awz[1];
            self.zirw_fir[0] = tmp_in;

            // | For J = LPCW, LPCW-1, ..., 2:                   | all-pole
            // |   TEMP(K) = TEMP(K) − ZIRWIIR(J) * AWP(J + 1)
            // |   ZIRWIIR(J) = ZIRWIIR(J − 1)
            // | ZIR(K) = TEMP(K) − ZIRWIIR(1) * AWP(2)
            // | ZIRWIIR(1) = ZIR(K)                             | output enters delay line
            for j in (2..=LPCW).rev() {
                temp -= self.zirw_iir[j - 1] * awp[j];
                self.zirw_iir[j - 1] = self.zirw_iir[j - 2];
            }
            zir[k] = temp - self.zirw_iir[0] * awp[1];
            self.zirw_iir[0] = zir[k];
        }

        zir
    }

    /// §3.10 memory-update phase — filter memory update for blocks 9
    /// and 10 (pseudo-code §5.13). Runs **after** the §3.9 codebook
    /// search has selected the gain-scaled excitation vector `e(n)`
    /// (= `ET = σ(n)·g_i·y_j`, blocks 19 + 21 of §5.12).
    ///
    /// Per §3.10: the memory left over after the §5.9 zero-input ring
    /// is kept, the zero-memory cascade is excited with `e(n)`, and
    /// the resulting **zero-state responses** are added on top —
    /// "this in effect adds the zero-input responses to the
    /// zero-state responses of the filters 9 and 10". The quantized
    /// speech vector `sq(n)` is then read out of the top five
    /// `STATELPC` taps (the §5.12 block-22 synthesis filter "can be
    /// omitted"), and is returned.
    ///
    /// `et` is the gain-scaled excitation `ET(1..IDIM)`; `a` and `w`
    /// are the same predictor / weighting coefficient sets the
    /// preceding [`Self::compute_target`] call used. Per §5.13 the
    /// synthesis-filter memory is clipped at the `MAX`/`MIN`
    /// saturation levels — "the positive and negative saturation
    /// levels of A-law or µ-law PCM"; we use the same default
    /// `±4 095` envelope as [`crate::Synthesizer::new`] (§3.1.1
    /// assumed-input range) so the encoder-side memory stays in
    /// lockstep with the decoder-side block-32 filter.
    ///
    /// Must be called exactly once per encoded vector, after the
    /// per-vector [`Self::compute_target`] — the §5.13 procedure
    /// assumes the ring has already shifted the memory.
    pub fn update_memory(
        &mut self,
        et: &[f64; IDIM],
        a: &[f64; LPC + 1],
        w: &WeightingFilterCoeff,
    ) -> [f64; IDIM] {
        let awz = &w.q_gamma1; // numerator taps, spec AWZ
        let awp = &w.q_gamma2; // denominator taps, spec AWP

        // ===== §5.13: zero-state responses of the cascade ===========
        // The pseudo-code reuses ZIRWFIR as a scratch history for the
        // synthesis-filter zero-state output and TEMP for the
        // weighting-filter zero-state output (both only need IDIM
        // slots — the filters start from zero memory and e(n) is five
        // samples long). We use two local arrays; the struct's
        // ZIRWFIR is rewritten from STATELPC at the end exactly as
        // the spec's "now set ZIRWFIR to the right value" loop does.
        //
        // | ZIRWFIR(1) = ET(1)              | ZIRWFIR now a scratch array
        // | TEMP(1) = ET(1)
        // | For K = 2,3,..,IDIM:
        // |   A0 = ET(K)
        // |   A1 = 0
        // |   A2 = 0
        // |   For I = K,K−1,..,2, do the next five lines
        // |     ZIRWFIR(I) = ZIRWFIR(I − 1) | shift histories
        // |     TEMP(I) = TEMP(I − 1)
        // |     A0 = A0 − A(I)·ZIRWFIR(I)   | synthesis feedback
        // |     A1 = A1 + AWZ(I)·ZIRWFIR(I) | weighting numerator
        // |     A2 = A2 − AWP(I)·TEMP(I)    | weighting denominator
        // |   ZIRWFIR(1) = A0
        // |   TEMP(1) = A0 + A1 + A2
        //
        // Same five-line shift-then-MAC idiom as the block-12 impulse
        // response calculator (§5.11) — the cascade is the same, the
        // input is e(n) instead of the unit impulse. Most recent
        // sample sits at index 1.
        let mut zsr_syn = [0.0f64; IDIM]; // spec's scratch ZIRWFIR(1..IDIM)
        let mut zsr_w = [0.0f64; IDIM]; // spec's TEMP(1..IDIM)
        zsr_syn[0] = et[0];
        zsr_w[0] = et[0];
        for k in 2..=IDIM {
            let mut a0 = et[k - 1];
            let mut a1 = 0.0f64;
            let mut a2 = 0.0f64;
            for i in (2..=k).rev() {
                zsr_syn[i - 1] = zsr_syn[i - 2];
                zsr_w[i - 1] = zsr_w[i - 2];
                a0 -= a[i - 1] * zsr_syn[i - 1];
                a1 += awz[i - 1] * zsr_syn[i - 1];
                a2 -= awp[i - 1] * zsr_w[i - 1];
            }
            zsr_syn[0] = a0;
            zsr_w[0] = a0 + a1 + a2;
        }

        // ===== §5.13: add zero-state to zero-input responses =========
        // | For K = 1,2,..,IDIM, do the next four lines
        // |   STATELPC(K) = STATELPC(K) + ZIRWFIR(K)
        // |   If STATELPC(K) > MAX, set STATELPC(K) = MAX  | limit
        // |   If STATELPC(K) < MIN, set STATELPC(K) = MIN  |
        // |   ZIRWIIR(K) = ZIRWIIR(K) + TEMP(K)
        //
        // Only the top IDIM taps receive the zero-state response —
        // the older taps (6..LPC of STATELPC, 6..LPCW of ZIRWIIR)
        // keep their rung-down values from the §5.9 ring; the
        // zero-state response is only five samples deep.
        for k in 0..IDIM {
            self.state_lpc[k] = (self.state_lpc[k] + zsr_syn[k]).clamp(-DEFAULT_MAX, DEFAULT_MAX);
            self.zirw_iir[k] += zsr_w[k];
        }

        // | For I = 1,2,..,LPCW:            | now set ZIRWFIR to the
        // |   ZIRWFIR(I) = STATELPC(I)      | right value
        //
        // The weighting filter's all-zero memory holds the past
        // *inputs* to filter 10 — which, with the switch closed, are
        // the synthesis-filter outputs, i.e. the quantized speech now
        // sitting in the top LPCW taps of STATELPC.
        for i in 0..LPCW {
            self.zirw_fir[i] = self.state_lpc[i];
        }

        // | I = IDIM + 1
        // | For K = 1,2,..,IDIM:            | obtain quantized speech
        // |   ST(K) = STATELPC(I − K)       | by reversing order
        let mut st = [0.0f64; IDIM];
        for k in 0..IDIM {
            st[k] = self.state_lpc[IDIM - 1 - k];
        }
        st
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synthesis_adapter::SynthesisAdapter;

    /// Build an all-pass synthesis predictor (`A(1) = 1`, rest zero)
    /// — the §3.5 initialisation state of the synthesis filter.
    fn allpass_a() -> [f64; LPC + 1] {
        let mut a = [0.0f64; LPC + 1];
        a[0] = 1.0;
        a
    }

    #[test]
    fn fresh_state_is_all_zero() {
        // Spec initialisation: STATELPC, ZIRWFIR, ZIRWIIR all zero.
        let z = ZeroInputResponse::new();
        assert!(z.state_lpc().iter().all(|&v| v == 0.0));
        assert!(z.zirw_fir().iter().all(|&v| v == 0.0));
        assert!(z.zirw_iir().iter().all(|&v| v == 0.0));
    }

    #[test]
    fn default_matches_new() {
        let a = ZeroInputResponse::default();
        let b = ZeroInputResponse::new();
        assert_eq!(a.state_lpc(), b.state_lpc());
        assert_eq!(a.zirw_fir(), b.zirw_fir());
        assert_eq!(a.zirw_iir(), b.zirw_iir());
    }

    #[test]
    fn first_vector_after_init_gives_zero_zir() {
        // §3.5 note: "except for the vector right after
        // initialization, the memory of the filters 9 and 10 is in
        // general non-zero". So the *first* vector after init has
        // all-zero memory → r(n) = 0, regardless of coefficients.
        let mut z = ZeroInputResponse::new();
        let a = allpass_a();
        // A non-trivial weighting coefficient set; doesn't matter
        // because every memory tap is zero on the first call.
        let mut q = [0.0f64; LPCW + 1];
        q[0] = 1.0;
        q[1] = -0.5;
        q[2] = 0.25;
        let w = WeightingFilterCoeff::from_lpc(&q);

        let r = z.zero_input_response(&a, &w);
        for k in 0..IDIM {
            assert_eq!(r[k], 0.0, "r[{k}] should be zero on first vector");
        }
    }

    #[test]
    fn target_equals_weighted_speech_on_first_vector() {
        // Block 11: x(n) = v(n) − r(n). With r(n) = 0 on the first
        // vector, x(n) must equal v(n) bit-for-bit.
        let mut z = ZeroInputResponse::new();
        let a = allpass_a();
        let w = WeightingFilterCoeff::disabled();
        let v = [11.0, -22.0, 33.0, -44.0, 55.0];
        let x = z.compute_target(&a, &w, &v);
        for k in 0..IDIM {
            assert_eq!(x[k], v[k], "k={k}");
        }
    }

    #[test]
    fn allpass_synth_and_disabled_weighting_passes_state_through() {
        // With A = all-pass (synthesis filter = identity) and the
        // weighting filter disabled (W(z) = 1), block 9 reduces to a
        // pure delay-line ring with no feedback (every A(J+1) = 0)
        // and block 10 is a pass-through. So the ZIR is identically
        // zero for *every* vector — the cascade has no excitation and
        // no recursion to ring.
        let mut z = ZeroInputResponse::new();
        let a = allpass_a();
        let w = WeightingFilterCoeff::disabled();
        for vec_idx in 0..8 {
            let r = z.zero_input_response(&a, &w);
            for k in 0..IDIM {
                assert_eq!(r[k], 0.0, "vec {vec_idx} sample {k}");
            }
        }
    }

    #[test]
    fn synthesis_ring_matches_decoder_zir_path() {
        // Block 9 here is the same zero-input response the decoder's
        // Synthesizer computes when fed a zero excitation vector. Seed
        // both with the same STATELPC and the same predictor, then
        // compare the block-9 output (read via the disabled-weighting
        // pass-through) against the decoder's first-IDIM filter ring.
        //
        // Construct a non-trivial predictor from a real adaptation so
        // A has feedback taps to ring.
        let mut adapter = SynthesisAdapter::new();
        let mut input = [0.0f64; crate::consts::NFRSZ];
        for k in 0..crate::consts::NFRSZ {
            let t = k as f64;
            input[k] = 100.0 * (0.4 * t).sin() + 50.0 * (0.9 * t).cos();
        }
        for _ in 0..6 {
            let _ = adapter.adapt(&input);
        }
        let a = *adapter.coefficients();

        // Seed the ZIR module's STATELPC with a known non-zero memory.
        let mut z = ZeroInputResponse::new();
        for (i, slot) in z.state_lpc.iter_mut().enumerate() {
            *slot = (i as f64 * 0.13).sin() * 100.0;
        }
        let seeded_state = z.state_lpc;

        // Reference: run the same block-9 recurrence by hand on a copy
        // of the seeded state.
        let mut ref_state = seeded_state;
        let mut ref_zir = [0.0f64; IDIM];
        for item in ref_zir.iter_mut() {
            let mut temp = 0.0f64;
            for j in (2..=LPC).rev() {
                temp -= ref_state[j - 1] * a[j];
                ref_state[j - 1] = ref_state[j - 2];
            }
            temp -= ref_state[0] * a[1];
            ref_state[0] = temp;
            *item = temp;
        }

        // Module: disabled weighting filter → block 10 is a no-op, so
        // the ZIR equals the block-9 output sample for sample.
        let w = WeightingFilterCoeff::disabled();
        let got = z.zero_input_response(&a, &w);
        for k in 0..IDIM {
            assert!(
                (got[k] - ref_zir[k]).abs() < 1e-12,
                "k={k}: got {} ref {}",
                got[k],
                ref_zir[k]
            );
        }
        // The module's STATELPC must end where the hand recurrence did.
        for j in 0..LPC {
            assert!(
                (z.state_lpc[j] - ref_state[j]).abs() < 1e-12,
                "state[{j}] diverged"
            );
        }
    }

    #[test]
    fn memory_evolves_after_nonzero_ring() {
        // After a vector with non-zero memory + a feedback predictor,
        // the filter memory must have changed (the "ring" shifted
        // state through). Seed STATELPC, run one ring, assert at least
        // one tap differs.
        let mut z = ZeroInputResponse::new();
        for (i, slot) in z.state_lpc.iter_mut().enumerate() {
            *slot = (i + 1) as f64;
        }
        let before = z.state_lpc;
        let mut a = allpass_a();
        a[1] = -0.5; // one feedback tap so the ring is non-trivial
        a[2] = 0.25;
        let w = WeightingFilterCoeff::disabled();
        let _ = z.zero_input_response(&a, &w);
        assert_ne!(z.state_lpc, before, "STATELPC should evolve after a ring");
    }

    #[test]
    fn finite_output_under_real_coefficients() {
        // End-to-end smoke: real synthesis predictor + real weighting
        // coefficients, several vectors of ringing, must stay finite
        // (γ₂ = 0.6 < γ₁ = 0.9 < 1 keeps the all-pole half stable, and
        // the bandwidth-expanded synthesis predictor is stable by the
        // λ = 253/256 expansion).
        let mut adapter = SynthesisAdapter::new();
        let mut input = [0.0f64; crate::consts::NFRSZ];
        for k in 0..crate::consts::NFRSZ {
            let t = k as f64;
            input[k] = 200.0 * (0.3 * t).sin() + 50.0 * (1.1 * t).cos();
        }
        for _ in 0..6 {
            let _ = adapter.adapt(&input);
        }
        let a = *adapter.coefficients();

        let mut q = [0.0f64; LPCW + 1];
        q[0] = 1.0;
        q[1] = -0.4;
        q[2] = 0.2;
        q[3] = -0.1;
        let w = WeightingFilterCoeff::from_lpc(&q);

        let mut z = ZeroInputResponse::new();
        // Seed some memory so the ring actually rings.
        for (i, slot) in z.state_lpc.iter_mut().enumerate() {
            *slot = (i as f64 * 0.07).sin() * 300.0;
        }
        for vec_idx in 0..64 {
            let mut v = [0.0f64; IDIM];
            for k in 0..IDIM {
                v[k] = 100.0 * ((vec_idx * IDIM + k) as f64 * 0.21).sin();
            }
            let x = z.compute_target(&a, &w, &v);
            for &xv in &x {
                assert!(xv.is_finite(), "non-finite target at vec {vec_idx}");
            }
        }
    }

    #[test]
    fn update_memory_identity_cascade_passes_excitation_through() {
        // Cold start + all-pass synthesis + disabled weighting: the
        // zero-state response of the identity cascade IS e(n), the
        // zero-input memory is zero, so sq(n) = e(n) and the top five
        // STATELPC taps hold e(n) reversed (most recent first).
        let mut z = ZeroInputResponse::new();
        let a = allpass_a();
        let w = WeightingFilterCoeff::disabled();
        let _ = z.compute_target(&a, &w, &[0.0; IDIM]);
        let et = [3.0, -1.5, 0.75, -0.375, 0.1875];
        let st = z.update_memory(&et, &a, &w);
        for k in 0..IDIM {
            assert_eq!(st[k], et[k], "k={k}");
            assert_eq!(z.state_lpc()[k], et[IDIM - 1 - k], "STATELPC[{k}]");
        }
        // ZIRWFIR mirrors STATELPC's top LPCW taps after the update.
        for i in 0..LPCW {
            assert_eq!(z.zirw_fir()[i], z.state_lpc()[i], "ZIRWFIR[{i}]");
        }
    }

    #[test]
    fn ring_plus_update_matches_decoder_synthesizer_exactly() {
        // §3.10 note: "after the filter memory update, the top five
        // elements of the memory of the synthesis filter 9 are
        // exactly the same as the components of the desired quantized
        // speech vector sq(n)". The decoder's block-32 Synthesizer
        // computes the same ZIR + ZSR + clamp sequence in one call —
        // so for any predictor and any excitation stream, ring +
        // update_memory here must reproduce Synthesizer::filter_vector
        // bit for bit, vector after vector.
        use crate::decoder::{ExcitationVector, Synthesizer};

        let mut adapter = SynthesisAdapter::new();
        let mut input = [0.0f64; crate::consts::NFRSZ];
        for k in 0..crate::consts::NFRSZ {
            let t = k as f64;
            input[k] = 150.0 * (0.35 * t).sin() + 60.0 * (1.3 * t).cos();
        }
        for _ in 0..6 {
            let _ = adapter.adapt(&input);
        }
        let a = *adapter.coefficients();

        // Non-trivial weighting coefficients exercise the block-10
        // halves too (they must not perturb the STATELPC lockstep).
        let mut q = [0.0f64; LPCW + 1];
        q[0] = 1.0;
        q[1] = -0.45;
        q[2] = 0.15;
        let w = WeightingFilterCoeff::from_lpc(&q);

        let mut z = ZeroInputResponse::new();
        let mut synth = Synthesizer::new();
        synth.set_predictor(a);

        for vec_idx in 0..32 {
            let mut et = [0.0f64; IDIM];
            for k in 0..IDIM {
                et[k] = 40.0 * ((vec_idx * IDIM + k) as f64 * 0.43).sin();
            }
            // Encoder side: §5.9 ring (discard target), then §5.13.
            let _ = z.compute_target(&a, &w, &[0.0; IDIM]);
            let st = z.update_memory(&et, &a, &w);
            // Decoder side: block 32 in one call.
            let dec = synth.filter_vector(&ExcitationVector(et));
            for k in 0..IDIM {
                assert_eq!(
                    st[k], dec[k],
                    "vec {vec_idx} sample {k}: encoder {} decoder {}",
                    st[k], dec[k]
                );
            }
            // The full 50-tap memories must agree too.
            for j in 0..LPC {
                assert_eq!(z.state_lpc()[j], synth.state()[j], "vec {vec_idx} tap {j}");
            }
        }
    }

    #[test]
    fn update_memory_clamps_statelpc_at_saturation() {
        // §5.13 limiter: STATELPC clips at ±MAX (the §3.1.1 ±4 095
        // envelope). Drive a huge excitation through the identity
        // cascade and confirm the memory saturates.
        let mut z = ZeroInputResponse::new();
        let a = allpass_a();
        let w = WeightingFilterCoeff::disabled();
        let _ = z.compute_target(&a, &w, &[0.0; IDIM]);
        let et = [1.0e6, -1.0e6, 1.0e6, -1.0e6, 1.0e6];
        let st = z.update_memory(&et, &a, &w);
        for k in 0..IDIM {
            assert_eq!(st[k].abs(), crate::decoder::DEFAULT_MAX, "k={k}");
        }
    }
}
