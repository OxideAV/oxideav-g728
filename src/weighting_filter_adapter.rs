//! Perceptual weighting filter adapter — encoder-side block 3 of
//! Figure 2/G.728 (the entry into the weighting-filter coefficient
//! chain of §3.3).
//!
//! The weighting filter `W(z) = Q̃(z/γ₁) / Q̃(z/γ₂)` (eq. 3-4a) is
//! computed from the order-10 LPC predictor `q_i` that the
//! coefficient calculator (block 38) consumes via
//! [`WeightingFilterCoeff::from_lpc`]. §3.3 of the Recommendation
//! lays the upstream chain that feeds block 38 with that `q_i`:
//!
//! 1. **Block 36 — hybrid window on input speech.** A 60-sample
//!    hybrid window (Annex A.3 `wnrw`, `NSBWSZ = LPCW + NONRW +
//!    NFRSZ = 10 + 30 + 20 = 60`) is applied to past input speech to
//!    produce 11 autocorrelation taps `R(1..LPCW+1)`. Same structural
//!    shape as blocks 43 and 49 — the [`crate::hybrid_window`] module
//!    transcribes the shared pseudocode and dispatches by parameter
//!    object.
//! 2. **Block 37 — Levinson-Durbin recursion.** The order-10
//!    Levinson recursion runs on the autocorrelations to produce the
//!    predictor `q_i` in canonical Levinson layout
//!    `[1.0, a_1, …, a_10]` (the polynomial `A(z) = 1 + Σ a_i z⁻ⁱ`).
//!
//! Block 38 (the substitutions `q_i ← q_i · γ₁ⁱ` / `q_i ← q_i · γ₂ⁱ`
//! of eqs. 3-4b / 3-4c) is already transcribed in
//! [`crate::weighting_filter_coeff`]; the adapter exposes a
//! [`WeightingFilterAdapter::predictor`] accessor that returns the
//! `q_i` vector this module produces so the encoder can call
//! [`crate::Encoder::set_weighting_filter_coeff_from_lpc`] without
//! reaching into private state.
//!
//! ## Cycle timing (§3.3)
//!
//! Per §3.3 the coefficient calculator (block 38) updates the
//! `W(z)` taps at the **third** vector of every 4-vector adaptation
//! cycle. The hybrid window + Levinson chain must therefore have
//! produced a fresh `q_i` by then, observing only input speech
//! through the **last vector of the previous adaptation cycle** —
//! the same 2-vector delay convention §3.7 spells out for the
//! synthesis-filter adapter (block 33), kept here so the encoder's
//! freeze-and-swap respects the spec's "do not use future samples"
//! invariant.
//!
//! In a one-shot software implementation we don't need to spread
//! the Levinson cost across vectors — but the *timing semantics*
//! still matter, because the adapter must observe the full previous
//! cycle's input speech before updating, not the current cycle's.
//! We therefore accept one buffered call per cycle (the
//! cycle-boundary delivery of `NFRSZ` samples) and produce one
//! `q_i` vector per cycle. This matches the
//! [`crate::SynthesisAdapter::adapt`] cadence exactly.
//!
//! ## Levinson failure semantics
//!
//! `levinson_durbin` flags three ill-conditioning exits: zero
//! signal (R(1) ≤ 0), unstable filter (|k_i| ≥ 1), and trailing
//! zero (R(LPCW+1) = 0). The base spec's §3.3 prose says nothing
//! explicit about block 37 failure, so we mirror block-33's policy
//! ([`crate::SynthesisAdapter::adapt`]): on Levinson failure keep
//! the previously cached `q_i`, propagate the error so the caller
//! can log/trace it, and let block 38 continue to operate on the
//! stale predictor — exactly the "keep using the previous cycle"
//! behaviour that the spec spells out for the order-50 synthesis
//! filter.

use crate::consts::{LPCW, NFRSZ, NONRW};
use crate::hybrid_window::{HybridWindow, HybridWindowState};
use crate::levinson::{levinson_durbin, LevinsonError};
use crate::tables::wnrw_f64;

/// Static (table-derived) parameters for the weighting-filter
/// hybrid window. Constructed once per adapter; passed by reference
/// to the per-cycle [`WeightingFilterAdapter::adapt`] call.
#[derive(Debug)]
struct WeightCfg {
    wnrw: [f64; 60],
}

impl WeightCfg {
    fn new() -> Self {
        Self { wnrw: wnrw_f64() }
    }
}

/// Perceptual weighting filter adapter (encoder-side block 3 of
/// Figure 2/G.728) state.
///
/// Carries the hybrid-window state (`SBW`, `REXPW` per Table 2)
/// and the most recently produced order-10 predictor `q_i` in
/// canonical Levinson layout: `predictor()[0] = 1.0`,
/// `predictor()[i]` is `a_i` of `A(z) = 1 + Σ a_i z⁻ⁱ` for
/// `1 ≤ i ≤ LPCW`. The caller passes this vector unchanged into
/// [`crate::Encoder::set_weighting_filter_coeff_from_lpc`], which
/// applies the block-38 substitution.
#[derive(Debug)]
pub struct WeightingFilterAdapter {
    cfg: WeightCfg,
    hw_state: HybridWindowState,
    /// Most recently produced order-10 predictor `q_i` in canonical
    /// Levinson `A` layout. Defaults to the all-pass predictor so
    /// the adapter is safe to read from before any speech has been
    /// pushed — block 38 with `q_1..q_10 = 0` collapses to
    /// `W(z) = 1`, which matches the §3.4 / §3.4.1 cold-start
    /// convention the existing [`crate::WeightingFilterCoeff::disabled`]
    /// constructor realises.
    last_predictor: [f64; LPCW + 1],
}

impl Default for WeightingFilterAdapter {
    fn default() -> Self {
        Self::new()
    }
}

impl WeightingFilterAdapter {
    /// Construct a fresh adapter with all internal buffers
    /// initialised to the Table 2/G.728 values (`SBW`, `REXPW`
    /// zeroed; predictor `q_0 = 1`, `q_i = 0` for `1 ≤ i ≤ LPCW`).
    pub fn new() -> Self {
        let cfg = WeightCfg::new();
        let hw = Self::window_descriptor(&cfg.wnrw);
        let hw_state = HybridWindowState::new(&hw);

        let mut last_predictor = [0.0f64; LPCW + 1];
        last_predictor[0] = 1.0;

        Self {
            cfg,
            hw_state,
            last_predictor,
        }
    }

    fn window_descriptor(wnrw: &[f64]) -> HybridWindow<'_> {
        HybridWindow {
            m: LPCW,
            l: NFRSZ,
            n: NONRW,
            window: wnrw,
            // Block 36's recursive decay is 1/2 (base spec §5.3
            // pseudo-code: "REXPW(I) = (1/2)·REXPW(I) + TMP"), unlike
            // the 3/4 of blocks 43/49.
            decay: 0.5,
        }
    }

    /// Most recently produced order-10 predictor `q_i` in canonical
    /// Levinson layout (`[1.0, a_1, …, a_10]`). Borrowed read-only
    /// so the encoder can pass it straight into
    /// [`crate::Encoder::set_weighting_filter_coeff_from_lpc`].
    pub fn predictor(&self) -> &[f64; LPCW + 1] {
        &self.last_predictor
    }

    /// Process one adaptation cycle of input-speech samples
    /// (`NFRSZ = 20` samples) through block 36 and block 37, and
    /// update the cached predictor.
    ///
    /// Returns the up-to-date predictor on success. On Levinson
    /// failure (any of the three ill-conditioning exits) the cached
    /// predictor from the previous cycle is left in place per the
    /// block-33 mirror policy (§3.3 is silent so we honour the
    /// spec's "keep previous predictor" rule from the order-50
    /// synthesis filter), and the error is propagated so the
    /// caller can log/trace it.
    pub fn adapt(&mut self, speech: &[f64]) -> Result<&[f64; LPCW + 1], LevinsonError> {
        assert_eq!(
            speech.len(),
            NFRSZ,
            "weighting-filter adapter input must be NFRSZ samples"
        );

        let hw = Self::window_descriptor(&self.cfg.wnrw);

        // ----- Block 36: hybrid window → R(1..LPCW+1) ----------------
        let mut rtmp = [0.0f64; LPCW + 1];
        self.hw_state.run(&hw, speech, &mut rtmp);

        // ----- Block 37: Levinson-Durbin → q_i -----------------------
        // The recursion writes the predictor into a fresh buffer; on
        // failure we keep the previously cached q_i. So we run into
        // a scratch buffer and only commit on Ok — exactly the same
        // pattern as block 33 ([`SynthesisAdapter::adapt`]).
        let mut qtmp = [0.0f64; LPCW + 1];
        levinson_durbin(&rtmp, &mut qtmp, LPCW)?;

        // Commit qtmp into the cached predictor (block 38 will pick
        // it up at the third vector of the cycle via
        // `Encoder::set_weighting_filter_coeff_from_lpc`; the
        // adapter itself does not gate on ICOUNT — that gating is
        // the encoder's job).
        self.last_predictor[..=LPCW].copy_from_slice(&qtmp[..=LPCW]);

        Ok(&self.last_predictor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fresh_adapter_has_allpass_predictor() {
        // Table 2 initial: q_0 = 1, q_i = 0 for 1 ≤ i ≤ LPCW. Block
        // 38 on this predictor collapses W(z) to 1 (the §3.4 /
        // §3.4.1 cold-start state), which is what the encoder
        // expects before the first adaptation cycle has run.
        let adapter = WeightingFilterAdapter::new();
        let q = adapter.predictor();
        assert_eq!(q[0], 1.0);
        assert!(q[1..].iter().all(|&v| v == 0.0));
    }

    #[test]
    fn default_matches_new() {
        // Default and explicit constructor must produce structurally
        // identical state: same all-pass predictor, same zeroed
        // hybrid-window buffers. We compare the public surface.
        let a = WeightingFilterAdapter::default();
        let b = WeightingFilterAdapter::new();
        assert_eq!(a.predictor(), b.predictor());
    }

    #[test]
    fn zero_input_fails_levinson_and_keeps_allpass() {
        // Driving zero speech vectors through the adapter produces a
        // zero autocorrelation. Block 37's R(1) ≤ 0 guard fires; per
        // the block-33 mirror policy we keep the all-pass predictor
        // and propagate the error.
        let mut adapter = WeightingFilterAdapter::new();
        let input = [0.0f64; NFRSZ];
        let result = adapter.adapt(&input);
        assert!(
            matches!(
                result,
                Err(LevinsonError::ZeroSignal | LevinsonError::TrailingZero)
            ),
            "expected Levinson to flag zero signal, got {:?}",
            result
        );
        let q = adapter.predictor();
        assert_eq!(q[0], 1.0);
        assert!(q[1..].iter().all(|&v| v == 0.0));
    }

    #[test]
    fn nonzero_input_produces_q1_eq_unity_and_nonzero_taps() {
        // A frequency-rich speech-like signal should drive the
        // autocorrelation away from zero and let Levinson produce a
        // non-trivial predictor. The q_0 = 1 invariant must always
        // hold; with 10 free taps and a non-white input at least one
        // q_i must be non-zero.
        let mut adapter = WeightingFilterAdapter::new();
        let mut input = [0.0f64; NFRSZ];
        for k in 0..NFRSZ {
            let t = k as f64;
            input[k] = 100.0 * (0.4 * t).sin() + 50.0 * (0.9 * t).cos() + (t * 13.0).fract();
        }
        // Run enough cycles to populate the recursive part of the
        // hybrid window; the synthesis adapter's test uses the same
        // "6 cycles" budget for the same reason.
        let mut last = None;
        for _ in 0..6 {
            if let Ok(q) = adapter.adapt(&input) {
                last = Some(*q);
            }
        }
        let q = last.expect("at least one adapt() should succeed on speech-like input");
        assert_eq!(q[0], 1.0, "q_0 = 1 invariant must hold");
        assert!(
            q[1..].iter().any(|&v| v != 0.0),
            "expected at least one non-zero predictor tap after adaptation"
        );
    }

    #[test]
    fn predictor_has_lpcw_plus_one_taps() {
        // §3.3 / eq. 3-3a: q_i for i = 0..=LPCW. The accessor must
        // expose all 11 (10 + leading unity) so block 38 can consume
        // every tap. Sanity guard against a future refactor narrowing
        // the layout.
        let adapter = WeightingFilterAdapter::new();
        assert_eq!(adapter.predictor().len(), LPCW + 1);
    }

    #[test]
    fn predictor_feeds_weighting_filter_coeff_round_trip() {
        // Direct integration: the q_i this adapter produces must
        // round-trip through WeightingFilterCoeff::from_lpc to a
        // valid (q_gamma1, q_gamma2) pair where both leading taps
        // are 1.0 (the implicit q_0 · γ_k^0 = 1 of eq. 3-4b /
        // 3-4c). This is the contract block 38 expects from
        // block 37.
        use crate::weighting_filter_coeff::WeightingFilterCoeff;

        let mut adapter = WeightingFilterAdapter::new();
        let mut input = [0.0f64; NFRSZ];
        for k in 0..NFRSZ {
            input[k] = 100.0 * (0.4 * k as f64).sin() + 30.0 * (0.9 * k as f64).cos();
        }
        for _ in 0..6 {
            let _ = adapter.adapt(&input);
        }
        let q = adapter.predictor();
        let coeff = WeightingFilterCoeff::from_lpc(q);
        // The implicit q_0 · γ_k^0 = 1 leading taps survive the
        // block-38 substitution by construction.
        assert_eq!(coeff.q_gamma1[0], 1.0);
        assert_eq!(coeff.q_gamma2[0], 1.0);
    }

    #[test]
    fn predictor_drives_encoder_set_lpc_path() {
        // Wire the adapter's predictor straight into the encoder's
        // block-38 setter. This is the spec's block 37 → block 38
        // handshake. The encoder must end up with non-disabled
        // (i.e. non-default) coefficients after the wiring.
        use crate::encoder::Encoder;
        use crate::weighting_filter_coeff::WeightingFilterCoeff;

        let mut adapter = WeightingFilterAdapter::new();
        let mut input = [0.0f64; NFRSZ];
        for k in 0..NFRSZ {
            let t = k as f64;
            input[k] = 200.0 * (0.3 * t).sin() + 50.0 * (1.1 * t).cos() + 7.0;
        }
        let mut got_predictor = None;
        for _ in 0..6 {
            if let Ok(q) = adapter.adapt(&input) {
                got_predictor = Some(*q);
            }
        }
        let q = got_predictor.expect("non-trivial input should land Levinson");

        let mut enc = Encoder::new();
        enc.set_weighting_filter_coeff_from_lpc(&q);

        let coeff = enc.weighting_filter_coeff();
        // After a real adaptation the coefficients can only equal
        // the disabled (all-pass) shape if q_1..q_10 happen to be
        // exactly zero — which the previous test already rules out
        // for this input.
        assert_ne!(*coeff, WeightingFilterCoeff::disabled());
    }

    #[test]
    fn predictor_q0_is_unity_after_any_number_of_cycles() {
        // The q_0 = 1 invariant (eq. 3-3a) must survive every cycle,
        // regardless of whether Levinson succeeds or fails. Drive a
        // mixed-success sequence: a few zero-vector cycles (Levinson
        // failures, predictor untouched) interleaved with non-trivial
        // cycles.
        let mut adapter = WeightingFilterAdapter::new();
        let zero = [0.0f64; NFRSZ];
        let mut sig = [0.0f64; NFRSZ];
        for k in 0..NFRSZ {
            sig[k] = 50.0 * (0.5 * k as f64).sin();
        }
        for cycle in 0..8 {
            let input = if cycle % 2 == 0 { &zero } else { &sig };
            let _ = adapter.adapt(input);
            assert_eq!(
                adapter.predictor()[0],
                1.0,
                "q_0 = 1 broken after cycle {}",
                cycle
            );
        }
    }

    #[test]
    fn hybrid_window_uses_lpcw_order_and_nonrw_tail() {
        // Block 36's hybrid-window dimensions are fixed by the
        // §3.3 spec: m = LPCW = 10, l = NFRSZ = 20, n = NONRW = 30,
        // total n3 = 60 (= len(wnrw)). Indirectly verify by
        // confirming the wnrw_f64 table is exactly 60 entries —
        // any divergence between the table size and the descriptor
        // would panic in HybridWindowState::run on the first adapt
        // call.
        let wnrw = wnrw_f64();
        assert_eq!(wnrw.len(), LPCW + NFRSZ + NONRW);
        assert_eq!(wnrw.len(), 60);
    }
}
