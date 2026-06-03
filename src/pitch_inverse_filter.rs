//! Pitch-extractor LPC inverse filter — block 81 of Figure 7/G.728
//! (§4.7 pitch-extraction front end).
//!
//! Implements the 10th-order LPC inverse filter `Ã(z)` of equation 4-6:
//!
//! ```text
//!     Ã(z)  =  1  −  Σ_{i=1..10} ã_i · z^{-i}
//! ```
//!
//! applied to the decoded-speech stream `sd(k)`. Per the decode-trace
//! §7.1 "pitch extraction" bullet list this stage is the **first** of
//! the four §4.7 steps:
//!
//! * **Block 81 (this module).** Run `Ã(z)` over `sd(k)` to obtain the
//!   LPC residual `d(k) = sd(k) − Σ ã_i · sd(k-i)`.
//! * **Block 82.** Lowpass + decimate the residual, correlate over
//!   lags 5..35, refine over `4τ−3..4τ+3`, fundamental-vs-multiple
//!   resolution against the previous frame's pitch.
//! * **Block 83.** Single-tap pitch predictor weight `β` over the
//!   decoded-speech buffer.
//! * **Block 84.** Map `β` → `(g_l, b)` per eq. 4-13 / 4-14; clamp `p`
//!   into `[KPMIN, KPMAX]`.
//!
//! The order-10 by-product `ã_1..ã_10` is supplied by the synthesis-
//! filter adapter (see
//! [`crate::synthesis_adapter::SynthesisAdapter::order10_predictor`]),
//! exactly the same source that already feeds the short-term postfilter
//! (block 72). Per §7.2 it refreshes at the first vector of each
//! adaptation cycle; the pitch extractor, per §7.1, uses the result for
//! the rest of that cycle's residual.
//!
//! ## Sign convention
//!
//! Our [`crate::levinson`] solver returns the polynomial in the form
//! `A(z) = 1 + a_1·z^{-1} + … + a_M·z^{-M}`, so the spec's `ã_i` (the
//! predictor coefficients of the underlying `s(n) ≈ Σ ã_i · s(n-i)`
//! model) are `ã_i = −a_i`. This module performs the sign flip
//! internally so callers may pass the raw `order10_predictor()` slice
//! straight through.
//!
//! ## Cold-start
//!
//! Before the first adaptation cycle commits a non-trivial predictor,
//! the order-10 by-product is the all-pass polynomial `(1, 0, …, 0)`
//! ([`crate::synthesis_adapter::SynthesisAdapter::new`] initial state).
//! Under that input the inverse filter reduces to `d(k) = sd(k)` —
//! the residual passes the decoded speech through unchanged, matching
//! the spec's "no predictor information available yet" cold start.
//!
//! ## State carry-over
//!
//! `Ã(z)` is a 10-tap FIR filter; its per-sample memory is the last
//! 10 input samples `sd(k-1)..sd(k-10)`, initialised to zero per the
//! Table 2/G.728 cold-start convention that applies uniformly to every
//! internal filter in the codec. The memory carries across vector
//! boundaries.

use crate::consts::IDIM;
use crate::synthesis_adapter::PF_LPC_ORDER;

/// Pitch-extractor LPC inverse filter (block 81) state.
///
/// One instance per decoder. Carries:
///
/// * The bandwidth-flipped 10th-order predictor `atilde_i = -a_i` for
///   `1 ≤ i ≤ 10` (1-based; index `[0]` is unused — the implicit
///   leading `1` of `Ã(z)` is realised in the difference equation).
/// * A 10-sample memory `mem` of the last decoded-speech inputs
///   `sd(k-1)..sd(k-10)`. `mem[0]` is the most recent past input,
///   `mem[9]` the oldest.
#[derive(Debug, Clone)]
pub struct PitchInverseFilter {
    /// Spec-convention predictor coefficients `ã_i = -a_i` for
    /// `1 ≤ i ≤ 10`. `[0]` is unused (implicit leading 1 in `Ã(z)`).
    atilde: [f64; PF_LPC_ORDER + 1],
    /// Per-sample input memory. `mem[i]` is `sd(k - (i+1))`.
    mem: [f64; PF_LPC_ORDER],
}

impl Default for PitchInverseFilter {
    fn default() -> Self {
        Self::new()
    }
}

impl PitchInverseFilter {
    /// Construct a fresh inverse filter on the all-pass `Ã(z) = 1`
    /// initial state — every coefficient zero, the input memory zero.
    /// Under this state the residual equals the input: `d(k) = sd(k)`.
    pub fn new() -> Self {
        Self {
            atilde: [0.0; PF_LPC_ORDER + 1],
            mem: [0.0; PF_LPC_ORDER],
        }
    }

    /// Refresh the inverse-filter coefficients from the order-10 by-
    /// product Levinson output of the synthesis-filter adapter.
    ///
    /// * `a10` — `[1.0, a_1, …, a_10]` in our Levinson `A` layout
    ///   (see [`crate::levinson`]). The sign flip `ã_i = −a_i` is
    ///   applied internally per the spec's predictor convention.
    ///
    /// Per §7.1 the pitch extractor should call this once per
    /// adaptation cycle (the same first-vector boundary at which the
    /// short-term postfilter, block 72, refreshes its coefficients).
    pub fn set_from_synthesis_byproduct(&mut self, a10: &[f64; PF_LPC_ORDER + 1]) {
        // Spec convention: ã_i = -a_i. `a10[0]` is the leading `1.0`
        // of our `A(z)` and is unused here (the leading `1` of `Ã(z)`
        // is realised by the `d(k) = sd(k) − Σ ã_i · sd(k-i)` form of
        // the difference equation).
        for i in 1..=PF_LPC_ORDER {
            self.atilde[i] = -a10[i];
        }
        self.atilde[0] = 0.0;
    }

    /// Read-only access to the spec-convention predictor coefficients
    /// `ã_i` for tests and audit. The useful entries are `[1..=10]`;
    /// `[0]` is always zero (placeholder for the implicit leading 1).
    pub fn coefficients(&self) -> &[f64; PF_LPC_ORDER + 1] {
        &self.atilde
    }

    /// Read-only access to the input memory for tests and audit.
    /// `mem[i]` is the input sample `i+1` samples in the past.
    pub fn memory(&self) -> &[f64; PF_LPC_ORDER] {
        &self.mem
    }

    /// Reset coefficients and memory to the cold-start defaults.
    /// Useful for tests.
    pub fn reset(&mut self) {
        self.atilde = [0.0; PF_LPC_ORDER + 1];
        self.mem = [0.0; PF_LPC_ORDER];
    }

    /// Filter one IDIM-sample decoded-speech vector `sd` through the
    /// 10th-order LPC inverse filter and return the residual vector
    /// `d(k)`.
    ///
    /// Per-sample dataflow (eq. 4-6):
    ///
    /// 1. `d(k) = sd(k) − Σ_{i=1..10} ã_i · mem[i-1]`, where `mem[i-1]`
    ///    is `sd(k-i)`.
    /// 2. Shift the memory one position deeper: `mem[i] ← mem[i-1]`
    ///    for `i = 9..1`, then `mem[0] ← sd(k)`.
    ///
    /// The memory advances one sample per output sample, preserving
    /// the 10-tap FIR state across vector boundaries.
    pub fn filter_vector(&mut self, sd: &[f64; IDIM]) -> [f64; IDIM] {
        let mut out = [0.0f64; IDIM];
        for k in 0..IDIM {
            // ----- d(k) = sd(k) − Σ ã_i · sd(k-i) ---------------------
            let mut d = sd[k];
            for i in 1..=PF_LPC_ORDER {
                d -= self.atilde[i] * self.mem[i - 1];
            }
            out[k] = d;

            // ----- Shift memory and push sd(k) at the front ----------
            // mem[i] = sd(k-(i+1)), so after the shift mem[0] holds
            // sd(k) (the new "k-1" for the next call) and mem[9] is
            // discarded.
            for i in (1..PF_LPC_ORDER).rev() {
                self.mem[i] = self.mem[i - 1];
            }
            self.mem[0] = sd[k];
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---------- Cold-start invariants ---------------------------------

    #[test]
    fn fresh_filter_passes_input_through_unchanged() {
        // All coefficients zero → `Ã(z) = 1` → d(k) = sd(k). Confirm
        // the first IDIM-sample input emerges unchanged.
        let mut pf = PitchInverseFilter::new();
        let sd = [100.0, -50.0, 0.0, 25.0, 80.0];
        let d = pf.filter_vector(&sd);
        for k in 0..IDIM {
            assert_eq!(d[k], sd[k], "k={k}");
        }
    }

    #[test]
    fn fresh_coefficients_are_zero() {
        let pf = PitchInverseFilter::new();
        assert!(pf.coefficients().iter().all(|&v| v == 0.0));
        assert!(pf.memory().iter().all(|&v| v == 0.0));
    }

    #[test]
    fn allpass_a10_drives_identity_output_across_many_vectors() {
        // Explicitly load the all-pass predictor `(1, 0, …, 0)` — the
        // same shape the synthesis adapter starts at — and confirm the
        // inverse filter is the identity across many vectors.
        let mut pf = PitchInverseFilter::new();
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        pf.set_from_synthesis_byproduct(&a10);
        for v in 0..32 {
            let sd = [
                v as f64,
                (v * 2 - 1) as f64,
                v as f64 + 7.0,
                -(v as f64 * 0.5),
                v as f64 * 3.0,
            ];
            let d = pf.filter_vector(&sd);
            for k in 0..IDIM {
                assert_eq!(d[k], sd[k], "vector {v} sample {k}");
            }
        }
    }

    // ---------- Sign-convention math ----------------------------------

    #[test]
    fn set_from_byproduct_applies_sign_flip() {
        // Our levinson returns A(z); the spec works with `ã_i = -a_i`.
        // Confirm the setter performs the flip.
        let mut pf = PitchInverseFilter::new();
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        for i in 1..=PF_LPC_ORDER {
            a10[i] = 0.1 * i as f64;
        }
        pf.set_from_synthesis_byproduct(&a10);
        for i in 1..=PF_LPC_ORDER {
            assert!(
                (pf.coefficients()[i] + a10[i]).abs() < 1e-15,
                "ã_{i} should be −a_{i}: got {} vs −{}",
                pf.coefficients()[i],
                a10[i]
            );
        }
        // The leading slot stays zeroed (placeholder; the implicit
        // leading 1 of Ã(z) is realised in the difference equation).
        assert_eq!(pf.coefficients()[0], 0.0);
    }

    // ---------- Eq. 4-6 against a controlled predictor ----------------

    #[test]
    fn impulse_response_matches_atilde_sign_flipped() {
        // With Ã(z) = 1 − Σ ã_i z^{-i} loaded with a known pattern,
        // feed a single impulse at sample 0 and read the impulse
        // response. The output should be:
        //   d(0) = 1
        //   d(i) = -ã_i for 1 ≤ i ≤ 10
        //   d(i) = 0 for i > 10
        let mut pf = PitchInverseFilter::new();
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        // Choose a non-trivial pattern; the loader flips signs so
        // the spec-convention atilde will be [0, -0.1, -0.2, …].
        for i in 1..=PF_LPC_ORDER {
            a10[i] = 0.1 * i as f64;
        }
        pf.set_from_synthesis_byproduct(&a10);

        // 15 samples total (3 vectors of IDIM = 5). The first sample
        // is a unit impulse, the rest are zero.
        let mut outs = [0.0f64; 15];
        for v in 0..3 {
            let mut sd = [0.0f64; IDIM];
            if v == 0 {
                sd[0] = 1.0;
            }
            let d = pf.filter_vector(&sd);
            for k in 0..IDIM {
                outs[v * IDIM + k] = d[k];
            }
        }
        // d(0) = 1 (impulse passes through Ã's leading 1).
        assert!((outs[0] - 1.0).abs() < 1e-15);
        // d(i) for 1 ≤ i ≤ 10: d(i) = −ã_i · sd(0) = +a_i (since the
        // setter stored ã_i = −a_i, the per-sample sum has sign:
        //   d(i) = 0 − ã_i · sd(0) = −ã_i = a_i = 0.1·i.
        for i in 1..=PF_LPC_ORDER {
            let expected = a10[i];
            assert!(
                (outs[i] - expected).abs() < 1e-14,
                "d({i}) expected {expected}, got {}",
                outs[i]
            );
        }
        // d(i) for i > 10: the impulse has fallen off the 10-tap
        // memory, so the residual is zero.
        for i in (PF_LPC_ORDER + 1)..15 {
            assert!(outs[i].abs() < 1e-15, "d({i}) expected 0, got {}", outs[i]);
        }
    }

    // ---------- Roundtrip vs synthesis ---------------------------------

    #[test]
    fn inverse_undoes_a_synthesis_step_on_a_single_tap_model() {
        // For a single-tap AR(1) synthesis model `s(n) = α · s(n-1) +
        // e(n)`, the matching inverse filter is `Ã(z) = 1 − α·z^{-1}`,
        // i.e. ã_1 = α and ã_i = 0 for i > 1. Equivalently, our `A(z)`
        // layout has a_1 = -α, a_2..10 = 0.
        //
        // Feed a known excitation `e(n) = 1, 0, 0, …, 0`, synthesize
        // `s(n) = α^n`, and confirm the inverse filter exactly
        // recovers `e(n)`.
        let alpha = 0.5_f64;
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        a10[1] = -alpha;

        let mut pf = PitchInverseFilter::new();
        pf.set_from_synthesis_byproduct(&a10);

        // Synthesize s(n) = α^n for n = 0..14 (3 vectors of 5).
        let mut s = [0.0f64; 15];
        s[0] = 1.0;
        for n in 1..15 {
            s[n] = alpha * s[n - 1];
        }

        // Push s through the inverse filter in IDIM-sized vectors.
        let mut recovered = [0.0f64; 15];
        for v in 0..3 {
            let mut sd = [0.0f64; IDIM];
            for k in 0..IDIM {
                sd[k] = s[v * IDIM + k];
            }
            let d = pf.filter_vector(&sd);
            for k in 0..IDIM {
                recovered[v * IDIM + k] = d[k];
            }
        }

        // First sample: e(0) = 1 (the impulse).
        assert!((recovered[0] - 1.0).abs() < 1e-12);
        // All subsequent samples should be zero (the inverse filter
        // perfectly cancels the AR(1) recursion).
        for n in 1..15 {
            assert!(
                recovered[n].abs() < 1e-12,
                "n={n}: residual should be 0, got {}",
                recovered[n]
            );
        }
    }

    // ---------- Memory carry-over across vectors ----------------------

    #[test]
    fn memory_carries_across_vector_boundaries() {
        // Set ã_1 = 1 (one-sample-back tap). With this, d(k) = sd(k) -
        // sd(k-1). Confirm the tap reaches BACK into the previous
        // vector — i.e. the residual at sample 0 of vector 1 depends
        // on sample IDIM-1 of vector 0.
        let mut pf = PitchInverseFilter::new();
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        a10[1] = -1.0; // ã_1 = +1 after sign flip.

        pf.set_from_synthesis_byproduct(&a10);
        let v0 = [10.0, 20.0, 30.0, 40.0, 50.0];
        let _ = pf.filter_vector(&v0);
        let v1 = [60.0, 70.0, 80.0, 90.0, 100.0];
        let d1 = pf.filter_vector(&v1);
        // d1[0] = sd[5] - sd[4] = 60 - 50 = 10
        assert!((d1[0] - 10.0).abs() < 1e-12);
        // d1[1] = sd[6] - sd[5] = 70 - 60 = 10
        assert!((d1[1] - 10.0).abs() < 1e-12);
        // d1[2..4] = 10 each (constant difference of 10).
        for k in 2..IDIM {
            assert!((d1[k] - 10.0).abs() < 1e-12, "k={k}: got {}", d1[k]);
        }
    }

    // ---------- Reset --------------------------------------------------

    #[test]
    fn reset_returns_to_identity() {
        let mut pf = PitchInverseFilter::new();
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        a10[1] = -0.7;
        pf.set_from_synthesis_byproduct(&a10);
        let _ = pf.filter_vector(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        pf.reset();
        assert!(pf.coefficients().iter().all(|&v| v == 0.0));
        assert!(pf.memory().iter().all(|&v| v == 0.0));
        // After reset the filter is the identity again.
        let sd = [42.0, -17.0, 99.0, 0.0, 1.0];
        let d = pf.filter_vector(&sd);
        for k in 0..IDIM {
            assert_eq!(d[k], sd[k]);
        }
    }

    // ---------- Finiteness smoke -------------------------------------

    #[test]
    fn output_finite_under_sinusoidal_drive() {
        // Drive a non-trivial predictor with bounded sinusoidal input
        // and confirm the residual stays finite.
        let mut pf = PitchInverseFilter::new();
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        for i in 1..=PF_LPC_ORDER {
            a10[i] = 0.1 / i as f64;
        }
        pf.set_from_synthesis_byproduct(&a10);
        for v in 0..256 {
            let mut sd = [0.0f64; IDIM];
            for k in 0..IDIM {
                let t = (v * IDIM + k) as f64;
                sd[k] = 2000.0 * (0.03 * t).sin();
            }
            let d = pf.filter_vector(&sd);
            for &x in &d {
                assert!(x.is_finite());
            }
        }
    }
}
