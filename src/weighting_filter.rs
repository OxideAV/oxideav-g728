//! Perceptual weighting filter (block 4 of Figure 2/G.728), §3.4.
//!
//! ## Scope
//!
//! Round 248 transcribed the **coefficient calculator** (block 38,
//! §3.3) — given the order-10 LPC predictor `q_i` from block 37, it
//! produces the two bandwidth-broadened sequences `qγ₁_i = q_i · γ₁ⁱ`
//! and `qγ₂_i = q_i · γ₂ⁱ` (eqs. 3-4b / 3-4c) defining the filter's
//! transfer function
//!
//! ```text
//!         Q̃(z/γ₁)     1 + Σᵢ qγ₁_i · z⁻ⁱ
//! W(z) = ────────── = ─────────────────────                       (3-4a)
//!         Q̃(z/γ₂)     1 + Σᵢ qγ₂_i · z⁻ⁱ
//! ```
//!
//! Round 249 lands the **application** of the filter (block 4):
//! given the current input speech vector `s(n)`, push every sample
//! through `W(z)` and emit the corresponding weighted speech vector
//! `v(n)`. §3.4 of the Recommendation is one paragraph:
//!
//! > In Figure 2/G.728, the current input speech vector `s(n)` is
//! > passed through the perceptual weighting filter (block 4),
//! > resulting in the weighted speech vector `v(n)`. Note that
//! > except during initialization, the filter memory (i.e. internal
//! > state variables, or the values held in the delay units of the
//! > filter) should not be reset to zero at any time.
//!
//! Block 10 (a second instance of `W(z)` used to compute the
//! zero-input response in cascade with the synthesis filter, §3.5)
//! has a different state-handling rule and is deliberately not part of
//! this round — it requires the cascade `F(z)·W(z)` with the synthesis
//! filter's pre-/post-save dance described in §3.10.
//!
//! ## Time-domain realisation
//!
//! Eq. 3-4a in difference form is the standard direct-form I pole-zero
//! filter:
//!
//! ```text
//! v(n) = s(n) + Σ_{i=1..10} qγ₁_i · s(n−i)
//!             − Σ_{i=1..10} qγ₂_i · v(n−i)
//! ```
//!
//! derived by cross-multiplying `V(z) · (1 + Σ qγ₂_i z⁻ⁱ) =
//! S(z) · (1 + Σ qγ₁_i z⁻ⁱ)` and solving for `v(n)`. The "leading 1"
//! taps in both numerator and denominator come from the implicit
//! `q_0 · γ_k^0 = 1` of `Q̃(z)` (eq. 3-3a) — they are not stored as
//! coefficients but appear explicitly as the `+ s(n)` and the standalone
//! `v(n)` on the two sides of the difference equation.
//!
//! ## State handling (§3.4)
//!
//! The spec is explicit: the per-sample memory of input samples
//! `s(n-1..n-10)` and output samples `v(n-1..n-10)` is **not** reset
//! between adaptation cycles, between calls, or at any cycle boundary
//! — only at construction.
//!
//! The §3.4.1 non-speech "disable" path sets `γ₁ = γ₂ = 0`, which the
//! coefficient calculator already realises by zeroing every non-unity
//! tap of both `q_gamma1` and `q_gamma2`. With every coefficient zero
//! the difference equation collapses to `v(n) = s(n)` and the filter
//! is a sample-rate pass-through — but the spec's "do not reset
//! memory" rule still applies. This module therefore never zeroes
//! `s_mem` / `v_mem` outside [`PerceptualWeightingFilter::new`].
//!
//! ## Coefficient update timing (§3.3 / §3.4)
//!
//! Per §3.3 the coefficient calculator (block 38) updates the
//! `W(z)` taps at the third vector of every 4-vector adaptation
//! cycle. While the new coefficients are loaded the filter memory
//! carries continuously — this is the standard "freeze + swap"
//! convention of every backward-adaptive predictor in the codec and is
//! what makes the filter's perceptual response evolve smoothly across
//! frames. The caller (the encoder) drives the coefficient swap; this
//! module just exposes [`set_coefficients`] and leaves
//! [`filter_vector`] to consume the current state.

use crate::consts::{IDIM, LPCW};
use crate::weighting_filter_coeff::WeightingFilterCoeff;

/// Per-sample perceptual weighting filter (block 4 of Figure 2/G.728).
///
/// Carries the current `W(z) = Q̃(z/γ₁) / Q̃(z/γ₂)` coefficient set
/// (output of block 38) plus the two 10-tap delay lines for past input
/// `s(n-1..n-10)` and past output `v(n-1..n-10)`. Per §3.4 the delay
/// lines are initialised to zero at construction and never reset
/// thereafter.
#[derive(Debug, Clone)]
pub struct PerceptualWeightingFilter {
    /// Current `W(z)` coefficients — `q_gamma1` is the numerator's
    /// `(qγ₁_i)` sequence with the implicit unity tap at index 0,
    /// `q_gamma2` is the denominator's `(qγ₂_i)` sequence with the same
    /// layout. Indices `[1..=LPCW]` are the broadened taps used by
    /// [`filter_vector`]; `[0]` is the implicit `q_0·γ_k^0 = 1` of
    /// `Q̃(z)` and is consumed implicitly by the difference equation,
    /// not by the inner sums.
    coeff: WeightingFilterCoeff,
    /// Memory of past input samples `s(n-1..n-10)`. Index 0 holds the
    /// most-recently-pushed sample (`s(n-1)`), index 9 holds the
    /// oldest (`s(n-10)`). Matches the layout used by
    /// [`crate::ShortTermPostfilter`] so the per-sample shift idiom is
    /// the same across the crate.
    s_mem: [f64; LPCW],
    /// Memory of past output samples `v(n-1..n-10)`. Same layout as
    /// `s_mem`.
    v_mem: [f64; LPCW],
}

impl Default for PerceptualWeightingFilter {
    /// Fresh filter — coefficients at the §3.4 / §3.4.1 all-pass
    /// `W(z) = 1` state (`q_gamma1 = q_gamma2 = [1, 0, ..., 0]`) and
    /// both delay lines at zero. This is the spec's "initialization"
    /// state per §3.4; the encoder constructs the filter once and
    /// thereafter only loads new coefficients via [`set_coefficients`].
    fn default() -> Self {
        Self::new()
    }
}

impl PerceptualWeightingFilter {
    /// Construct a fresh filter at the §3.4 initialisation state:
    /// `W(z) = 1` coefficients, both delay lines zero.
    ///
    /// The all-pass coefficient set is identical to
    /// [`WeightingFilterCoeff::disabled`] — the §3.4.1 "non-speech
    /// mode" value — so a freshly-constructed filter passes input
    /// straight through until the first call to [`set_coefficients`].
    pub fn new() -> Self {
        Self {
            coeff: WeightingFilterCoeff::disabled(),
            s_mem: [0.0; LPCW],
            v_mem: [0.0; LPCW],
        }
    }

    /// Read-only view of the current coefficient set.
    pub fn coefficients(&self) -> &WeightingFilterCoeff {
        &self.coeff
    }

    /// Read-only view of the input-side delay line. Index 0 is the
    /// most-recent past input `s(n-1)`. Useful for tests and audit.
    pub fn input_memory(&self) -> &[f64; LPCW] {
        &self.s_mem
    }

    /// Read-only view of the output-side delay line. Index 0 is the
    /// most-recent past output `v(n-1)`. Useful for tests and audit.
    pub fn output_memory(&self) -> &[f64; LPCW] {
        &self.v_mem
    }

    /// Load a new coefficient set produced by the block-38 calculator
    /// (output of [`WeightingFilterCoeff::from_lpc`]). Per §3.3 / §3.4
    /// the encoder calls this once per 4-vector adaptation cycle (at
    /// the third vector); the per-sample memory of the filter is left
    /// untouched, which is the standard "freeze + swap" convention.
    pub fn set_coefficients(&mut self, coeff: WeightingFilterCoeff) {
        self.coeff = coeff;
    }

    /// Filter one `IDIM`-sample input speech vector `s(n)` through
    /// `W(z)` and return the corresponding weighted vector `v(n)`.
    ///
    /// Per-sample dataflow (direct-form I pole-zero):
    ///
    /// 1. Compute the all-zero numerator branch
    ///    `u(n) = s(n) + Σ qγ₁_i · s(n-i)` using `s_mem[0..LPCW]` for
    ///    the past inputs.
    /// 2. Compute the all-pole denominator branch
    ///    `v(n) = u(n) − Σ qγ₂_i · v(n-i)` using `v_mem[0..LPCW]` for
    ///    the past outputs.
    /// 3. Shift both delay lines: insert `s(n)` at index 0 of `s_mem`
    ///    and `v(n)` at index 0 of `v_mem`, pushing older samples
    ///    one slot deeper.
    ///
    /// Both delay lines advance one slot per emitted output sample;
    /// the spec's "do not reset memory" rule of §3.4 is honoured by
    /// the fact that this routine never touches the arrays outside
    /// the per-sample shift.
    pub fn filter_vector(&mut self, s: &[f64; IDIM]) -> [f64; IDIM] {
        let mut out = [0.0f64; IDIM];
        let qg1 = &self.coeff.q_gamma1;
        let qg2 = &self.coeff.q_gamma2;
        for k in 0..IDIM {
            // All-zero numerator: u(n) = s(n) + Σ qγ₁_i · s(n-i).
            // The "+ s(n)" term is the implicit `q_0 · γ₁^0 = 1` of
            // `Q̃(z)` (eq. 3-3a). qg1[0] is also 1.0 by construction
            // and is *not* re-applied inside the loop — only the
            // broadened taps qg1[1..=LPCW] interact with past inputs.
            let mut u = s[k];
            for i in 1..=LPCW {
                u += qg1[i] * self.s_mem[i - 1];
            }
            // All-pole denominator: v(n) = u(n) − Σ qγ₂_i · v(n-i).
            // The "− Σ" sign comes from solving the cross-product
            // identity `V(z)·Q̃(z/γ₂) = S(z)·Q̃(z/γ₁)` for `v(n)` —
            // the denominator's past-sample contributions move to the
            // right-hand side with a sign flip. qg2[0] = 1 is again
            // the implicit leading tap and is not in the sum.
            let mut v = u;
            for i in 1..=LPCW {
                v -= qg2[i] * self.v_mem[i - 1];
            }
            // Advance the per-sample delay lines: push s(n) / v(n) at
            // index 0, shift older samples one slot deeper. The
            // spec's §3.4 "do not reset memory" rule means this is the
            // only path that ever modifies these arrays outside of
            // construction; the shift idiom matches the one in
            // crate::short_term_postfilter::ShortTermPostfilter so the
            // per-sample dataflow looks identical.
            for i in (1..LPCW).rev() {
                self.s_mem[i] = self.s_mem[i - 1];
                self.v_mem[i] = self.v_mem[i - 1];
            }
            self.s_mem[0] = s[k];
            self.v_mem[0] = v;
            out[k] = v;
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::{WPCF, WZCF};

    // ---------- §3.4 / §3.4.1 cold-start invariants -------------------

    #[test]
    fn fresh_filter_is_all_pass() {
        // §3.4 initialisation + §3.4.1 disabled-mode state: with the
        // disabled-mode coefficient set, every non-unity tap is zero
        // and the difference equation collapses to v(n) = s(n).
        let mut wf = PerceptualWeightingFilter::new();
        let s = [10.0, -20.0, 0.0, 30.0, -5.0];
        let v = wf.filter_vector(&s);
        for k in 0..IDIM {
            assert_eq!(v[k], s[k], "k={k}");
        }
    }

    #[test]
    fn fresh_memory_is_zero() {
        // §3.4 spec note: "the filter memory ... should not be reset
        // to zero at any time except during initialization." This is
        // exactly that initialization state.
        let wf = PerceptualWeightingFilter::new();
        for i in 0..LPCW {
            assert_eq!(wf.input_memory()[i], 0.0);
            assert_eq!(wf.output_memory()[i], 0.0);
        }
    }

    #[test]
    fn default_matches_new() {
        // Both construction paths produce the §3.4 initialisation
        // state — same coefficient set, same zeroed memory.
        let a = PerceptualWeightingFilter::default();
        let b = PerceptualWeightingFilter::new();
        assert_eq!(a.input_memory(), b.input_memory());
        assert_eq!(a.output_memory(), b.output_memory());
        assert_eq!(a.coefficients(), b.coefficients());
    }

    // ---------- Memory carry-over across vector boundaries ------------

    #[test]
    fn memory_advances_one_sample_per_filtered_output() {
        // Smoke: after one filter_vector call the input delay line
        // should carry the last IDIM input samples in reverse-time
        // order (index 0 = most recent).
        let mut wf = PerceptualWeightingFilter::new();
        let s = [1.0, 2.0, 3.0, 4.0, 5.0];
        let _ = wf.filter_vector(&s);
        // After processing s[0..5], the delay line should be:
        //   s_mem[0] = s[4] = 5
        //   s_mem[1] = s[3] = 4
        //   s_mem[2] = s[2] = 3
        //   s_mem[3] = s[1] = 2
        //   s_mem[4] = s[0] = 1
        //   s_mem[5..10] = 0 (never written this call)
        assert_eq!(wf.input_memory()[0], 5.0);
        assert_eq!(wf.input_memory()[1], 4.0);
        assert_eq!(wf.input_memory()[2], 3.0);
        assert_eq!(wf.input_memory()[3], 2.0);
        assert_eq!(wf.input_memory()[4], 1.0);
        for i in 5..LPCW {
            assert_eq!(wf.input_memory()[i], 0.0);
        }
    }

    #[test]
    fn memory_persists_across_vector_calls() {
        // §3.4 invariant: per-sample memory is not reset between
        // calls. After two IDIM-sample calls the input delay line
        // should hold the two newest samples of the *second* call at
        // indices 0..2 and the two oldest of the second call /
        // newest of the first at indices 3..5.
        let mut wf = PerceptualWeightingFilter::new();
        let s1 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let s2 = [6.0, 7.0, 8.0, 9.0, 10.0];
        let _ = wf.filter_vector(&s1);
        let _ = wf.filter_vector(&s2);
        // After both: most-recent at index 0 is s2[4] = 10, then
        // s2[3] = 9, s2[2] = 8, s2[1] = 7, s2[0] = 6, then the second
        // call's "past" which is the first call's: s1[4] = 5, s1[3] = 4,
        // s1[2] = 3, s1[1] = 2, s1[0] = 1.
        let expected = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
        for i in 0..LPCW {
            assert_eq!(wf.input_memory()[i], expected[i], "i={i}");
        }
    }

    // ---------- Coefficient swap path ---------------------------------

    #[test]
    fn set_coefficients_replaces_taps_without_touching_memory() {
        // §3.3 / §3.4: at the third vector of every adaptation cycle
        // the encoder reloads the W(z) coefficients via block 38, but
        // the per-sample memory must continue across the swap.
        let mut wf = PerceptualWeightingFilter::new();
        let s = [11.0, 22.0, 33.0, 44.0, 55.0];
        let _ = wf.filter_vector(&s);
        let mem_before = *wf.input_memory();

        // Swap in a non-trivial coefficient set.
        let mut q = [0.0f64; LPCW + 1];
        q[0] = 1.0;
        q[1] = -0.5;
        q[2] = 0.25;
        let new_coeff = WeightingFilterCoeff::from_lpc(&q);
        wf.set_coefficients(new_coeff);

        // Memory must be untouched by the swap.
        assert_eq!(*wf.input_memory(), mem_before);
        assert_eq!(*wf.coefficients(), new_coeff);
    }

    // ---------- Difference-equation correctness ------------------------

    #[test]
    fn first_sample_after_init_equals_input_sample_for_any_coefficients() {
        // At cold start every past s(n-i) and v(n-i) is zero, so the
        // very first emitted sample after a coefficient swap is
        // exactly s[0] regardless of the W(z) shape: every term in
        // both sums multiplies a zero memory entry. This pins the
        // difference equation's "+ s(n)" / standalone-v(n) sides.
        let mut q = [0.0f64; LPCW + 1];
        q[0] = 1.0;
        q[1] = -0.7;
        q[2] = 0.42;
        q[3] = -0.15;
        let coeff = WeightingFilterCoeff::from_lpc(&q);

        let mut wf = PerceptualWeightingFilter::new();
        wf.set_coefficients(coeff);
        let s = [123.456, 0.0, 0.0, 0.0, 0.0];
        let v = wf.filter_vector(&s);
        assert!(
            (v[0] - s[0]).abs() < 1e-12,
            "v[0] = {} vs s[0] = {}",
            v[0],
            s[0]
        );
    }

    #[test]
    fn impulse_response_matches_difference_equation_termwise() {
        // Drive a δ(n) input through a known W(z) and step through
        // the first several output samples by hand to pin the sign
        // convention of both sums.
        //
        // q[0] = 1, q[1] = a1, rest = 0  →
        //   qγ₁_1 = a1·γ₁,  qγ₂_1 = a1·γ₂,  higher taps = 0.
        // δ-input means s[0] = 1, s[k>0] = 0.
        //
        // Cold-start memory:
        //   v(0) = s[0] = 1
        //   s_mem after step 0: [1, 0, ..., 0]
        //   v_mem after step 0: [1, 0, ..., 0]
        //
        //   v(1) = 0 + qγ₁_1·1 − qγ₂_1·1 = a1·(γ₁ − γ₂)
        //   s_mem: [0, 1, 0, ...], v_mem: [a1·(γ₁−γ₂), 1, 0, ...]
        //
        //   v(2) = 0 + qγ₁_1·0 + qγ₁_2·1 − qγ₂_1·v(1) − qγ₂_2·1
        //        = 0 − qγ₂_1·a1·(γ₁−γ₂)    (since qγ₁_2 = qγ₂_2 = 0)
        //        = −a1·γ₂ · a1·(γ₁ − γ₂)
        //        = −a1² · γ₂ · (γ₁ − γ₂)
        let a1 = 0.5_f64;
        let mut q = [0.0f64; LPCW + 1];
        q[0] = 1.0;
        q[1] = a1;
        let coeff = WeightingFilterCoeff::from_lpc(&q);
        let mut wf = PerceptualWeightingFilter::new();
        wf.set_coefficients(coeff);

        let mut delta = [0.0_f64; IDIM];
        delta[0] = 1.0;
        let v = wf.filter_vector(&delta);

        let expected_v0 = 1.0;
        let expected_v1 = a1 * (WZCF - WPCF);
        let expected_v2 = -a1 * a1 * WPCF * (WZCF - WPCF);
        let expected_v3 = -a1 * WPCF * expected_v2;
        let expected_v4 = -a1 * WPCF * expected_v3;

        assert!((v[0] - expected_v0).abs() < 1e-12, "v[0] = {}", v[0]);
        assert!((v[1] - expected_v1).abs() < 1e-12, "v[1] = {}", v[1]);
        assert!((v[2] - expected_v2).abs() < 1e-12, "v[2] = {}", v[2]);
        assert!((v[3] - expected_v3).abs() < 1e-12, "v[3] = {}", v[3]);
        assert!((v[4] - expected_v4).abs() < 1e-12, "v[4] = {}", v[4]);
    }

    // ---------- Special-case: γ₁ = γ₂ degenerates to identity --------

    #[test]
    fn equal_gammas_collapse_filter_to_identity() {
        // When γ₁ = γ₂ the numerator and denominator polynomials are
        // identical and the transfer function reduces to W(z) = 1,
        // i.e. the filter is a pass-through *regardless of q_i*. This
        // is the same algebraic identity that backs the §3.4.1
        // γ₁ = γ₂ = 0 disable shortcut, but here we test it with
        // non-zero coefficients to exercise the full per-sample path
        // (every loop iteration runs at least once on non-zero terms).
        let q = [
            1.0_f64,
            0.5,
            -0.25,
            0.125,
            -0.0625,
            0.03125,
            -0.015_625,
            0.007_812_5,
            -0.003_906_25,
            0.001_953_125,
            -0.000_976_562_5,
        ];
        let coeff = WeightingFilterCoeff::from_lpc_with_gammas(&q, 0.8, 0.8);
        let mut wf = PerceptualWeightingFilter::new();
        wf.set_coefficients(coeff);

        // Drive several vectors of an arbitrary signal through and
        // verify the output exactly equals the input each time. The
        // per-sample feedback path means even rounding error would
        // accumulate noticeably across 64 vectors.
        for v_idx in 0..64 {
            let mut s = [0.0f64; IDIM];
            for k in 0..IDIM {
                s[k] = ((v_idx * IDIM + k) as f64 * 0.37).sin() * 500.0;
            }
            let out = wf.filter_vector(&s);
            for k in 0..IDIM {
                assert!(
                    (out[k] - s[k]).abs() < 1e-9,
                    "vec {v_idx} sample {k}: out={} s={}",
                    out[k],
                    s[k]
                );
            }
        }
    }

    // ---------- Stability smoke ---------------------------------------

    #[test]
    fn finite_output_under_long_sinusoidal_drive() {
        // §3.3 / §3.4 invariants γ₁ = 0.9 / γ₂ = 0.6 < 1 keep the
        // denominator's poles strictly inside the unit circle for any
        // q_i with |q_i| bounded — the standard bandwidth-expansion
        // robustness argument. Driving 1024 vectors of a synthetic
        // sinusoid through must produce a finite output (no Inf, no
        // NaN) for a generic non-trivial coefficient set.
        let mut q = [0.0f64; LPCW + 1];
        q[0] = 1.0;
        q[1] = -0.3;
        q[2] = 0.1;
        q[3] = -0.05;
        let coeff = WeightingFilterCoeff::from_lpc(&q);
        let mut wf = PerceptualWeightingFilter::new();
        wf.set_coefficients(coeff);
        for v_idx in 0..1024 {
            let mut s = [0.0f64; IDIM];
            for k in 0..IDIM {
                let t = (v_idx * IDIM + k) as f64;
                s[k] = 1500.0 * (0.31 * t).sin();
            }
            let out = wf.filter_vector(&s);
            for &x in &out {
                assert!(x.is_finite(), "non-finite output at vec {v_idx}");
            }
        }
    }
}
