//! Short-term (spectral) postfilter — block 72 of Figure 7/G.728
//! (postfilter sub-stage).
//!
//! Implements the §4.6 short-term postfilter transfer function
//! (equations 4-2 .. 4-5):
//!
//! ```text
//!                  1 − Σ b̄_i · z^{-i}
//!     H_s(z)  =  ─────────────────────── · (1 + µ · z^{-1})
//!                  1 − Σ ā_i · z^{-i}
//! ```
//!
//! with
//!
//! * `b̄_i = ã_i · SPFZCF^i`  (zero-control, `SPFZCF = 0.65`, eq. 4-3),
//! * `ā_i = ã_i · SPFPCF^i`  (pole-control, `SPFPCF = 0.75`, eq. 4-4),
//! * `µ   = TILTF · k1`      (tilt-compensation, `TILTF = 0.15`,
//!   eq. 4-5),
//!
//! where `ã_1..ã_10` are the order-10 LPC by-products of the synthesis-
//! filter Levinson recursion and `k1` is its first reflection
//! coefficient — both surfaced by
//! [`crate::synthesis_adapter::SynthesisAdapter::order10_predictor`]
//! and [`crate::synthesis_adapter::SynthesisAdapter::k1`].
//!
//! ## Coefficient sign convention
//!
//! Our Levinson-Durbin output uses the polynomial form
//! `A(z) = 1 + a_1 · z^{-1} + … + a_M · z^{-M}` (see
//! [`crate::levinson`]). The spec's predictor convention writes
//! `s(n) ≈ Σ ã_i · s(n-i)`, i.e. `ã_i = -a_i` of our Levinson output.
//! [`ShortTermPostfilter::set_from_synthesis_byproduct`] applies the
//! sign flip while it computes `b̄_i` and `ā_i`.
//!
//! ## Update timing
//!
//! Per §7.2 the coefficients update at the **first** vector of each
//! adaptation cycle (`ICOUNT = 1`), reading the `ã_i` / `k1` produced
//! by the most recent Levinson run (which lands at `ICOUNT = 3` of the
//! *previous* cycle in real-time hardware). In our batched cycle
//! model the producer commits its outputs at the cycle boundary and
//! the consumer reads them at the next vector — equivalent timing.
//!
//! ## State carry-over
//!
//! The pole-zero filter and the tilt all-zero filter both carry
//! per-tap memory across vector boundaries (the spec speaks of "the
//! filter memory") — initialised to zero per the Table 2/G.728
//! cold-start convention applied uniformly to every internal filter
//! in the codec.
//!
//! ## §4.6.1 non-speech "postfilter off" path
//!
//! When the postfilter is disabled the AGC sees `sf = sd` and
//! settles to unity (see [`crate::Agc`] docs). This module's
//! [`ShortTermPostfilter::filter_vector`] keeps producing real output
//! even with all coefficients zeroed — the cascade reduces to the
//! identity `sf = sd` (numerator = denominator = 1, µ = 0).

use crate::consts::{IDIM, SPFPCF, SPFZCF, TILTF};
use crate::synthesis_adapter::PF_LPC_ORDER;

/// Short-term postfilter (block 72) state.
///
/// One instance per decoder. Carries:
///
/// * The bandwidth-expanded zero / pole coefficient arrays `b̄_i` and
///   `ā_i` for `1 ≤ i ≤ 10`, plus the tilt scalar `µ`. Indexed 1-based
///   in the spec; arrays are sized `[f64; PF_LPC_ORDER + 1]` and
///   `[1]` is the first coefficient (`[0]` is unused for symmetry
///   with the spec layout — the leading "1" of `(1 - Σ ā_i z^{-i})`
///   is implicit in the difference equation).
/// * Per-sample memory for the pole-zero cascade `zero_mem` /
///   `pole_mem` (10 taps each) and the tilt filter `tilt_mem` (1 tap).
#[derive(Debug, Clone)]
pub struct ShortTermPostfilter {
    /// All-zero numerator coefficients `b̄_i` (1 ≤ i ≤ 10).
    /// `b_bar[0]` is unused (placeholder for the implicit leading 1).
    b_bar: [f64; PF_LPC_ORDER + 1],
    /// All-pole denominator coefficients `ā_i` (1 ≤ i ≤ 10).
    /// `a_bar[0]` is unused.
    a_bar: [f64; PF_LPC_ORDER + 1],
    /// Tilt-compensation factor `µ = TILTF · k1`.
    mu: f64,
    /// Memory of the all-zero numerator filter: `sd(n-1)..sd(n-10)`.
    zero_mem: [f64; PF_LPC_ORDER],
    /// Memory of the all-pole denominator filter: previous outputs of
    /// the pole-zero stage (i.e. tilt-filter inputs).
    pole_mem: [f64; PF_LPC_ORDER],
    /// One-tap memory of the tilt all-zero filter: previous pole-zero
    /// output (the tilt filter's input from the previous sample).
    tilt_mem: f64,
}

impl Default for ShortTermPostfilter {
    fn default() -> Self {
        Self::new()
    }
}

impl ShortTermPostfilter {
    /// Construct a fresh short-term postfilter with all coefficients
    /// zeroed and all memories at the Table-2 cold-start value of
    /// `0.0`.
    ///
    /// With all coefficients zero the filter reduces to the identity
    /// (`sf = sd`); the first call to [`Self::filter_vector`] before
    /// any [`Self::set_from_synthesis_byproduct`] therefore passes
    /// `sd` through unchanged, matching the spec's "postfilter off"
    /// §4.6.1 cold-start behaviour.
    pub fn new() -> Self {
        Self {
            b_bar: [0.0; PF_LPC_ORDER + 1],
            a_bar: [0.0; PF_LPC_ORDER + 1],
            mu: 0.0,
            zero_mem: [0.0; PF_LPC_ORDER],
            pole_mem: [0.0; PF_LPC_ORDER],
            tilt_mem: 0.0,
        }
    }

    /// Read-only access to the bandwidth-expanded numerator
    /// coefficients (for tests and audit). `[0]` is unused; the
    /// useful entries are `[1..=10]`.
    pub fn numerator(&self) -> &[f64; PF_LPC_ORDER + 1] {
        &self.b_bar
    }

    /// Read-only access to the bandwidth-expanded denominator
    /// coefficients. Same layout as [`Self::numerator`].
    pub fn denominator(&self) -> &[f64; PF_LPC_ORDER + 1] {
        &self.a_bar
    }

    /// Current tilt-compensation factor `µ = TILTF · k1`.
    pub fn tilt(&self) -> f64 {
        self.mu
    }

    /// Reset all coefficient values and all per-sample memories back
    /// to the cold-start defaults. Useful for tests.
    pub fn reset(&mut self) {
        self.b_bar = [0.0; PF_LPC_ORDER + 1];
        self.a_bar = [0.0; PF_LPC_ORDER + 1];
        self.mu = 0.0;
        self.zero_mem = [0.0; PF_LPC_ORDER];
        self.pole_mem = [0.0; PF_LPC_ORDER];
        self.tilt_mem = 0.0;
    }

    /// Refresh the postfilter coefficients from the order-10 by-
    /// product Levinson output of the synthesis-filter adapter.
    ///
    /// * `a10` — `[1.0, a_1, …, a_10]` in our Levinson `A` layout
    ///   (see [`crate::levinson`]).
    /// * `k1`  — first reflection coefficient `−R(2) / R(1)` from
    ///   the same RTMP.
    ///
    /// Computes:
    ///
    /// * `b̄_i = ã_i · SPFZCF^i`  with  `ã_i = -a_i`  (eq. 4-3),
    /// * `ā_i = ã_i · SPFPCF^i`  with the same sign flip (eq. 4-4),
    /// * `µ   = TILTF · k1`                              (eq. 4-5).
    ///
    /// Per §7.2 this should be called once at the first vector of
    /// each adaptation cycle.
    pub fn set_from_synthesis_byproduct(&mut self, a10: &[f64; PF_LPC_ORDER + 1], k1: f64) {
        // λ_z^i and λ_p^i — compute iteratively (one multiply per
        // step) instead of `.powi(i)`, both for speed and to keep
        // the floating-point path identical to the spec's tabulated
        // approach (SPFPCFV / SPFZCFV in Q14 are exactly λ_z^{i-1}
        // / λ_p^{i-1}; we want λ_z^i / λ_p^i, hence one extra
        // multiply per index).
        let mut lz_pow = SPFZCF;
        let mut lp_pow = SPFPCF;
        for i in 1..=PF_LPC_ORDER {
            // Spec sign convention: ã_i = -a_i (our Levinson A).
            let atilde_i = -a10[i];
            self.b_bar[i] = atilde_i * lz_pow;
            self.a_bar[i] = atilde_i * lp_pow;
            lz_pow *= SPFZCF;
            lp_pow *= SPFPCF;
        }
        // b_bar[0] and a_bar[0] are unused (implicit leading 1 in
        // both polynomials); leave them at zero.
        self.b_bar[0] = 0.0;
        self.a_bar[0] = 0.0;
        self.mu = TILTF * k1;
    }

    /// Filter one IDIM-sample synthesis-output vector `sd` through
    /// the cascade `H_s(z) · (1 + µ · z^{-1})` and return the post-
    /// filtered vector `sf`.
    ///
    /// Per-sample dataflow:
    ///
    /// 1. **All-zero numerator.** `w(n) = sd(n) - Σ b̄_i · sd(n-i)`,
    ///    where `sd(n-1..n-10)` come from `zero_mem`.
    /// 2. **All-pole denominator.** `y(n) = w(n) + Σ ā_i · y(n-i)`,
    ///    where `y(n-1..n-10)` come from `pole_mem`.
    /// 3. **Tilt compensation.** `sf(n) = y(n) + µ · y(n-1)`, with
    ///    `y(n-1)` from `tilt_mem`.
    ///
    /// Both filter memories advance one sample per output sample.
    pub fn filter_vector(&mut self, sd: &[f64; IDIM]) -> [f64; IDIM] {
        let mut out = [0.0f64; IDIM];
        for k in 0..IDIM {
            // ----- All-zero numerator: w(n) = sd(n) - Σ b̄_i sd(n-i) -
            let mut w = sd[k];
            for i in 1..=PF_LPC_ORDER {
                // zero_mem[i-1] = sd(n-i); the most recent past sample
                // is at index 0, oldest at index 9.
                w -= self.b_bar[i] * self.zero_mem[i - 1];
            }

            // ----- All-pole denominator: y(n) = w + Σ ā_i y(n-i) -----
            let mut y = w;
            for i in 1..=PF_LPC_ORDER {
                y += self.a_bar[i] * self.pole_mem[i - 1];
            }

            // ----- Tilt compensation: sf(n) = y(n) + µ · y(n-1) ------
            let sf = y + self.mu * self.tilt_mem;

            // ----- Advance per-sample state ------------------------
            // Push sd(n) into zero_mem, shifting older samples one
            // slot deeper. Likewise push y(n) into pole_mem and
            // record y(n) as the new tilt input for the next sample.
            for i in (1..PF_LPC_ORDER).rev() {
                self.zero_mem[i] = self.zero_mem[i - 1];
                self.pole_mem[i] = self.pole_mem[i - 1];
            }
            self.zero_mem[0] = sd[k];
            self.pole_mem[0] = y;
            self.tilt_mem = y;

            out[k] = sf;
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---------- Cold-start invariants ----------------------------------

    #[test]
    fn fresh_filter_is_identity() {
        // All coefficients zero → numerator = denominator = 1, µ = 0.
        // First IDIM-sample input passes through unchanged.
        let mut pf = ShortTermPostfilter::new();
        let sd = [100.0, -50.0, 0.0, 25.0, 80.0];
        let sf = pf.filter_vector(&sd);
        for k in 0..IDIM {
            assert_eq!(sf[k], sd[k], "k={k}");
        }
    }

    #[test]
    fn fresh_coefficients_are_zero() {
        let pf = ShortTermPostfilter::new();
        assert!(pf.numerator().iter().all(|&v| v == 0.0));
        assert!(pf.denominator().iter().all(|&v| v == 0.0));
        assert_eq!(pf.tilt(), 0.0);
    }

    // ---------- Bandwidth expansion math -------------------------------

    #[test]
    fn coefficients_apply_correct_lambda_powers() {
        // Drive a known order-10 predictor through set_from_synthesis
        // and verify each coefficient is `ã_i · λ^i` with the spec's
        // SPFZCF = 0.65 and SPFPCF = 0.75.
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        // Use a non-trivial predictor, no special structure.
        for i in 1..=PF_LPC_ORDER {
            a10[i] = 0.1 * i as f64; // arbitrary
        }
        let k1 = 0.42;
        let mut pf = ShortTermPostfilter::new();
        pf.set_from_synthesis_byproduct(&a10, k1);
        // Spec ã_i = -a_i; expect b̄_i = -a_i · 0.65^i; ā_i = -a_i · 0.75^i.
        let mut expected_lz = SPFZCF;
        let mut expected_lp = SPFPCF;
        for i in 1..=PF_LPC_ORDER {
            let expected_b = -a10[i] * expected_lz;
            let expected_a = -a10[i] * expected_lp;
            assert!(
                (pf.numerator()[i] - expected_b).abs() < 1e-12,
                "b_bar[{i}] = {} vs expected {}",
                pf.numerator()[i],
                expected_b
            );
            assert!(
                (pf.denominator()[i] - expected_a).abs() < 1e-12,
                "a_bar[{i}] = {} vs expected {}",
                pf.denominator()[i],
                expected_a
            );
            expected_lz *= SPFZCF;
            expected_lp *= SPFPCF;
        }
        assert!((pf.tilt() - TILTF * k1).abs() < 1e-12);
    }

    #[test]
    fn tilt_factor_uses_tiltf_constant() {
        // Direct invariant: µ = TILTF · k1 exactly.
        let a10 = [1.0; PF_LPC_ORDER + 1];
        let mut pf = ShortTermPostfilter::new();
        for &k1 in &[0.0_f64, 0.5, -0.7, 0.999] {
            pf.set_from_synthesis_byproduct(&a10, k1);
            assert!((pf.tilt() - TILTF * k1).abs() < 1e-12);
        }
    }

    // ---------- Filter behaviour on non-trivial coefficients ----------

    #[test]
    fn nonzero_coefficients_change_output() {
        // Once real coefficients are set, the output is no longer the
        // identity of the input.
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        a10[1] = -0.5;
        a10[2] = 0.3;
        let k1 = -0.5;
        let mut pf = ShortTermPostfilter::new();
        pf.set_from_synthesis_byproduct(&a10, k1);
        let sd = [200.0, 100.0, -50.0, 80.0, -20.0];
        let sf = pf.filter_vector(&sd);
        // At least one sample should differ from the identity output.
        // The very first sample is `sd[0] + b̄_1 · 0 + ā_1 · 0 + µ · 0 = sd[0]`
        // because all memories are zero at cold start; subsequent samples
        // see non-zero memory and diverge.
        assert!(sf.iter().zip(&sd).any(|(&a, &b)| (a - b).abs() > 1e-12));
    }

    #[test]
    fn filter_advances_memory_state() {
        // After one filter_vector call the memories should be the
        // most-recent input/output samples — not all-zero anymore.
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        a10[1] = -0.5;
        let k1 = 0.0;
        let mut pf = ShortTermPostfilter::new();
        pf.set_from_synthesis_byproduct(&a10, k1);
        let sd = [10.0, 20.0, 30.0, 40.0, 50.0];
        let _ = pf.filter_vector(&sd);
        // zero_mem[0] should now hold the last sample of sd (50.0).
        assert_eq!(pf.zero_mem[0], 50.0);
        // pole_mem[0] should be nonzero (driven by the input).
        assert_ne!(pf.pole_mem[0], 0.0);
    }

    #[test]
    fn finite_output_under_realistic_run() {
        // Smoke: drive 64 vectors of synthetic speech through and check
        // the output is always finite and bounded.
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        a10[1] = -0.5;
        a10[2] = 0.2;
        let k1 = -0.5;
        let mut pf = ShortTermPostfilter::new();
        pf.set_from_synthesis_byproduct(&a10, k1);
        for v in 0..64 {
            let mut sd = [0.0; IDIM];
            for k in 0..IDIM {
                let t = (v * IDIM + k) as f64;
                sd[k] = 1000.0 * (0.4 * t).sin();
            }
            let sf = pf.filter_vector(&sd);
            for &s in &sf {
                assert!(s.is_finite());
            }
        }
    }

    #[test]
    fn reset_returns_filter_to_identity_state() {
        // After driving the filter into a non-trivial state, reset
        // should make it pass the next vector through unchanged.
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        a10[1] = -0.5;
        let k1 = 0.5;
        let mut pf = ShortTermPostfilter::new();
        pf.set_from_synthesis_byproduct(&a10, k1);
        let sd = [10.0, 20.0, 30.0, 40.0, 50.0];
        let _ = pf.filter_vector(&sd);
        pf.reset();
        let sf = pf.filter_vector(&sd);
        for k in 0..IDIM {
            assert_eq!(sf[k], sd[k]);
        }
    }

    // ---------- Steady-state DC properties -----------------------------

    #[test]
    fn dc_input_produces_finite_dc_output() {
        // A constant DC input keeps the IIR taps finite (the spec's
        // bandwidth expansion ensures the all-pole denominator stays
        // strictly inside the unit circle). After many vectors the
        // output should still be finite.
        let mut a10 = [0.0; PF_LPC_ORDER + 1];
        a10[0] = 1.0;
        a10[1] = -0.3;
        a10[3] = 0.1;
        let mut pf = ShortTermPostfilter::new();
        pf.set_from_synthesis_byproduct(&a10, 0.0);
        let sd = [1.0_f64; IDIM];
        for _ in 0..2048 {
            let sf = pf.filter_vector(&sd);
            for &s in &sf {
                assert!(s.is_finite());
            }
        }
    }

    #[test]
    fn zero_coefficients_pass_through_for_all_inputs() {
        // §4.6.1 postfilter-off: with all coefficients zero the filter
        // is the identity for every input vector.
        let mut pf = ShortTermPostfilter::new();
        for v in 0..16 {
            let sd = [
                v as f64,
                v as f64 * 2.0,
                v as f64 * -3.0,
                v as f64 + 7.0,
                v as f64 - 1.0,
            ];
            let sf = pf.filter_vector(&sd);
            for k in 0..IDIM {
                assert_eq!(sf[k], sd[k]);
            }
        }
    }
}
