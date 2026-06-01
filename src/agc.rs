//! Output gain control (AGC) — blocks 73 / 74 / 75 / 76 / 77 of
//! Figure 7/G.728 (postfilter sub-stage).
//!
//! Per Recommendation §4.6 / §4.6.1 the AGC sits at the tail of the
//! postfilter chain and forces the post-filtered speech `sf(n)` to
//! carry roughly the same vector-rate power as the raw decoded speech
//! `sd(n)`. Without it, the long-term plus short-term postfilter
//! stages can produce occasional large gain excursions that would
//! audibly modulate level.
//!
//! ## Block dataflow (one IDIM-sample vector per call)
//!
//! ```text
//!   sd(1..IDIM) ─► [73: Σ|sd|] ─┐
//!                                 ├─► [75: SUMSD / SUMSF] ─► GAINSF (per-vector)
//!   sf(1..IDIM) ─► [74: Σ|sf|] ─┘                                │
//!                                                                 ▼
//!                                                          [76: 1st-order LP]
//!                                                                 │
//!                                                                 ▼ (per-sample SCALEFIL)
//!   sf(1..IDIM) ──────────────────────────────────────► [77: × SCALEFIL] ─► spf(1..IDIM)
//! ```
//!
//! ## Block-76 filter form (§4.6, paragraph after eq. 4-5)
//!
//! > "This scaling factor is then filtered by a first-order lowpass
//! > filter 76 to get a separate scaling factor for each of the five
//! > components of sf(n). The first-order lowpass filter 76 has a
//! > transfer function of `0.01 / (1 − 0.99·z⁻¹)`."
//!
//! In difference-equation form (one update per output sample):
//!
//! ```text
//!     SCALEFIL(n) = AGCFAC · SCALEFIL(n−1) + (1 − AGCFAC) · GAINSF
//! ```
//!
//! where `AGCFAC = 0.99` is the `Table 1/G.728` constant and
//! `1 − AGCFAC = 0.01` is the numerator gain that yields a unity
//! steady-state response (`H(1) = 0.01 / (1 − 0.99) = 1`). GAINSF is
//! held constant across the five samples of the current vector.
//!
//! ## Non-speech / postfilter-off path (§4.6.1)
//!
//! > "For some non-speech signals the performance of the coder is
//! > improved when the adaptive postfilter is turned off. Since the
//! > input to the adaptive postfilter is the output of the synthesis
//! > filter, this signal is always available."
//!
//! When the postfilter is disabled the API passes `sf = sd` to
//! [`Agc::apply`]; the ratio is then `1.0`, the IIR converges to
//! `1.0`, and the per-sample scaling becomes a no-op. Calling
//! [`Agc::apply`] every vector therefore stays correct even before
//! the long-term and short-term postfilter stages land.

use crate::consts::{AGCFAC, IDIM};

/// AGC state — one [`Agc`] instance per decoder.
///
/// Carries the block-76 lowpass memory across vector boundaries; the
/// filter has a single-tap state (`SCALEFIL(n−1)`). Initialised to
/// `1.0` so the first vector emerges at unit gain before the lowpass
/// has converged — see [`Agc::new`] for the rationale.
#[derive(Debug, Clone, Copy)]
pub struct Agc {
    /// Block-76 single-tap memory: `SCALEFIL(n−1)`. Updated once per
    /// output sample inside [`Agc::apply`].
    scalefil: f64,
}

impl Default for Agc {
    fn default() -> Self {
        Self::new()
    }
}

impl Agc {
    /// Construct a fresh AGC with the lowpass memory initialised to
    /// `1.0`.
    ///
    /// The spec does not list a Table-2 initial value for this single-
    /// tap state because, on a real DSP at cold start, the first
    /// adaptation cycle of decoded speech is silence (the synthesis
    /// filter starts with zeroed memory, §3.7 / Table 2). With both
    /// `sf(n)` and `sd(n)` near zero the block-73 / block-74 sums
    /// would underflow and block 75 would divide near zero. Pinning
    /// `SCALEFIL(−1) = 1.0` keeps the first-vector output equal to the
    /// raw `sf` vector while the IIR converges; this matches the spec's
    /// stated "force the post-filtered speech to have roughly the same
    /// power as the unfiltered speech" objective in the cold-start
    /// regime where there is no audible signal to gain-control yet.
    pub fn new() -> Self {
        Self { scalefil: 1.0 }
    }

    /// Read the current block-76 memory (for tests and audit).
    pub fn scalefil(&self) -> f64 {
        self.scalefil
    }

    /// Apply the AGC to one `IDIM`-sample vector.
    ///
    /// * `sd` — current decoded-speech vector (synthesis-filter output,
    ///   block 32). This is the reference power level.
    /// * `sf` — current post-filtered speech vector (output of the
    ///   short-term postfilter, block 72). When the long-term / short-
    ///   term postfilter stages are not yet wired up, callers may pass
    ///   `sd` here unchanged — see §4.6.1 and the module docs.
    ///
    /// Returns the gain-controlled output vector `spf(1..IDIM)`. The
    /// AGC state advances by one full IDIM-sample frame.
    ///
    /// Implements blocks 73 → 74 → 75 → 76 → 77 in spec order.
    pub fn apply(&mut self, sd: &[f64; IDIM], sf: &[f64; IDIM]) -> [f64; IDIM] {
        // ----- Block 73: Σ|sd| -----------------------------------------
        // | For K = 1, 2, ..., IDIM: SUMSD = SUMSD + |SD(K)|
        let mut sumsd = 0.0f64;
        for k in 0..IDIM {
            sumsd += sd[k].abs();
        }

        // ----- Block 74: Σ|sf| -----------------------------------------
        let mut sumsf = 0.0f64;
        for k in 0..IDIM {
            sumsf += sf[k].abs();
        }

        // ----- Block 75: GAINSF = SUMSD / SUMSF -------------------------
        // §4.6 leaves the SUMSF = 0 case (silence at the postfilter
        // output but non-zero raw sd, which can only happen on a wildly
        // mis-designed postfilter) unspecified. The pragmatic safe
        // value is `GAINSF = 1.0`: it asserts "leave sf untouched at
        // this vector" and lets the lowpass converge to whatever the
        // next vectors call for. Matches the cold-start initialisation
        // chosen in `Agc::new`.
        let gainsf = if sumsf > 0.0 { sumsd / sumsf } else { 1.0 };

        // ----- Blocks 76 + 77: per-sample lowpass + multiply -----------
        // | SCALEFIL(n) = AGCFAC * SCALEFIL(n-1) + (1 - AGCFAC) * GAINSF
        // | SPF(K)      = SCALEFIL(K) * SF(K)
        //
        // The lowpass is a one-tap IIR; the (1 - AGCFAC) numerator gain
        // yields a unity DC response (H(1) = (1-AGCFAC) / (1-AGCFAC) =
        // 1), so a constant input of value `g` settles to `SCALEFIL = g`
        // — the property the AGC depends on.
        let one_minus_fac = 1.0 - AGCFAC;
        let mut out = [0.0f64; IDIM];
        for k in 0..IDIM {
            self.scalefil = AGCFAC * self.scalefil + one_minus_fac * gainsf;
            out[k] = self.scalefil * sf[k];
        }

        out
    }

    /// Reset the AGC state to its construction-time initial value.
    /// Useful for tests and for callers that drive a fresh stream
    /// through an existing decoder without rebuilding it.
    pub fn reset(&mut self) {
        self.scalefil = 1.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::AGCFAC;

    // ---------- Cold-start invariants ----------------------------------

    #[test]
    fn fresh_agc_holds_unity_scalefil() {
        // The cold-start convention pins SCALEFIL = 1.0 so the very
        // first output sample is `1.0 * sf[0]` (modulo one IIR step).
        let a = Agc::new();
        assert_eq!(a.scalefil(), 1.0);
    }

    #[test]
    fn default_is_equivalent_to_new() {
        let a = Agc::new();
        let b = Agc::default();
        assert_eq!(a.scalefil(), b.scalefil());
    }

    // ---------- Block-76 lowpass mathematical properties ----------------

    #[test]
    fn lowpass_dc_gain_is_unity() {
        // With AGCFAC = 0.99 the numerator gain is 1 - AGCFAC = 0.01;
        // H(z=1) = 0.01 / (1 - 0.99) = 1. After many vectors of a
        // constant `gainsf = g`, the IIR memory converges to `g` and
        // every output sample is `g * sf[k]`.
        let mut agc = Agc::new();
        // Drive a constant ratio of 2.0 (sd = 2 * sf) for many cycles.
        let sf = [1.0f64; IDIM];
        let sd = [2.0f64; IDIM];
        for _ in 0..1024 {
            agc.apply(&sd, &sf);
        }
        // Steady-state SCALEFIL → 2.0; allow a small numerical margin
        // for the IIR's finite settling time.
        assert!(
            (agc.scalefil() - 2.0).abs() < 1e-6,
            "expected SCALEFIL → 2.0, got {}",
            agc.scalefil()
        );
    }

    #[test]
    fn lowpass_settling_obeys_agcfac_geometric_decay() {
        // Drive a non-zero sf with zero sd: GAINSF = sumsd / sumsf =
        // 0 / 5 = 0, so the IIR input is identically 0 and the memory
        // decays geometrically by AGCFAC per sample. After ONE IDIM-
        // sample vector the memory should be exactly AGCFAC^IDIM
        // (no input contribution to the recursion).
        let mut agc = Agc::new();
        let sf = [1.0f64; IDIM]; // sumsf = 5
        let sd = [0.0f64; IDIM]; // sumsd = 0 → gainsf = 0
        agc.apply(&sd, &sf);
        let expected = AGCFAC.powi(IDIM as i32);
        assert!(
            (agc.scalefil() - expected).abs() < 1e-12,
            "after 1 vector with gainsf=0, SCALEFIL = AGCFAC^IDIM = {}, got {}",
            expected,
            agc.scalefil()
        );
    }

    // ---------- Block-75 ratio behaviour --------------------------------

    #[test]
    fn passthrough_when_sf_eq_sd() {
        // The §4.6.1 "postfilter off" path: when sf == sd, GAINSF = 1
        // and SCALEFIL stays at its cold-start value of 1.0; output
        // equals sf unchanged.
        let mut agc = Agc::new();
        let sd = [100.0, -200.0, 50.0, 0.0, 75.0];
        let sf = sd;
        let out = agc.apply(&sd, &sf);
        for k in 0..IDIM {
            assert!(
                (out[k] - sf[k]).abs() < 1e-12,
                "passthrough should be exact, got out[{}] = {} vs sf = {}",
                k,
                out[k],
                sf[k]
            );
        }
        // And SCALEFIL stays at 1.0 because the IIR input is 1.0 and
        // the initial state is 1.0 — every IIR step yields 1.0 again.
        assert_eq!(agc.scalefil(), 1.0);
    }

    #[test]
    fn zero_sf_falls_back_to_safe_gainsf_unity() {
        // SUMSF = 0 is the spec-unspecified divide-by-zero case; our
        // policy is GAINSF = 1.0 (= "leave sf untouched"). Output
        // should be zero everywhere because sf is zero.
        let mut agc = Agc::new();
        let sd = [1.0, 2.0, 3.0, 4.0, 5.0];
        let sf = [0.0f64; IDIM];
        let out = agc.apply(&sd, &sf);
        for &v in &out {
            assert_eq!(v, 0.0);
        }
    }

    // ---------- Output structural properties ----------------------------

    #[test]
    fn output_is_finite_for_finite_inputs() {
        let mut agc = Agc::new();
        let sd = [123.4, -567.8, 0.0, 99.99, -1.0];
        let sf = [12.3, 56.7, 89.0, -45.6, 7.8];
        for _ in 0..32 {
            let out = agc.apply(&sd, &sf);
            for &v in &out {
                assert!(v.is_finite());
            }
        }
    }

    #[test]
    fn reset_restores_unity_scalefil() {
        // After driving the AGC into a non-trivial steady state, reset
        // must take SCALEFIL back to the cold-start value.
        let mut agc = Agc::new();
        let sd = [10.0f64; IDIM];
        let sf = [1.0f64; IDIM];
        for _ in 0..512 {
            agc.apply(&sd, &sf);
        }
        assert!(
            (agc.scalefil() - 10.0).abs() < 1e-6,
            "pre-reset SCALEFIL should have settled near 10.0"
        );
        agc.reset();
        assert_eq!(agc.scalefil(), 1.0);
    }

    #[test]
    fn applied_gain_matches_per_sample_iir_state() {
        // Direct invariant: out[k] = scalefil_after_step_k * sf[k].
        // We can reconstruct the IIR trajectory and check sample-by
        // sample.
        let mut agc = Agc::new();
        let sf = [2.0f64; IDIM];
        let sd = [6.0f64; IDIM]; // gainsf = 30 / 10 = 3.0
        let mut expected_state = agc.scalefil();
        let one_minus_fac = 1.0 - AGCFAC;
        let gainsf = 3.0f64;
        let out = agc.apply(&sd, &sf);
        for k in 0..IDIM {
            expected_state = AGCFAC * expected_state + one_minus_fac * gainsf;
            let expected_out = expected_state * sf[k];
            assert!(
                (out[k] - expected_out).abs() < 1e-12,
                "sample {}: expected {}, got {}",
                k,
                expected_out,
                out[k]
            );
        }
    }

    #[test]
    fn many_vectors_with_constant_ratio_do_not_drift() {
        // The IIR should converge to GAINSF and stay there indefinitely
        // without oscillation or drift.
        let mut agc = Agc::new();
        let sd = [3.0f64; IDIM];
        let sf = [1.0f64; IDIM]; // ratio 15/5 = 3
        for _ in 0..4096 {
            agc.apply(&sd, &sf);
        }
        let s1 = agc.scalefil();
        for _ in 0..4096 {
            agc.apply(&sd, &sf);
        }
        let s2 = agc.scalefil();
        assert!(
            (s1 - s2).abs() < 1e-12,
            "steady-state should not drift: {} vs {}",
            s1,
            s2
        );
        assert!((s2 - 3.0).abs() < 1e-9);
    }
}
