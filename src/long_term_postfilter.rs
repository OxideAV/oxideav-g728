//! Long-term (pitch) postfilter — block 71 of Figure 7/G.728.
//!
//! Implements the §4.6 long-term postfilter transfer function
//! (equation 4-1):
//!
//! ```text
//!     H_l(z)  =  g_l · (1 + b · z^{-p})
//! ```
//!
//! where:
//!
//! * `p`   — fundamental pitch period in samples, `KPMIN ≤ p ≤ KPMAX`
//!   (Table 1/G.728: `KPMIN = 20`, `KPMAX = 140`),
//! * `b`   — pitch-predictor scaled tap coefficient,
//! * `g_l` — overall comb-filter gain.
//!
//! In the spec the three values `(g_l, b, p)` are produced by §4.7
//! blocks 81..84 — the 10th-order LPC inverse filter (block 81), the
//! pitch period extractor (block 82), the pitch predictor tap
//! calculator (block 83) and the long-term postfilter coefficient
//! calculator (block 84) — and refreshed at the **third** vector of
//! every adaptation cycle. This module owns just the comb filter
//! itself; callers supply `(g_l, b, p)` via
//! [`LongTermPostfilter::set_coefficients`].
//!
//! ## §4.6.1 "postfilter off" behaviour
//!
//! When the §4.7 pitch detector reports an unvoiced or weakly-voiced
//! frame (`β < PPFTH = 0.6`, decode trace §7.1; equations 4-13/4-14),
//! the spec turns the long-term stage off by setting `b = 0` and
//! `g_l = 1`. With `b = 0` the `b · z^{-p}` term vanishes and the gain
//! is unity, so `H_l(z) = 1` and the comb filter is the identity.
//! [`LongTermPostfilter::new`] starts in exactly this state; until
//! [`LongTermPostfilter::set_coefficients`] is called with a
//! non-trivial `b ≠ 0`, the cascade is `sf = sd`.
//!
//! Both the §4.7 pitch-extraction front end (blocks 81..84) and the
//! §4.7 coefficient calculator (block 84) land in later rounds. Until
//! they do, the [`Decoder`](crate::Decoder) keeps `(b, g_l) = (0, 1)`
//! and `p = KPMIN` so the long-term stage is a provable identity
//! filter regardless of any decoded speech.
//!
//! ## State carry-over
//!
//! The comb filter `(1 + b · z^{-p})` is purely all-zero (FIR): every
//! output sample depends on the input sample `p` samples ago. The
//! filter therefore carries a single delay line of length `KPMAX`
//! samples (the maximum pitch period). On every input sample we read
//! one tap at index `p` and push the new input. The line is
//! initialised to zero per the Table 2/G.728 cold-start convention
//! applied uniformly to every internal filter in the codec.

use crate::consts::{IDIM, KPMAX, KPMIN};

/// Long-term (pitch) postfilter (block 71) state.
///
/// One instance per decoder. Carries:
///
/// * The comb-filter coefficients `(g_l, b, p)` — set externally via
///   [`Self::set_coefficients`]. Cold-start values are
///   `g_l = 1.0, b = 0.0, p = KPMIN` (the spec's §4.6.1 "postfilter
///   off" passthrough).
/// * A `KPMAX`-sample circular delay line `delay` holding the most
///   recent `KPMAX` input samples (`sd(n-1)..sd(n-KPMAX)`), indexed
///   by `delay_idx` (cursor of the next slot to write).
#[derive(Debug, Clone)]
pub struct LongTermPostfilter {
    /// Overall comb-filter gain `g_l` (eq. 4-1).
    g_l: f64,
    /// Pitch-predictor scaled tap `b` (eq. 4-1).
    b: f64,
    /// Fundamental pitch period `p`, clamped to `[KPMIN, KPMAX]`.
    p: usize,
    /// Circular delay line carrying the last `KPMAX` input samples.
    /// `delay[(delay_idx + KPMAX - k) % KPMAX]` is the input sample
    /// `k` samples in the past (`k = 1..=KPMAX`).
    delay: [f64; KPMAX],
    /// Cursor: the next slot in `delay` to write the *current* input
    /// sample. After writing, advance by one (mod `KPMAX`).
    delay_idx: usize,
}

impl Default for LongTermPostfilter {
    fn default() -> Self {
        Self::new()
    }
}

impl LongTermPostfilter {
    /// Construct a fresh long-term postfilter on the §4.6.1
    /// "postfilter off" passthrough setting:
    ///
    /// * `g_l = 1.0`,
    /// * `b   = 0.0`,
    /// * `p   = KPMIN`,
    ///
    /// with an all-zero delay line. Until [`Self::set_coefficients`]
    /// is called with a non-zero `b`, the filter is the identity
    /// (`sf = sd`).
    pub fn new() -> Self {
        Self {
            g_l: 1.0,
            b: 0.0,
            p: KPMIN,
            delay: [0.0; KPMAX],
            delay_idx: 0,
        }
    }

    /// Read-only access to `g_l` (for tests and audit).
    pub fn g_l(&self) -> f64 {
        self.g_l
    }

    /// Read-only access to `b` (for tests and audit).
    pub fn b(&self) -> f64 {
        self.b
    }

    /// Read-only access to the pitch period `p` (for tests and audit).
    pub fn p(&self) -> usize {
        self.p
    }

    /// Reset all coefficient values back to the §4.6.1 "off"
    /// passthrough and clear the delay line. Useful for tests.
    pub fn reset(&mut self) {
        self.g_l = 1.0;
        self.b = 0.0;
        self.p = KPMIN;
        self.delay = [0.0; KPMAX];
        self.delay_idx = 0;
    }

    /// Set the comb-filter coefficients `(g_l, b, p)`.
    ///
    /// In the full §4.7 pipeline this is called once per adaptation
    /// cycle, at the third vector (see decode-trace §7.1 and
    /// equations 4-12 / 4-13 / 4-14), with `(g_l, b, p)` computed by
    /// block 84.
    ///
    /// Until that pipeline lands, callers should keep the §4.6.1
    /// passthrough by leaving the filter at its cold-start defaults,
    /// or by calling `set_coefficients(1.0, 0.0, KPMIN)` explicitly.
    ///
    /// `p` is clamped to `[KPMIN, KPMAX]` per Table 1/G.728. Out-of-
    /// range values are saturated to the nearest endpoint rather than
    /// rejected — the spec's pitch extractor (block 82) never emits a
    /// lag outside this band, so saturating is a safe transcription
    /// guard rather than a behavioural change.
    pub fn set_coefficients(&mut self, g_l: f64, b: f64, p: usize) {
        self.g_l = g_l;
        self.b = b;
        self.p = p.clamp(KPMIN, KPMAX);
    }

    /// Filter one IDIM-sample synthesis-output vector `sd` through the
    /// comb cascade `g_l · (1 + b · z^{-p})` and return the
    /// long-term-postfiltered vector.
    ///
    /// Per-sample dataflow:
    ///
    /// 1. Read the past sample `sd(n - p)` from the circular delay
    ///    line. Indexing convention: the delay line stores the most
    ///    recent `KPMAX` input samples `sd(n-1)..sd(n-KPMAX)`; the
    ///    sample `p` samples in the past lives at offset `p` back
    ///    from the write cursor.
    /// 2. Compute `y(n) = g_l · (sd(n) + b · sd(n - p))`.
    /// 3. Push `sd(n)` into the delay line at `delay_idx`, then
    ///    advance `delay_idx` by one (mod `KPMAX`).
    ///
    /// The delay line advances one sample per output sample,
    /// preserving the FIR memory across vector boundaries.
    pub fn filter_vector(&mut self, sd: &[f64; IDIM]) -> [f64; IDIM] {
        let mut out = [0.0f64; IDIM];
        for k in 0..IDIM {
            // ----- Read sd(n - p) from the circular delay line ---------
            // The slot we are about to write (`delay_idx`) holds the
            // OLDEST sample currently in memory: it will be overwritten
            // with the current input. The sample `p` samples in the
            // past sits at `(delay_idx + KPMAX - p) % KPMAX`, because:
            //   - the slot just before `delay_idx` (offset 1) holds
            //     sd(n - 1),
            //   - the slot `p` slots before `delay_idx` holds sd(n - p).
            let past_idx = (self.delay_idx + KPMAX - self.p) % KPMAX;
            let sd_past = self.delay[past_idx];

            // ----- Eq. 4-1: y(n) = g_l · (sd(n) + b · sd(n - p)) ------
            out[k] = self.g_l * (sd[k] + self.b * sd_past);

            // ----- Push sd(n) into the delay line ---------------------
            self.delay[self.delay_idx] = sd[k];
            self.delay_idx = (self.delay_idx + 1) % KPMAX;
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
        // §4.6.1 cold start: g_l = 1, b = 0 → H_l(z) = 1 (identity).
        let mut pf = LongTermPostfilter::new();
        let sd = [100.0, -50.0, 0.0, 25.0, 80.0];
        let sf = pf.filter_vector(&sd);
        for k in 0..IDIM {
            assert_eq!(sf[k], sd[k], "k={k}");
        }
    }

    #[test]
    fn fresh_coefficients_are_passthrough() {
        let pf = LongTermPostfilter::new();
        assert_eq!(pf.g_l(), 1.0);
        assert_eq!(pf.b(), 0.0);
        assert_eq!(pf.p(), KPMIN);
    }

    #[test]
    fn identity_holds_for_many_vectors() {
        // Drive 64 vectors through the cold-start filter and confirm
        // the output matches the input exactly — even though the
        // delay line is being populated, `b = 0` zeroes its
        // contribution.
        let mut pf = LongTermPostfilter::new();
        for v in 0..64 {
            let sd = [
                v as f64,
                v as f64 * 2.0,
                v as f64 - 1.0,
                v as f64 * -0.5,
                v as f64 + 7.0,
            ];
            let sf = pf.filter_vector(&sd);
            for k in 0..IDIM {
                assert_eq!(sf[k], sd[k], "vector {v} sample {k}");
            }
        }
    }

    // ---------- Coefficient setter -------------------------------------

    #[test]
    fn set_coefficients_clamps_pitch_into_range() {
        let mut pf = LongTermPostfilter::new();
        // Below KPMIN saturates up.
        pf.set_coefficients(0.5, 0.1, 0);
        assert_eq!(pf.p(), KPMIN);
        // Above KPMAX saturates down.
        pf.set_coefficients(0.5, 0.1, KPMAX + 100);
        assert_eq!(pf.p(), KPMAX);
        // Inside the band passes through unchanged.
        pf.set_coefficients(0.5, 0.1, 50);
        assert_eq!(pf.p(), 50);
    }

    #[test]
    fn set_coefficients_stores_gain_and_tap() {
        let mut pf = LongTermPostfilter::new();
        pf.set_coefficients(0.875, 0.13, 42);
        assert_eq!(pf.g_l(), 0.875);
        assert_eq!(pf.b(), 0.13);
        assert_eq!(pf.p(), 42);
    }

    // ---------- Eq. 4-1 verification on a synthetic impulse trail -----

    #[test]
    fn comb_filter_adds_delayed_input_scaled_by_b_and_g_l() {
        // Set p = KPMIN = 20 and b = 0.5, g_l = 1.0. Feed a 25-sample
        // unit-impulse-at-sample-0 sequence (5 vectors of IDIM = 5
        // samples). The output is:
        //   out[0]  = g_l · (sd[0] + b · sd[-20])  = 1 · (1 + 0.5·0) = 1
        //   out[1..19] = 0
        //   out[20] = g_l · (sd[20] + b · sd[0])   = 1 · (0 + 0.5·1) = 0.5
        //   out[21..24] = 0
        // The impulse re-emerges p = 20 samples later, scaled by b.
        let mut pf = LongTermPostfilter::new();
        pf.set_coefficients(1.0, 0.5, KPMIN);
        // 5 vectors × 5 samples = 25 samples total.
        let mut outs = [0.0f64; 25];
        for v in 0..5 {
            let mut sd = [0.0f64; IDIM];
            if v == 0 {
                sd[0] = 1.0;
            }
            let sf = pf.filter_vector(&sd);
            for k in 0..IDIM {
                outs[v * IDIM + k] = sf[k];
            }
        }
        assert_eq!(outs[0], 1.0, "impulse passes through immediately");
        for i in 1..KPMIN {
            assert_eq!(outs[i], 0.0, "no echo before lag p");
        }
        assert_eq!(outs[KPMIN], 0.5, "echo at lag p, scaled by b");
        for i in (KPMIN + 1)..25 {
            assert_eq!(outs[i], 0.0, "no further echo (FIR)");
        }
    }

    #[test]
    fn gain_g_l_scales_every_output_sample_uniformly() {
        // With b = 0 the filter reduces to a pure gain `y(n) = g_l · sd(n)`;
        // confirm it scales every sample by exactly `g_l`.
        let mut pf = LongTermPostfilter::new();
        let g_l = 0.625;
        pf.set_coefficients(g_l, 0.0, KPMIN);
        let sd = [10.0, -20.0, 30.5, 0.0, -7.5];
        let sf = pf.filter_vector(&sd);
        for k in 0..IDIM {
            assert!((sf[k] - g_l * sd[k]).abs() < 1e-12);
        }
    }

    #[test]
    fn b_equals_zero_is_pure_gain_regardless_of_p() {
        // With b = 0 the choice of p must not affect the output —
        // the second term in eq. 4-1 vanishes.
        let mut a = LongTermPostfilter::new();
        let mut b = LongTermPostfilter::new();
        a.set_coefficients(0.5, 0.0, KPMIN);
        b.set_coefficients(0.5, 0.0, KPMAX);
        for v in 0..16 {
            let sd = [
                v as f64,
                v as f64 + 1.0,
                v as f64 - 2.0,
                v as f64 * 0.5,
                -(v as f64),
            ];
            let sf_a = a.filter_vector(&sd);
            let sf_b = b.filter_vector(&sd);
            for k in 0..IDIM {
                assert!((sf_a[k] - sf_b[k]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn comb_filter_at_kpmax_echoes_at_correct_lag() {
        // Same shape as the KPMIN test but with p = KPMAX (the longest
        // pitch period in the spec's band): the echo must appear
        // exactly KPMAX samples after the impulse, never sooner and
        // never later.
        let mut pf = LongTermPostfilter::new();
        let b = 0.3_f64;
        pf.set_coefficients(1.0, b, KPMAX);
        // We need KPMAX + 1 samples to see the echo; do 30 vectors of
        // 5 = 150 samples (> KPMAX = 140).
        let total = 30 * IDIM;
        let mut outs = vec![0.0f64; total];
        for v in 0..30 {
            let mut sd = [0.0f64; IDIM];
            if v == 0 {
                sd[0] = 1.0;
            }
            let sf = pf.filter_vector(&sd);
            for k in 0..IDIM {
                outs[v * IDIM + k] = sf[k];
            }
        }
        assert_eq!(outs[0], 1.0);
        for i in 1..KPMAX {
            assert_eq!(outs[i], 0.0, "echo arrived too soon at i={i}");
        }
        assert!(
            (outs[KPMAX] - b).abs() < 1e-12,
            "echo at lag KPMAX = {KPMAX} should be b={b}, got {}",
            outs[KPMAX]
        );
    }

    #[test]
    fn comb_filter_is_fir_no_recursion() {
        // The long-term postfilter is FIR by construction (b · z^{-p}
        // is all-zero). Set b = 1 (perfect echo) and verify a single
        // impulse produces exactly TWO non-zero output samples in any
        // 2·p-sample window — no second-generation recursion.
        let mut pf = LongTermPostfilter::new();
        pf.set_coefficients(1.0, 1.0, KPMIN);
        let total = 3 * KPMIN; // covers two would-be recursive hops.
        let n_vecs = total.div_ceil(IDIM);
        let mut outs = vec![0.0f64; n_vecs * IDIM];
        for v in 0..n_vecs {
            let mut sd = [0.0f64; IDIM];
            if v == 0 {
                sd[0] = 1.0;
            }
            let sf = pf.filter_vector(&sd);
            for k in 0..IDIM {
                outs[v * IDIM + k] = sf[k];
            }
        }
        let nz = outs.iter().filter(|&&v| v != 0.0).count();
        assert_eq!(nz, 2, "FIR filter must emit exactly 2 nonzero samples");
        // And the two non-zero samples must sit at lags 0 and KPMIN.
        assert_eq!(outs[0], 1.0);
        assert_eq!(outs[KPMIN], 1.0);
    }

    // ---------- Reset ---------------------------------------------------

    #[test]
    fn reset_restores_passthrough() {
        // Drive into a non-trivial state, reset, and confirm the next
        // vector emerges unchanged.
        let mut pf = LongTermPostfilter::new();
        pf.set_coefficients(0.5, 0.3, 50);
        let sd = [10.0, 20.0, 30.0, 40.0, 50.0];
        let _ = pf.filter_vector(&sd);
        pf.reset();
        let sf = pf.filter_vector(&sd);
        for k in 0..IDIM {
            assert_eq!(sf[k], sd[k]);
        }
        assert_eq!(pf.g_l(), 1.0);
        assert_eq!(pf.b(), 0.0);
        assert_eq!(pf.p(), KPMIN);
    }

    // ---------- Stability / finiteness ---------------------------------

    #[test]
    fn output_finite_under_sinusoidal_drive() {
        // FIR comb filter must keep finite output for any bounded
        // input — sanity smoke at a realistic operating point.
        let mut pf = LongTermPostfilter::new();
        pf.set_coefficients(0.9, 0.15, 80);
        for v in 0..256 {
            let mut sd = [0.0f64; IDIM];
            for k in 0..IDIM {
                let t = (v * IDIM + k) as f64;
                sd[k] = 2000.0 * (0.05 * t).sin();
            }
            let sf = pf.filter_vector(&sd);
            for &s in &sf {
                assert!(s.is_finite());
                // Conservative bound: |g_l| · (|sd| + |b| · |sd|) ≤
                //                     1.0  · 2000 · (1 + 0.15) = 2300.
                assert!(s.abs() < 2400.0);
            }
        }
    }

    // ---------- Delay-line wrap-around ---------------------------------

    #[test]
    fn delay_line_wraps_correctly_past_kpmax() {
        // Run more than KPMAX samples through the filter with a
        // single impulse at the start and verify the circular delay
        // line correctly wraps without emitting a phantom second echo.
        let mut pf = LongTermPostfilter::new();
        let p = 30; // arbitrary in-range pitch period.
        pf.set_coefficients(1.0, 0.5, p);
        // 2 × KPMAX samples = ~280 samples = 56 vectors of IDIM = 5.
        let n_vecs = (2 * KPMAX).div_ceil(IDIM);
        let mut outs = vec![0.0f64; n_vecs * IDIM];
        for v in 0..n_vecs {
            let mut sd = [0.0f64; IDIM];
            if v == 0 {
                sd[0] = 1.0;
            }
            let sf = pf.filter_vector(&sd);
            for k in 0..IDIM {
                outs[v * IDIM + k] = sf[k];
            }
        }
        // Exactly one echo at lag p, no spurious echo at lag p + KPMAX
        // (which would indicate a circular-buffer aliasing bug).
        let nz_idxs: Vec<usize> = outs
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v != 0.0 { Some(i) } else { None })
            .collect();
        assert_eq!(nz_idxs, vec![0, p]);
    }
}
