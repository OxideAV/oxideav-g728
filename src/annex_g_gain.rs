//! Annex G §G.2.1 — reformulated backward vector gain adapter (block 20).
//!
//! ITU-T G.728 Annex G (1994-11) replaces the once-per-vector logarithm
//! of the floating-point gain adapter (blocks 39 / 40 in Figure 6/G.728:
//! square-RMS of the gain-scaled excitation `e(n)` followed by a runtime
//! `10·log10`) with a *mathematically equivalent* method that reads the
//! gain and shape codebook indices chosen for the previous vector and
//! looks their dB contributions up in two precomputed tables (Figure
//! G.1, blocks 93 / 94). The two methods produce the same
//! offset-removed log-gain `δ(n−1)`; the reformulation just avoids the
//! per-vector logarithm and a division, which are expensive on a
//! fixed-point DSP (§G.2.1 advantage (a)).
//!
//! ## Derivation (§G.2.1, equations G-5 … G-14)
//!
//! With gain index `i` and shape index `j` chosen for vector `n`, the
//! emitted excitation is `e(n) = σ(n)·g_i·y_j` (eq. G-5). Its power is
//!
//! ```text
//! P[e(n)] = (1/5)·Σ_k e_k(n)²
//!         = σ²(n)·g_i²·( (1/5)·Σ_k y_{jk}² )    (eq. G-10 … G-11)
//!         = σ²(n)·g_i²·P[y_j]                    (eq. G-12)
//! ```
//!
//! where `P[x]` is the *power* of a vector — its energy divided by the
//! vector dimension `IDIM = 5`. Substituting eq. G-12 into the
//! floating-point definition `δ(n) = 10·log10·P[e(n)] − 32` (eq. G-8)
//! and using `σ(n) = 10^((δ̂(n)+32)/20)` (eq. G-4) gives the
//! reformulation:
//!
//! ```text
//! δ(n) = δ̂(n) + 20·log10|g_i| + 10·log10·P[y_j]   (eq. G-14, corrected)
//! ```
//!
//! i.e. the predicted log-gain `δ̂(n)` plus two "correction terms":
//!
//! 1. `20·log10|g_i|` — the dB value of the chosen gain level
//!    ([`GAIN_LOG_DB`], block 93). Only the four distinct magnitudes
//!    `|GQ[i]|` exist, so this is an 8-entry table indexed by the gain
//!    codebook index (the two signs share a magnitude).
//! 2. `10·log10·P[y_j]` — the dB value of the chosen shape codevector's
//!    power ([`shape_log_db`], block 94). 128 values.
//!
//! ### Typo note (eq. G-14)
//!
//! The printed equation G-14 in the 11/94 PDF shows the second
//! correction term as `+ log10·P[y_j]` (coefficient 1), but every other
//! source in the same subclause — the defining `δ = 10·log10·P − 32`
//! (eq. G-8), the substituted eq. G-13 (`… + 10·log10·P[y_j]`), and the
//! prose list item (2) which spells out "`10·log10·P[y_j]`, the dB value
//! of the power of the best shape codevector" — carry the factor 10.
//! The missing `10` in the printed G-14 is a typesetting slip; this
//! module follows the prose-consistent factor 10. The per-test
//! [`tests::reformulation_matches_direct_log`] proves the factor-10 form
//! reproduces the direct `10·log10·P[e(n)]−32` value bit-for-bit (to
//! float tolerance), which the factor-1 form does not.
//!
//! ## Adder 96 / limiter 97 — `δ(n−1)` reconstruction
//!
//! Figure G.1's adder 96 sums the two table look-ups to form the
//! unclipped offset-removed log-gain of the *previous* vector; limiter
//! 97 then enforces the eq. G-9 floor `δ(n−1) ≥ −32` (which follows from
//! the `P[e(n)] ≥ 1` clip of eq. G-7). The output of limiter 97 is
//! mathematically the output of adder 42 in Figure 6/G.728 — the value
//! that feeds the hybrid-window / Levinson chain (blocks 43 / 44 / 45)
//! and the 1-sample delay 95 of the log-gain predictor. This module
//! provides exactly that reconstruction in [`offset_removed_log_gain`].
//!
//! The downstream limiter 98 (`δ̂(n)` clamped to `[−32, 28]`, eq. G-3),
//! adder 99 (`+32`), and the inverse-logarithm block 48 are unchanged
//! from the floating-point adapter (they live in [`crate::gain_adapter`]
//! as the `[0, 60]`-dB limiter and `10^(·/20)`); this module supplies
//! only the §G.2.1 reformulated *measurement* of the realised log-gain.
//!
//! Floating-point only: the dB tables use `f64::log10` so they cannot be
//! `const`-evaluated (mirroring the runtime accessors in
//! [`crate::tables`]); the Q-format fixed-point packing of these tables
//! (§G.5 / §G.6) is deferred behind the floating-point build.

use crate::consts::{DIMINV, GOFF, NCWD, NG};
use crate::tables::{GQ, Y_ENERGY};

/// Floor for the offset-removed log-gain `δ(n−1)`, in dB — the limiter
/// 97 bound of Figure G.1/G.728 (`δ(n−1) ≥ −32`, eq. G-9). It is the dB
/// image of the `P[e(n)] ≥ 1` power clip (eq. G-7): `10·log10(1) − 32 =
/// −32`. Equal to `−GOFF`.
pub const DELTA_FLOOR_DB: f64 = -GOFF;

/// Block 93 — gain-codebook log-gain table `20·log10|g_i|`, in dB,
/// indexed by the 3-bit gain codebook index `i` (0 ≤ i < NG = 8).
///
/// Per §G.2.1 there are only four distinct magnitudes `|GQ[i]|` (the
/// gain levels and their sign mirrors), so entries `i` and `i+4` are
/// equal. The table is derived from the transcribed [`GQ`] levels via
/// the eq. G-13 / G-14 correction term `20·log10|g_i|`; the per-test
/// [`tests::gain_log_db_is_sign_symmetric`] proves the magnitude
/// symmetry and [`tests::gain_log_db_adjacent_levels_differ_by_ratio`]
/// cross-checks the `7/4` inter-level ratio against the spec footnote.
#[must_use]
pub fn gain_log_db() -> [f64; NG] {
    let mut out = [0.0f64; NG];
    for (i, &g) in GQ.iter().enumerate() {
        // 20·log10|g_i| — the dB value of the chosen gain level.
        out[i] = 20.0 * g.abs().log10();
    }
    out
}

/// Block 94 — shape-codebook log-gain table `10·log10·P[y_j]`, in dB,
/// indexed by the 7-bit shape codebook index `j` (0 ≤ j < NCWD = 128).
///
/// `P[y_j]` is the *power* of shape codevector `j`: its energy divided
/// by the vector dimension (`P[y_j] = E_j / IDIM = Y_ENERGY[j]·DIMINV`,
/// per the §G.2.1 definition of `P[x]`). Derived from the precomputed
/// [`Y_ENERGY`] shape-energy table; the per-test
/// [`tests::shape_log_db_matches_power_of_codevector`] cross-checks each
/// entry against a direct `10·log10((1/IDIM)·Σ_k y_j(k)²)` computation.
///
/// No shape codevector is the zero vector (the codebook has unit-ish
/// energy rows), so `P[y_j] > 0` and the logarithm is always finite; the
/// per-test [`tests::shape_log_db_all_finite`] enforces this.
#[must_use]
pub fn shape_log_db() -> [f64; NCWD] {
    let mut out = [0.0f64; NCWD];
    for (j, &energy) in Y_ENERGY.iter().enumerate() {
        // P[y_j] = E_j / IDIM = E_j · DIMINV; 10·log10·P[y_j].
        out[j] = 10.0 * (energy * DIMINV).log10();
    }
    out
}

/// Adders 96 + limiter 97 of Figure G.1/G.728 — reconstruct the
/// offset-removed log-gain `δ(n−1)` of the *previous* vector from its
/// chosen gain index `i` and shape index `j`, via the §G.2.1
/// reformulation (eq. G-14, corrected).
///
/// ```text
/// δ(n−1) = δ̂(n−1) + 20·log10|g_i| + 10·log10·P[y_j]   (adder 96)
///        → max(δ(n−1), −32)                            (limiter 97)
/// ```
///
/// `predicted_log_gain` is the predicted (and possibly already range-
/// limited) `δ̂(n−1)` held by the 1-sample delay 95 — i.e. the output of
/// the log-gain linear predictor (block 46) for the previous vector,
/// *before* the `+GOFF` offset is added back (so it lives in the
/// offset-removed, `−32`-floored domain, not the `[0, 60]` σ-dB domain).
///
/// `gain_index` / `shape_index` are the previous vector's codebook
/// indices (the outputs of the 1-index delays 91 / 92).
///
/// Returns the limiter-97 output: the offset-removed log-gain that is
/// mathematically the adder-42 output of Figure 6/G.728. This value
/// feeds the hybrid-window / Levinson-Durbin chain and the log-gain
/// predictor memory.
#[must_use]
pub fn offset_removed_log_gain(
    predicted_log_gain: f64,
    gain_index: usize,
    shape_index: usize,
) -> f64 {
    debug_assert!(gain_index < NG, "gain index out of range");
    debug_assert!(shape_index < NCWD, "shape index out of range");

    // Adder 96: δ̂(n−1) + 20·log10|g_i| + 10·log10·P[y_j].
    let gain_db = 20.0 * GQ[gain_index].abs().log10();
    let shape_db = 10.0 * (Y_ENERGY[shape_index] * DIMINV).log10();
    let delta = predicted_log_gain + gain_db + shape_db;

    // Limiter 97: enforce eq. G-9 floor δ(n−1) ≥ −32.
    delta.max(DELTA_FLOOR_DB)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::IDIM;
    use crate::tables::{y_f64, GSQ};

    // --- Block 93: gain-codebook log-gain table -------------------

    #[test]
    fn gain_log_db_is_sign_symmetric() {
        // §G.2.1: only four distinct |g_i| exist; entries i and i+4
        // (a level and its sign mirror) share a magnitude, hence the
        // same 20·log10|g_i|.
        let t = gain_log_db();
        for i in 0..(NG / 2) {
            assert!(
                (t[i] - t[i + 4]).abs() < 1e-12,
                "gain dB asymmetric at {i}: {} vs {}",
                t[i],
                t[i + 4]
            );
        }
    }

    #[test]
    fn gain_log_db_matches_direct_20log10() {
        // Cross-check each entry against a direct 20·log10|GQ[i]|.
        let t = gain_log_db();
        for i in 0..NG {
            let want = 20.0 * GQ[i].abs().log10();
            assert!((t[i] - want).abs() < 1e-12, "entry {i}");
        }
    }

    #[test]
    fn gain_log_db_adjacent_levels_differ_by_ratio() {
        // The spec footnote gives GQ(i) = (7/4)·GQ(i−1); in dB the
        // adjacent levels therefore differ by exactly 20·log10(7/4).
        let t = gain_log_db();
        let step = 20.0 * (7.0f64 / 4.0).log10();
        for i in 1..(NG / 2) {
            assert!(
                (t[i] - t[i - 1] - step).abs() < 1e-9,
                "level {i} step {} != {step}",
                t[i] - t[i - 1]
            );
        }
    }

    #[test]
    fn gain_log_db_first_level_value() {
        // GQ(1) = 33/64 = 0.515625 ⇒ 20·log10(0.515625) ≈ −5.7533 dB.
        let t = gain_log_db();
        assert!((t[0] - (-5.753321)).abs() < 1e-5, "got {}", t[0]);
    }

    // --- Block 94: shape-codebook log-gain table ------------------

    #[test]
    fn shape_log_db_matches_power_of_codevector() {
        // Each entry must equal 10·log10((1/IDIM)·Σ_k y_j(k)²) computed
        // directly from the float-view shape codebook.
        let t = shape_log_db();
        let y = y_f64();
        for j in 0..NCWD {
            let mut energy = 0.0f64;
            for k in 0..IDIM {
                energy += y[j][k] * y[j][k];
            }
            let want = 10.0 * (energy / IDIM as f64).log10();
            assert!(
                (t[j] - want).abs() < 1e-9,
                "shape dB[{j}] = {} vs direct {want}",
                t[j]
            );
        }
    }

    #[test]
    fn shape_log_db_all_finite() {
        // No zero shape codevector ⇒ every power is positive ⇒ every dB
        // value is finite (the logarithm never sees a zero argument).
        for v in shape_log_db() {
            assert!(v.is_finite(), "non-finite shape dB {v}");
        }
    }

    // --- Reformulation equivalence (eq. G-14 vs eq. G-8) ----------

    #[test]
    fn reformulation_matches_direct_log() {
        // The whole point of §G.2.1: the reformulated δ(n) must equal
        // the direct floating-point δ(n) = 10·log10·P[e(n)] − 32 for the
        // same chosen σ, gain index and shape index. We exercise a
        // spread of σ, i, j and assert agreement to float tolerance.
        let gain_db = gain_log_db();
        let shape_db = shape_log_db();
        let y = y_f64();

        for &sigma in &[1.0f64, 3.5, 12.0, 200.0, 900.0] {
            // δ̂(n) such that σ = 10^((δ̂+32)/20) (eq. G-4 inverse).
            let delta_hat = 20.0 * sigma.log10() - 32.0;
            for ig in [0usize, 1, 3, 5, 7] {
                for js in [0usize, 1, 42, 99, 127] {
                    // Reformulated δ(n) (eq. G-14, corrected factor 10).
                    let reform = delta_hat + gain_db[ig] + shape_db[js];

                    // Direct δ(n): build e(n) = σ·g_i·y_j, take its
                    // power, clip to ≥ 1, then 10·log10·P − 32 (eq. G-8).
                    let mut p = 0.0f64;
                    for k in 0..IDIM {
                        let e = sigma * GQ[ig] * y[js][k];
                        p += e * e;
                    }
                    p /= IDIM as f64;
                    let p_clipped = p.max(1.0);
                    let direct = 10.0 * p_clipped.log10() - 32.0;

                    // When P[e(n)] ≥ 1 (the common case here) the two
                    // forms agree exactly bar float round-off.
                    if p >= 1.0 {
                        assert!(
                            (reform - direct).abs() < 1e-6,
                            "σ={sigma} i={ig} j={js}: reform {reform} vs direct {direct}"
                        );
                    } else {
                        // Below the clip the reformulation (which has no
                        // clip in this term) must sit at or below the
                        // clipped direct value; limiter 97 then floors it.
                        assert!(reform <= direct + 1e-6);
                    }
                }
            }
        }
    }

    #[test]
    fn factor_one_form_would_diverge() {
        // Guard against the eq. G-14 printed typo: using coefficient 1
        // on the shape term (instead of 10) breaks the equivalence, so
        // the factor-10 form is demonstrably the correct one.
        let y = y_f64();
        let sigma = 50.0f64;
        let delta_hat = 20.0 * sigma.log10() - 32.0;
        let ig = 2usize;
        let js = 17usize;

        let mut p = 0.0f64;
        for k in 0..IDIM {
            let e = sigma * GQ[ig] * y[js][k];
            p += e * e;
        }
        p /= IDIM as f64;
        let direct = 10.0 * p.max(1.0).log10() - 32.0;

        let gain_db = 20.0 * GQ[ig].abs().log10();
        let shape_pow = Y_ENERGY[js] * DIMINV;
        let reform10 = delta_hat + gain_db + 10.0 * shape_pow.log10();
        let reform1 = delta_hat + gain_db + shape_pow.log10();

        assert!((reform10 - direct).abs() < 1e-6, "factor-10 must match");
        assert!(
            (reform1 - direct).abs() > 1e-3,
            "factor-1 should diverge, but matched"
        );
    }

    // --- Adder 96 / limiter 97 ------------------------------------

    #[test]
    fn offset_removed_log_gain_applies_floor() {
        // A very negative predicted log-gain plus the quietest codebook
        // entries must be floored to −32 by limiter 97 (eq. G-9).
        let d = offset_removed_log_gain(-100.0, 0, 0);
        assert!((d - DELTA_FLOOR_DB).abs() < 1e-12, "got {d}");
    }

    #[test]
    fn offset_removed_log_gain_matches_table_sum_above_floor() {
        // Above the floor the result is exactly the table sum.
        let gain_db = gain_log_db();
        let shape_db = shape_log_db();
        let predicted = 5.0f64;
        let ig = 3usize;
        let js = 64usize;
        let want = predicted + gain_db[ig] + shape_db[js];
        // Confirm the test case is above the floor before asserting.
        assert!(want > DELTA_FLOOR_DB);
        let got = offset_removed_log_gain(predicted, ig, js);
        assert!((got - want).abs() < 1e-12, "got {got} want {want}");
    }

    #[test]
    fn delta_floor_is_negative_goff() {
        // Limiter 97's floor is −GOFF = −32 dB (the dB image of the
        // P[e(n)] ≥ 1 clip).
        assert!((DELTA_FLOOR_DB - (-GOFF)).abs() < 1e-12);
        assert_eq!(DELTA_FLOOR_DB, -32.0);
    }

    #[test]
    fn gsq_consistency_sanity() {
        // GSQ[i] = GQ[i]² (already a table); the gain dB term equals
        // 10·log10(GSQ[i]) as well as 20·log10|GQ[i]|. Cross-check the
        // two routes agree, anchoring the dB derivation to the existing
        // GSQ table.
        let t = gain_log_db();
        for i in 0..NG {
            let via_gsq = 10.0 * GSQ[i].log10();
            assert!((t[i] - via_gsq).abs() < 1e-9, "entry {i}");
        }
    }
}
