//! Annex G §G.3 — fixed-point short-term adaptive-postfilter coefficient
//! calculator (block 85).
//!
//! The §4.6 adaptive postfilter's short-term section is a 10th-order
//! pole-zero filter cascaded with a first-order spectral-tilt compensator
//! (eq. 4-2):
//!
//! ```text
//! H_s(z) = (1 − Σ b̄_i z⁻ⁱ) / (1 − Σ ā_i z⁻ⁱ) · (1 + µ·z⁻¹)
//! ```
//!
//! with (eq. 4-3 / 4-4 / 4-5)
//!
//! * `b̄_i = ã_i · SPFZCF^i`  (zero-control, `SPFZCF = 0.65`),
//! * `ā_i = ã_i · SPFPCF^i`  (pole-control, `SPFPCF = 0.75`),
//! * `µ   = TILTF · k1`      (tilt, `TILTF = 0.15`),
//!
//! where `ã_i = −a_i` is the order-10 synthesis-filter LPC by-product (the
//! coefficients the §G.2.2 Levinson recursion produces at the `MINC0 = 10`
//! decoder restart) and `k1` is its first reflection coefficient.
//!
//! Block 85 is exactly two bandwidth expansions of the same order-10
//! predictor — the postfilter analogue of the weighting-filter block 38
//! ([`crate::annex_g_weight_adapter`]) and of blocks 45 / 51 — plus the
//! one-multiply tilt term. This module computes the fixed-point
//! [`crate::annex_g_postfilter::ShortTermCoeff`] (`AZ` / `AP` / `TILTZ`,
//! all Q14) that the §G.3.20–§G.3.23 postfilter
//! ([`crate::annex_g_postfilter::PostfilterFixed`]) consumes.
//!
//! ## Q-format provenance
//!
//! `AZ` / `AP` reuse the *already-specified* Q14 postfilter coefficient
//! format (Table G.2/G.728) and the shared §G.3.15 / §G.3.19 bandwidth-
//! expansion core [`bandwidth_expand_q14`] with the staged Q14 broadening
//! tables `SPFZCFV_Q14` / `SPFPCFV_Q14`. The first reflection coefficient
//! `k1` is the §G.2.2 recursion's `RC1`, a Q15 value with magnitude `< 1`;
//! the tilt term `TILTZ = TILTF·k1` is a Q14 quantity formed from the
//! Q15·Q15 = Q30 product `TILTF_Q15 · k1` whose rounded high word is Q14.
//! The multiplier constant is `TILTF_Q15 = 4915` (`§G.3.28`'s inline
//! "`TILTF = 4915 in Q15`", `0.15·2¹⁵ = 4915.2` truncated), **not** the
//! `2458` a Q14 reading of `TILTF` would give — see the `N2` note in
//! `docs/audio/g728/g728-errata.md`; the Annex-G conformance vectors are
//! bit-exact only with `4915`.
//! No new dB-domain Q-format is introduced — the log-gain adapter's
//! §G.3.12–§G.3.16 dB scaling remains a documented gap, unrelated to this
//! block.
//!
//! The floating-point reference for the same calculator lives in
//! [`crate::short_term_postfilter::ShortTermPostfilter::set_from_synthesis_byproduct`];
//! the per-test cross-checks below run this fixed-point path against it on
//! identical order-10 predictors.

use crate::annex_g_postfilter::{ShortTermCoeff, PF_ORDER};
use crate::annex_g_synth_adapter::bandwidth_expand_q14;
use crate::tables::{SPFPCFV_Q14, SPFZCFV_Q14};

/// `TILTF = 0.15` in Q15 — the §G.3.28 fixed-point pseudo-code's inline
/// constant ("TILTF = 4915 in Q15"; `0.15 · 2¹⁵ = 4915.2`). The spectral-
/// tilt controlling factor of Table 1/G.728 (`TILTF`), the multiplier for
/// the `TILTZ = RND(TILTF·RC1)` term of eq. 4-5 (`Q15 · Q15 = Q30`, whose
/// `RND` high word is the Q14 `TILTZ`).
pub const TILTF_Q15: i32 = 4915;

/// Block 85 — derive the fixed-point short-term postfilter coefficients
/// [`ShortTermCoeff`] from the order-10 synthesis-filter LPC by-product.
///
/// * `atil` is the order-10 predictor `ATMP(1..=PF_ORDER+1)` in the
///   block-floating `Q13/Q14/Q15` format signalled by `nlsatmp` (the
///   §G.2.2 Levinson output, `atil[0]` the implicit unity tap, `atil[i]`
///   the `ã_i = −a_i` taps). This is the order-10 by-product the decoder
///   captures at the `MINC0 = 10` recursion restart.
/// * `k1_q15` is the first reflection coefficient `RC1` (Q15, `|k1| < 1`).
/// * `nlsatmp` is the predictor's `Q` signal (`13` / `14` / `15`).
///
/// Returns `Some(ShortTermCoeff)` on success, or `None` when either
/// bandwidth expansion overflows Q14 / `nlsatmp` is out of range — in which
/// case the caller keeps the previous cycle's postfilter coefficients (the
/// §4.6 "frozen until the next valid update" behaviour).
#[must_use]
pub fn short_term_coeff_fixed(
    atil: &[i16; PF_ORDER + 1],
    k1_q15: i16,
    nlsatmp: i32,
) -> Option<ShortTermCoeff> {
    // Eq. 4-3: AZ(i+1) = ã_i · SPFZCF^i  (numerator, zero-control).
    let az = bandwidth_expand_q14(&atil[..=PF_ORDER], &SPFZCFV_Q14, nlsatmp)?;
    // Eq. 4-4: AP(i+1) = ã_i · SPFPCF^i  (denominator, pole-control).
    let ap = bandwidth_expand_q14(&atil[..=PF_ORDER], &SPFPCFV_Q14, nlsatmp)?;

    let mut az_q14 = [0i32; PF_ORDER + 1];
    let mut ap_q14 = [0i32; PF_ORDER + 1];
    for i in 0..=PF_ORDER {
        az_q14[i] = az[i] as i32;
        ap_q14[i] = ap[i] as i32;
    }

    // Eq. 4-5 / §G.3.28: AA0 = TILTF·RC1 (Q15·Q15 = Q30);
    // TILTZ = RND(AA0) — the rounded high word is the Q14 TILTZ.
    let prod = TILTF_Q15 * i32::from(k1_q15);
    let tiltz_q14 = i32::from(crate::annex_g_arith::rnd(prod));

    Some(ShortTermCoeff {
        az_q14,
        ap_q14,
        tiltz_q14,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::{SPFPCF, SPFZCF, TILTF};
    use crate::short_term_postfilter::ShortTermPostfilter;

    /// Q14 unity head builder for a small known predictor in Q14
    /// (`nlsatmp = 14`).
    fn q14_predictor(taps: &[f64]) -> [i16; PF_ORDER + 1] {
        let mut a = [0i16; PF_ORDER + 1];
        a[0] = 1 << 14;
        for (i, &t) in taps.iter().enumerate() {
            a[i + 1] = (t * 16384.0).round() as i16;
        }
        a
    }

    #[test]
    fn tiltf_q15_matches_constant() {
        // TILTF_Q15 must encode TILTF = 0.15 in Q15 (§G.3.28's inline
        // "TILTF = 4915 in Q15"; 0.15·2¹⁵ = 4915.2 truncates — NOT
        // rounds — to the annex constant 4915, and the Annex-G
        // conformance vectors are bit-exact only with 4915).
        assert_eq!(TILTF_Q15, (TILTF * 32768.0).floor() as i32);
        assert_eq!(TILTF_Q15, 4915);
    }

    #[test]
    fn unity_heads_and_orders() {
        // A modest predictor expands to AZ/AP with Q14 unity heads and the
        // full PF_ORDER+1 length.
        let atil = q14_predictor(&[-0.5, 0.25, -0.1, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let sc = short_term_coeff_fixed(&atil, 0, 14).expect("well-formed expands");
        assert_eq!(sc.az_q14.len(), PF_ORDER + 1);
        assert_eq!(sc.ap_q14.len(), PF_ORDER + 1);
        assert_eq!(sc.az_q14[0], 1 << 14, "AZ(1) Q14 unity");
        assert_eq!(sc.ap_q14[0], 1 << 14, "AP(1) Q14 unity");
    }

    #[test]
    fn az_decays_faster_than_ap() {
        // SPFZCF = 0.65 < SPFPCF = 0.75 ⇒ for the same ã_i the all-zero
        // (numerator) coefficients are pulled radially inward faster than
        // the all-pole (denominator) ones, so |AZ(i+1)| <= |AP(i+1)|.
        const _: () = assert!(SPFZCF < SPFPCF);
        let atil = q14_predictor(&[0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]);
        let sc = short_term_coeff_fixed(&atil, 0, 14).unwrap();
        for i in 1..=PF_ORDER {
            assert!(
                sc.az_q14[i].abs() <= sc.ap_q14[i].abs() + 1,
                "tap {i}: |AZ| {} should not exceed |AP| {}",
                sc.az_q14[i],
                sc.ap_q14[i]
            );
        }
    }

    #[test]
    fn tilt_term_matches_float() {
        // TILTZ = TILTF·k1. Sweep k1 across [−1, 1) Q15 and confirm the
        // fixed-point Q14 result tracks the float 0.15·k1 within a half-LSB.
        let atil = q14_predictor(&[0.0; PF_ORDER]);
        for &k1f in &[-0.9f64, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9] {
            let k1_q15 = (k1f * 32768.0).round().clamp(-32768.0, 32767.0) as i16;
            let sc = short_term_coeff_fixed(&atil, k1_q15, 14).unwrap();
            let float = TILTF * (k1_q15 as f64 / 32768.0);
            let fixed = sc.tiltz_q14 as f64 / 16384.0;
            assert!(
                (fixed - float).abs() < 1.0 / 16384.0 + 1e-9,
                "k1={k1f}: TILTZ fixed {fixed} vs float {float}"
            );
        }
    }

    #[test]
    fn tracks_floating_point_block85() {
        // The fixed-point AZ/AP/TILTZ must track the floating-point block-85
        // calculator (set_from_synthesis_byproduct) for the same order-10
        // predictor. The float path takes a10 (the synthesis predictor
        // a_i, NOT ã_i) and flips the sign internally (ã_i = −a_i); the
        // fixed path takes ATMP which already stores ã_i = −a_i, so we feed
        // the float path the negation of our ã_i taps.
        let atil_taps = [-0.6f64, 0.42, -0.3, 0.18, -0.1, 0.05, -0.02, 0.0, 0.0, 0.0];
        let atil = q14_predictor(&atil_taps);
        let k1f = -0.45f64;
        let k1_q15 = (k1f * 32768.0).round() as i16;

        let sc = short_term_coeff_fixed(&atil, k1_q15, 14).unwrap();

        // Float reference: a10[i] = a_i = −ã_i (so a10 = −atil taps).
        let mut a10 = [0.0f64; PF_ORDER + 1];
        a10[0] = 1.0;
        for i in 1..=PF_ORDER {
            a10[i] = -atil_taps[i - 1];
        }
        let mut pf = ShortTermPostfilter::new();
        pf.set_from_synthesis_byproduct(&a10, k1f);

        // Compare AZ(i+1) = b̄_i and AP(i+1) = ā_i for i = 1..=PF_ORDER.
        for i in 1..=PF_ORDER {
            let az_fixed = sc.az_q14[i] as f64 / 16384.0;
            let ap_fixed = sc.ap_q14[i] as f64 / 16384.0;
            let az_float = pf.numerator()[i];
            let ap_float = pf.denominator()[i];
            assert!(
                (az_fixed - az_float).abs() < 2.0 / 16384.0 + 1e-9,
                "AZ({}) fixed {az_fixed} vs float {az_float}",
                i + 1
            );
            assert!(
                (ap_fixed - ap_float).abs() < 2.0 / 16384.0 + 1e-9,
                "AP({}) fixed {ap_fixed} vs float {ap_float}",
                i + 1
            );
        }
        // Tilt term.
        let tilt_fixed = sc.tiltz_q14 as f64 / 16384.0;
        assert!((tilt_fixed - pf.tilt()).abs() < 2.0 / 16384.0 + 1e-9);
    }

    #[test]
    fn declines_on_bad_nlsatmp() {
        let atil = q14_predictor(&[0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(short_term_coeff_fixed(&atil, 0, 12).is_none());
        assert!(short_term_coeff_fixed(&atil, 0, 16).is_none());
    }
}
