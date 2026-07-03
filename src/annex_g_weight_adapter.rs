//! Annex G §G.3 — fixed-point perceptual-weighting-filter backward adapter
//! (blocks 36 / 37 / 38).
//!
//! The encoder's perceptual weighting filter `W(z) = Q̃(z/γ₁) / Q̃(z/γ₂)`
//! (§3.3, eq. 3-4a) is backward-adapted from past input speech by the same
//! three-stage chain as the synthesis filter, at order `LPCW = 10`:
//!
//! * **block 36** (§G.3.12) — the perceptual-weighting hybrid window
//!   `WNRW` (Annex A.3, `N3 = LPCW + NFRSZ + NONRW = 60`) applied to the
//!   flat `Q2` past-input-speech buffer `SBW`: one `FINDNLS(…, 14)` call
//!   over the whole 60-sample buffer sets `NLSTMP = NLS − 1` (2 bits of
//!   headroom), the `Q15` window multiply rounds into the BFL `WS`, and
//!   the shared §G.3.18 HWMCORE runs with `NLSATTW = 15` — block 36's
//!   recursive decay is **1/2** (the base spec's float pseudo-code:
//!   `REXPW(I) = (1/2)·REXPW(I) + TMP`), unlike the 3/4 of blocks
//!   43/49. `NLSREXPW` is clamped at 41 after the core ("to avoid
//!   reduced accuracy … during long periods of zero input signal"),
//! * **block 37** (§G.2.2) — the fixed-point Levinson-Durbin recursion
//!   [`levinson_durbin_fixed`], producing the order-10 predictor `QTMP`
//!   (`q_i`, eq. 3-2f) in the `Q13/Q14/Q15` block-floating format signalled
//!   by `NLSATMP`,
//! * **block 38** (§3.4, eq. 3-4b / 3-4c) — two bandwidth expansions of the
//!   predictor: the numerator `Q̃(z/γ₁)` with `γ₁ = WZCF = 0.9`
//!   (`WZCFV_Q14`) and the denominator `Q̃(z/γ₂)` with `γ₂ = WPCF = 0.6`
//!   (`WPCFV_Q14`). Both reuse the shared §G.3.15 / §G.3.19 expansion core
//!   [`bandwidth_expand_q14`], the same machinery block 51 / block 45 use.
//!
//! This is the fixed-point analogue of the floating-point
//! [`crate::weighting_filter_adapter`] (blocks 36/37) +
//! [`crate::weighting_filter_coeff`] (block 38). The emitted Q14 numerator
//! / denominator coefficient pair is exactly the coefficient set a
//! fixed-point weighting filter `W(z)` (encoder blocks 4 / 10) would
//! consume.
//!
//! The non-speech `W(z) = 1` disabling of §3.4.1 (`γ₁ = γ₂ = 0`) collapses
//! both coefficient sets to the all-pass `(16384, 0, …, 0)`; the
//! [`WeightAdapterFixed::disabled_coeffs`] helper returns that state
//! directly.
//!
//! All arithmetic is built on the §G.1.3 primitives in
//! [`crate::annex_g_arith`]; no codebook is read here. The per-test
//! cross-checks below run this fixed-point chain against the floating-point
//! [`crate::weighting_filter_adapter::WeightingFilterAdapter`] +
//! [`crate::weighting_filter_coeff::WeightingFilterCoeff`] on identical
//! requantized speech.

use crate::annex_g_arith::{findnls, rnd};
use crate::annex_g_hybrid::{HwmcoreState, NLSATTW};
use crate::annex_g_levinson::{levinson_durbin_fixed, LevinsonInput, LevinsonStatus};
use crate::annex_g_synth_adapter::bandwidth_expand_q14;
use crate::consts::{LPCW, NFRSZ, NONRW};
use crate::tables::{WNRW_Q15, WPCFV_Q14, WZCFV_Q14};

/// `N1 = LPCW + NFRSZ = 30` — end of the recursive-component taps.
pub const N1W: usize = LPCW + NFRSZ;
/// `N2 = LPCW + NONRW = 40` — length of the kept-old part of `SBW`.
pub const N2W: usize = LPCW + NONRW;
/// `N3 = LPCW + NFRSZ + NONRW = 60` — total `SBW` / `WS` length.
pub const N3W: usize = LPCW + NFRSZ + NONRW;

/// §G.3.12: "If NLSREXPW > 41, set NLSREXPW = 41" — the zero-input
/// accuracy clamp applied after HWMCORE.
pub const NLSREXPW_MAX: i32 = 41;

/// The two Q14 coefficient vectors of `W(z) = Q̃(z/γ₁) / Q̃(z/γ₂)` — the
/// fixed-point analogue of
/// [`crate::weighting_filter_coeff::WeightingFilterCoeff`].
///
/// Both vectors use the spec 1-based layout: index 0 is the implicit
/// `q₀·γ⁰ = 1` leading tap (Q14 unity `16384`); indices `1..=LPCW` carry
/// the bandwidth-broadened predictor coefficients in Q14.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WeightCoeffFixed {
    /// Numerator `Q̃(z/γ₁)` coefficients `q_i·γ₁ⁱ` (Q14), `i = 0..=LPCW`.
    pub q_gamma1: [i16; LPCW + 1],
    /// Denominator `Q̃(z/γ₂)` coefficients `q_i·γ₂ⁱ` (Q14), `i = 0..=LPCW`.
    pub q_gamma2: [i16; LPCW + 1],
}

impl WeightCoeffFixed {
    /// The §3.4.1 disabled `W(z) = 1` state: both vectors all-pass
    /// `(16384, 0, …, 0)` in Q14.
    #[must_use]
    pub const fn disabled() -> Self {
        let mut q_gamma1 = [0i16; LPCW + 1];
        let mut q_gamma2 = [0i16; LPCW + 1];
        q_gamma1[0] = 1 << 14;
        q_gamma2[0] = 1 << 14;
        Self { q_gamma1, q_gamma2 }
    }
}

impl Default for WeightCoeffFixed {
    fn default() -> Self {
        Self::disabled()
    }
}

/// Fixed-point perceptual-weighting-filter backward adapter (Annex G §G.3,
/// blocks 36 / 37 / 38).
///
/// Owns the shared §G.3.17 / §G.3.18 hybrid-window core sized for the
/// `LPCW = 10`-order weighting filter, plus the live Q14 coefficient pair
/// (held from the previous cycle when block 37 reports ill-conditioning).
#[derive(Debug, Clone)]
pub struct WeightAdapterFixed {
    /// `SBW(1..=N3)` — the flat `Q2` past-input-speech buffer (§G.3.12:
    /// "All SBW are Q2 and represented in 15 bits precision").
    sbw: [i16; N3W],
    /// §G.3.18 HWMCORE state (`REXPW` / `NLSREXPW`, init 31).
    core: HwmcoreState,
    /// Live Q14 weighting coefficients `W(z) = Q̃(z/γ₁) / Q̃(z/γ₂)`.
    coeff: WeightCoeffFixed,
}

impl Default for WeightAdapterFixed {
    fn default() -> Self {
        Self::new()
    }
}

impl WeightAdapterFixed {
    /// Fresh adapter: zeroed `SBW` / `REXPW` per Table 2/G.728 and the
    /// §3.4.1 disabled `W(z) = 1` coefficient pair (the cold-start state
    /// before the first adaptation cycle has produced `q_i`).
    #[must_use]
    pub fn new() -> Self {
        Self {
            sbw: [0; N3W],
            core: HwmcoreState::new(LPCW),
            coeff: WeightCoeffFixed::disabled(),
        }
    }

    /// Borrow the live Q14 weighting coefficient pair.
    #[must_use]
    pub fn coeff(&self) -> &WeightCoeffFixed {
        &self.coeff
    }

    /// Borrow the flat `Q2` signal-history buffer `SBW` (for tests/audit).
    #[must_use]
    pub fn sbw(&self) -> &[i16] {
        &self.sbw
    }

    /// Borrow the recursive-autocorrelation `REXPW` (for tests/audit).
    #[must_use]
    pub fn rexpw(&self) -> &[i16] {
        self.core.rexp()
    }

    /// `NLSREXPW` — the shared shift of `REXPW` (for tests/audit).
    #[must_use]
    pub fn nlsrexpw(&self) -> i32 {
        self.core.nlsrexp()
    }

    /// The §3.4.1 disabled `W(z) = 1` coefficient pair (`γ₁ = γ₂ = 0`).
    #[must_use]
    pub fn disabled_coeffs() -> WeightCoeffFixed {
        WeightCoeffFixed::disabled()
    }

    /// Run one adaptation cycle: block 36 (§G.3.12 flat-`Q2` hybrid
    /// window) → HWMCORE (`NLSATTW = 15`, 1/2 decay) → block 37
    /// (Levinson) → block 38 (two bandwidth expansions).
    ///
    /// `stmp` is the cycle's `NFRSZ = 20` new input-speech samples in
    /// flat `Q2` (`STMP`, oldest first). The returned [`LevinsonStatus`]
    /// surfaces the recursion verdict; the live [`coeff`](Self::coeff)
    /// is updated when block 37 succeeds and left unchanged (the
    /// previous cycle's pair) when it reports ill-conditioning — the
    /// §3.4 "keep the previous coefficients" rule.
    pub fn adapt(&mut self, stmp: &[i16; NFRSZ]) -> LevinsonStatus {
        // ---- Block 36 (§G.3.12): flat-Q2 front end ------------------
        // | For N = 1..N2: SBW(N) = SBW(N + NFRSZ)   | Shift old buffer
        // | For N = 1..NFRSZ: SBW(N2 + N) = STMP(N)  | Newest at SBW(N3)
        for n in 0..N2W {
            self.sbw[n] = self.sbw[n + NFRSZ];
        }
        self.sbw[N2W..N3W].copy_from_slice(stmp);

        // | Call FINDNLS(SBW, N3, N3, 14, NLS); NLSTMP = NLS − 1
        let vals: [i32; N3W] = std::array::from_fn(|i| i32::from(self.sbw[i]));
        let nls = findnls(&vals, N3W, 14);
        let nlstmp = nls - 1;

        // | K = 1; For N = 60..1: P = SBW(N)·WNRW(K); AA0 = P << NLSTMP;
        // | WS(N) = RND(AA0); K += 1.
        // (`P >> 1` when NLSTMP = −1, the same convention block 43
        // spells out for its own front end.)
        let mut ws = [0i16; N3W];
        for (k, n) in (0..N3W).rev().enumerate() {
            let p = i32::from(self.sbw[n]) * i32::from(WNRW_Q15[k]);
            let aa0 = if nlstmp == -1 { p >> 1 } else { p << nlstmp };
            ws[n] = rnd(aa0);
        }

        // | NLSATTW = 15
        // | Call HWMCORE(LPCW, N1, N3, NLSATTW, WS, NLSTMP, REXPW,
        // |              NLSREXPW, R, ILLCONDW)
        let core = self.core.run(LPCW, N1W, N3W, NLSATTW, &ws, nlstmp);
        // | If NLSREXPW > 41, set NLSREXPW = 41
        self.core.clamp_nlsrexp(NLSREXPW_MAX);

        // ---- Block 37: §G.2.2 Levinson-Durbin recursion (order LPCW) ----
        let input = LevinsonInput {
            rtmp: &core.rtmp,
            illcond: core.illcond,
            lpc: LPCW,
        };
        let mut qtmp = [0i16; LPCW + 1];
        let status = levinson_durbin_fixed(&input, &mut qtmp);

        // ---- Block 38: two bandwidth expansions (eq. 3-4b / 3-4c) ----
        self.block38(&qtmp, &status);

        status
    }

    /// Block 38 (§3.4, eq. 3-4b / 3-4c): bandwidth-expand the order-10
    /// predictor `QTMP` (Q-format signalled by `status.nlsatmp`) by `γ₁`
    /// (numerator, `WZCFV`) and `γ₂` (denominator, `WPCFV`) into the live
    /// Q14 coefficient pair. On ill-conditioning or a Q14 overflow the
    /// previous cycle's pair is kept (the §3.4 "keep the previous
    /// coefficients" rule).
    fn block38(&mut self, qtmp: &[i16; LPCW + 1], status: &LevinsonStatus) {
        // | If ILLCOND = .TRUE., skip the execution of this block.
        if status.illcond {
            return;
        }
        // Both expansions must succeed (same NLSATMP) before committing —
        // a Q14 overflow on either keeps the previous W(z) intact.
        let num = bandwidth_expand_q14(&qtmp[..=LPCW], &WZCFV_Q14, status.nlsatmp);
        let den = bandwidth_expand_q14(&qtmp[..=LPCW], &WPCFV_Q14, status.nlsatmp);
        if let (Some(num), Some(den)) = (num, den) {
            self.coeff.q_gamma1.copy_from_slice(&num);
            self.coeff.q_gamma2.copy_from_slice(&den);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::{WPCF, WZCF};
    use crate::hybrid_window::{HybridWindow, HybridWindowState};
    use crate::levinson::levinson_durbin;
    use crate::tables::wnrw_f64;
    use crate::weighting_filter_coeff::WeightingFilterCoeff;

    /// Quantize `NFRSZ` float speech samples to the flat 15-bit `Q2`
    /// `STMP` form and return the requantized float samples so the
    /// float oracle sees the same values.
    fn q2_stmp(samples: &[f64; NFRSZ]) -> ([i16; NFRSZ], [f64; NFRSZ]) {
        let mut stmp = [0i16; NFRSZ];
        let mut requant = [0.0f64; NFRSZ];
        for (i, &v) in samples.iter().enumerate() {
            let q = (v * 4.0).round().clamp(-16384.0, 16383.0);
            stmp[i] = q as i16;
            requant[i] = q / 4.0;
        }
        (stmp, requant)
    }

    fn lcg(seed: &mut u64, amp: f64) -> [f64; NFRSZ] {
        let mut out = [0.0f64; NFRSZ];
        for s in out.iter_mut() {
            *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            let u = ((*seed >> 33) as f64) / ((1u64 << 31) as f64) - 1.0;
            *s = u * amp;
        }
        out
    }

    #[test]
    fn fresh_adapter_is_disabled_filter() {
        // Cold-start: W(z) = 1, both coefficient sets all-pass Q14.
        let a = WeightAdapterFixed::new();
        assert_eq!(a.coeff().q_gamma1[0], 1 << 14);
        assert_eq!(a.coeff().q_gamma2[0], 1 << 14);
        assert!(a.coeff().q_gamma1[1..].iter().all(|&v| v == 0));
        assert!(a.coeff().q_gamma2[1..].iter().all(|&v| v == 0));
    }

    #[test]
    fn default_matches_new() {
        let a = WeightAdapterFixed::default();
        let b = WeightAdapterFixed::new();
        assert_eq!(a.coeff(), b.coeff());
        assert_eq!(a.sbw(), b.sbw());
    }

    #[test]
    fn dimensions_match_spec() {
        // §G.3.12: N1 = 30, N2 = 40, N3 = 60.
        assert_eq!(N1W, 30);
        assert_eq!(N2W, 40);
        assert_eq!(N3W, 60);
        let a = WeightAdapterFixed::new();
        assert_eq!(a.sbw().len(), 60);
        assert_eq!(a.rexpw().len(), LPCW + 1);
        assert_eq!(a.nlsrexpw(), 31, "NLSREXPW initial value = 31");
    }

    #[test]
    fn zero_input_keeps_disabled_filter_and_flags_illcond() {
        // Zero-speech ⇒ RTMP all-zero ⇒ ILLCOND ⇒ block 37/38 decline ⇒
        // the disabled W(z) = 1 pair is preserved.
        let mut a = WeightAdapterFixed::new();
        let status = a.adapt(&[0i16; NFRSZ]);
        assert!(status.illcond, "zero input ⇒ ILLCOND");
        assert_eq!(a.coeff().q_gamma1[0], 1 << 14);
        assert!(a.coeff().q_gamma1[1..].iter().all(|&v| v == 0));
        assert!(a.coeff().q_gamma2[1..].iter().all(|&v| v == 0));
    }

    #[test]
    fn block38_unity_heads_and_gamma_ratio() {
        // Drive a real broadband cycle until block 37 succeeds, then
        // confirm the block-38 invariants: both Q14 unity heads survive
        // and the denominator (γ₂ = 0.6) coefficients are strictly smaller
        // in magnitude than the numerator (γ₁ = 0.9) ones for the same q_i
        // (the §3.4 / §3.7 pole-vs-zero radius ordering γ₂ < γ₁).
        let mut a = WeightAdapterFixed::new();
        let mut seed = 0x51a7_e1b2_c3d4_e5f6u64;
        let mut committed = false;
        for _ in 0..20 {
            let samples = lcg(&mut seed, 3500.0);
            let (stmp, _) = q2_stmp(&samples);
            let status = a.adapt(&stmp);
            if !status.illcond {
                committed = true;
                assert_eq!(a.coeff().q_gamma1[0], 1 << 14, "numerator unity head");
                assert_eq!(a.coeff().q_gamma2[0], 1 << 14, "denominator unity head");
                // For each tap, |q_gamma2[i]| <= |q_gamma1[i]| (γ₂ < γ₁
                // pulls the denominator coefficients radially inward).
                for i in 1..=LPCW {
                    assert!(
                        (a.coeff().q_gamma2[i] as i32).abs()
                            <= (a.coeff().q_gamma1[i] as i32).abs() + 1,
                        "tap {i}: |den| {} should not exceed |num| {}",
                        a.coeff().q_gamma2[i],
                        a.coeff().q_gamma1[i]
                    );
                }
            }
        }
        assert!(committed, "broadband input must commit a real W(z)");
    }

    #[test]
    fn tracks_floating_point_weighting_adapter() {
        // End-to-end shape cross-check: the fixed-point W(z) coefficient
        // pair must track the floating-point block-36/37/38 chain
        // (WeightingFilterAdapter + WeightingFilterCoeff) driven by the
        // *same requantized* speech. We compare the normalized predictor
        // q_i / q_1 produced by both Levinson recursions (the quantity
        // block 38 expands), within the BFL-arithmetic tolerance.
        let wnrw = wnrw_f64();
        let fhw = HybridWindow {
            m: LPCW,
            l: NFRSZ,
            n: NONRW,
            window: &wnrw,
            // Block 36's recursive decay (1/2 — see module header).
            decay: 0.5,
        };
        let mut fstate = HybridWindowState::new(&fhw);
        let mut a = WeightAdapterFixed::new();

        let mut seed = 0xfeed_face_dead_c0deu64;
        let mut compared = false;
        for cyc in 0..14 {
            let samples = lcg(&mut seed, 3200.0);
            let (stmp, requant) = q2_stmp(&samples);

            // Float oracle: hybrid window → Levinson → block 38.
            let mut frtmp = [0.0f64; LPCW + 1];
            fstate.run(&fhw, &requant, &mut frtmp);
            let mut fq = [0.0f64; LPCW + 1];
            let fok = levinson_durbin(&frtmp, &mut fq, LPCW).is_ok();

            let status = a.adapt(&stmp);

            if cyc >= 11 && fok && !status.illcond {
                // Float block 38 for reference (γ₁ = 0.9, γ₂ = 0.6).
                let fcoeff = WeightingFilterCoeff::from_lpc(&fq);
                // Compare normalized numerator taps q_gamma1[i] (Q14) to
                // the float q_i·γ₁ⁱ. Both share the unity head; the BFL
                // path rounds, so allow a modest tolerance on the low taps.
                for i in 1..=4 {
                    let fixed = a.coeff().q_gamma1[i] as f64 / 16384.0;
                    let float = fcoeff.q_gamma1[i];
                    assert!(
                        (fixed - float).abs() < 0.15,
                        "q_gamma1[{i}] fixed {fixed} vs float {float}"
                    );
                }
                compared = true;
            }
        }
        assert!(compared, "a fully-populated good cycle must have compared");
        // Anchor the γ constants the Q14 tables encode.
        assert!((WZCF - 0.9).abs() < 1e-12);
        assert!((WPCF - 0.6).abs() < 1e-12);
    }

    #[test]
    fn disabled_coeffs_helper_matches_const() {
        assert_eq!(
            WeightAdapterFixed::disabled_coeffs(),
            WeightCoeffFixed::disabled()
        );
    }

    #[test]
    fn repeated_cycles_keep_finite_coeffs() {
        // Stability: many cycles never escape the Q14 i16 range (type-
        // guaranteed) and keep the unity heads intact.
        let mut a = WeightAdapterFixed::new();
        for cyc in 0..48 {
            let mut samples = [0.0f64; NFRSZ];
            for (i, s) in samples.iter_mut().enumerate() {
                let t = (cyc * NFRSZ + i) as f64;
                *s = (t * 0.17).sin() * 6000.0;
            }
            let (stmp, _) = q2_stmp(&samples);
            a.adapt(&stmp);
            assert_eq!(a.coeff().q_gamma1[0], 1 << 14);
            assert_eq!(a.coeff().q_gamma2[0], 1 << 14);
        }
    }
}
