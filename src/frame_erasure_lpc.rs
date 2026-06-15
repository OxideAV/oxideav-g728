//! Annex I §I.4.2 — frame-erasure LPC filter "softening".
//!
//! During an erased frame the G.728 decoder has no channel index, so it
//! cannot run the normal backward synthesis-filter adapter (blocks
//! 49/50/51) to produce a fresh predictor. Annex I instead reuses the
//! **last good** 50th-order LPC predictor, but "softens" it by an
//! *extra* bandwidth expansion so the synthesis filter's spectral peaks
//! are flattened — this masks the otherwise audible artefacts of feeding
//! an extrapolated excitation (§I.4.1) through an unchanged sharp filter.
//!
//! The softening is the same operation as block 51's bandwidth expansion
//! (`5.6/G.728`), except the expansion factor [`FEFAC`] is `0.97` rather
//! than the normal `FAC = 253/256 ≈ 0.9883` (§I.4.2):
//!
//! > Let `a_i` be the i-th LPC coefficients used in the last adaptation
//! > cycle of the last good frame. Then, the new set of softened LPC
//! > coefficients is obtained as `a′_i = (0.97)^i · a_i`, i = 1, …, 50.
//!
//! and the expansion is *progressive*: the farther into the erasure, the
//! more it is softened. The first erased frame uses `(0.97)^i` (step
//! `k = 1`); if the erasure lasts beyond 10 ms the predictor is
//! re-softened to `(0.97)^{2i}` (step `k = 2`), and in general after
//! `k * 10 ms` of erasure the factor is `(0.97)^{k·i}` (§I.4.2):
//!
//! > the bandwidth expansion factor is `(0.97)^k`. The farther away from
//! > the last good frame, the more bandwidth expansion we apply.
//!
//! One adaptation cycle is `NFRSZ = 20` samples = 2.5 ms, so 10 ms is
//! [`FE_LPC_CYCLES_PER_STEP`] `= 4` adaptation cycles; the spec pins the
//! cadence as "updated again at the third vector in the 5th adaptation
//! cycle" — i.e. the second softening step lands once four adaptation
//! cycles (one of which is the partial onset cycle) have elapsed.
//!
//! This module realises the pure transform and the step bookkeeping. The
//! base predictor `a_i` it is handed is the last-good order-50 predictor
//! in the crate's spec `A` layout `[1.0, a_1, …, a_50]` (the same layout
//! [`crate::synthesis_adapter::SynthesisAdapter::coefficients`]
//! produces). Element 0 (the implicit `a_0 = 1`) is never scaled, exactly
//! as block 51 leaves `FACV(1) = 1.0` untouched.
//!
//! Clean-room: every value and rule here is transcribed from the §I.4.2
//! prose of `T-REC-G.728-199905-AnnI.pdf`; no reference C / external
//! implementation was consulted.

use crate::consts::{FEFAC, FE_LPC_CYCLES_PER_STEP, LPC};

/// Softened (frame-erasure bandwidth-expanded) order-50 LPC predictor.
///
/// `softened[0] == base[0]` (the implicit `a_0 = 1` is never scaled) and
/// `softened[i] == FEFAC.powi((step·i) as i32) · base[i]` for
/// `i = 1..=LPC`.
///
/// `step` is the §I.4.2 expansion step `k` (≥ 1): `k = 1` for the first
/// erased frame, incremented by one per additional 10 ms of erasure.
/// `step = 0` returns `base` unchanged (no softening), which is the
/// last-good / non-erased state.
#[must_use]
pub fn soften_predictor(base: &[f64; LPC + 1], step: usize) -> [f64; LPC + 1] {
    // Block-51-style bandwidth expansion: a′_i = factor^i · a_i, where
    // the per-erasure factor is FEFAC^step (§I.4.2). The implicit a_0 = 1
    // tap is left untouched (block 51 keeps FACV(1) = 1.0).
    let mut out = *base;
    if step == 0 {
        return out;
    }
    for (i, slot) in out.iter_mut().enumerate().take(LPC + 1).skip(1) {
        // factor^i = (FEFAC^step)^i = FEFAC^(step·i); compute directly
        // as one powi so the geometric progression is bit-identical to
        // the closed form the spec writes ((0.97)^{k·i}).
        let exponent = (step * i) as i32;
        *slot *= FEFAC.powi(exponent);
    }
    out
}

/// Tracks the §I.4.2 softening **step** `k` across a run of erased
/// adaptation cycles.
///
/// The decoder drives one [`Self::on_erased_cycle`] per erased
/// adaptation cycle (2.5 ms) and one [`Self::on_good_cycle`] per good
/// cycle. The current step is exposed via [`Self::step`] and fed to
/// [`soften_predictor`]; it is `0` outside an erasure, `1` for the first
/// 10 ms of erasure, `2` for the next 10 ms, etc. — re-incrementing
/// every [`FE_LPC_CYCLES_PER_STEP`] erased cycles, and resetting to `0`
/// as soon as a good frame returns ("next time there is a bad frame
/// again, the process starts from a bandwidth expansion factor of 0.97
/// again").
#[derive(Debug, Clone, Default)]
pub struct FrameErasureLpc {
    /// Number of consecutively erased adaptation cycles so far in the
    /// current erasure (0 when not in an erasure).
    erased_cycles: usize,
}

impl FrameErasureLpc {
    /// Fresh tracker, not in an erasure (`step == 0`).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Advance one **erased** adaptation cycle and return the new step.
    ///
    /// The first erased cycle moves the step to `1`; the step then
    /// increments by one once every `FE_LPC_CYCLES_PER_STEP` erased
    /// cycles thereafter (§I.4.2's "every additional 10 ms").
    pub fn on_erased_cycle(&mut self) -> usize {
        self.erased_cycles += 1;
        self.step()
    }

    /// A good (non-erased) adaptation cycle returns the tracker to the
    /// non-erased state so the next erasure restarts at step 1
    /// (`FEFAC^1`), per §I.4.2.
    pub fn on_good_cycle(&mut self) {
        self.erased_cycles = 0;
    }

    /// Current softening step `k` (§I.4.2): `0` outside an erasure,
    /// otherwise `ceil(erased_cycles / FE_LPC_CYCLES_PER_STEP)` so the
    /// first 10 ms (cycles 1..=4) is step 1, the next 10 ms is step 2,
    /// and so on.
    #[must_use]
    pub fn step(&self) -> usize {
        if self.erased_cycles == 0 {
            0
        } else {
            // ceil division: cycles 1..=4 → 1, 5..=8 → 2, …
            self.erased_cycles.div_ceil(FE_LPC_CYCLES_PER_STEP)
        }
    }

    /// Number of consecutively erased adaptation cycles in the current
    /// erasure (0 outside an erasure). Exposed for tests / audit.
    #[must_use]
    pub fn erased_cycles(&self) -> usize {
        self.erased_cycles
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_predictor() -> [f64; LPC + 1] {
        // A hand-built non-trivial predictor: A(z) layout with a_0 = 1
        // and decaying taps so the geometric softening is observable at
        // every order. The exact values don't matter — only that the
        // scaling is applied tap-for-tap.
        let mut a = [0.0f64; LPC + 1];
        a[0] = 1.0;
        for (i, slot) in a.iter_mut().enumerate().take(LPC + 1).skip(1) {
            *slot = 0.5 * (-0.9f64).powi(i as i32);
        }
        a
    }

    #[test]
    fn step_zero_is_identity() {
        // step = 0 is the non-erased state — no softening at all.
        let base = base_predictor();
        let out = soften_predictor(&base, 0);
        assert_eq!(out, base);
    }

    #[test]
    fn step_one_matches_spec_closed_form() {
        // §I.4.2: a′_i = (0.97)^i · a_i for the first erased frame.
        let base = base_predictor();
        let out = soften_predictor(&base, 1);
        assert_eq!(out[0], base[0], "a_0 = 1 is never scaled (FACV(1)=1)");
        for i in 1..=LPC {
            let expected = FEFAC.powi(i as i32) * base[i];
            assert!(
                (out[i] - expected).abs() < 1e-15,
                "tap {i}: {} vs {}",
                out[i],
                expected
            );
        }
    }

    #[test]
    fn step_k_is_factor_to_the_k() {
        // §I.4.2: after k·10 ms of erasure the factor is (0.97)^k, so
        // a′_i = (0.97)^{k·i}·a_i. Cross-check step 3 against the
        // "compose the step-1 transform three times" reading: applying
        // FEFAC^i three times to the *base* gives FEFAC^{3i}.
        let base = base_predictor();
        let out = soften_predictor(&base, 3);
        for i in 1..=LPC {
            // (FEFAC^3)^i, written the way the spec phrases it.
            let per_erasure_factor = FEFAC.powi(3);
            let expected = per_erasure_factor.powi(i as i32) * base[i];
            assert!(
                (out[i] - expected).abs() < 1e-14,
                "tap {i}: {} vs {}",
                out[i],
                expected
            );
        }
    }

    #[test]
    fn softening_strictly_flattens_high_order_taps() {
        // Sanity: bandwidth expansion shrinks every non-leading tap
        // toward zero, and shrinks higher-order taps proportionally more
        // (factor^i with factor < 1). |a′_i| < |a_i| and the ratio
        // a′_i/a_i = factor^i decreases with i.
        let base = base_predictor();
        let out = soften_predictor(&base, 1);
        let mut prev_ratio = 1.0;
        for i in 1..=LPC {
            assert!(out[i].abs() <= base[i].abs());
            let ratio = FEFAC.powi(i as i32);
            assert!(ratio < prev_ratio, "ratio must shrink with order");
            prev_ratio = ratio;
        }
    }

    #[test]
    fn tracker_starts_non_erased() {
        let t = FrameErasureLpc::new();
        assert_eq!(t.step(), 0);
        assert_eq!(t.erased_cycles(), 0);
    }

    #[test]
    fn tracker_first_erased_cycle_is_step_one() {
        // §I.4.2: the first erased frame softens by (0.97)^i — step 1.
        let mut t = FrameErasureLpc::new();
        assert_eq!(t.on_erased_cycle(), 1);
    }

    #[test]
    fn tracker_step_increments_every_four_cycles() {
        // 10 ms = FE_LPC_CYCLES_PER_STEP = 4 adaptation cycles. Cycles
        // 1..=4 stay at step 1 ("first 10 ms"); cycle 5 ("third vector
        // in the 5th adaptation cycle") moves to step 2; cycles 5..=8 at
        // step 2; cycle 9 → step 3.
        let mut t = FrameErasureLpc::new();
        let steps: Vec<usize> = (0..12).map(|_| t.on_erased_cycle()).collect();
        assert_eq!(steps, vec![1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]);
    }

    #[test]
    fn tracker_good_cycle_resets_to_step_zero() {
        // "next time there is a bad frame again, the process starts from
        // a bandwidth expansion factor of 0.97 again" — a good cycle
        // resets, so the next erasure restarts at step 1.
        let mut t = FrameErasureLpc::new();
        t.on_erased_cycle();
        t.on_erased_cycle();
        t.on_erased_cycle();
        t.on_erased_cycle();
        t.on_erased_cycle(); // now step 2
        assert_eq!(t.step(), 2);
        t.on_good_cycle();
        assert_eq!(t.step(), 0);
        assert_eq!(t.erased_cycles(), 0);
        assert_eq!(t.on_erased_cycle(), 1, "next erasure restarts at step 1");
    }

    #[test]
    fn fefac_matches_spec_paragraph_value() {
        // §I.4.2 spells out 0.97 rather than 253/256; guard the literal.
        assert!((FEFAC - 0.97).abs() < 1e-15);
        // and it is genuinely softer (smaller) than the normal block-51
        // factor, so the erased-frame filter is flatter than the good one.
        const _: () = assert!(FEFAC < crate::consts::FAC);
    }
}
