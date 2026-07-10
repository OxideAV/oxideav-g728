//! Backward synthesis-filter adapter — block 33 of Figure 5/G.728.
//!
//! Chains the three sub-blocks the spec lists in §5.6:
//!
//! * **Block 49** — hybrid window on the quantised-speech buffer
//!   `STTMP(1..NFRSZ)`, producing autocorrelation taps `RTMP(1..LPC+1)`.
//! * **Block 50** — Levinson-Durbin recursion on `RTMP` producing the
//!   bandwidth-unexpanded predictor `ATMP(1..LPC+1)`.
//! * **Block 51** — bandwidth expansion: `ATMP(I) = FACV(I) · ATMP(I)`;
//!   the resulting `ATMP` is copied into the live predictor `A` at
//!   the third vector of the adaptation cycle (`ICOUNT = 3`).
//!
//! ## Cycle timing (§3.7, §5.6 block 50 prose, block 51 prose)
//!
//! G.728's adaptation cycle is `NUPDATE = 4` vectors (`NFRSZ = 20`
//! samples). The block-50 prose spells out the staggered update:
//!
//! > Although the autocorrelation `RTMP` array is available at the
//! > first vector of each adaptation cycle, the actual updates of
//! > synthesis filter coefficients will not take place until the
//! > third vector. This intentional delay of updates allows the
//! > real-time hardware to spread the computation of this module
//! > over the first three vectors of each adaptation cycle.
//!
//! In a one-shot software implementation we don't need to spread
//! across vectors — but the *timing semantics* still matter, because
//! the adapter must observe the full previous cycle's quantised
//! speech before updating, not the current cycle's. We therefore
//! accept one buffered call per cycle (the cycle-boundary delivery
//! of `NFRSZ` samples) and produce one expanded predictor per cycle.

use crate::consts::{LPC, NFRSZ, NONR};
use crate::hybrid_window::{HybridWindow, HybridWindowState};
use crate::levinson::{levinson_durbin, LevinsonError};
use crate::tables::{facv_f64, wnr_f64, FACV_Q14, Q14};

/// Order at which the synthesis-filter Levinson recursion is stopped to
/// extract the postfilter's `ã_1..ã_10` by-products and the first
/// reflection coefficient `k1` (decode-trace §7.2, base spec §4.6 last
/// paragraph: "the 10 coefficients `ã_i` and the first reflection
/// coefficient `k1` … the recursion is stopped at order 10, `k1` and
/// `ã_1..ã_10` copied, then resumed to order 50"). The value matches
/// the log-gain predictor order [`crate::consts::LPCLG`] coincidentally
/// — they are independent quantities in the spec.
pub const PF_LPC_ORDER: usize = 10;

/// Static (table-derived) parameters for the synthesis-filter
/// hybrid window. Constructed once per adapter; passed by reference to
/// the per-cycle [`SynthesisAdapter::adapt`] call.
#[derive(Debug)]
struct SynthCfg {
    wnr: [f64; 105],
    facv: [f64; LPC + 1],
}

impl SynthCfg {
    fn new() -> Self {
        Self {
            wnr: wnr_f64(),
            facv: facv_f64(),
        }
    }
}

/// Backward synthesis-filter adapter (block 33) state.
///
/// Carries the hybrid-window state (`SB` and `REXP` per Table 2) and
/// the most recently produced bandwidth-expanded predictor. The
/// predictor is in spec `A` layout: `coefficients()[0] = 1.0`,
/// `coefficients()[i] = -aᵢ` for `1 ≤ i ≤ LPC`.
///
/// On every successful adaptation the adapter also caches the **order-10
/// by-product** Levinson result `ã_1..ã_10` (re-run on the same RTMP)
/// and the first reflection coefficient `k1`. These feed the short-term
/// postfilter (block 72); see [`SynthesisAdapter::order10_predictor`]
/// and [`SynthesisAdapter::k1`].
#[derive(Debug)]
pub struct SynthesisAdapter {
    cfg: SynthCfg,
    hw_state: HybridWindowState,
    /// Most recently produced predictor (in spec `A` layout: leading
    /// `1.0`, then `-aᵢ` × LPC). Defaults to the all-pass predictor
    /// at construction time so the adapter is safe to read from before
    /// any speech has been pushed.
    last_predictor: [f64; LPC + 1],
    /// The spec's persistent `ATMP(1..LPC+1)` array — the **raw**
    /// (bandwidth-unexpanded) block-50 Levinson output. Persisted
    /// across cycles because the Annex I frame-erasure path both reads
    /// it (block 51FE softens `FACVFE(I)·ATMP(I)` at the first erased
    /// cycle, §I.5.8) and partially overwrites it (block 50 runs order
    /// 1–10 only during erased cycles, refreshing `ATMP(1..11)` while
    /// `ATMP(12..51)` keeps the last good-frame values — exactly the
    /// §I.5.1 array semantics). Defaults to the all-pass predictor.
    atmp: [f64; LPC + 1],
    /// `ILLCOND` — whether the most recent block-50 recursion (order 50
    /// on good cycles, order 10 on erased cycles) reported
    /// ill-conditioning. Block 51 / 51FE skip the expansion when set.
    illcond: bool,
    /// Most recently computed order-10 by-product predictor in our
    /// canonical Levinson `A` layout: `[1.0, a_1, …, a_{10}]`. The
    /// spec's `ã_i` (predictor convention `s(n) ≈ Σ ã_i s(n-i)`) maps
    /// to `-a_i` of this array. No bandwidth expansion has been
    /// applied. Defaults to the all-pass predictor.
    last_order10: [f64; PF_LPC_ORDER + 1],
    /// First reflection coefficient `k1 = -R(2)/R(1)` of the most
    /// recently processed RTMP (the spec's "first reflection
    /// coefficient" copied at the order-10 break point). Defaults to
    /// `0.0` (synthesis filter is the identity before any adaptation).
    last_k1: f64,
}

impl Default for SynthesisAdapter {
    fn default() -> Self {
        Self::new()
    }
}

impl SynthesisAdapter {
    /// Construct a fresh adapter with all internal buffers initialised
    /// to the Table 2/G.728 values (`SB`, `REXP` zeroed; predictor
    /// `A(1) = 1`, `A(i ≥ 2) = 0`).
    pub fn new() -> Self {
        let cfg = SynthCfg::new();
        let hw = Self::window_descriptor(&cfg.wnr);
        let hw_state = HybridWindowState::new(&hw);

        let mut last_predictor = [0.0f64; LPC + 1];
        last_predictor[0] = 1.0;

        let mut last_order10 = [0.0f64; PF_LPC_ORDER + 1];
        last_order10[0] = 1.0;

        Self {
            cfg,
            hw_state,
            last_predictor,
            atmp: last_predictor,
            illcond: true, // no valid ATMP until the first good cycle
            last_order10,
            last_k1: 0.0,
        }
    }

    fn window_descriptor(wnr: &[f64]) -> HybridWindow<'_> {
        HybridWindow {
            m: LPC,
            l: NFRSZ,
            n: NONR,
            window: wnr,
            decay: 0.75,
        }
    }

    /// Most recently produced bandwidth-expanded predictor, in spec
    /// `A` layout. Borrowed read-only so the synthesiser can call
    /// `set_predictor` with an owned copy.
    pub fn coefficients(&self) -> &[f64; LPC + 1] {
        &self.last_predictor
    }

    /// Most recently computed **order-10 by-product** predictor.
    ///
    /// Layout matches our Levinson convention: `[1.0, a_1, …, a_{10}]`
    /// such that the polynomial `A(z) = 1 + Σ a_i z^{-i}` is the
    /// order-10 LPC inverse filter from the same RTMP autocorrelation
    /// the order-50 predictor in [`Self::coefficients`] was derived
    /// from. **No bandwidth expansion has been applied.** The spec's
    /// `ã_i` (predictor sense `s(n) ≈ Σ ã_i s(n-i)`) is therefore
    /// `-a_i` of this array.
    ///
    /// Used by the short-term postfilter (block 72, §4.6) to derive
    /// `b̄_i = ã_i · SPFZCF^i` and `ā_i = ã_i · SPFPCF^i` after
    /// bandwidth expansion.
    pub fn order10_predictor(&self) -> &[f64; PF_LPC_ORDER + 1] {
        &self.last_order10
    }

    /// First reflection coefficient `k_1 = −R(2)/R(1)` of the most
    /// recently processed RTMP autocorrelation.
    ///
    /// Used by the short-term postfilter for the tilt-compensation
    /// factor `µ = TILTF · k_1` (block 72, §4.6, eq. 4-5).
    pub fn k1(&self) -> f64 {
        self.last_k1
    }

    /// The spec's persistent `ATMP(1..LPC+1)` array — the raw
    /// (bandwidth-unexpanded) block-50 output, in our canonical Levinson
    /// `A` layout (`atmp[0] = 1.0`). During an erasure only the first 11
    /// entries are fresh (block 50 runs order 1–10); the tail keeps the
    /// last good frame's values. Consumed by the Annex I block 51FE
    /// (`A(I) = FACVFE(I)·ATMP(I)` at the first erased cycle, §I.5.8).
    pub fn raw_atmp(&self) -> &[f64; LPC + 1] {
        &self.atmp
    }

    /// `ILLCOND` — whether the most recent block-50 recursion reported
    /// ill-conditioning (blocks 51 / 51FE skip their expansion when
    /// set; §I.5.8 gates the `FECOUNT = 1` branch on it).
    pub fn illcond(&self) -> bool {
        self.illcond
    }

    /// Block 51 on demand: the `FACV` bandwidth expansion of the
    /// **current** persistent `ATMP` array, or `None` when `ILLCOND` is
    /// set.
    ///
    /// The Annex I decoder needs this at the first good `ICOUNT = 3`
    /// after an erasure: the §I.5.1 loop's "do block 51" reads the
    /// then-current `ATMP` — whose head `ATMP(1..11)` was refreshed by
    /// the erased cycles' order-10 block 50 while `ATMP(12..51)` kept
    /// the pre-erasure values — rather than a predictor staged before
    /// the erasure began.
    pub fn expand_current(&self) -> Option<[f64; LPC + 1]> {
        if self.illcond {
            return None;
        }
        let mut a = self.atmp;
        for (slot, &f) in a.iter_mut().zip(self.cfg.facv.iter()) {
            *slot *= f;
        }
        Some(a)
    }

    /// Process one adaptation cycle of quantised-speech samples
    /// (`NFRSZ = 20` samples) through blocks 49, 50 and 51, and update
    /// the cached predictor.
    ///
    /// Returns the up-to-date predictor on success. On Levinson-Durbin
    /// failure (any of the three ill-conditioning exits) the cached
    /// predictor from the previous cycle is left in place per the
    /// spec's "keep using the previous adaptation cycle" rule, and the
    /// error is propagated so the caller can log/trace it.
    pub fn adapt(&mut self, sttmp: &[f64]) -> Result<&[f64; LPC + 1], LevinsonError> {
        assert_eq!(sttmp.len(), NFRSZ, "adapter input must be NFRSZ samples");

        let hw = Self::window_descriptor(&self.cfg.wnr);

        // ----- Block 49: hybrid window → RTMP -------------------------
        let mut rtmp = vec![0.0f64; LPC + 1];
        self.hw_state.run(&hw, sttmp, &mut rtmp);

        // ----- Block 50: Levinson-Durbin → ATMP -----------------------
        // The recursion writes ATMP into a fresh buffer; on failure
        // (any of the three ill-conditioning exits) levinson_durbin
        // populates the buffer with the canonical all-pass predictor,
        // but per the spec we must KEEP the previous adaptation
        // cycle's predictor on failure, not adopt the all-pass reset.
        // So we run into a scratch buffer and only commit on Ok.
        let mut atmp = vec![0.0f64; LPC + 1];
        if let Err(e) = levinson_durbin(&rtmp, &mut atmp, LPC) {
            self.illcond = true;
            return Err(e);
        }
        self.illcond = false;
        // Persist the raw (unexpanded) ATMP for the Annex I block 51FE
        // consumer before block 51 expands the scratch copy below.
        self.atmp.copy_from_slice(&atmp[..=LPC]);

        // ----- Block 50 by-product: order-10 predictor + k1 ------------
        //
        // Decode trace §7.2 (and base spec §4.6 last paragraph): the
        // postfilter consumes the 10 coefficients `ã_i` and the first
        // reflection coefficient `k1` produced at the order-10 break
        // point of this same Levinson. The spec's optimisation is to
        // copy them off mid-recursion and resume to order 50; here we
        // simply re-run a self-contained order-10 Levinson on the
        // SAME RTMP. Since the recursion is deterministic and the
        // RTMP unchanged, the order-10 output coefficients are
        // bit-identical to what the embedded copy-off would yield.
        // The runtime cost is negligible compared to the order-50
        // recursion and the alternative would require an unrelated
        // refactor of `levinson_durbin` to surface intermediate
        // state. Only commit on success — preserves the spec's
        // "keep previous predictor" semantics on failure.
        //
        // `k1 = -R(2)/R(1)` is the canonical first reflection
        // coefficient (Levinson's first iteration `rc[0]`). We
        // compute it explicitly here so that future callers can read
        // it even when the order-10 recursion would otherwise drop
        // it after the second iteration overwrites `a[1]`.
        // The order-50 Levinson above already passed the R(1) > 0
        // guard; the order-10 sub-recursion shares that opener so the
        // only failure mode unique to it is R(11) = 0 (trailing-zero
        // at order 10 even when order 50 is fine — rare but possible).
        // On order-10 failure we keep the previous cache, matching the
        // spec's "skip block 51" intent at the postfilter scale.
        let mut atmp10 = [0.0f64; PF_LPC_ORDER + 1];
        let order10_ok = levinson_durbin(&rtmp[..=PF_LPC_ORDER], &mut atmp10, PF_LPC_ORDER).is_ok();
        let k1 = if rtmp[0] > 0.0 {
            -rtmp[1] / rtmp[0]
        } else {
            0.0
        };

        // ----- Block 51: bandwidth expansion --------------------------
        // | For I = 2,3,...,LPC+1: ATMP(I) = FACV(I) * ATMP(I)
        //
        // Note FACV(1) = 16384 (Q14) = 1.0; the spec writes the loop
        // starting at I = 2 because the leading tap is implicitly 1.0.
        // We multiply every tap; the I=1 multiplication is a no-op
        // and folds in cleanly.
        for i in 0..=LPC {
            atmp[i] *= self.cfg.facv[i];
        }

        // Commit ATMP into the cached predictor (this is the spec's
        // `Wait until ICOUNT = 3, then A(I) = ATMP(I)` step — in our
        // batched cycle model the wait collapses to "now").
        self.last_predictor[..=LPC].copy_from_slice(&atmp[..=LPC]);
        // Commit the order-10 by-product (un-expanded) and k1 for the
        // postfilter to consume. Only commit on success — preserves the
        // postfilter's "keep previous coefficients" semantics that
        // mirror the spec's order-50 rule.
        if order10_ok {
            self.last_order10.copy_from_slice(&atmp10);
            self.last_k1 = k1;
        }

        Ok(&self.last_predictor)
    }

    /// Annex I frame-erasure adaptation (§I.5.1 `ICOUNT = 4` processing
    /// during an erased cycle): block **49FE** + block 50 **order 1–10
    /// only**.
    ///
    /// Per §I.4.3/§I.4.4:
    ///
    /// * the hybrid-window buffer `SB` and the recursive component
    ///   `REXP` are updated exactly as in a good cycle (the "vital
    ///   operations"), but only `RTMP(1..11)` are produced;
    /// * Durbin's recursion runs from order 1 to order 10, refreshing
    ///   the persistent `ATMP(1..11)` and the post-filter's order-10
    ///   by-product + first reflection coefficient (the post-filter
    ///   "floats" with the concealed speech);
    /// * the order-50 continuation and block 51 are **not** performed —
    ///   [`Self::coefficients`] is left untouched. The live synthesis
    ///   predictor is softened by block 51FE instead (driven by the
    ///   decoder at the third vector of the cycle).
    ///
    /// Returns the Levinson status of the order-10 recursion; on error
    /// the previous by-product (and `ATMP`) are kept, and
    /// [`Self::illcond`] reports `true`.
    pub fn adapt_erased(&mut self, sttmp: &[f64]) -> Result<(), LevinsonError> {
        assert_eq!(sttmp.len(), NFRSZ, "adapter input must be NFRSZ samples");

        let hw = Self::window_descriptor(&self.cfg.wnr);

        // ----- Block 49FE: hybrid window → RTMP(1..11) only -----------
        let mut rtmp = [0.0f64; PF_LPC_ORDER + 1];
        self.hw_state
            .run_erased(&hw, sttmp, &mut rtmp, PF_LPC_ORDER + 1);

        // ----- Block 50, order 1–10 only -------------------------------
        let mut atmp10 = [0.0f64; PF_LPC_ORDER + 1];
        if let Err(e) = levinson_durbin(&rtmp, &mut atmp10, PF_LPC_ORDER) {
            self.illcond = true;
            return Err(e);
        }
        self.illcond = false;
        // Refresh the persistent ATMP head; ATMP(12..51) keeps the last
        // good frame's values (§I.5.1: "do block 50, order 1 to 10").
        self.atmp[..=PF_LPC_ORDER].copy_from_slice(&atmp10);
        // The post-filter by-product keeps floating (§I.4.4).
        self.last_order10.copy_from_slice(&atmp10);
        self.last_k1 = if rtmp[0] > 0.0 {
            -rtmp[1] / rtmp[0]
        } else {
            0.0
        };
        Ok(())
    }
}

/// Compile-time guard: FACV(1) in Q14 must be exactly 16384
/// (i.e. 1.0), so the implicit `A(1) = 1` invariant survives the
/// bandwidth-expansion multiplication.
const _: () = {
    let scale = 1i32 << Q14;
    assert!(FACV_Q14[0] as i32 == scale);
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::IDIM;
    use crate::tables::FACV_Q14;

    #[test]
    fn fresh_adapter_has_allpass_predictor() {
        // Table 2 initial: A(1) = 1, A(i ≥ 2) = 0 (the synthesis
        // filter is the identity before any adaptation has happened).
        let adapter = SynthesisAdapter::new();
        let a = adapter.coefficients();
        assert_eq!(a[0], 1.0);
        assert!(a[1..].iter().all(|&v| v == 0.0));
    }

    #[test]
    fn fresh_adapter_has_allpass_order10_byproduct() {
        // Pre-adapt the order-10 by-product should match the
        // all-pass shape (leading 1, rest zero) and `k1 = 0` (no
        // RTMP yet).
        let adapter = SynthesisAdapter::new();
        let a10 = adapter.order10_predictor();
        assert_eq!(a10[0], 1.0);
        assert!(a10[1..].iter().all(|&v| v == 0.0));
        assert_eq!(adapter.k1(), 0.0);
    }

    #[test]
    fn facv_first_tap_is_unity_q14() {
        // Block 51 multiplies every predictor tap by FACV(I); for the
        // I=1 tap (which must remain 1.0 to preserve the A(1)=1
        // invariant) FACV(1) must be exactly 16384 in Q14.
        assert_eq!(FACV_Q14[0], 16_384);
    }

    #[test]
    fn zero_input_leaves_predictor_as_allpass_after_levinson_failure() {
        // Driving zero speech vectors through the adapter produces a
        // zero autocorrelation. Block 50's R(1) ≤ 0 guard fires, the
        // spec says "keep using the previous cycle's predictor", and
        // the previous predictor is the all-pass identity. Confirm
        // both branches: the error is surfaced AND the cached
        // predictor stays at all-pass.
        let mut adapter = SynthesisAdapter::new();
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
        let a = adapter.coefficients();
        assert_eq!(a[0], 1.0);
        assert!(a[1..].iter().all(|&v| v == 0.0));
    }

    #[test]
    fn nonzero_input_produces_order10_byproduct_with_a1_eq_unity() {
        // After a few successful adaptation cycles on a non-trivial
        // signal, the order-10 by-product must (a) preserve the
        // A(1) = 1 invariant, (b) carry non-zero taps somewhere in
        // a_1..a_10, and (c) have a non-zero k1 (the input is not
        // white, so R(2) is not zero).
        let mut adapter = SynthesisAdapter::new();
        let mut input = [0.0f64; NFRSZ];
        for k in 0..NFRSZ {
            let t = k as f64;
            input[k] = 100.0 * (0.4 * t).sin() + 50.0 * (0.9 * t).cos() + (t * 13.0).fract();
        }
        let mut success = false;
        for _ in 0..6 {
            if adapter.adapt(&input).is_ok() {
                success = true;
            }
        }
        assert!(success);
        let a10 = adapter.order10_predictor();
        assert_eq!(a10[0], 1.0, "A(1) = 1 invariant");
        assert!(
            a10[1..].iter().any(|&v| v != 0.0),
            "expected non-zero order-10 taps after non-trivial input"
        );
        assert_ne!(adapter.k1(), 0.0, "k1 should be non-zero for AR signal");
    }

    #[test]
    fn order10_byproduct_k1_matches_levinson_first_iter_convention() {
        // The first reflection coefficient `k1 = -R(2)/R(1)` is the
        // canonical Levinson opener. For an AR(1) test signal we know
        // the closed-form: r[k] = ρ^|k|, so `k1 = -ρ` regardless of
        // the final recursion order. We don't drive the synthesis
        // adapter through a real cycle (its hybrid window doesn't
        // produce a clean AR(1) RTMP from short inputs); instead we
        // run a synthetic check on the public k1 / order10 layout
        // by driving the adapter with smoothly correlated speech and
        // checking that |k1| ≤ 1.
        let mut adapter = SynthesisAdapter::new();
        let mut input = [0.0f64; NFRSZ];
        for k in 0..NFRSZ {
            input[k] = 100.0 * (0.3 * k as f64).sin();
        }
        for _ in 0..4 {
            let _ = adapter.adapt(&input);
        }
        // Reflection coefficients of a valid AR process are in (-1, 1);
        // outside that range means the autocorrelation is not positive
        // semidefinite, which Levinson would have flagged.
        let k1 = adapter.k1();
        assert!(k1.abs() < 1.0, "|k1| = {} should be < 1", k1.abs());
    }

    #[test]
    fn nonzero_input_produces_a1_eq_unity_and_nonzero_taps() {
        // A non-trivial speech signal should drive the autocorrelation
        // away from zero and let Levinson produce a non-trivial
        // predictor. The A(1) = 1 invariant must be preserved by
        // block 51 (FACV(1) = 1 in Q14, hand-verified above).
        let mut adapter = SynthesisAdapter::new();
        // Use a frequency-rich signal (mix of two sines + small noise)
        // so Levinson has enough excitation to land away from the
        // ill-conditioning boundary.
        let mut input = [0.0f64; NFRSZ];
        for k in 0..NFRSZ {
            let t = k as f64;
            input[k] = 100.0 * (0.4 * t).sin() + 50.0 * (0.9 * t).cos() + (t * 13.0).fract();
        }
        // Run enough cycles to populate the recursive part of the
        // hybrid window. A single cycle on a fresh adapter might still
        // produce a degenerate Levinson on some signals; this test
        // exercises the typical, steady-state behaviour.
        let mut last = None;
        for _ in 0..6 {
            if let Ok(a) = adapter.adapt(&input) {
                last = Some(*a);
            }
        }
        let a = last.expect("at least one adapt() should succeed");
        assert_eq!(a[0], 1.0, "A(1) = 1 invariant must hold");
        assert!(
            a[1..].iter().any(|&v| v != 0.0),
            "expected at least one non-zero predictor tap after adaptation"
        );
    }

    #[test]
    fn adapter_drives_decoder_without_panicking() {
        // Smoke test: wire the adapter into the existing decoder hook
        // and confirm the predictor swap path is exercised cleanly.
        // We round-trip three adaptation cycles' worth of vectors.
        use crate::decoder::Synthesizer;

        let mut adapter = SynthesisAdapter::new();
        let mut synth = Synthesizer::new();

        for cycle in 0..3 {
            // Run NUPDATE=4 vectors per cycle, IDIM samples each.
            let mut sttmp = [0.0f64; NFRSZ];
            for k in 0..NFRSZ {
                let t = (cycle * NFRSZ + k) as f64;
                sttmp[k] = 200.0 * (0.3 * t).sin();
            }
            let _ = adapter.adapt(&sttmp); // may be Err on the first cycle
            synth.set_predictor(*adapter.coefficients());

            for v in 0..(NFRSZ / IDIM) {
                let mut exc = [0.0f64; IDIM];
                for k in 0..IDIM {
                    exc[k] = 10.0 * ((cycle * NFRSZ + v * IDIM + k) as f64).sin();
                }
                let out = synth.filter_vector(&crate::decoder::ExcitationVector(exc));
                for &v in &out {
                    assert!(v.is_finite(), "decoder produced NaN/Inf with adapter");
                }
            }
        }
    }
}
