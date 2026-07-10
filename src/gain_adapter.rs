//! Backward vector gain adapter — block 30 of Figure 6/G.728.
//!
//! Per §5.7 the block-30 chain reads the previous gain-scaled
//! excitation vector `ET(1..IDIM)` and produces the current vector's
//! predicted excitation gain `σ(n)` (= `GAIN` after block 48).
//!
//! ## Per-vector dataflow
//!
//! Every speech vector runs blocks 67/39/40/42/46/47/48; blocks
//! 43/44/45 only run when `ICOUNT = 2`.
//!
//! ```text
//! ET(1..IDIM) ─► [67: 1-vec delay]
//!              ─► [39: RMS]      = sqrt(Σ ET²) / sqrt(IDIM)
//!              ─► [40: log10·20] = 10·log10(Σ ET² / IDIM)   (clamped ≥ 1)
//!              ─► [42: -GOFF]    = ETRMS - 32
//!                                  → GSTATE(1)
//!
//! Every 4th vector (ICOUNT=2):
//! GTMP(1..4) ─► [43: hybrid window] → R(1..LPCLG+1)
//!            ─► [44: Levinson]     → GPTMP(1..LPCLG+1)
//!            ─► [45: BW expand]    → GP(I) = FACGPV(I) · GPTMP(I)
//!
//! Every vector:
//! GP, GSTATE ─► [46: log-gain LP]  → GAIN = -Σ GP(I+1)·GSTATE(I)
//!             ─► [+GOFF]            → GAIN += 32
//!             ─► [47: limiter]      → clamp to [0, 60]
//!             ─► [48: 10^(·/20)]    → GAIN = 10^(GAIN/20)   (= σ(n))
//! ```
//!
//! ## GTMP / GSTATE buffering convention
//!
//! Per §5.6 block 43 prose: "The GTMP array below contains four
//! offset-removed log-gain values, starting from the log-gain of the
//! second vector of the previous adaptation cycle to the log-gain of
//! the first vector of the current adaptation cycle, which is
//! GTMP(1). GTMP(4) is the offset-removed log-gain value from the
//! first vector of the current adaptation cycle, the newest value."
//!
//! Wait — the spec actually says GTMP(1) is the *second vector of
//! the previous cycle* (oldest) and GTMP(4) is the *first vector of
//! the current cycle* (newest). We mirror that ordering exactly: a
//! 4-deep FIFO that pushes the newest offset-removed log-gain at
//! `gtmp[3]` and shifts older values down.

use crate::consts::{DIMINV, GOFF, IDIM, LPCLG, NONRLG, NUPDATE};
use crate::hybrid_window::{HybridWindow, HybridWindowState};
use crate::levinson::levinson_durbin;
use crate::tables::{facgpv_f64, wnrg_f64, FACGPV_Q14, Q14};

/// Spec-stated initial value for every entry of the GSTATE and GTMP
/// log-gain buffers (Table 2/G.728). The choice of `-32` (the
/// negative of `GOFF`) makes the first decoded vector's predicted
/// gain start at the all-pass `10^0 = 1.0` once GOFF is added back.
pub const GSTATE_INIT_DB: f64 = -GOFF;

/// Backward vector gain adapter (block 30) state.
#[derive(Debug)]
pub struct GainAdapter {
    /// Static window data — wnrg (length 34).
    wnrg: [f64; 34],
    /// Static bandwidth-broadening vector — FACGPV (length LPCLG+1).
    facgpv: [f64; LPCLG + 1],
    /// Hybrid-window state for block 43 (SBLG buffer + REXPLG).
    hw_state: HybridWindowState,
    /// 4-deep buffer of offset-removed log-gain values seen during
    /// the current and previous adaptation cycle's first half.
    /// `gtmp[0]` is the oldest, `gtmp[3]` is the newest. Initialised
    /// to GSTATE_INIT_DB = -32 dB per Table 2/G.728 (GTMP initial
    /// `-32,-32,-32,-32`).
    gtmp: [f64; NUPDATE],
    /// Log-gain predictor memory `GSTATE(1..LPCLG)`. `gstate[0]` is
    /// the spec's `GSTATE(1)` — the most recent offset-removed
    /// log-gain. Initialised to GSTATE_INIT_DB = -32 dB per Table 2
    /// (GSTATE initial `-32,-32,-32,...,-32`).
    gstate: [f64; LPCLG],
    /// Most recent bandwidth-expanded log-gain predictor `GP(1..LPCLG+1)`.
    /// `gp[0]` = `GP(1)` is not used by block 46; we still store it
    /// to keep indexing 1:1 with the spec. Defaults to the all-pass
    /// predictor `(1, 0, 0, ..., 0)` so the first predicted log-gain
    /// is simply `GOFF` (= unity linear gain).
    gp: [f64; LPCLG + 1],
    /// Vector index within the current adaptation cycle (1..=NUPDATE,
    /// 1-based to match the spec's `ICOUNT`). Starts at 1.
    icount: usize,
}

impl Default for GainAdapter {
    fn default() -> Self {
        Self::new()
    }
}

impl GainAdapter {
    /// Construct a fresh gain adapter with Table 2/G.728 initial
    /// values for every internal buffer.
    pub fn new() -> Self {
        let wnrg = wnrg_f64();
        let facgpv = facgpv_f64();
        let hw = HybridWindow {
            m: LPCLG,
            l: NUPDATE,
            n: NONRLG,
            window: &wnrg,
            decay: 0.75,
        };
        let hw_state = HybridWindowState::new(&hw);

        // Table 2/G.728 initial predictor: GP = (1, −1, 0, …, 0)
        // ("Initial value 1,–1,0,0,..."). The −1 second tap makes the
        // initial prediction δ̂(n) = δ(n−1); with GSTATE pre-loaded at
        // −32 dB the first block-46 evaluation is −(−1)·(−32) = −32 dB,
        // +GOFF cancels to 0 dB, so the first predicted σ(n) is
        // exactly 1.0.
        let mut gp = [0.0f64; LPCLG + 1];
        gp[0] = 1.0;
        gp[1] = -1.0;

        Self {
            wnrg,
            facgpv,
            hw_state,
            gtmp: [GSTATE_INIT_DB; NUPDATE],
            gstate: [GSTATE_INIT_DB; LPCLG],
            gp,
            icount: 1,
        }
    }

    /// Borrow the current bandwidth-expanded log-gain predictor (for
    /// tests and audit).
    pub fn gp(&self) -> &[f64; LPCLG + 1] {
        &self.gp
    }

    /// Borrow the log-gain predictor memory (for tests).
    pub fn gstate(&self) -> &[f64; LPCLG] {
        &self.gstate
    }

    /// Current `ICOUNT` (1-based, in `1..=NUPDATE`).
    pub fn icount(&self) -> usize {
        self.icount
    }

    /// Run block-30 for the next decoded vector.
    ///
    /// `et_prev` is the **previous** gain-scaled excitation vector
    /// `ET(1..IDIM)` (the block-67 1-vector delay output). The
    /// function returns the predicted linear gain `σ(n)` for the
    /// vector that is about to be decoded.
    ///
    /// Advances `ICOUNT` modulo NUPDATE after each call.
    pub fn predict_next(&mut self, et_prev: &[f64; IDIM]) -> f64 {
        // With OGAINDB/AFTERFE at their inert values the Annex I block
        // 47AF reduces exactly to the base block-47 limiter, so the
        // non-erasure path is unchanged.
        self.predict_next_limited(et_prev, 0.0, 0, false).0
    }

    /// Blocks 67 + 39 + 40 (§5.7) — energy of the previous excitation
    /// vector in dB (`ETRMS`, the value block 42 offsets into
    /// `GSTATE(1)` and Annex I's block 97FE assigns to `OGAINDB`).
    fn etrms_db(et_prev: &[f64; IDIM]) -> f64 {
        // | ETRMS = ET(1) * ET(1)
        // | For K = 2,3,...,IDIM: ETRMS = ETRMS + ET(K) * ET(K)
        // | ETRMS = ETRMS * DIMINV                  | divide by IDIM
        // | If ETRMS < 1, set ETRMS = 1             | clip to avoid log overflow
        // | ETRMS = 10 * log10(ETRMS)
        let mut etrms = 0.0f64;
        for k in 0..IDIM {
            etrms += et_prev[k] * et_prev[k];
        }
        etrms *= DIMINV;
        if etrms < 1.0 {
            etrms = 1.0;
        }
        10.0 * etrms.log10()
    }

    /// Block 42 + the GTMP FIFO push — insert the previous vector's
    /// offset-removed log-gain at the head of the predictor memory.
    fn push_log_gain(&mut self, etrms: f64) {
        // ----- Block 42: log-gain offset subtractor ------------------
        // | GSTATE(1) = ETRMS - GOFF
        //
        // The spec writes `GSTATE(1) = ETRMS - GOFF` here, but
        // GSTATE(1) is also the first tap of the predictor's memory
        // — the new value must enter at the head and the older
        // values must shift back. We shift first so the new entry
        // landing at gstate[0] doesn't overwrite an unsaved value.
        for i in (1..LPCLG).rev() {
            self.gstate[i] = self.gstate[i - 1];
        }
        let new_gstate1 = etrms - GOFF;
        self.gstate[0] = new_gstate1;

        // Also push the new value into the GTMP FIFO for block 43.
        // GTMP(4) is the newest, GTMP(1) is the oldest. Shift-then-write.
        for i in 0..(NUPDATE - 1) {
            self.gtmp[i] = self.gtmp[i + 1];
        }
        self.gtmp[NUPDATE - 1] = new_gstate1;
    }

    /// Annex I erased-vector advance (§I.4.3 "vital operations"): run
    /// blocks 67/39/40/42 on the previous (extrapolated) excitation —
    /// this **is** the float block 97FE, which the annex describes as
    /// "the same code as Blocks 67, 39, 40 and 42" — and, at
    /// `ICOUNT = 2`, run block **43FE** (update the `SBLG` buffer and
    /// the recursive component `REXPLG`, skip the autocorrelation
    /// output and blocks 44/45 whose results would be discarded).
    ///
    /// No predicted gain is produced (block 31FE supplies the erased
    /// vector's excitation directly). Returns the previous vector's
    /// `ETRMS` in dB — the §I.5.1 `OGAINDB = output of block 97FE`
    /// value the caller must latch while the erasure lasts.
    ///
    /// Advances `ICOUNT` modulo NUPDATE, exactly like
    /// [`Self::predict_next`], so good/erased vectors interleave.
    pub fn advance_erased(&mut self, et_prev: &[f64; IDIM]) -> f64 {
        let etrms = Self::etrms_db(et_prev);
        self.push_log_gain(etrms);

        // Block 43FE at ICOUNT = 2: buffer + recursive-component update
        // only ("If FERROR = .TRUE., skip all the lines below").
        if self.icount == 2 {
            let hw = HybridWindow {
                m: LPCLG,
                l: NUPDATE,
                n: NONRLG,
                window: &self.wnrg,
                decay: 0.75,
            };
            let gtmp = self.gtmp;
            self.hw_state.advance_erased(&hw, &gtmp);
        }

        self.icount = (self.icount % NUPDATE) + 1;
        etrms
    }

    /// Run block-30 for the next decoded vector with the Annex I block
    /// **47AF** post-erasure gain-growth limiter in place of the base
    /// block 47 (§I.5.6).
    ///
    /// * `ogaindb` — the last predicted gain in dB (§I.5.1 `OGAINDB`);
    /// * `afterfe` — the remaining gain-clamp cycles (0 ⇒ clamp off,
    ///   making this identical to [`Self::predict_next`]);
    /// * `prev_erased` — whether the previous vector was concealed; if
    ///   so, `et_prev` is the extrapolated excitation and its `ETRMS`
    ///   (the block-97FE output) supersedes `ogaindb` per §I.5.1.
    ///
    /// Returns `(σ(n), limited log-gain in dB)`; the caller stores the
    /// latter as the next `OGAINDB`.
    pub fn predict_next_limited(
        &mut self,
        et_prev: &[f64; IDIM],
        ogaindb: f64,
        afterfe: usize,
        prev_erased: bool,
    ) -> (f64, f64) {
        let etrms = Self::etrms_db(et_prev);
        let ogaindb = if prev_erased { etrms } else { ogaindb };
        self.push_log_gain(etrms);

        // ----- ICOUNT = 2: re-derive log-gain predictor --------------
        // | This block is executed only when ICOUNT = 2, after
        // | block 43 is executed.
        if self.icount == 2 {
            self.update_log_gain_predictor();
        }

        // ----- Block 46: log-gain linear predictor -------------------
        // | GAIN = 0
        // | For I = LGLPC, LGLPC-1, ..., 3, 2: GAIN = GAIN - GP(I+1)·GSTATE(I); GSTATE(I) = GSTATE(I-1)
        // | GAIN = GAIN - GP(2) * GSTATE(1)
        //
        // The spec walks GSTATE(I) backwards from LPCLG..2 and
        // shifts GSTATE(I) ← GSTATE(I-1) as it goes. We already did
        // the shift up at "Block 42" above to make room for the new
        // GSTATE(1) entry, so here we ONLY compute the dot product
        // — re-shifting would corrupt the predictor memory. The
        // spec's combined shift+dot pattern is equivalent to a
        // shift-then-dot when laid out across the full vector loop.
        let mut gain = 0.0f64;
        for i_spec in (2..=LPCLG).rev() {
            // GP(I+1) → gp[i_spec]; GSTATE(I) → gstate[i_spec - 1]
            gain -= self.gp[i_spec] * self.gstate[i_spec - 1];
        }
        gain -= self.gp[1] * self.gstate[0];

        // ----- "Log-gain offset adder" (between blocks 46 and 47) ---
        // | GAIN = GAIN + GOFF
        gain += GOFF;

        // ----- Block 47AF: log-gain limiter (§I.5.6) -----------------
        // | If GAIN < 0., set GAIN = 0.    | linear gain 1
        // | If GAIN > 60., set GAIN = 60. | linear gain 1000
        // | TMP = GAIN - OGAINDB
        // | If AFTERFE > 0 and TMP > FEGAINMAX, GAIN = OGAINDB + FEGAINMAX
        //
        // With `afterfe == 0` this is exactly the base block 47 (§5.7):
        // the [0, 60] range clamp and nothing else.
        gain = crate::gain_growth_limiter::limit_log_gain(gain, ogaindb, afterfe);
        let gain_db = gain;

        // ----- Block 48: inverse logarithm ---------------------------
        // | GAIN = 10^(GAIN/20)
        gain = 10.0f64.powf(gain / 20.0);

        // Advance ICOUNT modulo NUPDATE.
        self.icount = (self.icount % NUPDATE) + 1;

        (gain, gain_db)
    }

    /// Inner helper: run blocks 43 (hybrid window), 44 (Levinson) and
    /// 45 (bandwidth expansion) on the current GTMP buffer to update
    /// the bandwidth-expanded log-gain predictor `gp`. Called from
    /// `predict_next` when `ICOUNT = 2`.
    ///
    /// If Levinson-Durbin reports ill-conditioning, the spec says
    /// "skip blocks 44 and 45 without updating the log-gain predictor
    /// coefficients (...we keep using the old log-gain predictor
    /// coefficients determined in the previous adaptation cycle)" —
    /// so we leave `self.gp` unchanged on error.
    fn update_log_gain_predictor(&mut self) {
        let hw = HybridWindow {
            m: LPCLG,
            l: NUPDATE,
            n: NONRLG,
            window: &self.wnrg,
            decay: 0.75,
        };

        // ----- Block 43: hybrid window → R ---------------------------
        let mut r = vec![0.0f64; LPCLG + 1];
        self.hw_state.run(&hw, &self.gtmp, &mut r);

        // ----- Block 44: Levinson-Durbin → GPTMP ---------------------
        let mut gptmp = vec![0.0f64; LPCLG + 1];
        if levinson_durbin(&r, &mut gptmp, LPCLG).is_err() {
            // Spec: keep previous-cycle gp on Levinson failure.
            return;
        }

        // ----- Block 45: bandwidth expansion -------------------------
        // | For I = 2,3,...,LPCLG+1: GP(I) = FACGPV(I) * GPTMP(I)
        //
        // FACGPV(1) = 16384 (Q14) = 1.0, so multiplying the I=1 tap
        // is a no-op; we run the loop from 0 for symmetry with block
        // 51 (same compile-time guard as `synthesis_adapter`).
        for i in 0..=LPCLG {
            self.gp[i] = self.facgpv[i] * gptmp[i];
        }
    }
}

/// Compile-time guard: FACGPV(1) in Q14 must be exactly 16384, i.e.
/// 1.0. (Mirrors the synthesis-adapter check on FACV.)
const _: () = {
    let scale = 1i32 << Q14;
    assert!(FACGPV_Q14[0] as i32 == scale);
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fresh_adapter_has_table2_initial_state() {
        // Table 2/G.728 lists `-32,...,-32` for GSTATE and GTMP, and
        // `1,-1,0,0,...` for GP (the "predict the previous log-gain"
        // seed predictor).
        let a = GainAdapter::new();
        for &v in a.gstate() {
            assert_eq!(v, GSTATE_INIT_DB);
        }
        assert_eq!(a.gp()[0], 1.0);
        assert_eq!(a.gp()[1], -1.0);
        assert!(a.gp()[2..].iter().all(|&v| v == 0.0));
        assert_eq!(a.icount(), 1);
    }

    #[test]
    fn first_predict_with_zero_et_returns_unity_gain() {
        // With ET(1..IDIM) = 0, ETRMS is clipped to 1 → log10(1)·10 = 0
        // → GSTATE(1) = 0 - GOFF = -32 dB. Table 2's GP = (1, −1, 0…)
        // predictor then evaluates block 46 as −GP(2)·GSTATE(1) =
        // −(−1)·(−32) = −32 dB; adding GOFF back gives 0 dB, so the
        // first predicted σ(n) is exactly 1.0 — the spec-correct
        // "first vector" behaviour (the fixed-point Annex G adapter
        // produces the same unity gain from its GP/GSTATE inits).
        let mut a = GainAdapter::new();
        let et_zero = [0.0f64; IDIM];
        let sigma = a.predict_next(&et_zero);
        assert!(
            (sigma - 1.0).abs() < 1e-9,
            "expected first-vector sigma = 1.0, got {}",
            sigma
        );
    }

    #[test]
    fn predict_advances_icount_modulo_nupdate() {
        // ICOUNT walks 1 → 2 → 3 → 4 → 1 → 2 → ... per spec §3.7.
        let mut a = GainAdapter::new();
        let et_zero = [0.0f64; IDIM];
        let expected_seq = [2usize, 3, 4, 1, 2, 3, 4, 1];
        for &want in &expected_seq {
            a.predict_next(&et_zero);
            assert_eq!(a.icount(), want);
        }
    }

    #[test]
    fn predict_clips_to_limiter_range() {
        // The limiter clamps GAIN ∈ [0, 60] dB ⇒ σ ∈ [1, 1000].
        // Drive massive excitation magnitudes through the adapter
        // and confirm sigma never escapes [1, 1000].
        let mut a = GainAdapter::new();
        // Pump in maximum-magnitude vectors for many cycles; the
        // adapter should clamp without panicking.
        let et_loud = [4_095.0f64; IDIM];
        for _ in 0..32 {
            let sigma = a.predict_next(&et_loud);
            assert!(
                (1.0..=1_000.0).contains(&sigma),
                "sigma {} escaped limiter range [1, 1000]",
                sigma
            );
        }
        // And silence: should bottom-clamp to >= 1.0 (linear gain 1).
        let et_silent = [0.0f64; IDIM];
        for _ in 0..32 {
            let sigma = a.predict_next(&et_silent);
            assert!(sigma >= 1.0);
            assert!(sigma <= 1_000.0);
        }
    }

    #[test]
    fn facgpv_first_tap_is_unity_q14() {
        // Sister test of synthesis_adapter::facv_first_tap_is_unity_q14.
        // Block 45's I=2 starting index relies on FACGPV(1) = 1.0
        // (in Q14: 16384) so the leading tap survives a no-op
        // multiply.
        assert_eq!(FACGPV_Q14[0], 16_384);
    }

    #[test]
    fn predictor_updates_only_at_icount_eq_2() {
        // Drive vectors and confirm the GP predictor changes only
        // when ICOUNT was 2 (immediately after the third predict_next
        // call: 1→2, then *during* the 2→3 call the spec runs
        // blocks 43/44/45 because the call enters with icount = 2).
        let mut a = GainAdapter::new();
        // Use a frequency-rich ET to give Levinson enough excitation.
        let mut et = [0.0f64; IDIM];
        for k in 0..IDIM {
            et[k] = 500.0 * (k as f64 * 0.7).sin();
        }
        let gp_initial = *a.gp();
        // vector 1 (icount=1) → gp unchanged
        a.predict_next(&et);
        assert_eq!(*a.gp(), gp_initial);
        // vector 2 enters with icount=2, runs blocks 43/44/45; but
        // with only 4 GTMP slots of which 3 are still the initial -32
        // and 1 was just written, the autocorrelation may degenerate
        // and Levinson may keep the old predictor. We don't require
        // the predictor to change YET — only that the path runs
        // without panic, exiting with icount=3.
        a.predict_next(&et);
        assert_eq!(a.icount(), 3);
        // vector 3 (icount=3) → gp unchanged
        let gp_at_icount3 = *a.gp();
        a.predict_next(&et);
        assert_eq!(*a.gp(), gp_at_icount3);
        // vector 4 (icount=4) → gp unchanged
        let gp_at_icount4 = *a.gp();
        a.predict_next(&et);
        assert_eq!(*a.gp(), gp_at_icount4);
    }
}
