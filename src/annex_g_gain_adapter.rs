//! Annex G — the fixed-point backward vector gain adapter (Figure
//! G.1/G.728, blocks 43 – 48 and 93 – 99).
//!
//! ITU-T G.728 Annex G (1994-11) reformulates the backward vector gain
//! adapter (block 20 of the base Recommendation) for fixed point:
//! instead of computing the RMS of the previous gain-scaled excitation
//! and taking its logarithm (blocks 39/40/42 of Figure 6/G.728), the
//! offset-removed log-gain `δ(n−1)` is reconstructed from the *chosen
//! codebook indices* via two dB look-up tables (§G.2.1, eq. G-14):
//!
//! ```text
//! δ(n−1) = δ̂(n−1) + 20·log10|g_i| + 10·log10 P[y_j]
//! ```
//!
//! where `δ̂(n−1)` is the previous predicted log-gain (block 95), the
//! gain-table term comes from block 93 ([`crate::tables::GCBLG_Q11`])
//! and the shape-table term from block 94
//! ([`crate::tables::SHAPELG_Q11`]). The adder 96 / limiter 97 write
//! the result into `GSTATE(1)` (§G.3.16 fixed-point pseudo-code).
//!
//! The once-per-cycle predictor update is the block 43 hybrid window
//! (§G.3.14, [`LogGainWindowFixed`]) → block 44 Levinson-Durbin
//! (§G.2.2, [`crate::annex_g_levinson`]) → block 45 bandwidth
//! expansion (§G.3.15,
//! [`crate::annex_g_synth_adapter::log_gain_bandwidth_expand`]); the
//! per-vector prediction is block 46 (log-gain linear prediction) →
//! block 98 (log-gain limiter, −32 … +28 dB) → block 99 (offset adder,
//! `GOFF = 32 dB = 16384` Q9) → block 48 (inverse-logarithm, §G.3.16)
//! producing the scalar-floating-point excitation gain
//! `σ(n) = GAIN · 2^(−NLSGAIN)` consumed by blocks 16 and 21 of the
//! §G.3 coder ([`crate::annex_g_codebook`]).
//!
//! Everything in this module is transcribed from the §G.3.14 / §G.3.16
//! fixed-point pseudo-code and the §G.2.1 prose; the dB tables are the
//! Table G.3 / Table G.4 integer columns.

use crate::annex_g_arith::{findnls, rnd};
use crate::annex_g_hybrid::{HwmcoreOut, HwmcoreState};
use crate::consts::{LPCLG, NONRLG, NUPDATE};
use crate::tables::{GCBLG_Q11, SHAPELG_Q11, WNRG_Q15};

/// `N1 = LPCLG + NUPDATE = 14` — end of the recursive-component taps.
pub const N1LG: usize = LPCLG + NUPDATE;
/// `N2 = LPCLG + NONRLG = 30` — length of the kept-old part of `SBLG`.
pub const N2LG: usize = LPCLG + NONRLG;
/// `N3 = LPCLG + NUPDATE + NONRLG = 34` — total `SBLG` / `WS` length.
pub const N3LG: usize = LPCLG + NUPDATE + NONRLG;

/// §G.3.14 block 43 — the fixed-point hybrid windowing module of the
/// log-gain predictor.
///
/// Unlike blocks 36/49 (whose speech input arrives in segmented
/// block-floating-point form), block 43's input is the flat `Q9`
/// log-gain history `SBLG(1..=34)`: four new `GTMP` scalars (the
/// offset-removed log-gains of the previous adaptation cycle, oldest
/// first) are shifted in per cycle, the whole buffer is normalized with
/// one `FINDNLS(…, 14)` call (`NLSTMP = NLS − 1` leaves 2 bits of
/// headroom), multiplied by the `Q15` window `WNRLG`, and handed to the
/// shared §G.3.18 HWMCORE with `NLSATTLG = 14` (the 3/4 recursive
/// decay).
#[derive(Debug, Clone)]
pub struct LogGainWindowFixed {
    /// `SBLG(1..=N3)` — flat `Q9` log-gain history, 0-based. Table
    /// 2/G.728 lists `SBLG` all-zero at reset.
    sblg: [i16; N3LG],
    /// §G.3.18 HWMCORE state (`REXPLG` / `NLSREXPLG`, the latter with
    /// its §G.3.14 initial value 31).
    core: HwmcoreState,
}

impl Default for LogGainWindowFixed {
    fn default() -> Self {
        Self::new()
    }
}

impl LogGainWindowFixed {
    /// Fresh state: `SBLG` all zero (Table 2/G.728), `REXPLG` zero with
    /// `NLSREXPLG = 31` (§G.3.14: "NLSREXPLG … is initialized with a
    /// value of 31").
    #[must_use]
    pub fn new() -> Self {
        Self {
            sblg: [0; N3LG],
            core: HwmcoreState::new(LPCLG),
        }
    }

    /// Borrow the `Q9` log-gain history `SBLG` (for tests/audit).
    #[must_use]
    pub fn sblg(&self) -> &[i16; N3LG] {
        &self.sblg
    }

    /// Borrow the recursive-autocorrelation `REXPLG` (for tests/audit).
    #[must_use]
    pub fn rexplg(&self) -> &[i16] {
        self.core.rexp()
    }

    /// `NLSREXPLG` — the shared shift of `REXPLG` (for tests/audit).
    #[must_use]
    pub fn nlsrexplg(&self) -> i32 {
        self.core.nlsrexp()
    }

    /// Run one §G.3.14 block-43 cycle on `gtmp` (the `NUPDATE = 4`
    /// offset-removed log-gains of the previous adaptation cycle in
    /// `Q9`, oldest first: `GTMP(1) = GSTATE(4) … GTMP(4) = GSTATE(1)`)
    /// and return the autocorrelation `R(1..=LPCLG+1)` plus the
    /// `ILLCONDG` verdict.
    pub fn run(&mut self, gtmp: &[i16; NUPDATE]) -> HwmcoreOut {
        // | For N = 1, …, N2: SBLG(N) = SBLG(N + NUPDATE)  | Shift old buffer
        for n in 0..N2LG {
            self.sblg[n] = self.sblg[n + NUPDATE];
        }
        // | For N = 1, …, NUPDATE: SBLG(N2 + N) = GTMP(N)
        // | SBLG(N3) is the newest sample; all SBLG are Q9, 16 bits.
        self.sblg[N2LG..N3LG].copy_from_slice(gtmp);

        // | Call FINDNLS(SBLG, N3, N3, 14, NLS)  | left shifts for the
        // | NLSTMP = NLS − 1                     | next loop, 2 bits headroom
        let vals: [i32; N3LG] = std::array::from_fn(|i| i32::from(self.sblg[i]));
        let nls = findnls(&vals, N3LG, 14);
        let nlstmp = nls - 1;

        // | K = 1
        // | For N = 34, 33, …, 1:
        // |   P = SBLG(N) * WNRLG(K)          | WNRLG is Q15
        // |   If NLSTMP = −1, AA0 = P >> 1
        // |   If NLSTMP > −1, AA0 = P << NLSTMP
        // |   WS(N) = RND(AA0)                | WS(N) is 14 bits or less
        // |   K = K + 1
        let mut ws = [0i16; N3LG];
        for (k, n) in (0..N3LG).rev().enumerate() {
            let p = i32::from(self.sblg[n]) * i32::from(WNRG_Q15[k]);
            let aa0 = if nlstmp == -1 { p >> 1 } else { p << nlstmp };
            ws[n] = rnd(aa0);
        }

        // | NLSATTLG = 14
        // | Call HWMCORE(LPCLG, N1, N3, NLSATTLG, WS, NLSTMP,
        // |              REXPLG, NLSREXPLG, R, ILLCONDG)
        self.core.run(
            LPCLG,
            N1LG,
            N3LG,
            crate::annex_g_hybrid::NLSATT50,
            &ws,
            nlstmp,
        )
    }
}

// ---------------------------------------------------------------------
// §G.3.16 — the per-vector blocks: 46 (log-gain linear prediction),
// 98 (log-gain limiter), 99 (offset adder), 48 (inverse logarithm) and
// 93/94/96/97 (the GSTATE(1) update from the chosen codebook indices).
// ---------------------------------------------------------------------

/// `GSTATE` / `GTMP` initial value — Table G.2/G.728: `−16384` in Q9
/// (= −32 dB, the offset-removed log-gain floor).
pub const GSTATE_INIT_Q9: i16 = -16384;

/// Block 98's log-gain limits in Q9: `−32 dB ≤ δ̂(n) ≤ +28 dB`
/// (§G.3.16: "If LOGGAIN > 14336 … If LOGGAIN < −16384").
pub const LOGGAIN_MAX_Q9: i32 = 14336;
/// See [`LOGGAIN_MAX_Q9`].
pub const LOGGAIN_MIN_Q9: i32 = -16384;

/// `GOFF = 32 dB` in Q9 (Table G.1/G.728: 16384, "Q9") — the log-gain
/// offset added back by block 99 before the inverse logarithm.
pub const GOFF_Q9: i32 = 16384;

/// §G.3.16 block 46 — fixed-point log-gain linear prediction.
///
/// Computes `LOGGAIN = −Σ_{i=1..LPCLG} GP(i+1)·GSTATE(i)` (`GP` Q14,
/// `GSTATE` Q9 → `AA0 >> 14` back to Q9) while shifting the `GSTATE`
/// delay line down one position (the annex folds the shift into the
/// prediction loop). Returns the *unlimited* Q9 log-gain; block 98
/// ([`limit_log_gain_q9`]) is applied next.
pub fn log_gain_predict(gp: &[i16; LPCLG + 1], gstate: &mut [i16; LPCLG]) -> i32 {
    // | AA0 = 0
    // | For I = LPCLG, LPCLG − 1, …, 3, 2:
    // |   P = GP(I + 1) * GSTATE(I); AA0 = AA0 − P
    // |   GSTATE(I) = GSTATE(I − 1)
    let mut aa0: i64 = 0;
    for i in (2..=LPCLG).rev() {
        let p = i32::from(gp[i]) * i32::from(gstate[i - 1]);
        aa0 -= i64::from(p);
        gstate[i - 1] = gstate[i - 2];
    }
    // | P = GP(2) * GSTATE(1); AA0 = AA0 − P
    let p = i32::from(gp[1]) * i32::from(gstate[0]);
    aa0 -= i64::from(p);
    // | AA0 = AA0 >> 14; LOGGAIN = AA0
    // (saturate the 32-bit accumulator the annex's DSP would have used
    // before the shift; block 98 clamps to the Q9 dB range right after).
    let aa0 = aa0.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
    aa0 >> 14
}

/// §G.3.16 block 98 — the log-gain limiter: clip the predicted Q9
/// log-gain to `[−32 dB, +28 dB]` (`[−16384, 14336]`).
#[must_use]
pub fn limit_log_gain_q9(loggain: i32) -> i16 {
    loggain.clamp(LOGGAIN_MIN_Q9, LOGGAIN_MAX_Q9) as i16
}

/// §G.3.16 blocks 99 + 48 — the log-gain offset adder and the inverse
/// logarithm calculator.
///
/// Adds `GOFF = 16384` (32 dB, Q9) to the limited log-gain and converts
/// to the linear domain: `σ = 10^(Z/20) = 2^(0.1660964·Z)` with
/// `0.1660964` split as `10·2⁻⁶ + 20649·2⁻²¹` for precision, the
/// integer part `[X]` becoming the exponent and the fractional part
/// `x ∈ [0, 1)` evaluated by the §G.3.16 Q15 Taylor polynomial
/// `2^x = ((c₄x + c₃)x + c₂)x² + c₁x + c₀`.
///
/// Returns `(GAIN, NLSGAIN)`: the Q14 mantissa of `2^x` and
/// `NLSGAIN = 14 − [X]`, i.e. `σ = GAIN · 2^(−NLSGAIN)` — the scalar
/// floating-point excitation gain blocks 16 / 21 consume.
#[must_use]
pub fn inverse_log_gain(loggain: i16) -> (i16, i32) {
    // Block 99: Z = LOGGAIN + GOFF (Q9, 0 … 30720 = 0 … 60 dB).
    let z = i32::from(loggain) + GOFF_Q9;

    // | AA0 = 10 * Z        | Z is Q9, 10 is Q6, so AA0 is Q15
    // | AA1 = 20649 * Z     | 20649 is Q21, so AA1 is Q30
    // | AA1 = AA1 << 1      | Make AA1 Q31
    // | AA1 = RND(AA1)      | Round AA1 to Q15 in low word
    // | AA0 = AA0 + AA1     | AA0 = [X] + x in Q15
    let aa0 = 10 * z;
    let aa1 = 20649 * z;
    let aa1 = i32::from(rnd(aa1 << 1));
    let aa0 = aa0 + aa1;

    // | AA1 = AA0 >> 15; NLS = AA1   | NLS = [X]
    // | AA1 = AA1 << 15; x = AA0 − AA1
    let nls = aa0 >> 15;
    let x = aa0 - (nls << 15);

    // Taylor coefficients (§G.3.16): c4..c1 in Q15, c0 in Q14.
    const C4: i32 = 323; // 0.0098571
    const C3: i32 = 1874; // 0.0571899
    const C2: i32 = 7866; // 0.2400512
    const C1: i32 = 22702; // 0.6928100
    const C0: i32 = 16384; // 1.0 in Q14

    // Horner stages: each computes AA0 = TMP·x (Q30) << 1 (Q31), adds
    // the next coefficient promoted to the high word (cᵢ << 16, Q31)
    // and rounds back to a Q15 TMP.
    let tmp = i32::from(rnd(((C4 * x) << 1) + (C3 << 16)));
    let tmp = i32::from(rnd(((tmp * x) << 1) + (C2 << 16)));
    let tmp = i32::from(rnd(((tmp * x) << 1) + (C1 << 16)));
    // | AA0 = TMP * x       | Q30 — no left shift this time!
    // | AA1 = c0 << 16      | Q30
    // | GAIN = RND(AA0 + AA1)  | Q14, contains 2^x
    let gain = rnd((tmp * x) + (C0 << 16));

    // | NLSGAIN = 14 − NLS
    (gain, 14 - nls)
}

/// §G.3.16 blocks 93 / 94 / 96 / 97 — the `GSTATE(1)` update.
///
/// Reconstructs the offset-removed log-gain of the *previous* vector
/// from the predicted log-gain (block 95's `LOGGAIN`, post-limiter)
/// plus the two Q11 dB table look-ups for the chosen gain (`ig`) and
/// shape (`is`) indices (both 1-based), aligning at the Q16 boundary
/// (`LOGGAIN << 7`, tables `<< 5`), shifting back to Q9 and flooring at
/// `−16384` (= −32 dB, limiter 97 / eq. G-9).
#[must_use]
pub fn gstate1_update(loggain: i16, ig: usize, is: usize) -> i16 {
    // | AA0 = LOGGAIN << 7                    | Align decimal points
    // | AA0 = AA0 + (GCBLG(IG) << 5)
    // | AA0 = AA0 + (SHAPELG(IS) << 5)
    let mut aa0 = i32::from(loggain) << 7;
    aa0 += i32::from(GCBLG_Q11[ig - 1]) << 5;
    aa0 += i32::from(SHAPELG_Q11[is - 1]) << 5;
    // | AA0 = AA0 >> 7                        | Back to Q9
    aa0 >>= 7;
    // | IF AA0 < −16384, set AA0 = −16384     | Lower limit
    if aa0 < LOGGAIN_MIN_Q9 {
        aa0 = LOGGAIN_MIN_Q9;
    }
    // | GSTATE(1) = AA0                       | Lower 16-bit word saved
    aa0 as i16
}

// ---------------------------------------------------------------------
// The complete Figure G.1 fixed-point backward vector gain adapter.
// ---------------------------------------------------------------------

/// The full fixed-point backward vector gain adapter (Figure G.1/G.728
/// blocks 43 – 48 + 91 – 99) — the fixed-point analogue of the
/// floating-point [`crate::gain_adapter::GainAdapter`], with the §G.2.1
/// index-table reformulation replacing the RMS/logarithm path.
///
/// Drive it per the §G.7 main-program order (`icount` is the 1-based
/// vector counter within the adaptation cycle):
///
/// 1. [`predict`](Self::predict)`(icount)` before the codebook search /
///    excitation decode — applies block 45 at `ICOUNT = 2`, then runs
///    blocks 46/98/99/48 and returns the scalar-floating excitation
///    gain `(GAIN, NLSGAIN)`.
/// 2. [`update`](Self::update)`(icount, ig, is)` after the indices are
///    known — blocks 93/94/96/97 write `GSTATE(1)`, and at `ICOUNT = 1`
///    the per-cycle blocks 43 (hybrid window) + 44 (Levinson) refresh
///    the pending predictor `GPTMP`.
#[derive(Debug, Clone)]
pub struct GainAdapterFixed {
    /// `GSTATE(1..=LPCLG)` — Q9 offset-removed log-gain memory,
    /// Table G.2 initial value −16384.
    gstate: [i16; LPCLG],
    /// `GP(1..=LPCLG+1)` — live Q14 log-gain predictor, Table G.2
    /// initial value (16384, −16384, 0, …, 0).
    gp: [i16; LPCLG + 1],
    /// `GPTMP` — the block-44 output pending block 45's `ICOUNT = 2`
    /// commit, in the Q format signalled by `nlsgptmp`.
    gptmp: [i16; LPCLG + 1],
    /// `NLSGPTMP` — block 44's precision signal (13/14/15).
    nlsgptmp: i32,
    /// `ILLCONDG` — the log-gain ill-conditioning flag (block 43's
    /// zero-`R(LPCLG+1)` verdict or a block-44 failure).
    illcondg: bool,
    /// Block 95 — the previous predicted (and limited) log-gain, Q9.
    loggain: i16,
    /// Block 43 state (`SBLG` / `REXPLG` / `NLSREXPLG`).
    window: LogGainWindowFixed,
}

impl Default for GainAdapterFixed {
    fn default() -> Self {
        Self::new()
    }
}

impl GainAdapterFixed {
    /// Fresh adapter with the Table G.2/G.728 initial state.
    #[must_use]
    pub fn new() -> Self {
        let mut gp = [0i16; LPCLG + 1];
        gp[0] = 16384; // GP(1) = 1.0 Q14
        gp[1] = -16384; // GP(2) = −1.0 Q14 (δ̂(n) = δ(n−1) at reset)
        Self {
            gstate: [GSTATE_INIT_Q9; LPCLG],
            gp,
            gptmp: [0; LPCLG + 1],
            nlsgptmp: 0,
            illcondg: false,
            loggain: GSTATE_INIT_Q9,
            window: LogGainWindowFixed::new(),
        }
    }

    /// Borrow the live Q14 log-gain predictor `GP` (for tests/audit).
    #[must_use]
    pub fn gp(&self) -> &[i16; LPCLG + 1] {
        &self.gp
    }

    /// Borrow the Q9 log-gain memory `GSTATE` (for tests/audit).
    #[must_use]
    pub fn gstate(&self) -> &[i16; LPCLG] {
        &self.gstate
    }

    /// The previous predicted+limited Q9 log-gain (block 95's content).
    #[must_use]
    pub fn loggain(&self) -> i16 {
        self.loggain
    }

    /// `ILLCONDG` — the current log-gain ill-conditioning verdict.
    #[must_use]
    pub fn illcondg(&self) -> bool {
        self.illcondg
    }

    /// Per-vector prediction: block 45 at `ICOUNT = 2` (gated by
    /// `ILLCONDG` / Q14 overflow — keep the old `GP` on decline), then
    /// blocks 46 → 98 → 99 → 48. Returns `(GAIN, NLSGAIN)` — the Q14
    /// mantissa of `2^x` and its shift, i.e. `σ(n) = GAIN·2^(−NLSGAIN)`.
    pub fn predict(&mut self, icount: usize) -> (i16, i32) {
        // OGAINDB/AFTERFE at their inert values make block 98AF the
        // plain block 98, so the non-erasure path is unchanged.
        self.predict_limited(icount, 0, 0)
    }

    /// [`Self::predict`] with the Annex I block **98AF** (§I.5.10)
    /// post-erasure growth limiter in place of block 98: after the
    /// `[−16384, 14336]` range clamp, while `afterfe > 0` the Q9
    /// log-gain may exceed `ogaindb_q9` (the §I.5.1 `OGAINDB`) by at
    /// most `FEGAINMAX = 1024` (+2 dB). Identical to [`Self::predict`]
    /// when `afterfe == 0`.
    pub fn predict_limited(
        &mut self,
        icount: usize,
        ogaindb_q9: i32,
        afterfe: usize,
    ) -> (i16, i32) {
        debug_assert!((1..=NUPDATE).contains(&icount));
        // | If ICOUNT = 2 and ILLCONDG = .FALSE., then do block 45
        if icount == 2 {
            if let Some(gp) = crate::annex_g_synth_adapter::log_gain_bandwidth_expand(
                &self.gptmp,
                self.nlsgptmp,
                self.illcondg,
            ) {
                self.gp = gp;
            }
        }
        // | do blocks 46, 98AF, 99, and 48
        let lg = log_gain_predict(&self.gp, &mut self.gstate);
        self.loggain = crate::annex_i_fixed::limit_log_gain_after_fe_q9(lg, ogaindb_q9, afterfe);
        inverse_log_gain(self.loggain)
    }

    /// Post-search update: blocks 93/94/96/97 write `GSTATE(1)` from
    /// the chosen 1-based gain / shape indices, then — at `ICOUNT = 1`
    /// — the per-cycle blocks 43 + 44 refresh `GPTMP` / `ILLCONDG`
    /// (block 45 commits it at the next `ICOUNT = 2` [`predict`]).
    /// Annex I erased-vector update (§I.4.3 "vital operations"):
    /// `GSTATE(1)` receives the block-97FE log-gain of the extrapolated
    /// excitation instead of the blocks-93/94/96/97 index reconstruction,
    /// and — at `ICOUNT = 1` — block **43FE** runs the hybrid window
    /// (buffer + recursive-component update; the fixed-point 43FE
    /// computes its full `R(1..11)` output anyway, which is simply
    /// discarded) while block 44 is skipped ("If FERROR = .FALSE., do
    /// block 44", §I.5.1). `ILLCONDG` is refreshed from 43FE's verdict
    /// per the §I.5.1 `output ill-condition flag = ILLCONDG` line.
    pub fn update_erased(&mut self, icount: usize, gstate1_97fe_q9: i16) {
        debug_assert!((1..=NUPDATE).contains(&icount));
        // | do block 97FE; GSTATE(1) = output of block 97FE
        self.gstate[0] = gstate1_97fe_q9;

        // | If ICOUNT = 1: GTMP in one shot; do block 43FE (no block 44).
        if icount == 1 {
            let gtmp = [
                self.gstate[3],
                self.gstate[2],
                self.gstate[1],
                self.gstate[0],
            ];
            let core = self.window.run(&gtmp);
            self.illcondg = core.illcond;
        }
    }

    pub fn update(&mut self, icount: usize, ig: usize, is: usize) {
        debug_assert!((1..=NUPDATE).contains(&icount));
        // | do blocks 93, 94, 96, and 97; GSTATE(1) = output of block 97
        self.gstate[0] = gstate1_update(self.loggain, ig, is);

        // | If ICOUNT = 1: GTMP(1) = GSTATE(4) … GTMP(4) = GSTATE(1);
        // | do block 43; do block 44.
        if icount == 1 {
            let gtmp = [
                self.gstate[3],
                self.gstate[2],
                self.gstate[1],
                self.gstate[0],
            ];
            let core = self.window.run(&gtmp);
            let mut gptmp = [0i16; LPCLG + 1];
            let status = crate::annex_g_levinson::levinson_durbin_fixed(
                &crate::annex_g_levinson::LevinsonInput {
                    rtmp: &core.rtmp,
                    illcond: core.illcond,
                    lpc: LPCLG,
                },
                &mut gptmp,
            );
            self.illcondg = status.illcond;
            if !status.illcond {
                self.gptmp = gptmp;
                self.nlsgptmp = status.nlsatmp;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hybrid_window::{HybridWindow, HybridWindowState};
    use crate::tables::wnrg_f64;

    #[test]
    fn adapter_fresh_state_predicts_unity() {
        // First vector: GP = (1, −1, 0…), GSTATE = −16384 ⇒ δ̂ = −32 dB
        // ⇒ σ = 1 exactly.
        let mut a = GainAdapterFixed::new();
        let (gain, nlsgain) = a.predict(1);
        assert_eq!((gain, nlsgain), (16384, 14));
        assert_eq!(a.loggain(), -16384);
    }

    #[test]
    fn adapter_tracks_floating_point_gain_adapter() {
        // Drive the fixed-point Figure G.1 adapter and the floating-point
        // block-20/30 adapter with the SAME (IG, IS) index stream. The
        // float takes the previous gain-scaled excitation ET = σ·gq·y —
        // eq. G-14 proves the index-table path is mathematically
        // equivalent — so feed it ET built from ITS OWN σ, and compare
        // the two σ trajectories in dB. Differences come only from Q9
        // log-gains, the Q11 dB tables, the Q14 predictor and the Q15
        // antilog polynomial.
        use crate::gain_adapter::GainAdapter;
        use crate::tables::{y_f64, GQ};

        let mut fixed = GainAdapterFixed::new();
        let mut float = GainAdapter::new();
        let y = y_f64();
        let gc_db = crate::annex_g_gain::gain_log_db();
        let sh_db = crate::annex_g_gain::shape_log_db();

        let mut et_prev = [0.0f64; crate::consts::IDIM];
        let mut worst_db = 0.0f64;
        let mut sum_db = 0.0f64;
        let mut count = 0usize;

        for n in 0..400 {
            let icount = (n % NUPDATE) + 1;

            // Fixed side: predict σ, float side: ingest previous ET and
            // predict σ (predict_next both updates and predicts).
            let (gain, nlsgain) = fixed.predict(icount);
            let sigma_x = f64::from(gain) * (2.0f64).powi(-nlsgain);
            let sigma_f = float.predict_next(&et_prev);

            // Choose the next indices the way a real encoder does: track
            // a smooth log-gain target (a slow ±10 dB swell) by picking
            // the (IG, IS) whose dB contribution lands δ(n) closest to
            // it. Uncorrelated random indices would instead exercise the
            // backward feedback loop as a chaos amplifier, where the two
            // arithmetics legitimately diverge.
            let delta_hat_f = 20.0 * sigma_f.log10() - 32.0;
            let target_db = 10.0 * ((n as f64) * 0.045).sin() - 6.0;
            let mut best = (1usize, 1usize, f64::INFINITY);
            for ig in 1..=4usize {
                for is in 1..=128usize {
                    let d = (target_db - (delta_hat_f + gc_db[ig - 1] + sh_db[is - 1])).abs();
                    if d < best.2 {
                        best = (ig, is, d);
                    }
                }
            }
            let (ig, is, _) = best;

            fixed.update(icount, ig, is);
            for k in 0..crate::consts::IDIM {
                et_prev[k] = sigma_f * GQ[ig - 1] * y[is - 1][k];
            }

            if n >= 40 {
                let db = 20.0 * (sigma_x / sigma_f).log10();
                worst_db = worst_db.max(db.abs());
                sum_db += db.abs();
                count += 1;
            }
        }
        let mean_db = sum_db / count as f64;
        assert!(
            worst_db < 1.0,
            "worst σ divergence {worst_db} dB (mean {mean_db})"
        );
        assert!(mean_db < 0.2, "mean σ divergence {mean_db} dB");
    }

    #[test]
    fn inverse_log_gain_tracks_exact_antilog() {
        // Block 48 vs the exact 10^(Z/20) over the whole legal Q9 range
        // (Z = LOGGAIN + 32 dB ∈ [0, 60] dB). The Q15 Taylor polynomial
        // is accurate to a few 1e-4 relative; assert 1e-3.
        let mut lg = LOGGAIN_MIN_Q9;
        while lg <= LOGGAIN_MAX_Q9 {
            let (gain, nlsgain) = inverse_log_gain(lg as i16);
            let sigma = f64::from(gain) * (2.0f64).powi(-nlsgain);
            let z_db = (f64::from(lg) + f64::from(GOFF_Q9)) / 512.0;
            let want = 10.0f64.powf(z_db / 20.0);
            let rel = (sigma - want).abs() / want;
            assert!(
                rel < 1e-3,
                "LOGGAIN {lg}: σ fixed {sigma} vs exact {want} (rel {rel})"
            );
            lg += 37; // dense-but-fast sweep
        }
    }

    #[test]
    fn inverse_log_gain_unity_at_floor() {
        // At the −32 dB floor, Z = 0 dB ⇒ σ = 1 exactly: [X] = 0,
        // x = 0, the polynomial collapses to c0 ⇒ GAIN = 16384 Q14.
        let (gain, nlsgain) = inverse_log_gain(LOGGAIN_MIN_Q9 as i16);
        assert_eq!(gain, 16384);
        assert_eq!(nlsgain, 14);
    }

    #[test]
    fn predict_on_fresh_state_returns_floor() {
        // Table G.2 initial state: GP = (1, −1, 0, …) Q14 and
        // GSTATE = −16384 Q9 everywhere ⇒ the prediction is
        // −GP(2)·GSTATE(1) >> 14 = −16384 (−32 dB) ⇒ unity σ.
        let mut gp = [0i16; LPCLG + 1];
        gp[0] = 16384;
        gp[1] = -16384;
        let mut gstate = [GSTATE_INIT_Q9; LPCLG];
        let lg = log_gain_predict(&gp, &mut gstate);
        assert_eq!(lg, -16384);
        let lg = limit_log_gain_q9(lg);
        let (gain, nlsgain) = inverse_log_gain(lg);
        assert_eq!((gain, nlsgain), (16384, 14));
    }

    #[test]
    fn predict_shifts_gstate_delay_line() {
        let mut gp = [0i16; LPCLG + 1];
        gp[0] = 16384;
        let mut gstate: [i16; LPCLG] = std::array::from_fn(|i| (i as i16 + 1) * 100);
        let _ = log_gain_predict(&gp, &mut gstate);
        // GSTATE(i) = GSTATE(i−1): everything moves one down, GSTATE(1)
        // untouched (the driver overwrites it via blocks 96/97).
        let want: [i16; LPCLG] = std::array::from_fn(|i| if i == 0 { 100 } else { i as i16 * 100 });
        assert_eq!(gstate, want);
    }

    #[test]
    fn limiter_clips_both_ends() {
        assert_eq!(limit_log_gain_q9(20000), LOGGAIN_MAX_Q9 as i16);
        assert_eq!(limit_log_gain_q9(-20000), LOGGAIN_MIN_Q9 as i16);
        assert_eq!(limit_log_gain_q9(1234), 1234);
    }

    #[test]
    fn gstate1_update_matches_float_formula() {
        // Blocks 93/94/96/97 vs the float eq. G-14 on the Q11 tables:
        // δ(n−1) = δ̂(n−1) + GCBLG + SHAPELG, floored at −32 dB. The Q16
        // alignment + >> 7 truncation loses at most 1 Q9 quantum.
        for &(lg, ig, is) in &[
            (0i16, 1usize, 1usize),
            (5120, 4, 21),   // large positive contributions
            (-8000, 1, 96),  // strong negative shape term
            (-16384, 1, 96), // floor engaged
            (14336, 8, 21),  // sign-mirrored gain level
        ] {
            let got = gstate1_update(lg, ig, is);
            let f = f64::from(lg) / 512.0
                + f64::from(crate::tables::GCBLG_Q11[ig - 1]) / 2048.0
                + f64::from(crate::tables::SHAPELG_Q11[is - 1]) / 2048.0;
            let want = (f.max(-32.0) * 512.0).floor();
            assert!(
                (f64::from(got) - want).abs() <= 1.0,
                "gstate1({lg}, {ig}, {is}) = {got} vs float {want}"
            );
        }
    }

    #[test]
    fn dimensions_match_spec() {
        // §G.3.14: N1 = 14, N2 = 30, N3 = 34.
        assert_eq!(N1LG, 14);
        assert_eq!(N2LG, 30);
        assert_eq!(N3LG, 34);
        assert_eq!(WNRG_Q15.len(), N3LG);
    }

    #[test]
    fn fresh_state_matches_table_g2() {
        let w = LogGainWindowFixed::new();
        assert!(w.sblg().iter().all(|&v| v == 0), "SBLG all zero at reset");
        assert!(w.rexplg().iter().all(|&v| v == 0), "REXPLG zero at reset");
        assert_eq!(w.nlsrexplg(), 31, "NLSREXPLG initial value = 31");
        assert_eq!(w.rexplg().len(), LPCLG + 1);
    }

    #[test]
    fn newest_gain_lands_at_buffer_tail() {
        // §G.3.14: SBLG(N3) is the newest sample; GTMP(4) = GSTATE(1)
        // is the newest log-gain.
        let mut w = LogGainWindowFixed::new();
        let gtmp: [i16; NUPDATE] = [-100, -200, -300, -400];
        w.run(&gtmp);
        assert_eq!(&w.sblg()[N2LG..], &gtmp[..]);
        // A second cycle shifts them 4 places down.
        let gtmp2: [i16; NUPDATE] = [11, 22, 33, 44];
        w.run(&gtmp2);
        assert_eq!(&w.sblg()[N2LG - NUPDATE..N2LG], &gtmp[..]);
        assert_eq!(&w.sblg()[N2LG..], &gtmp2[..]);
    }

    #[test]
    fn zero_history_flags_illcondg() {
        // An all-zero log-gain history windows to zero ⇒ R(LPCLG+1)'s
        // 32-bit accumulator is zero ⇒ ILLCONDG.
        let mut w = LogGainWindowFixed::new();
        let out = w.run(&[0; NUPDATE]);
        assert!(out.illcond);
        assert!(out.rtmp.iter().all(|&v| v == 0));
    }

    #[test]
    fn autocorrelation_tracks_floating_point_window() {
        // Drive the fixed-point block 43 and the floating-point hybrid
        // window (same LPCLG/NUPDATE/NONRLG dimensions, same requantized
        // Q9 log-gains) side by side; the normalized autocorrelation
        // sequences must agree closely once the recursive history fills.
        let fwin = wnrg_f64();
        let fhw = HybridWindow {
            m: LPCLG,
            l: NUPDATE,
            n: NONRLG,
            window: &fwin,
            decay: 0.75,
        };
        let mut fstate = HybridWindowState::new(&fhw);
        let mut xstate = LogGainWindowFixed::new();

        // A slowly-meandering log-gain trajectory in dB (Q9), similar in
        // scale to real GSTATE contents (−32 … +28 dB).
        let mut seed = 0x1234_5678_9abc_def0u64;
        let mut compared = false;
        for cyc in 0..24 {
            let mut gtmp_q9 = [0i16; NUPDATE];
            let mut gtmp_f = [0.0f64; NUPDATE];
            for k in 0..NUPDATE {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                let u = ((seed >> 33) as f64) / ((1u64 << 31) as f64) - 1.0;
                let db = u * 20.0 + 4.0 * ((cyc as f64) * 0.7).sin();
                let q9 = (db * 512.0).round().clamp(-16384.0, 14336.0);
                gtmp_q9[k] = q9 as i16;
                gtmp_f[k] = q9 / 512.0; // the value the Q9 mantissa represents
            }

            let mut frtmp = vec![0.0f64; LPCLG + 1];
            fstate.run(&fhw, &gtmp_f, &mut frtmp);
            let core = xstate.run(&gtmp_q9);

            if cyc >= 20 && frtmp[0].abs() > 1e-9 && core.rtmp[0] != 0 {
                compared = true;
                for k in 0..=LPCLG {
                    let f = frtmp[k] / frtmp[0];
                    let x = f64::from(core.rtmp[k]) / f64::from(core.rtmp[0]);
                    assert!(
                        (f - x).abs() < 0.05,
                        "normalized R[{k}]: fixed {x} vs float {f}"
                    );
                }
                assert!(!core.illcond, "well-conditioned input clears ILLCONDG");
            }
        }
        assert!(compared, "comparison cycles must have run");
    }
}
