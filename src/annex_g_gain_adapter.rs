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
use crate::tables::WNRG_Q15;

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
        self.core.run(LPCLG, N1LG, N3LG, &ws, nlstmp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hybrid_window::{HybridWindow, HybridWindowState};
    use crate::tables::wnrg_f64;

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
