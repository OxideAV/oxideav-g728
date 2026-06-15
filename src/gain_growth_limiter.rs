//! Annex I §I.4.5 — limitation of gain growth after a frame erasure.
//!
//! After a long frame erasure (or one that straddles a low-amplitude →
//! high-amplitude transition), the log-gain buffer `SBG()` inside the
//! hybrid window holds very low gains at the end of the erasure. When the
//! first good frames arrive with much higher gains, the backward LPC
//! analysis produces a gain predictor that catches up *too* fast,
//! overshoots, and produces an audible "pop" (§I.4.5):
//!
//! > the predicted (i.e. backward adapted) gain therefore grew too fast
//! > and too much, causing a gain overshoot and the resulting "pop".
//!
//! Annex I fixes this by clamping the *growth* of the backward-adapted
//! log-gain for the first few good frames after an erasure. This module
//! realises **Block 47AF** (§I.5.6, the floating-point replacement for
//! the normal Block 47 log-gain limiter) together with the **AFTERFE /
//! FECOUNT** bookkeeping from the §I.5.1 decoder main loop that decides
//! *how long* the clamp stays active.
//!
//! ## Block 47AF (floating-point, §I.5.6)
//!
//! ```text
//! If GAIN < 0., set GAIN = 0.            | corresponds to linear gain 1
//! If GAIN > 60., set GAIN = 60.          | corresponds to linear gain 1000
//! TMP = GAIN - OGAINDB
//! If AFTERFE > 0 and TMP > FEGAINMAX, do GAIN = OGAINDB + FEGAINMAX
//! ```
//!
//! i.e. the normal block-47 clamp to the dB range `[0, 60]` (already done
//! for non-erased operation in [`crate::gain_adapter`]) followed by the
//! Annex I addition: while the clamp is active (`AFTERFE > 0`), the new
//! log-gain is not allowed to exceed the last log-gain `OGAINDB` by more
//! than [`FEGAINMAX`] dB. `OGAINDB` is the dB value of the gain computed
//! for the previous vector ("old prediction gain in dB", Table I.2).
//!
//! ## AFTERFE bookkeeping (§I.5.1)
//!
//! The clamp duration is driven once per 2.5 ms adaptation cycle, at the
//! top of the cycle (`ICOUNT = 4` in the spec's main loop):
//!
//! ```text
//! If FERROR = .FALSE. and AFTERFE > 0, do AFTERFE = AFTERFE - 1
//! If <start of a new frame> :
//!     read FERROR
//!     If FERROR = .FALSE. and FECOUNT > 0 :          | first good frame after FE
//!         AFTERFE = AFTERFE + FECOUNT
//!         If AFTERFE > AFTERFEMAX, AFTERFE = AFTERFEMAX
//!         FECOUNT = 0
//!     If FERROR = .TRUE. :                            | still erased
//!         FECOUNT = FECOUNT + 1
//! ```
//!
//! `FECOUNT` counts the consecutively-erased adaptation cycles; at the
//! first good frame after the erasure, `AFTERFE` is loaded with that
//! length (so the clamp lasts as long as the erasure) but never beyond
//! [`AFTERFEMAX`] cycles (= 40 ms). Each subsequent good cycle decrements
//! `AFTERFE` until it reaches zero, after which the gain growth is no
//! longer limited.
//!
//! Clean-room: every value and rule here is transcribed from the §I.4.5
//! prose, the §I.5.6 Block 47AF floating-point pseudo-code, the §I.5.1
//! main-loop pseudo-code, and Table I.1/I.2 of
//! `T-REC-G.728-199905-AnnI.pdf`. No reference C / external
//! implementation was consulted.

use crate::consts::{AFTERFEMAX, FEGAINMAX, OGAINDB_INIT};

/// Block 47AF — the Annex I §I.5.6 log-gain limiter (floating-point).
///
/// `gain` is the log-gain in dB just produced by the backward gain
/// predictor (Block 46 + the offset adder), *before* the inverse
/// logarithm (Block 48). `ogaindb` is the dB log-gain of the previous
/// vector (Table I.2's `OGAINDB`). `afterfe` is the §I.5.1 gain-clamping
/// counter: when it is non-zero the post-erasure growth clamp is active.
///
/// Returns the limited log-gain in dB. This is the value that is both fed
/// to the inverse-logarithm block *and* saved as the next `OGAINDB`.
///
/// When `afterfe == 0` this reduces exactly to the normal Block 47 clamp
/// to `[0, 60]` dB, so it is safe to route every vector through this
/// function regardless of erasure state.
#[must_use]
pub fn limit_log_gain(gain: f64, ogaindb: f64, afterfe: usize) -> f64 {
    // Normal block-47 range clamp: GAIN ∈ [0, 60] dB.
    //   If GAIN < 0., set GAIN = 0.   | linear gain 1
    //   If GAIN > 60., set GAIN = 60. | linear gain 1000
    let mut gain = gain.clamp(0.0, 60.0);

    // Annex I addition: while the post-erasure clamp is active, limit the
    // growth relative to the previous log-gain to FEGAINMAX dB/vector.
    //   TMP = GAIN - OGAINDB
    //   If AFTERFE > 0 and TMP > FEGAINMAX, do GAIN = OGAINDB + FEGAINMAX
    if afterfe > 0 {
        let tmp = gain - ogaindb;
        if tmp > FEGAINMAX {
            gain = ogaindb + FEGAINMAX;
        }
    }
    gain
}

/// Drives the §I.5.1 `AFTERFE` / `FECOUNT` / `OGAINDB` bookkeeping that
/// the Annex I decoder uses to scope the [`limit_log_gain`] clamp.
///
/// One [`Self::on_cycle`] call is made per 2.5 ms adaptation cycle, at
/// the top of the cycle, carrying the frame-erasure flag for the frame
/// that cycle belongs to. The decoder reads [`Self::afterfe`] (and the
/// running [`Self::ogaindb`]) when it limits each vector's log-gain, and
/// reports each limited log-gain back via [`Self::set_ogaindb`].
#[derive(Debug, Clone)]
pub struct GainGrowthLimiter {
    /// §I.5.1 `AFTERFE`: number of 2.5 ms adaptation cycles after an
    /// erasure during which the gain growth is still clamped. Counts
    /// down; `0` means the clamp is inactive.
    afterfe: usize,
    /// §I.5.1 `FECOUNT`: number of consecutively erased adaptation cycles
    /// in the current erasure (0 outside an erasure).
    fecount: usize,
    /// §I.5.1 `OGAINDB`: the last (limited) predicted log-gain, in dB.
    /// Initialised to [`OGAINDB_INIT`] (= -32 dB, the log-gain floor).
    ogaindb: f64,
}

impl Default for GainGrowthLimiter {
    fn default() -> Self {
        Self {
            afterfe: 0,
            fecount: 0,
            ogaindb: OGAINDB_INIT,
        }
    }
}

impl GainGrowthLimiter {
    /// Fresh limiter in the non-erased state (`AFTERFE = FECOUNT = 0`,
    /// `OGAINDB = -32`).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Advance one **adaptation cycle**, given whether the frame that
    /// cycle belongs to is erased (`erased == FERROR`).
    ///
    /// This realises the §I.5.1 main-loop bookkeeping that runs at the
    /// top of each adaptation cycle:
    ///
    /// * a good cycle decrements an active `AFTERFE`;
    /// * the first good cycle after an erasure loads `AFTERFE` with the
    ///   erasure length `FECOUNT` (saturated at [`AFTERFEMAX`]) and
    ///   resets `FECOUNT`;
    /// * an erased cycle increments `FECOUNT`.
    ///
    /// Returns the updated [`Self::afterfe`].
    pub fn on_cycle(&mut self, erased: bool) -> usize {
        // If FERROR = .FALSE. and AFTERFE > 0, do AFTERFE = AFTERFE - 1
        //   | the last adaptation cycle was not erased, and the gain was
        //   | clamped: decrease the number of cycles left to clamp.
        if !erased && self.afterfe > 0 {
            self.afterfe -= 1;
        }

        if erased {
            // If FERROR = .TRUE., do FECOUNT = FECOUNT + 1 | length of FE.
            self.fecount += 1;
        } else if self.fecount > 0 {
            // First good frame after an erasure: set AFTERFE so the gain
            // is clamped for the next few good frames.
            //   AFTERFE = AFTERFE + FECOUNT
            //   If AFTERFE > AFTERFEMAX, AFTERFE = AFTERFEMAX
            //   FECOUNT = 0
            self.afterfe += self.fecount;
            if self.afterfe > AFTERFEMAX {
                self.afterfe = AFTERFEMAX;
            }
            self.fecount = 0;
        }

        self.afterfe
    }

    /// Limit one vector's log-gain (Block 47AF) against the current
    /// `OGAINDB` and `AFTERFE`, **and** save the result as the next
    /// `OGAINDB`. Convenience wrapper around [`limit_log_gain`] for the
    /// common decoder path.
    pub fn limit(&mut self, gain: f64) -> f64 {
        let limited = limit_log_gain(gain, self.ogaindb, self.afterfe);
        self.ogaindb = limited;
        limited
    }

    /// Current §I.5.1 `AFTERFE` (gain-clamp cycle counter).
    #[must_use]
    pub fn afterfe(&self) -> usize {
        self.afterfe
    }

    /// Current §I.5.1 `FECOUNT` (consecutively erased cycles).
    #[must_use]
    pub fn fecount(&self) -> usize {
        self.fecount
    }

    /// Current §I.5.1 `OGAINDB` (last limited log-gain, dB).
    #[must_use]
    pub fn ogaindb(&self) -> f64 {
        self.ogaindb
    }

    /// Override the stored `OGAINDB` (Table I.2's "old prediction gain in
    /// dB"). The §I.5.1 main loop sets `OGAINDB` to the output of Block
    /// 47AF for good frames and to the output of Block 97FE (the log-gain
    /// of the extrapolated excitation) for erased frames; the decoder
    /// drives that via this setter when it does not route the value
    /// through [`Self::limit`].
    pub fn set_ogaindb(&mut self, ogaindb: f64) {
        self.ogaindb = ogaindb;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_clamp_when_afterfe_zero_is_plain_range_clamp() {
        // With AFTERFE = 0 the function is exactly the normal block-47
        // [0, 60] dB clamp: OGAINDB is irrelevant.
        assert_eq!(limit_log_gain(-5.0, 100.0, 0), 0.0);
        assert_eq!(limit_log_gain(75.0, -100.0, 0), 60.0);
        assert_eq!(limit_log_gain(12.5, 0.0, 0), 12.5);
    }

    #[test]
    fn range_clamp_runs_before_growth_clamp() {
        // §I.5.6 applies the [0, 60] clamp first, then the growth clamp.
        // A gain of -3 with an active clamp clamps to 0 first; 0 - OGAINDB
        // can still trip the growth clamp if OGAINDB is below 0 - FEGAINMAX.
        let out = limit_log_gain(-3.0, -10.0, 4);
        // GAIN -> 0; TMP = 0 - (-10) = 10 > 2 -> GAIN = -10 + 2 = -8.
        assert!((out - (-8.0)).abs() < 1e-12);
    }

    #[test]
    fn growth_clamp_limits_to_ogaindb_plus_fegainmax() {
        // §I.4.5: at most +2 dB/vector over the last log-gain while the
        // clamp is active. OGAINDB = 10, GAIN = 20 -> TMP = 10 > 2 ->
        // GAIN = 10 + 2 = 12.
        let out = limit_log_gain(20.0, 10.0, 1);
        assert!((out - 12.0).abs() < 1e-12);
    }

    #[test]
    fn growth_clamp_passes_through_when_within_budget() {
        // Growth within FEGAINMAX is untouched even when the clamp is
        // active: OGAINDB = 10, GAIN = 11 -> TMP = 1 <= 2 -> GAIN = 11.
        let out = limit_log_gain(11.0, 10.0, 1);
        assert!((out - 11.0).abs() < 1e-12);
    }

    #[test]
    fn growth_clamp_never_limits_a_gain_decrease() {
        // A drop in gain has TMP < 0 < FEGAINMAX, so it is never clamped
        // (the limiter only caps growth, not decay).
        let out = limit_log_gain(3.0, 10.0, 8);
        assert!((out - 3.0).abs() < 1e-12);
    }

    #[test]
    fn tracker_starts_non_erased() {
        let t = GainGrowthLimiter::new();
        assert_eq!(t.afterfe(), 0);
        assert_eq!(t.fecount(), 0);
        assert!((t.ogaindb() - OGAINDB_INIT).abs() < 1e-15);
    }

    #[test]
    fn tracker_counts_erased_cycles_in_fecount() {
        // Each erased cycle bumps FECOUNT and leaves AFTERFE at 0 (the
        // clamp activates only once the erasure ends).
        let mut t = GainGrowthLimiter::new();
        for n in 1..=5 {
            assert_eq!(t.on_cycle(true), 0, "AFTERFE stays 0 during erasure");
            assert_eq!(t.fecount(), n);
        }
    }

    #[test]
    fn tracker_loads_afterfe_at_first_good_cycle() {
        // §I.5.1: at the first good frame after an erasure of FECOUNT
        // cycles, AFTERFE = FECOUNT (clamp lasts as long as the erasure).
        let mut t = GainGrowthLimiter::new();
        for _ in 0..3 {
            t.on_cycle(true);
        }
        assert_eq!(t.fecount(), 3);
        // First good cycle: AFTERFE = 0 + 3 = 3, FECOUNT reset.
        assert_eq!(t.on_cycle(false), 3);
        assert_eq!(t.fecount(), 0);
    }

    #[test]
    fn tracker_afterfe_decrements_each_good_cycle() {
        // After loading AFTERFE, each subsequent good cycle decrements it
        // until it reaches 0 and the clamp deactivates.
        let mut t = GainGrowthLimiter::new();
        for _ in 0..3 {
            t.on_cycle(true);
        }
        assert_eq!(t.on_cycle(false), 3); // first good: load 3
        assert_eq!(t.on_cycle(false), 2);
        assert_eq!(t.on_cycle(false), 1);
        assert_eq!(t.on_cycle(false), 0);
        assert_eq!(t.on_cycle(false), 0, "stays at 0");
    }

    #[test]
    fn tracker_afterfe_saturates_at_afterfemax() {
        // §I.4.5: "If the frame erasure is longer than 40 ms ... the gain
        // clamping is limited to the first 40 ms" => AFTERFE <= AFTERFEMAX.
        let mut t = GainGrowthLimiter::new();
        for _ in 0..(AFTERFEMAX + 10) {
            t.on_cycle(true);
        }
        assert_eq!(t.fecount(), AFTERFEMAX + 10);
        assert_eq!(t.on_cycle(false), AFTERFEMAX, "AFTERFE saturates");
    }

    #[test]
    fn tracker_limit_updates_ogaindb() {
        // limit() applies Block 47AF and stores the result as OGAINDB so
        // the next vector's growth is measured from this one.
        let mut t = GainGrowthLimiter::new();
        for _ in 0..2 {
            t.on_cycle(true);
        }
        t.on_cycle(false); // AFTERFE = 2, OGAINDB still -32
                           // First good vector: GAIN = 0 (silence floor). TMP = 0 - (-32) =
                           // 32 > 2 -> clamped to -32 + 2 = -30; OGAINDB becomes -30.
        let g1 = t.limit(0.0);
        assert!((g1 - (-30.0)).abs() < 1e-12);
        assert!((t.ogaindb() - (-30.0)).abs() < 1e-12);
        // Next vector still within the clamp window: +2 dB/vector growth.
        let g2 = t.limit(0.0);
        assert!((g2 - (-28.0)).abs() < 1e-12);
    }

    #[test]
    fn fegainmax_and_afterfemax_match_table_i1() {
        // Guard the Table I.1 literals against a future typo.
        assert!((FEGAINMAX - 2.0).abs() < 1e-15);
        assert_eq!(AFTERFEMAX, 16);
        // 16 cycles × 2.5 ms = 40 ms, the spec's stated clamp ceiling.
        assert_eq!(AFTERFEMAX as f64 * 2.5, 40.0);
    }
}
