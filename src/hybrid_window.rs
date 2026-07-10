//! Hybrid windowing — generic transcription of block 49 (synthesis
//! filter), block 36 (perceptual weighting filter) and block 43
//! (log-gain predictor).
//!
//! The three blocks share identical pseudocode shape; only the order
//! (`m`), the per-call input size (`l`), the non-recursive tail length
//! (`n`), the window data and the recursive decay factor differ. This
//! module exposes the structure as a single [`HybridWindow`] parameter
//! object plus a [`HybridWindowState`] carrier with [`HybridWindowState::run`]
//! transcribing the spec pseudocode line-for-line.
//!
//! ## Spec mapping
//!
//! Recommendation §5.6 block 49 ("very similar to block 36"):
//!
//! ```text
//! N1 = LPC + NFRSZ
//! N2 = LPC + NONR
//! N3 = LPC + NFRSZ + NONR
//! For N = 1,...,N2: SB(N) = SB(N + NFRSZ)             | shift old buffer
//! For N = 1,...,NFRSZ: SB(N2 + N) = STTMP(N)          | shift in new
//! K = 1
//! For N = N3, N3-1, ..., 1: WS(N) = SB(N) * WNR(K); K = K + 1
//! For I = 1,...,LPC+1:
//!     TMP = 0
//!     For N = LPC+1, ..., N1: TMP = TMP + WS(N)*WS(N+1-I)
//!     REXP(I) = (3/4) * REXP(I) + TMP                 | recursive part
//! For I = 1,...,LPC+1:
//!     RTMP(I) = REXP(I)
//!     For N = N1+1, ..., N3: RTMP(I) = RTMP(I) + WS(N)*WS(N+1-I)
//! RTMP(1) = RTMP(1) * WNCF                            | white-noise corr
//! ```
//!
//! Block 36 (weighting filter) is the same structure with
//! `(LPC, NFRSZ, NONR, WNR, REXP)` swapped for
//! `(LPCW, NFRSZ, NONRW, WNRW, REXPW)`. Block 43 (log-gain) is the
//! same structure with `(LPCLG, NUPDATE, NONRLG, WNRLG, REXPLG)` — note
//! the per-call input length drops from 20 (NFRSZ) to 4 (NUPDATE)
//! because the log-gain predictor only sees four offset-removed
//! log-gain values per adaptation cycle, not 20 speech samples (§5.6
//! block 43 explanatory prose: "only four (rather than 20) gain samples
//! are fed to this block each time").
//!
//! ## Recursive decay
//!
//! The recursive component is decayed by `3/4` per cycle in all three
//! blocks (the prose constant `(3/4)^(1/40)` in §3.7 is the per-sample
//! decay; the per-cycle factor is `(3/4)^(N/40)` for an N-sample
//! shift). The spec's pseudo-code hard-codes the per-cycle factor as
//! `3/4` for blocks 36/43/49, which is what we transcribe here.

use crate::consts::WNCF;

/// Pure-data description of a hybrid window's parameters and table.
///
/// All three G.728 hybrid-window blocks share this shape: a window
/// table sized `n3 = m + l + n`, an LPC order `m`, an input length `l`
/// (samples per call), a non-recursive tail length `n` and a recursive
/// decay factor: `3/4` for the synthesis-filter (block 49) and log-gain
/// (block 43) windows, but `1/2` for the perceptual-weighting-filter
/// window (block 36) — the base Recommendation's block-36 pseudo-code
/// spells `REXPW(I) = (1/2)·REXPW(I) + TMP`.
#[derive(Debug, Clone, Copy)]
pub struct HybridWindow<'a> {
    /// LPC order — `m` in the spec (LPC = 50, LPCW = 10, LPCLG = 10).
    pub m: usize,
    /// Per-call input length — `l` in the spec (NFRSZ = 20 for blocks
    /// 36/49, NUPDATE = 4 for block 43).
    pub l: usize,
    /// Non-recursive tail length — `n` in the spec (NONR = 35,
    /// NONRW = 30, NONRLG = 20).
    pub n: usize,
    /// Per-sample window weights, length `m + l + n` (`n3` in the
    /// spec). Element 0 corresponds to spec `WNR(1)` (the rightmost /
    /// newest sample's weight in the right-to-left walk).
    pub window: &'a [f64],
    /// Recursive-component decay: `0.75` for blocks 43/49, `0.5` for
    /// block 36.
    pub decay: f64,
}

impl<'a> HybridWindow<'a> {
    /// Total window length `n3 = m + l + n`.
    pub const fn n3(&self) -> usize {
        self.m + self.l + self.n
    }

    /// Length of the "kept old samples" portion `n2 = m + n`.
    pub const fn n2(&self) -> usize {
        self.m + self.n
    }

    /// Length of the "kept history before this call's autocorrelation
    /// tap" portion `n1 = m + l`.
    pub const fn n1(&self) -> usize {
        self.m + self.l
    }
}

/// Mutable state of one hybrid-window adaptation chain (one of:
/// block 36, block 43, block 49). Owns the signal-history buffer
/// `SB`/`SBW`/`SBLG` of length `n3` and the recursive-autocorrelation
/// buffer `REXP`/`REXPW`/`REXPLG` of length `m + 1`.
#[derive(Debug, Clone)]
pub struct HybridWindowState {
    /// Order — duplicated so methods can size temporaries without
    /// borrowing the parameter object.
    m: usize,
    /// Per-call input length — duplicated for the same reason.
    l: usize,
    /// Non-recursive tail length — duplicated for the same reason.
    n: usize,
    /// Signal-history buffer `SB(1..n3)` in spec form, 0-based here.
    /// Initialised to zero per Table 2/G.728 (`SB`, `SBW`, `SBLG`).
    sb: Vec<f64>,
    /// Recursive-autocorrelation buffer `REXP(1..m+1)` in spec form,
    /// 0-based here. Initialised to zero per Table 2/G.728.
    rexp: Vec<f64>,
}

impl HybridWindowState {
    /// Fresh state with all buffers zeroed (the spec's Table 2 initial
    /// values for `SB`, `SBW`, `SBLG` and `REXP`, `REXPW`, `REXPLG`).
    pub fn new(window: &HybridWindow<'_>) -> Self {
        Self {
            m: window.m,
            l: window.l,
            n: window.n,
            sb: vec![0.0; window.n3()],
            rexp: vec![0.0; window.m + 1],
        }
    }

    /// Borrow the signal-history buffer (for tests).
    pub fn sb(&self) -> &[f64] {
        &self.sb
    }

    /// Borrow the recursive-autocorrelation buffer (for tests).
    pub fn rexp(&self) -> &[f64] {
        &self.rexp
    }

    /// Run one cycle of the hybrid window with `input.len() == l` new
    /// samples and write the `m + 1` autocorrelation coefficients into
    /// `rtmp_out`.
    ///
    /// `rtmp_out` must have length at least `m + 1`. The spec's
    /// `WNCF = 257/256` white-noise correction is applied to
    /// `rtmp_out[0]` (the spec's `R(1)`) for **all three** hybrid
    /// windows — block 36, block 43 and block 49 all spell the line
    /// `R(1) = R(1) * WNCF`.
    pub fn run(&mut self, window: &HybridWindow<'_>, input: &[f64], rtmp_out: &mut [f64]) {
        self.run_taps(window, input, rtmp_out, self.m + 1);
    }

    /// Annex I frame-erasure variant, block **43FE** shape (§I.5.5):
    /// shift the signal buffer, shift in the new input and update the
    /// recursive autocorrelation component `REXP*` — but skip the
    /// non-recursive output accumulation entirely ("If FERROR = .TRUE.,
    /// skip all the lines below", including the `WNCF` correction).
    ///
    /// The buffer and `REXP*` states after this call are identical to
    /// what [`Self::run`] would have left, so the adaptation chain
    /// resynchronises seamlessly at the next good frame (§I.4.3: the
    /// "vital operations" of the backward adapters continue during
    /// erased frames).
    pub fn advance_erased(&mut self, window: &HybridWindow<'_>, input: &[f64]) {
        self.run_taps(window, input, &mut [], 0);
    }

    /// Annex I frame-erasure variant, block **49FE** shape (§I.5.7):
    /// full buffer shift and recursive `REXP` update, but only the first
    /// `taps` non-recursive output coefficients are produced (the
    /// synthesis-filter window computes `RTMP(1..11)` during an erasure
    /// — enough for the post-filter's order-10 by-product — instead of
    /// `RTMP(1..51)`). The `WNCF` white-noise correction still applies
    /// to `R(1)` (the §I.5.7 float listing runs it unconditionally).
    pub fn run_erased(
        &mut self,
        window: &HybridWindow<'_>,
        input: &[f64],
        rtmp_out: &mut [f64],
        taps: usize,
    ) {
        assert!(taps >= 1 && taps <= self.m + 1, "taps out of range");
        self.run_taps(window, input, rtmp_out, taps);
    }

    /// Shared core: shift + window + recursive update over all `m + 1`
    /// lags, then produce the first `taps` output coefficients
    /// (`taps == 0` ⇒ pure state advance, no output, no WNCF).
    fn run_taps(
        &mut self,
        window: &HybridWindow<'_>,
        input: &[f64],
        rtmp_out: &mut [f64],
        taps: usize,
    ) {
        assert_eq!(input.len(), self.l, "hybrid-window input length mismatch");
        assert!(
            taps == 0 || rtmp_out.len() >= taps,
            "rtmp_out buffer too short for requested taps"
        );
        assert_eq!(window.window.len(), window.n3(), "window length mismatch");
        debug_assert_eq!(self.m, window.m);
        debug_assert_eq!(self.l, window.l);
        debug_assert_eq!(self.n, window.n);

        let m = self.m;
        let n_call = self.l;
        let n2 = window.n2();
        let n3 = window.n3();
        let n1 = window.n1();

        // | For N = 1,...,N2, do the next line          | Shift the old signal buffer
        // |     SB(N) = SB(N + NFRSZ)
        for ni in 0..n2 {
            self.sb[ni] = self.sb[ni + n_call];
        }
        // | For N = 1,...,NFRSZ, do the next line       | Shift in the new signal
        // |     SB(N2 + N) = STTMP(N)
        self.sb[n2..n2 + n_call].copy_from_slice(&input[..n_call]);

        // | K = 1
        // | For N = N3, N3-1, ..., 1: WS(N) = SB(N) * WNR(K); K = K + 1
        //
        // Spec is 1-based: SB(N) for N = N3..1 is sb[N-1] for N = N3..1.
        // The window WNR(K) for K = 1..N3 is window.window[K-1].
        // So pair sb[N-1] ↔ window.window[N3 - N] — i.e. sb[i] ↔
        // window.window[n3 - 1 - i] for 0-based i.
        let mut ws = vec![0.0f64; n3];
        for ni in 0..n3 {
            ws[ni] = self.sb[ni] * window.window[n3 - 1 - ni];
        }

        // | For I = 1,...,LPC+1, do the next four lines  | Update recursive part
        // |     TMP = 0
        // |     For N = LPC+1, LPC+2, ..., N1, do the next line
        // |         TMP = TMP + WS(N) * WS(N+1-I)
        // |     REXP(I) = (3/4) * REXP(I) + TMP
        //
        // Spec 1-based:
        //   WS(N)         for N = m+1..n1 is ws[N-1] for N = m+1..n1.
        //   WS(N + 1 - I) for I = 1..m+1 is ws[N - I] (1-based) = ws[N-1-(I-1)] in 0-based.
        // The 0-based equivalents: for `i_out in 0..=m`, for
        // `ni in m..n1`: tmp += ws[ni] * ws[ni - i_out].
        for i_out in 0..=m {
            let mut tmp = 0.0f64;
            for ni in m..n1 {
                tmp += ws[ni] * ws[ni - i_out];
            }
            self.rexp[i_out] = window.decay * self.rexp[i_out] + tmp;
        }

        // | For I = 1,...,LPC+1, do the next three lines | Add non-recursive part
        // |     RTMP(I) = REXP(I)
        // |     For N = N1+1, ..., N3, do the next line
        // |         RTMP(I) = RTMP(I) + WS(N) * WS(N+1-I)
        //
        // Annex I: during an erased frame block 43FE skips this part
        // entirely (`taps == 0`) and block 49FE truncates it to the
        // first 11 lags (`taps == 11`); see `advance_erased` /
        // `run_erased` above.
        for i_out in 0..taps {
            rtmp_out[i_out] = self.rexp[i_out];
            for ni in n1..n3 {
                rtmp_out[i_out] += ws[ni] * ws[ni - i_out];
            }
        }

        // | R(1) = R(1) * WNCF                           | White-noise correction
        if taps > 0 {
            rtmp_out[0] *= WNCF;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::{LPC, LPCLG, NFRSZ, NONR, NONRLG, NUPDATE};
    use crate::tables::{wnr_f64, wnrg_f64};

    #[test]
    fn synthesis_window_dimensions_match_table_2() {
        // Table 2/G.728 lists SB(1..105) and REXP(1..LPC+1=51). The
        // computed dimensions for the synthesis-filter hybrid window
        // must match.
        let wnr = wnr_f64();
        let hw = HybridWindow {
            m: LPC,
            l: NFRSZ,
            n: NONR,
            window: &wnr,
            decay: 0.75,
        };
        assert_eq!(hw.n3(), 105, "SB buffer length should be 105");
        assert_eq!(hw.n2(), LPC + NONR);
        assert_eq!(hw.n1(), LPC + NFRSZ);
        let st = HybridWindowState::new(&hw);
        assert_eq!(st.sb().len(), 105);
        assert_eq!(st.rexp().len(), LPC + 1);
    }

    #[test]
    fn loggain_window_dimensions_match_table_2() {
        // Table 2/G.728 lists SBLG(1..34) and REXPLG(1..LPCLG+1=11).
        let wnrg = wnrg_f64();
        let hw = HybridWindow {
            m: LPCLG,
            l: NUPDATE,
            n: NONRLG,
            window: &wnrg,
            decay: 0.75,
        };
        assert_eq!(hw.n3(), 34, "SBLG buffer length should be 34");
        let st = HybridWindowState::new(&hw);
        assert_eq!(st.sb().len(), 34);
        assert_eq!(st.rexp().len(), LPCLG + 1);
    }

    #[test]
    fn fresh_state_buffers_are_zero() {
        // Table 2/G.728 lists `0,0,...,0` as the initial value for
        // every SB / REXP buffer.
        let wnr = wnr_f64();
        let hw = HybridWindow {
            m: LPC,
            l: NFRSZ,
            n: NONR,
            window: &wnr,
            decay: 0.75,
        };
        let st = HybridWindowState::new(&hw);
        assert!(st.sb().iter().all(|&v| v == 0.0));
        assert!(st.rexp().iter().all(|&v| v == 0.0));
    }

    #[test]
    fn zero_input_yields_zero_autocorr() {
        // With zero history and a zero-valued input vector, every
        // autocorrelation tap is zero; the WNCF white-noise correction
        // multiplies zero by a non-zero scalar and stays zero.
        let wnr = wnr_f64();
        let hw = HybridWindow {
            m: LPC,
            l: NFRSZ,
            n: NONR,
            window: &wnr,
            decay: 0.75,
        };
        let mut st = HybridWindowState::new(&hw);
        let input = vec![0.0; NFRSZ];
        let mut rtmp = vec![0.0; LPC + 1];
        st.run(&hw, &input, &mut rtmp);
        assert!(rtmp.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn wncf_factor_is_applied_to_rtmp0() {
        // Drive in a single non-zero sample at the very newest slot,
        // and check that rtmp[0] picks up the WNCF=257/256 factor.
        //
        // For the loggain window (smallest dims, easiest to reason
        // about), the newest sample is multiplied by window[n3-1-(n3-1)]
        // = window[0]. With only one non-zero sample at sb[n3-1], the
        // entire WS vector has one non-zero entry: ws[n3-1] = x ·
        // window[0]. The R(1) autocorrelation tap = sum over N=N1+1..N3
        // of ws[N-1]·ws[N-1]; the only non-zero term is at N=N3, giving
        // R(1) = (x·window[0])². WNCF then multiplies this by 257/256.
        let wnrg = wnrg_f64();
        let hw = HybridWindow {
            m: LPCLG,
            l: NUPDATE,
            n: NONRLG,
            window: &wnrg,
            decay: 0.75,
        };
        let mut st = HybridWindowState::new(&hw);
        let mut input = vec![0.0; NUPDATE];
        input[NUPDATE - 1] = 1.0; // the newest sample
        let mut rtmp = vec![0.0; LPCLG + 1];
        st.run(&hw, &input, &mut rtmp);
        let w0 = wnrg[0];
        let expected = w0 * w0 * WNCF;
        assert!(
            (rtmp[0] - expected).abs() < 1e-12,
            "R(1) should pick up WNCF correction: got {} vs expected {}",
            rtmp[0],
            expected
        );
    }

    #[test]
    fn autocorrelation_r0_is_nonnegative_for_arbitrary_input() {
        // R(1) = Σ x[n]² · w[n]² · WNCF is a sum of squares times a
        // positive constant, so it must always be non-negative. Pin
        // that invariant for a random-looking input pattern.
        let wnrg = wnrg_f64();
        let hw = HybridWindow {
            m: LPCLG,
            l: NUPDATE,
            n: NONRLG,
            window: &wnrg,
            decay: 0.75,
        };
        let mut st = HybridWindowState::new(&hw);
        // Drive several cycles to populate the recursive part.
        for cycle in 0..5 {
            let mut input = [0.0; NUPDATE];
            for k in 0..NUPDATE {
                input[k] = ((cycle * NUPDATE + k) as f64).sin() * 10.0;
            }
            let mut rtmp = vec![0.0; LPCLG + 1];
            st.run(&hw, &input, &mut rtmp);
            assert!(
                rtmp[0] >= 0.0,
                "R(1) must be non-negative; got {} on cycle {}",
                rtmp[0],
                cycle
            );
        }
    }
}
