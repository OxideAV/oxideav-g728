//! Annex G §G.3.17 / §G.3.18 — the shared fixed-point hybrid-window core.
//!
//! ITU-T G.728 Annex G (1994-11) realises the three backward-adaptation
//! hybrid windows — block 49 (synthesis filter, §G.3.17), block 36
//! (perceptual weighting filter) and block 43 (log-gain predictor) — with
//! the *same* fixed-point structure. The floating-point reference already
//! factors this commonality into one [`crate::hybrid_window::HybridWindowState`]
//! parameter object; this module is the fixed-point analogue.
//!
//! All three blocks
//!
//! 1. buffer their windowed input signal in **segmented block
//!    floating-point** (SBFL) form: the signal history is a concatenation
//!    of `IDIM`-sample segments, each with its own number-of-left-shifts
//!    (`NLSSB(.)`). Block 49 / block 36 receive `NFRSZ = 20` new speech
//!    samples (4 segments) per cycle; block 43 receives `NUPDATE = 4`
//!    log-gain scalars (which this core treats as `NUPDATE / IDIM`… — see
//!    the per-block front-end notes),
//! 2. shift the buffer, align every segment to the common minimum shift
//!    `NLSTMP`, multiply by the Q15 window and round into the BFL windowed
//!    signal `WS` (block 49 / the §G.3.17 windowing step), and
//! 3. accumulate the recursive (3/4-decay) and non-recursive
//!    autocorrelation components, apply the white-noise correction and
//!    emit `RTMP(1..=order+1)` (BFL, 16-bit) plus the `ILLCOND` verdict
//!    (the full 32-bit `R(order+1)` accumulator being zero — the §G.2.2
//!    interoperability fix) (HWMCORE, §G.3.18).
//!
//! Only the *dimensions* differ between blocks: the LPC order `order`, the
//! per-cycle new-segment count, the non-recursive tail length and the
//! window table. [`HybridWindowFixed`] captures those as a parameter
//! object; [`HybridWindowFixedState`] carries the permanent SBFL buffer +
//! recursive-autocorrelation state and runs the two §G.3.17 / §G.3.18
//! steps. The recursive decay is `3/4` for every G.728 block, realised in
//! HWMCORE as the `−RREC << NLSATT50 + RREC << 16` combination with
//! [`NLSATT50`] `= 14` (`(−2¹⁴ + 2¹⁶)/2¹⁶ = 3/4`).
//!
//! The synthesis-filter adapter [`crate::annex_g_synth_adapter`] and the
//! weighting-filter adapter [`crate::annex_g_weight_adapter`] both drive
//! this core; the per-block bandwidth-expansion and Levinson wiring stay in
//! those modules.

use crate::annex_g_arith::{findnls, rnd, shl_sat, shr_sat};
use crate::consts::IDIM;

/// Initial `NLSSB` per-segment shift — Table G.2/G.728 lists `NLSSB`
/// (and `NLSSTATE` / `NLSSTTMP`) with "initial value = 16": a zeroed
/// signal-history buffer is held at the maximum useful left shift.
pub const NLSSB_INIT_G2: i32 = 16;

/// Initial `NLSREXP` / `NLSREXPW` / `NLSREXPLG` — Table G.2/G.728 lists
/// all three recursive-autocorrelation shifts with "initial value = 31".
pub const NLSREXP_INIT_G2: i32 = 31;

/// `NLSATT50 = 14` — the attenuation NLS realising the `3/4` recursive
/// decay of the synthesis-filter (block 49) and log-gain (block 43)
/// hybrid windows: HWMCORE forms `−RREC << 14 + RREC << 16` over a
/// `>> 16`, i.e. `(−2¹⁴ + 2¹⁶)/2¹⁶ = 49152/65536 = 0.75`.
pub const NLSATT50: i32 = 14;

/// `NLSATTW = 15` — the attenuation NLS realising block 36's `1/2`
/// recursive decay (§G.3.12; the base spec's float block 36 spells
/// `REXPW(I) = (1/2)·REXPW(I) + TMP`): `(−2¹⁵ + 2¹⁶)/2¹⁶ = 0.5`. Block 36
/// decays by `1/2`, **not** the `3/4` of sibling blocks 43/49 — the
/// `3/4` listed first in HWMCORE's shared margin comment does not apply
/// to block 36 (erratum E2 in `docs/audio/g728/g728-errata.md`).
pub const NLSATTW: i32 = 15;

/// `MLS` for the double-precision accumulator passed to the §G.1.3
/// `VSCALE` / `FINDNLS` inside HWMCORE (the 32-bit `AA0` / `AA1`
/// accumulators normalize against `MLS = 30`).
const MLS_ACC: i32 = 30;

/// One signal segment in 14-bit BFL form: `IDIM` mantissas sharing one
/// number-of-left-shifts `nls`. This is the per-segment unit the §G.3.11
/// block-32 synthesis filter produces (the `ST` / `NLSST` pair) and the
/// unit the windowing step buffers in its SBFL array.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BflSegment {
    /// `IDIM` 14-bit BFL mantissas, newest-last within the vector.
    pub samples: [i16; IDIM],
    /// Shared number-of-left-shifts for this segment.
    pub nls: i32,
}

impl BflSegment {
    /// A zero segment held at the given shift.
    #[must_use]
    pub const fn zero(nls: i32) -> Self {
        Self {
            samples: [0; IDIM],
            nls,
        }
    }
}

/// Pure-data description of one fixed-point hybrid window's dimensions and
/// table — the fixed-point analogue of
/// [`crate::hybrid_window::HybridWindow`].
#[derive(Debug, Clone, Copy)]
pub struct HybridWindowFixed<'a> {
    /// LPC order — `m` in the spec (LPC = 50, LPCW = 10, LPCLG = 10).
    pub order: usize,
    /// Per-call input length — `l` in the spec (NFRSZ = 20 for blocks
    /// 36/49). Must be a multiple of `IDIM`.
    pub l: usize,
    /// Non-recursive tail length — `n` in the spec (NONR = 35, NONRW = 30).
    pub n: usize,
    /// Q15 window weights, length `order + l + n` (`N3`). Element 0 is the
    /// spec's `WNR(1)` (the newest sample's weight in the right-to-left
    /// walk).
    pub window: &'a [i16],
    /// Attenuation NLS for the recursive decay ([`NLSATT50`] = 14 for
    /// blocks 43/49's 3/4, [`NLSATTW`] = 15 for block 36's 1/2).
    pub nlsatt: i32,
}

impl HybridWindowFixed<'_> {
    /// `N1 = order + l` — end of the recursive-component tap range.
    #[must_use]
    pub const fn n1(&self) -> usize {
        self.order + self.l
    }
    /// `N2 = order + n` — length of the kept-old portion of `SB`.
    #[must_use]
    pub const fn n2(&self) -> usize {
        self.order + self.n
    }
    /// `N3 = order + l + n` — total length of `SB` / `WS`.
    #[must_use]
    pub const fn n3(&self) -> usize {
        self.order + self.l + self.n
    }
    /// `N4 = N3 / IDIM` — number of SBFL segments in `SB`.
    #[must_use]
    pub const fn n4(&self) -> usize {
        self.n3() / IDIM
    }
    /// `N5 = l / IDIM` — number of new segments shifted in per cycle.
    #[must_use]
    pub const fn n5(&self) -> usize {
        self.l / IDIM
    }
    /// `N6 = N4 − N5` — number of `NLSSB` entries kept across the shift.
    #[must_use]
    pub const fn n6(&self) -> usize {
        self.n4() - self.n5()
    }
}

/// Result of one §G.3.18 HWMCORE call: the autocorrelation `RTMP` (BFL,
/// 16-bit mantissas, Q15-normalized so `RTMP(1)` dominates) plus the
/// `ILLCOND` verdict (`RTMP(order+1)`'s 32-bit accumulator was zero).
#[derive(Debug, Clone)]
pub struct HwmcoreOut {
    /// `R(1..=order+1)` BFL mantissas (`rtmp[0]` is the spec's `R(1)`).
    pub rtmp: Vec<i16>,
    /// `ILLCOND` — `true` when the 32-bit `R(order+1)` accumulator is zero.
    pub illcond: bool,
}

/// The permanent state of one §G.3.18 HWMCORE instance: the recursive-
/// autocorrelation `RREC`/`REXP` BFL mantissas plus their shared shift
/// `NLSREXP`. Every G.728 hybrid window (blocks 36 / 43 / 49) owns one;
/// the segmented-SBFL front end ([`HybridWindowFixedState`]) and the
/// flat-`Q9` log-gain front end (§G.3.14 block 43) both drive it.
#[derive(Debug, Clone)]
pub struct HwmcoreState {
    /// `REXP(1..=order+1)` recursive-autocorrelation BFL mantissas.
    rexp: Vec<i16>,
    /// `NLSREXP` — the shared shift of the BFL `REXP` block.
    nlsrexp: i32,
}

impl HwmcoreState {
    /// Fresh HWMCORE state for an `order`-tap window: zeroed `REXP`
    /// (Table G.2 lists `REXP` / `REXPW` / `REXPLG` all zero at reset)
    /// with `NLSREXP` at its Table G.2 initial value 31.
    #[must_use]
    pub fn new(order: usize) -> Self {
        Self {
            rexp: vec![0; order + 1],
            nlsrexp: NLSREXP_INIT_G2,
        }
    }

    /// Borrow the recursive-autocorrelation `REXP` (for tests/audit).
    #[must_use]
    pub fn rexp(&self) -> &[i16] {
        &self.rexp
    }

    /// `NLSREXP` — the shared shift of `REXP` (for tests/audit).
    #[must_use]
    pub fn nlsrexp(&self) -> i32 {
        self.nlsrexp
    }

    /// §G.3.12's post-core clamp ("If NLSREXPW > 41, set NLSREXPW =
    /// 41") — cap the recursive-autocorrelation shift to avoid reduced
    /// accuracy during long zero-input periods.
    pub fn clamp_nlsrexp(&mut self, max: i32) {
        if self.nlsrexp > max {
            self.nlsrexp = max;
        }
    }
}

/// Mutable state of one fixed-point hybrid-window adaptation chain: the
/// permanent SBFL signal-history buffer `SB` with its per-segment `NLSSB`,
/// and the recursive-autocorrelation core [`HwmcoreState`] (`REXP` /
/// `NLSREXP`). Sized to the [`HybridWindowFixed`] it is constructed for.
#[derive(Debug, Clone)]
pub struct HybridWindowFixedState {
    order: usize,
    l: usize,
    n: usize,
    /// `SB(1..=N3)` SBFL mantissas, 0-based. 14-bit precision.
    sb: Vec<i16>,
    /// `NLSSB(1..=N4)` per-segment shifts, 0-based.
    nlssb: Vec<i32>,
    /// Blocks §G.3.18 — the recursive-autocorrelation core.
    core: HwmcoreState,
}

impl HybridWindowFixedState {
    /// Fresh state: zeroed `SB` / `REXP` per Table G.2/G.728, with the
    /// Table G.2 initial shifts (`NLSSB = 16`, `NLSREXP = 31`).
    #[must_use]
    pub fn new(win: &HybridWindowFixed<'_>) -> Self {
        Self {
            order: win.order,
            l: win.l,
            n: win.n,
            sb: vec![0; win.n3()],
            nlssb: vec![NLSSB_INIT_G2; win.n4()],
            core: HwmcoreState::new(win.order),
        }
    }

    /// Borrow the SBFL signal-history buffer `SB` (for tests/audit).
    #[must_use]
    pub fn sb(&self) -> &[i16] {
        &self.sb
    }

    /// Borrow the per-segment shift array `NLSSB` (for tests/audit).
    #[must_use]
    pub fn nlssb(&self) -> &[i32] {
        &self.nlssb
    }

    /// Borrow the recursive-autocorrelation `REXP` (for tests/audit).
    #[must_use]
    pub fn rexp(&self) -> &[i16] {
        self.core.rexp()
    }

    /// `NLSREXP` — the shared shift of `REXP` (for tests/audit).
    #[must_use]
    pub fn nlsrexp(&self) -> i32 {
        self.core.nlsrexp()
    }

    /// Run one hybrid-window + HWMCORE cycle on `new_segs` (the `N5` new
    /// `IDIM`-sample segments of this cycle, newest-vector last) and return
    /// the autocorrelation `RTMP` plus the `ILLCOND` verdict.
    pub fn run(&mut self, win: &HybridWindowFixed<'_>, new_segs: &[BflSegment]) -> HwmcoreOut {
        debug_assert_eq!(self.order, win.order);
        debug_assert_eq!(self.l, win.l);
        debug_assert_eq!(self.n, win.n);
        debug_assert_eq!(new_segs.len(), win.n5(), "new-segment count must be N5");

        let (ws, nlstmp) = self.window_step(win, new_segs);
        self.core
            .run(win.order, win.n1(), win.n3(), win.nlsatt, &ws, nlstmp)
    }

    /// Annex I block **49FE** (§I.5.7): the same window step, but the
    /// core runs as `HWMCOREFE` with the non-recursive output truncated
    /// to `R(1..=lpfe+1)` (`lpfe = 10` during a frame erasure — enough
    /// for the post-filter's order-10 by-product). The signal buffer
    /// and the recursive `REXP` update are identical to [`Self::run`],
    /// so the adaptation chain resynchronises seamlessly at the next
    /// good frame.
    pub fn run_erased(
        &mut self,
        win: &HybridWindowFixed<'_>,
        new_segs: &[BflSegment],
        lpfe: usize,
    ) -> HwmcoreOut {
        debug_assert_eq!(new_segs.len(), win.n5(), "new-segment count must be N5");
        let (ws, nlstmp) = self.window_step(win, new_segs);
        self.core
            .run_fe(win.order, win.n1(), win.n3(), win.nlsatt, &ws, nlstmp, lpfe)
    }

    /// §G.3.17 windowing step (block 49 / block 36): shift the SBFL `SB`
    /// buffer, align every segment to the common minimum shift `NLSTMP`,
    /// multiply by the Q15 window, and return the BFL windowed signal `WS`
    /// with that shared shift.
    fn window_step(
        &mut self,
        win: &HybridWindowFixed<'_>,
        new_segs: &[BflSegment],
    ) -> (Vec<i16>, i32) {
        let n2 = win.n2();
        let n3 = win.n3();
        let n4 = win.n4();
        let n5 = win.n5();
        let n6 = win.n6();
        let l = win.l;

        // | For N = 1, ..., N2: SB(N) = SB(N + l)   | Shift old part of SB
        for ni in 0..n2 {
            self.sb[ni] = self.sb[ni + l];
        }
        // | For N = 1, ..., N6: NLSSB(N) = NLSSB(N + N5) | Shift old NLSSB
        for ni in 0..n6 {
            self.nlssb[ni] = self.nlssb[ni + n5];
        }
        // | For N = 1, ..., l: SB(N2 + N) = STTMP(N) | Shift in new part
        // | For N = 1, ..., N5: NLSSB(N6 + N) = NLSSTTMP(N)
        for (j, seg) in new_segs.iter().enumerate() {
            let base = n2 + j * IDIM;
            self.sb[base..base + IDIM].copy_from_slice(&seg.samples);
            self.nlssb[n6 + j] = seg.nls;
        }

        // | NLSTMP = Min{NLSSB(1), ..., NLSSB(N4)}   | determines NLSWS
        let nlstmp = *self.nlssb[..n4].iter().min().expect("N4 > 0");

        // | K = 1; N = N3
        // | For J = 1, ..., N4:
        // |   NRSH = NLSSB(J) − NLSTMP − 1           | −1 for Q15 mult
        // |   For M = 1, ..., IDIM:
        // |     P = SB(K) * WNR(N)                    | WNR is Q15
        // |     If NRSH = −1, AA0 = P << 1
        // |     If NRSH > −1, AA0 = P >> NRSH
        // |     WS(K) = RND(AA0); N = N − 1; K = K + 1
        let mut ws = vec![0i16; n3];
        let mut k = 0usize;
        for j in 0..n4 {
            let nrsh = self.nlssb[j] - nlstmp - 1;
            for _m in 0..IDIM {
                let p = self.sb[k] as i32 * win.window[n3 - 1 - k] as i32;
                let aa0 = if nrsh == -1 {
                    shl_sat(p, 1)
                } else {
                    shr_sat(p, nrsh as u32)
                };
                ws[k] = rnd(aa0);
                k += 1;
            }
        }

        (ws, nlstmp)
    }
}

impl HwmcoreState {
    /// §G.3.18 HWMCORE: accumulate the recursive (3/4-decay) and
    /// non-recursive autocorrelation components from the BFL windowed
    /// signal `ws` (shared shift `nlstmp`), thread the BFL shift through
    /// the `REXP` update, apply the white-noise correction, and emit the
    /// BFL `RTMP` plus the `ILLCOND` verdict. `order` / `n1` / `n3` are
    /// the caller's window dimensions (`LPC`+`NFRSZ`… for block 49,
    /// `LPCLG`+`NUPDATE`… for block 43, `LPCW`+`NFRSZ`… for block 36)
    /// and `nlsatt` selects the recursive decay (`NLSATT50` = 14 ⇒ 3/4,
    /// `NLSATTW` = 15 ⇒ 1/2).
    pub fn run(
        &mut self,
        order: usize,
        n1: usize,
        n3: usize,
        nlsatt: i32,
        ws: &[i16],
        nlstmp: i32,
    ) -> HwmcoreOut {
        self.run_fe(order, n1, n3, nlsatt, ws, nlstmp, order)
    }

    /// Annex I `HWMCOREFE` (§I.5.7) — identical to [`Self::run`] except
    /// that the **non-recursive** output is truncated to `R(1..lpfe+1)`
    /// ("in the event of a frame erasure, only the values of RTMP(1),
    /// …, RTMP(11) are calculated, although REXP(1), …, REXP(51) are
    /// updated"): the recursive `REXP` update always spans the full
    /// `order + 1` lags, while the returned `rtmp` has `lpfe + 1`
    /// entries and the `ILLCOND` verdict tests the 32-bit `R(lpfe+1)`
    /// accumulator (`LPFE = LPO` on good frames, `LPFE = 10` during an
    /// erasure — the §I.5.7 `If FERROR = .TRUE., set LPFE = 10` line).
    #[allow(clippy::too_many_arguments)]
    pub fn run_fe(
        &mut self,
        order: usize,
        n1: usize,
        n3: usize,
        nlsatt: i32,
        ws: &[i16],
        nlstmp: i32,
        lpfe: usize,
    ) -> HwmcoreOut {
        let lpo = order;
        debug_assert!(lpfe <= lpo);
        debug_assert_eq!(self.rexp.len(), lpo + 1);
        debug_assert!(ws.len() >= n3);

        // | NLSAA0 = 2 * NLSTMP    | scale of WS²
        let nls_aa0 = 2 * nlstmp;

        // Recursive-range energy AA0 = Σ_{N=LPO+1..N1} WS(N)·WS(N−i).
        let energy = |i: usize| -> i64 {
            let mut acc = 0i64;
            for n in (lpo + 1)..=n1 {
                acc += ws[n - 1] as i64 * ws[n - 1 - i] as i64;
            }
            acc
        };

        let nls_rrec = self.nlsrexp;
        let nlsre: i32;
        let mut new_rexp = vec![0i16; lpo + 1];

        // R(1) recursive part (I = 0) across the three NLS-alignment cases.
        let aa0_e0 = energy(0);
        if nls_rrec > nls_aa0 {
            // Case 1.
            let ir = nls_rrec - nls_aa0 + 1;
            let aa0 = aa0_e0 >> 1;
            let r1 = self.rexp[0] as i64;
            let aa1 = (-(r1 << nlsatt)) + (r1 << 16);
            let aa1 = aa1 >> ir;
            let combined = aa0 + aa1;
            nlsre = findnls(&[clamp_i32(combined)], 1, MLS_ACC);
            new_rexp[0] = rnd(clamp_i32(shl_acc(combined, nlsre)));
            for i in 1..=lpo {
                let mut aa0 = energy(i) >> 1;
                let ri = self.rexp[i] as i64;
                // The published §G.3.18 case-1 loop line reads
                // "AA1 = AA1 + (RREC(I+1) << 16)" — a sign typo (it
                // would scale RREC by 5/4). Its own margin comment says
                // "Scale RREC by 3/4", the R(1) line three lines above
                // and all four case-2/case-3 instances spell
                // "AA1 = −AA1 + (RREC << 16)", and the floating-point
                // recursion is REXP = (3/4)·REXP + TMP; we follow the
                // five consistent instances — erratum E1 in
                // `docs/audio/g728/g728-errata.md` (the Annex I republished
                // `HWMCOREFE` also prints the corrected minus form).
                let aa1 = ((-(ri << nlsatt)) + (ri << 16)) >> ir;
                aa0 += aa1;
                aa0 = shl_acc(aa0, nlsre);
                new_rexp[i] = rnd(clamp_i32(aa0));
            }
            self.nlsrexp = nls_aa0 - 1 + nlsre;
        } else if nls_rrec == nls_aa0 {
            // Case 2.
            let r1 = self.rexp[0] as i64;
            let aa1 = (-(r1 << nlsatt)) + (r1 << 16);
            let aa0 = aa0_e0 >> 1;
            let aa1 = aa1 >> 1;
            let combined = aa0 + aa1;
            nlsre = findnls(&[clamp_i32(combined)], 1, MLS_ACC);
            new_rexp[0] = rnd(clamp_i32(shl_acc(combined, nlsre)));
            for i in 1..=lpo {
                let mut aa0 = energy(i) >> 1;
                let ri = self.rexp[i] as i64;
                let aa1 = ((-(ri << nlsatt)) + (ri << 16)) >> 1;
                aa0 += aa1;
                aa0 = shl_acc(aa0, nlsre);
                new_rexp[i] = rnd(clamp_i32(aa0));
            }
            self.nlsrexp = nls_rrec - 1 + nlsre;
        } else {
            // Case 3: nls_rrec < nls_aa0.
            let ir = nls_aa0 - nls_rrec + 1;
            let aa0 = aa0_e0 >> ir;
            let r1 = self.rexp[0] as i64;
            let aa1 = ((-(r1 << nlsatt)) + (r1 << 16)) >> 1;
            let combined = aa0 + aa1;
            nlsre = findnls(&[clamp_i32(combined)], 1, MLS_ACC);
            new_rexp[0] = rnd(clamp_i32(shl_acc(combined, nlsre)));
            for i in 1..=lpo {
                let mut aa0 = energy(i) >> ir;
                let ri = self.rexp[i] as i64;
                let aa1 = ((-(ri << nlsatt)) + (ri << 16)) >> 1;
                aa0 += aa1;
                aa0 = shl_acc(aa0, nlsre);
                new_rexp[i] = rnd(clamp_i32(aa0));
            }
            self.nlsrexp = nls_rrec - 1 + nlsre;
        }
        self.rexp.copy_from_slice(&new_rexp);

        // Non-recursive component → R = RTMP, aligned with the just-updated
        // RREC (now at shift self.nlsrexp).
        let nls_rrec2 = self.nlsrexp;
        let nonrec = |i: usize| -> i64 {
            let mut acc = 0i64;
            for n in (n1 + 1)..=n3 {
                acc += ws[n - 1] as i64 * ws[n - 1 - i] as i64;
            }
            acc
        };

        // §I.5.7 HWMCOREFE: `LPFE = LPO` normally, 10 during a frame
        // erasure — only R(1..LPFE+1) are produced below.
        let mut rtmp = vec![0i16; lpfe + 1];
        let mut illcond = false;
        let nlsrr: i32;

        if nls_rrec2 > nls_aa0 {
            // Case 1.
            let ir = nls_rrec2 - nls_aa0 + 1;
            let r1 = self.rexp[0] as i64;
            let aa1 = (r1 << 16) >> ir;
            let aa0 = nonrec(0) >> 1;
            let mut aa1 = aa0 + aa1;
            let wn = aa1 >> 8;
            aa1 += wn;
            nlsrr = findnls(&[clamp_i32(aa1)], 1, MLS_ACC);
            rtmp[0] = rnd(clamp_i32(shl_acc(aa1, nlsrr)));
            for i in 1..=lpfe {
                let aa0 = nonrec(i) >> 1;
                let ri = self.rexp[i] as i64;
                let aa1 = (ri << 16) >> ir;
                let aa1 = aa0 + aa1;
                let aa1 = shl_acc(aa1, nlsrr);
                rtmp[i] = rnd(clamp_i32(aa1));
                if i == lpfe && aa1 == 0 {
                    illcond = true;
                }
            }
        } else if nls_rrec2 == nls_aa0 {
            // Case 2.
            let r1 = self.rexp[0] as i64;
            let aa0 = nonrec(0) >> 1;
            let aa1 = r1 << 15;
            let mut aa1 = aa0 + aa1;
            let wn = aa1 >> 8;
            aa1 += wn;
            nlsrr = findnls(&[clamp_i32(aa1)], 1, MLS_ACC);
            rtmp[0] = rnd(clamp_i32(shl_acc(aa1, nlsrr)));
            for i in 1..=lpfe {
                let aa0 = nonrec(i) >> 1;
                let ri = self.rexp[i] as i64;
                let aa1 = ri << 15;
                let aa1 = aa0 + aa1;
                let aa1 = shl_acc(aa1, nlsrr);
                rtmp[i] = rnd(clamp_i32(aa1));
                if i == lpfe && aa1 == 0 {
                    illcond = true;
                }
            }
        } else {
            // Case 3.
            let ir = nls_aa0 - nls_rrec2 + 1;
            let r1 = self.rexp[0] as i64;
            let aa0 = nonrec(0) >> ir;
            let aa1 = r1 << 15;
            let mut aa1 = aa0 + aa1;
            let wn = aa1 >> 8;
            aa1 += wn;
            nlsrr = findnls(&[clamp_i32(aa1)], 1, MLS_ACC);
            rtmp[0] = rnd(clamp_i32(shl_acc(aa1, nlsrr)));
            for i in 1..=lpfe {
                let aa0 = nonrec(i) >> ir;
                let ri = self.rexp[i] as i64;
                let aa1 = ri << 15;
                let aa1 = aa0 + aa1;
                let aa1 = shl_acc(aa1, nlsrr);
                rtmp[i] = rnd(clamp_i32(aa1));
                if i == lpfe && aa1 == 0 {
                    illcond = true;
                }
            }
        }

        HwmcoreOut { rtmp, illcond }
    }
}

/// Left-shift a 64-bit accumulator by a possibly-negative count: a
/// positive `k` shifts left, a negative `k` shifts (arithmetic) right.
/// Used for the HWMCORE `AA0 << NLSRE` line where `NLSRE` may be negative.
#[inline]
fn shl_acc(acc: i64, k: i32) -> i64 {
    if k >= 0 {
        acc << k
    } else {
        acc >> (-k)
    }
}

/// Saturate a 64-bit accumulator value into the 32-bit range so the
/// §G.1.3 primitives (`findnls` / `rnd`), which operate on `i32`, see the
/// same saturated value the annex's 32-bit accumulator holds.
#[inline]
fn clamp_i32(v: i64) -> i32 {
    v.clamp(i32::MIN as i64, i32::MAX as i64) as i32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::{LPC, LPCW, NFRSZ, NONR, NONRW};
    use crate::tables::{WNRW_Q15, WNR_Q15};

    // `lcg_samples` / `sbfl_segments_requant` are defined further down and
    // shared by the float-tracking + broadband tests.

    fn synth_win() -> HybridWindowFixed<'static> {
        HybridWindowFixed {
            order: LPC,
            l: NFRSZ,
            n: NONR,
            window: &WNR_Q15,
            nlsatt: NLSATT50,
        }
    }

    fn weight_win() -> HybridWindowFixed<'static> {
        HybridWindowFixed {
            order: LPCW,
            l: NFRSZ,
            n: NONRW,
            window: &WNRW_Q15,
            nlsatt: NLSATTW,
        }
    }

    #[test]
    fn synth_dimensions_match_spec() {
        let w = synth_win();
        assert_eq!(w.n1(), 70);
        assert_eq!(w.n2(), 85);
        assert_eq!(w.n3(), 105);
        assert_eq!(w.n4(), 21);
        assert_eq!(w.n5(), 4);
        assert_eq!(w.n6(), 17);
    }

    #[test]
    fn weight_dimensions_match_spec() {
        let w = weight_win();
        assert_eq!(w.n1(), 30);
        assert_eq!(w.n2(), 40);
        assert_eq!(w.n3(), 60);
        assert_eq!(w.n4(), 12);
        assert_eq!(w.n5(), 4);
        assert_eq!(w.n6(), 8);
        // Window table length must equal N3.
        assert_eq!(w.window.len(), w.n3());
    }

    #[test]
    fn fresh_state_buffers_are_zero() {
        let w = synth_win();
        let s = HybridWindowFixedState::new(&w);
        assert_eq!(s.sb().len(), w.n3());
        assert_eq!(s.nlssb().len(), w.n4());
        assert_eq!(s.rexp().len(), w.order + 1);
        assert!(s.sb().iter().all(|&v| v == 0));
        assert!(s.rexp().iter().all(|&v| v == 0));
    }

    #[test]
    fn zero_input_flags_illcond() {
        // A zero-speech cycle drives RTMP to all-zero ⇒ the upstream
        // ILLCOND verdict fires (R(order+1) == 0).
        for w in [synth_win(), weight_win()] {
            let mut s = HybridWindowFixedState::new(&w);
            let segs = vec![BflSegment::zero(15); w.n5()];
            let out = s.run(&w, &segs);
            assert!(out.illcond, "zero input ⇒ ILLCOND");
            assert!(out.rtmp.iter().all(|&v| v == 0));
            assert_eq!(out.rtmp.len(), w.order + 1);
        }
    }

    #[test]
    fn new_segments_land_at_buffer_tail() {
        let w = synth_win();
        let mut s = HybridWindowFixedState::new(&w);
        let mut segs = vec![BflSegment::zero(12); w.n5()];
        for (j, seg) in segs.iter_mut().enumerate() {
            for (m, v) in seg.samples.iter_mut().enumerate() {
                *v = (100 * (j as i16 + 1)) + m as i16;
            }
        }
        s.run(&w, &segs);
        // The newest segment lands at SB(N2 + (N5-1)*IDIM ..).
        let base = w.n2() + (w.n5() - 1) * IDIM;
        assert_eq!(&s.sb()[base..base + IDIM], &segs[w.n5() - 1].samples[..]);
        assert_eq!(s.nlssb()[w.n4() - 1], segs[w.n5() - 1].nls);
    }

    #[test]
    fn autocorrelation_tracks_floating_point_window() {
        // The fixed-point RTMP must track the floating-point hybrid
        // window's RTMP in *shape* (the normalized autocorrelation
        // sequence). Both are driven by the *same requantized* 14-bit
        // speech, so any divergence is purely the §G.3.18 BFL arithmetic
        // versus the exact floating-point recursion. We compare the
        // low-order normalized autocorrelation R(k)/R(1) — the taps that
        // drive the predictor — within a tolerance reflecting the 16-bit
        // RTMP rounding and the segmented-BFL alignment shifts. Runs for
        // both the LPC (block 49) and LPCW (block 36) windows.
        use crate::hybrid_window::{HybridWindow, HybridWindowState};
        use crate::tables::{wnr_f64, wnrw_f64};

        for (order, l, n, fwin, fxwin, decay, nlsatt) in [
            (
                LPC,
                NFRSZ,
                NONR,
                wnr_f64().to_vec(),
                WNR_Q15.to_vec(),
                0.75,
                NLSATT50,
            ),
            (
                LPCW,
                NFRSZ,
                NONRW,
                wnrw_f64().to_vec(),
                WNRW_Q15.to_vec(),
                0.5,
                NLSATTW,
            ),
        ] {
            let fhw = HybridWindow {
                m: order,
                l,
                n,
                window: &fwin,
                decay,
            };
            let xhw = HybridWindowFixed {
                order,
                l,
                n,
                window: &fxwin,
                nlsatt,
            };
            let mut fstate = HybridWindowState::new(&fhw);
            let mut xstate = HybridWindowFixedState::new(&xhw);

            let mut last_norm_fixed = vec![0.0f64; order + 1];
            let mut last_norm_float = vec![0.0f64; order + 1];
            let mut compared = false;

            let mut seed = 0x0bad_c0de_dead_beefu64;
            for cyc in 0..12 {
                let samples = lcg_samples(&mut seed, NFRSZ, 7000.0);
                let (segs, requant) = sbfl_segments_requant(&samples, xhw.n5());

                let mut frtmp = vec![0.0f64; order + 1];
                fstate.run(&fhw, &requant, &mut frtmp);
                let core = xstate.run(&xhw, &segs);

                if cyc >= 9 && frtmp[0].abs() > 1e-6 && core.rtmp[0] != 0 {
                    let fr1 = frtmp[0];
                    let xr1 = core.rtmp[0] as f64;
                    for k in 0..=order {
                        last_norm_float[k] = frtmp[k] / fr1;
                        last_norm_fixed[k] = core.rtmp[k] as f64 / xr1;
                    }
                    compared = true;
                }
            }

            assert!(compared, "comparison cycle must have produced finite R(1)");
            for k in 0..=8 {
                let diff = (last_norm_fixed[k] - last_norm_float[k]).abs();
                assert!(
                    diff < 0.12,
                    "order {order}: normalized RTMP[{k}] mismatch: fixed {} vs float {} (diff {diff})",
                    last_norm_fixed[k],
                    last_norm_float[k]
                );
            }
        }
    }

    /// LCG noise generator producing `len` samples in `[-amp, amp)`.
    fn lcg_samples(seed: &mut u64, len: usize, amp: f64) -> Vec<f64> {
        let mut out = vec![0.0f64; len];
        for s in out.iter_mut() {
            *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            let u = ((*seed >> 33) as f64) / ((1u64 << 31) as f64) - 1.0;
            *s = u * amp;
        }
        out
    }

    /// Build the `n5` SBFL segments AND return the float samples actually
    /// represented (re-expanded from the mantissas) so the float oracle
    /// sees the same quantized values the fixed-point core does.
    fn sbfl_segments_requant(samples: &[f64], n5: usize) -> (Vec<BflSegment>, Vec<f64>) {
        let mut segs = vec![BflSegment::zero(0); n5];
        let mut requant = vec![0.0f64; n5 * IDIM];
        for (j, seg) in segs.iter_mut().enumerate() {
            let vec: Vec<f64> = (0..IDIM).map(|m| samples[j * IDIM + m]).collect();
            let peak = vec.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
            let nls = if peak <= 0.0 {
                15
            } else {
                let mut n = 0;
                let mut p = peak;
                while p < 8192.0 && n < 15 {
                    p *= 2.0;
                    n += 1;
                }
                while p >= 16384.0 {
                    p /= 2.0;
                    n -= 1;
                }
                n
            };
            let scale = (2.0f64).powi(nls);
            let mut samp = [0i16; IDIM];
            for m in 0..IDIM {
                let q = (vec[m] * scale).round().clamp(-32768.0, 32767.0) as i16;
                samp[m] = q;
                requant[j * IDIM + m] = q as f64 / scale;
            }
            *seg = BflSegment { samples: samp, nls };
        }
        (segs, requant)
    }

    #[test]
    fn broadband_input_eventually_clears_illcond() {
        // A broadband (varied-sine) excitation is well conditioned: once
        // the recursive history fills, HWMCORE produces a non-zero
        // R(order+1) (ILLCOND clears).
        let w = synth_win();
        let mut s = HybridWindowFixedState::new(&w);
        let mut cleared = false;
        for cyc in 0..16 {
            let mut samples = [0.0f64; NFRSZ];
            for (i, v) in samples.iter_mut().enumerate() {
                let t = (cyc * NFRSZ + i) as f64;
                *v = (t * 0.21).sin() * 7000.0 + (t * 0.07).cos() * 3000.0;
            }
            let (segs, _) = sbfl_segments_requant(&samples, w.n5());
            let out = s.run(&w, &segs);
            if !out.illcond {
                cleared = true;
                assert!(out.rtmp[0] != 0, "R(1) must be non-zero on a good cycle");
            }
        }
        assert!(cleared, "well-conditioned input must clear ILLCOND");
    }
}
