//! Annex G §G.3.1 – §G.3.3 — the fixed-point encoder filter loop:
//! block 4 (perceptual weighting filter), "blockzir" (blocks 9 and 10
//! during zero-input-response computation) and blocks 9/10 during
//! filter memory update.
//!
//! Together with the §G.3.4 – §G.3.10 codebook search
//! ([`crate::annex_g_codebook`]) and the backward adapters
//! ([`crate::annex_g_synth_adapter`], [`crate::annex_g_weight_adapter`],
//! [`crate::annex_g_gain_adapter`]) these close the fixed-point
//! encoder's per-vector dataflow (Figure 2/G.728):
//!
//! ```text
//! S (Q2) ──block 4──► SW ──┐
//!                          ├─ block 11 target … blocks 17/18 search
//! blockzir ────────► ZIR ──┘
//! blocks 19/21 ► ET ──blocks 9/10 update──► ST (BFL) + refreshed
//!                                            STATELPC / ZIRW* memories
//! ```
//!
//! Fixed-point specifics (all from the §G.3.1 – §G.3.3 pseudo-code):
//!
//! * The input speech `S`, the weighting-filter memories `WFIR` /
//!   `WIIR` and the ZIR memories `ZIRWFIR` / `ZIRWIIR` are flat `Q2`
//!   (`NLSS = 2`); `A` / `AWZ` / `AWP` are `Q14`.
//! * `STATELPC` is **segmented block floating point**: ten
//!   `IDIM`-sample segments sharing `NLSSTATE(1..=10)` (newest segment
//!   = `NLSSTATE(10)`), plus the scratch slot `NLSSTATE(11)` holding
//!   the common minimum during the ZIR pass. Initial value 16.
//! * The memory update's LPC zero-state response carries the `AA1 =
//!   AA0 << 3` overflow probe with the `LABEL1` "shift `ET` right one
//!   bit and retry" loop, guaranteeing the response fits 15 bits
//!   before it is added to the 15-bit `STATELPC`.
//! * The updated `STATELPC` is clipped at `±4095` *at segment scale*
//!   (the `AA1 = 4095 << NLSSTATE(10)` alignment) exactly like the
//!   §G.3.11 decoder filter, keeping encoder and decoder quantized
//!   speech in lockstep.

use crate::consts::{IDIM, LPC, LPCW};

/// Number of `IDIM`-sample segments in `STATELPC` (`LPC / IDIM = 10`).
pub const NSEG_LPC: usize = LPC / IDIM;

/// Table G.2/G.728 initial value for every `NLSSTATE` slot.
pub const NLSSTATE_INIT: i32 = 16;

/// The quantized-speech clip level (±4095 at `Q0`, i.e. the §G.3.3
/// "4095 = STATELPC clipping level" aligned to the segment scale).
pub const STATELPC_CLIP: i32 = 4095;

/// One §G.3.3 memory-update result: the reversed quantized-speech
/// vector `ST` (BFL mantissas + shared `NLSST`) plus the possibly
/// re-scaled excitation (`ET` is right-shifted by the `LABEL1` overflow
/// probe when the zero-state response would exceed 15 bits).
#[derive(Debug, Clone, Copy)]
pub struct MemoryUpdateOut {
    /// Quantized speech `ST(1..=IDIM)` (0-based), BFL mantissas.
    pub st: [i16; IDIM],
    /// `NLSST` — the shared shift of `st` (= `NLSSTATE(10)`).
    pub nlsst: i32,
    /// The (possibly right-shifted) excitation actually filtered.
    pub et: [i16; IDIM],
    /// The matching `NLSET` after any `LABEL1` shifts.
    pub nlset: i32,
}

/// The fixed-point encoder filter memories (§G.3.1 – §G.3.3): the
/// segmented-BFL synthesis-filter state `STATELPC` / `NLSSTATE`, the
/// block-4 weighting memories `WFIR` / `WIIR` and the blockzir
/// weighting memories `ZIRWFIR` / `ZIRWIIR` (all `Q2`).
#[derive(Debug, Clone)]
pub struct EncoderFiltersFixed {
    /// `STATELPC(1..=LPC)`, 0-based (`statelpc[j-1]` = `STATELPC(J)`,
    /// `STATELPC(1)` the newest sample). 16-bit mantissas.
    statelpc: [i16; LPC],
    /// `NLSSTATE(1..=11)`, 0-based. Slot 10 (spec 11) is the scratch
    /// common-minimum slot used during the ZIR pass.
    nlsstate: [i32; NSEG_LPC + 1],
    /// Block-4 all-zero memory `WFIR(1..=LPCW)` (Q2).
    wfir: [i16; LPCW],
    /// Block-4 all-pole memory `WIIR(1..=LPCW)` (Q2).
    wiir: [i16; LPCW],
    /// Blockzir all-zero memory `ZIRWFIR(1..=LPCW)` (Q2). Doubles as
    /// the §G.3.3 zero-state-response scratch array.
    zirwfir: [i16; LPCW],
    /// Blockzir all-pole memory `ZIRWIIR(1..=LPCW)` (Q2).
    zirwiir: [i16; LPCW],
}

impl Default for EncoderFiltersFixed {
    fn default() -> Self {
        Self::new()
    }
}

impl EncoderFiltersFixed {
    /// Fresh state: all memories zero, `NLSSTATE` at its Table G.2
    /// initial value 16.
    #[must_use]
    pub fn new() -> Self {
        Self {
            statelpc: [0; LPC],
            nlsstate: [NLSSTATE_INIT; NSEG_LPC + 1],
            wfir: [0; LPCW],
            wiir: [0; LPCW],
            zirwfir: [0; LPCW],
            zirwiir: [0; LPCW],
        }
    }

    /// Borrow `STATELPC` (for tests/audit).
    #[must_use]
    pub fn statelpc(&self) -> &[i16; LPC] {
        &self.statelpc
    }

    /// Borrow `NLSSTATE` (for tests/audit).
    #[must_use]
    pub fn nlsstate(&self) -> &[i32; NSEG_LPC + 1] {
        &self.nlsstate
    }

    /// §G.3.1 block 4 — filter one `Q2` input-speech vector through the
    /// perceptual weighting filter `W(z)` (`AWZ` / `AWP` in `Q14`,
    /// 0-based with the unity head at index 0). Returns `SW` (`Q2`).
    pub fn weighting_filter(
        &mut self,
        s: &[i16; IDIM],
        awz: &[i16; LPCW + 1],
        awp: &[i16; LPCW + 1],
    ) -> [i16; IDIM] {
        let mut sw = [0i16; IDIM];
        for k in 0..IDIM {
            // | AA0 = S(K); AA0 = AA0 << 14
            let mut aa0 = i64::from(s[k]) << 14;
            // | For J = LPCW..2: AA0 += WFIR(J)·AWZ(J+1); WFIR(J) = WFIR(J−1)
            for j in (2..=LPCW).rev() {
                aa0 += i64::from(self.wfir[j - 1]) * i64::from(awz[j]);
                self.wfir[j - 1] = self.wfir[j - 2];
            }
            // | AA0 += WFIR(1)·AWZ(2); WFIR(1) = S(K)
            aa0 += i64::from(self.wfir[0]) * i64::from(awz[1]);
            self.wfir[0] = s[k];
            // | For J = LPCW..2: AA0 −= WIIR(J)·AWP(J+1); WIIR(J) = WIIR(J−1)
            for j in (2..=LPCW).rev() {
                aa0 -= i64::from(self.wiir[j - 1]) * i64::from(awp[j]);
                self.wiir[j - 1] = self.wiir[j - 2];
            }
            // | AA0 −= WIIR(1)·AWP(2); AA0 = AA0 >> 14; clip to 16 bits
            aa0 -= i64::from(self.wiir[0]) * i64::from(awp[1]);
            let out = clip_i16_64(aa0 >> 14);
            // | WIIR(1) = AA0; SW(K) = AA0     | SW is Q2
            self.wiir[0] = out;
            sw[k] = out;
        }
        sw
    }

    /// §G.3.2 "blockzir" — compute the zero-input response of the
    /// synthesis filter (block 9, on the segmented-BFL `STATELPC`) and
    /// of the weighting filter (block 10, on the `Q2` `ZIRW*`
    /// memories). `a` / `awz` / `awp` are the `Q14` coefficient arrays
    /// (0-based, unity heads). Returns `ZIR` (`Q2`).
    pub fn zero_input_response(
        &mut self,
        a: &[i16; LPC + 1],
        awz: &[i16; LPCW + 1],
        awp: &[i16; LPCW + 1],
    ) -> [i16; IDIM] {
        // | NLSSTATE(11) = min{NLSSTATE(1..10)}
        let min_nls = *self.nlsstate[..NSEG_LPC].iter().min().expect("NSEG > 0");
        self.nlsstate[NSEG_LPC] = min_nls;

        let mut temp = [0i16; IDIM];
        for k in 1..=IDIM {
            // | I = 1; L = 6 − K; J = LPC; AA0 = 0
            let mut j = LPC;
            let mut aa0: i64 = 0;
            let l = IDIM + 1 - k;
            // | For LL = 1..L: AA0 −= STATELPC(J)·A(J+1); shift; J −= 1
            for _ll in 0..l {
                aa0 -= i64::from(self.statelpc[j - 1]) * i64::from(a[j]);
                self.statelpc[j - 1] = self.statelpc[j - 2];
                j -= 1;
            }
            // | NLS = NLSSTATE(1) − NLSSTATE(11); AA1 = AA0 >> NLS
            let nls = self.nlsstate[0] - min_nls;
            let mut aa1 = aa0 >> nls;

            // | For I = 2..10: five-tap partial sums, aligned per segment.
            // `aa0` is zeroed **once** here, before the inner `LL` loop,
            // then accumulates (`aa0 −=`) across the five taps. The
            // §G.3.2 fixed-point `Blockzir` block-9 listing mis-prints the
            // inner MAC as `AA0 = 0 − …` (a reset each tap, keeping only
            // the last); the accumulating form matches the first/third MAC
            // loops and is required for conformance — erratum E6 in
            // `docs/audio/g728/g728-errata.md`.
            for i in 2..=NSEG_LPC {
                let mut aa0: i64 = 0;
                for _ll in 0..IDIM {
                    aa0 -= i64::from(self.statelpc[j - 1]) * i64::from(a[j]);
                    // | STATELPC(0) = garbage if J = 1; it is OK (the
                    // | slot is overwritten below) — skip the shift.
                    if j >= 2 {
                        self.statelpc[j - 1] = self.statelpc[j - 2];
                    }
                    j -= 1;
                }
                let nls = self.nlsstate[i - 1] - min_nls;
                aa1 += aa0 >> nls;
            }

            // | If K = 1, go to SHIFT2; else L = K − 1 more taps at the
            // | scratch scale (no alignment shift needed).
            if k != 1 {
                let mut aa0: i64 = 0;
                for _ll in 0..(k - 1) {
                    aa0 -= i64::from(self.statelpc[j - 1]) * i64::from(a[j]);
                    if j >= 2 {
                        self.statelpc[j - 1] = self.statelpc[j - 2];
                    }
                    j -= 1;
                }
                aa1 += aa0;
            }

            // SHIFT2:
            // | AA1 = AA1 >> 14; clip to 16 bits; STATELPC(1) = AA1
            let out = clip_i16_64(aa1 >> 14);
            self.statelpc[0] = out;
            // | IR = NLSSTATE(11) − 2; TEMP(K) = AA1 aligned to Q2
            let ir = min_nls - 2;
            temp[k - 1] = clip_i16_64(shr_signed(i64::from(out), ir));
        }

        // | Call VSCALE(STATELPC, IDIM, IDIM, 13, STATELPC, NLS)
        // | NLSSTATE(11) = NLSSTATE(11) + NLS
        let nls = self.vscale_top(13);
        self.nlsstate[NSEG_LPC] += nls;

        // | For L = 1..10: NLSSTATE(L) = NLSSTATE(L+1)
        for l in 0..NSEG_LPC {
            self.nlsstate[l] = self.nlsstate[l + 1];
        }

        // ---- Block 10 during ZIR (all Q2) ---------------------------
        let mut zir = [0i16; IDIM];
        for k in 0..IDIM {
            // | AA0 = TEMP(K) << 14
            let mut aa0 = i64::from(temp[k]) << 14;
            // | For J = LPCW..2: AA0 += ZIRWFIR(J)·AWZ(J+1); shift
            for j in (2..=LPCW).rev() {
                aa0 += i64::from(self.zirwfir[j - 1]) * i64::from(awz[j]);
                self.zirwfir[j - 1] = self.zirwfir[j - 2];
            }
            // | AA0 += ZIRWFIR(1)·AWZ(2); ZIRWFIR(1) = TEMP(K)
            aa0 += i64::from(self.zirwfir[0]) * i64::from(awz[1]);
            self.zirwfir[0] = temp[k];
            // | For J = LPCW..2: AA0 −= ZIRWIIR(J)·AWP(J+1); shift
            for j in (2..=LPCW).rev() {
                aa0 -= i64::from(self.zirwiir[j - 1]) * i64::from(awp[j]);
                self.zirwiir[j - 1] = self.zirwiir[j - 2];
            }
            // | AA0 −= ZIRWIIR(1)·AWP(2); AA0 >>= 14; clip
            aa0 -= i64::from(self.zirwiir[0]) * i64::from(awp[1]);
            let out = clip_i16_64(aa0 >> 14);
            // | ZIR(K) = AA0; ZIRWIIR(1) = AA0
            zir[k] = out;
            self.zirwiir[0] = out;
        }
        zir
    }

    /// §G.3.3 blocks 9/10 — filter-memory update: add the zero-state
    /// responses of the gain-scaled excitation `ET` (BFL mantissas +
    /// `nlset`) to the zero-input responses left in `STATELPC` /
    /// `ZIRWIIR` by [`zero_input_response`](Self::zero_input_response),
    /// clip the quantized speech at ±4095, re-normalize, and emit the
    /// reversed `ST` / `NLSST`.
    pub fn memory_update(
        &mut self,
        et: &[i16; IDIM],
        nlset: i32,
        a: &[i16; LPC + 1],
        awz: &[i16; LPCW + 1],
        awp: &[i16; LPCW + 1],
    ) -> MemoryUpdateOut {
        let mut et = *et;
        let mut nlset = nlset;

        // ---- LPC zero-state response into ZIRWFIR (scratch) ---------
        // LABEL1 retry: if the ZSR would exceed 15 bits, ET >>= 1 and
        // repeat (§G.3.3: "empirically the process repeats at most 3
        // times").
        'label1: loop {
            // | ZIRWFIR(1) = ET(1)
            self.zirwfir[0] = et[0];
            for k in 2..=IDIM {
                // | AA0 = ET(K) << 14        | A(1) = 16384
                let mut aa0 = i64::from(et[k - 1]) << 14;
                // | For I = K..2: ZIRWFIR(I) = ZIRWFIR(I−1);
                // |               AA0 −= A(I)·ZIRWFIR(I)
                for i in (2..=k).rev() {
                    self.zirwfir[i - 1] = self.zirwfir[i - 2];
                    aa0 -= i64::from(a[i - 1]) * i64::from(self.zirwfir[i - 1]);
                }
                // | AA1 = AA0 << 3; on overflow: ET >>= 1; NLSET −= 1;
                // | go to LABEL1
                let aa1 = aa0 << 3;
                if aa1 > i64::from(i32::MAX) || aa1 < i64::from(i32::MIN) {
                    for e in et.iter_mut() {
                        *e >>= 1;
                    }
                    nlset -= 1;
                    continue 'label1;
                }
                // | AA0 = AA0 >> 14; ZIRWFIR(1) = AA0 (lowest 16 bits)
                self.zirwfir[0] = clip_i16_64(aa0 >> 14);
            }
            break;
        }

        // ---- Weighting-filter zero-state response into TEMP ---------
        // | N = IDIM + 1; TEMP(1) = ZIRWFIR(IDIM)
        let mut temp = [0i16; IDIM];
        temp[0] = self.zirwfir[IDIM - 1];
        for k in 2..=IDIM {
            // | AA1 = ZIRWFIR(N − K) << 14    | AWZ(1) = 16384
            let mut aa1 = i64::from(self.zirwfir[IDIM - k]) << 14;
            // | M = IDIM − K
            let m = IDIM - k;
            // | For I = K..2: TEMP(I) = TEMP(I−1);
            // |   AA1 += AWZ(I)·ZIRWFIR(I+M); AA1 −= AWP(I)·TEMP(I)
            for i in (2..=k).rev() {
                temp[i - 1] = temp[i - 2];
                aa1 += i64::from(awz[i - 1]) * i64::from(self.zirwfir[i + m - 1]);
                aa1 -= i64::from(awp[i - 1]) * i64::from(temp[i - 1]);
            }
            // | AA1 >>= 14; clip; TEMP(1) = AA1
            temp[0] = clip_i16_64(aa1 >> 14);
        }

        // | IR = NLSET − 2; shift TEMP to Q2 like ZIRWIIR
        let ir = nlset - 2;
        for t in temp.iter_mut() {
            *t = clip_i16_64(shr_signed(i64::from(*t), ir));
        }

        // ---- Match the NLS of ZIRWFIR (ZSR) and STATELPC (ZIR) ------
        // | If NLSET = NLSSTATE(10), go to LABEL2
        // | If NLSET < NLSSTATE(10): STATELPC >>= NLSD; NLSSTATE(10) = NLSET
        // | Else: ZIRWFIR >>= NLSD
        let nls_top = self.nlsstate[NSEG_LPC - 1];
        match nlset.cmp(&nls_top) {
            std::cmp::Ordering::Less => {
                let nlsd = (nls_top - nlset) as u32;
                for k in 0..IDIM {
                    self.statelpc[k] >>= nlsd.min(15);
                }
                self.nlsstate[NSEG_LPC - 1] = nlset;
            }
            std::cmp::Ordering::Greater => {
                let nlsd = (nlset - nls_top) as u32;
                for k in 0..IDIM {
                    self.zirwfir[k] >>= nlsd.min(15);
                }
            }
            std::cmp::Ordering::Equal => {}
        }

        // LABEL2:
        // | AA1 = 4095 aligned to the STATELPC segment scale
        let nls_top = self.nlsstate[NSEG_LPC - 1];
        let clip = shr_signed(i64::from(STATELPC_CLIP), -nls_top);

        for k in 0..IDIM {
            // | AA0 = STATELPC(K) + ZIRWFIR(K); clip at ±AA1, then to
            // | 16 bits; STATELPC(K) = AA0
            let mut aa0 = i64::from(self.statelpc[k]) + i64::from(self.zirwfir[k]);
            aa0 = aa0.clamp(-clip, clip);
            self.statelpc[k] = clip_i16_64(aa0);
            // | AA0 = ZIRWIIR(K) + TEMP(K); clip to 16 bits
            let aa0 = i64::from(self.zirwiir[k]) + i64::from(temp[k]);
            self.zirwiir[k] = clip_i16_64(aa0);
        }

        // | Call VSCALE(STATELPC, IDIM, IDIM, 12, STATELPC, NLS)
        // | NLSSTATE(10) = NLSSTATE(10) + NLS
        let nls = self.vscale_top(12);
        self.nlsstate[NSEG_LPC - 1] += nls;

        // | Set ZIRWFIR(1..10) to the Q2 rendering of STATELPC(1..10)
        let ir = self.nlsstate[NSEG_LPC - 1] - 2;
        for i in 0..IDIM {
            self.zirwfir[i] = clip_i16_64(shr_signed(i64::from(self.statelpc[i]), ir));
        }
        let ir = self.nlsstate[NSEG_LPC - 2] - 2;
        for i in IDIM..LPCW {
            self.zirwfir[i] = clip_i16_64(shr_signed(i64::from(self.statelpc[i]), ir));
        }

        // | ST(K) = STATELPC(IDIM + 1 − K); NLSST = NLSSTATE(10)
        let mut st = [0i16; IDIM];
        for k in 0..IDIM {
            st[k] = self.statelpc[IDIM - 1 - k];
        }
        MemoryUpdateOut {
            st,
            nlsst: self.nlsstate[NSEG_LPC - 1],
            et,
            nlset,
        }
    }

    /// The §G.1.3 `VSCALE(STATELPC, IDIM, IDIM, mls, …)` calls of
    /// blockzir / memory update: block-normalize the top `IDIM`
    /// entries of `STATELPC` in place and return the applied `NLS`.
    fn vscale_top(&mut self, mls: i32) -> i32 {
        let vals: [i32; IDIM] = std::array::from_fn(|i| i32::from(self.statelpc[i]));
        let (scaled, nls) = crate::annex_g_arith::vscale(&vals, IDIM, mls);
        for (i, &v) in scaled.iter().enumerate() {
            self.statelpc[i] = clip_i16_64(i64::from(v));
        }
        nls
    }
}

// ---------------------------------------------------------------------
// Shared clip / shift helpers.
// ---------------------------------------------------------------------

/// Clip a 64-bit accumulator to the signed 16-bit range (the §G.3
/// "If AA0 > 32767 / < −32768" saturation lines).
#[inline]
fn clip_i16_64(v: i64) -> i16 {
    v.clamp(i64::from(i16::MIN), i64::from(i16::MAX)) as i16
}

/// Right shift by a possibly-negative count (negative ⇒ left shift) —
/// the §G.3 `If IR > 0 … If IR < 0` alignment idiom.
#[inline]
fn shr_signed(v: i64, ir: i32) -> i64 {
    if ir >= 0 {
        v >> ir
    } else {
        v << (-ir)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Q14 unity-head coefficient array builders.
    fn q14_lpc(taps: &[f64]) -> [i16; LPC + 1] {
        let mut a = [0i16; LPC + 1];
        a[0] = 1 << 14;
        for (i, &t) in taps.iter().enumerate() {
            a[i + 1] = (t * 16384.0).round() as i16;
        }
        a
    }
    fn q14_w(taps: &[f64]) -> [i16; LPCW + 1] {
        let mut a = [0i16; LPCW + 1];
        a[0] = 1 << 14;
        for (i, &t) in taps.iter().enumerate() {
            a[i + 1] = (t * 16384.0).round() as i16;
        }
        a
    }

    /// Floating-point transcription of the §G.3.1 block-4 pseudo-code
    /// (the annex prints it right above the fixed-point version).
    struct FloatWeighting {
        wfir: [f64; LPCW],
        wiir: [f64; LPCW],
    }
    impl FloatWeighting {
        fn new() -> Self {
            Self {
                wfir: [0.0; LPCW],
                wiir: [0.0; LPCW],
            }
        }
        fn run(
            &mut self,
            s: &[f64; IDIM],
            awz: &[f64; LPCW + 1],
            awp: &[f64; LPCW + 1],
        ) -> [f64; IDIM] {
            let mut sw = [0.0; IDIM];
            for k in 0..IDIM {
                let mut acc = s[k];
                for j in (2..=LPCW).rev() {
                    acc += self.wfir[j - 1] * awz[j];
                    self.wfir[j - 1] = self.wfir[j - 2];
                }
                acc += self.wfir[0] * awz[1];
                self.wfir[0] = s[k];
                for j in (2..=LPCW).rev() {
                    acc -= self.wiir[j - 1] * awp[j];
                    self.wiir[j - 1] = self.wiir[j - 2];
                }
                acc -= self.wiir[0] * awp[1];
                self.wiir[0] = acc;
                sw[k] = acc;
            }
            sw
        }
    }

    #[test]
    fn weighting_filter_tracks_float_pseudocode() {
        // Block 4 fixed vs the float pseudo-code on identical
        // requantized inputs / Q14 coefficients.
        let wz = [-0.9f64, 0.3, -0.1, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let wp = [-0.6f64, 0.15, -0.04, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let awz = q14_w(&wz);
        let awp = q14_w(&wp);
        let awz_f: [f64; LPCW + 1] = std::array::from_fn(|i| f64::from(awz[i]) / 16384.0);
        let awp_f: [f64; LPCW + 1] = std::array::from_fn(|i| f64::from(awp[i]) / 16384.0);

        let mut fixed = EncoderFiltersFixed::new();
        let mut float = FloatWeighting::new();
        let mut seed = 0x0123_4567_89ab_cdefu64;
        for _ in 0..20 {
            let mut s = [0i16; IDIM];
            let mut sf = [0.0f64; IDIM];
            for k in 0..IDIM {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                let v = ((seed >> 40) as i32 % 16384) - 8192; // Q2 of ±2048
                s[k] = v as i16;
                sf[k] = f64::from(v);
            }
            let sw = fixed.weighting_filter(&s, &awz, &awp);
            let swf = float.run(&sf, &awz_f, &awp_f);
            for k in 0..IDIM {
                assert!(
                    (f64::from(sw[k]) - swf[k]).abs() <= 4.0,
                    "SW[{k}] fixed {} vs float {}",
                    sw[k],
                    swf[k]
                );
            }
        }
    }

    #[test]
    fn identity_cascade_zir_is_zero_and_et_passes_through() {
        // With the all-pass cascade (A = AWZ = AWP = δ) and zero state,
        // the ZIR is zero and the memory update writes ET straight into
        // STATELPC — so ST is ET reversed at the same scale.
        let a = q14_lpc(&[]);
        let awz = q14_w(&[]);
        let awp = q14_w(&[]);
        let mut f = EncoderFiltersFixed::new();

        let zir = f.zero_input_response(&a, &awz, &awp);
        assert_eq!(zir, [0i16; IDIM]);

        let et = [1000i16, -2000, 3000, -4000, 5000];
        let nlset = 6;
        let out = f.memory_update(&et, nlset, &a, &awz, &awp);
        // ST is the reversed excitation, possibly re-normalized: compare
        // as real values.
        for k in 0..IDIM {
            let st_val = f64::from(out.st[k]) / (2.0f64).powi(out.nlsst);
            let et_val = f64::from(et[k]) / (2.0f64).powi(nlset);
            assert!(
                (st_val - et_val).abs() < 0.51,
                "ST[{k}] {st_val} vs ET {et_val}"
            );
        }
    }

    /// Floating-point reference for the whole ZIR + update cycle: a
    /// plain direct-form synthesis filter 1/A(z) and weighting filter
    /// W(z) run on the same excitation (the float §G.3.2 / §G.3.3
    /// pseudo-code collapses to exactly this).
    struct FloatLoop {
        lpc_mem: [f64; LPC], // y(n−1) … y(n−LPC)
        wfir: [f64; LPCW],   // W(z) all-zero memory (past y)
        wiir: [f64; LPCW],   // W(z) all-pole memory
    }
    impl FloatLoop {
        fn new() -> Self {
            Self {
                lpc_mem: [0.0; LPC],
                wfir: [0.0; LPCW],
                wiir: [0.0; LPCW],
            }
        }
        /// One vector: returns (zir_weighted, st) for excitation `et`.
        fn run(
            &mut self,
            et: &[f64; IDIM],
            a: &[f64; LPC + 1],
            awz: &[f64; LPCW + 1],
            awp: &[f64; LPCW + 1],
        ) -> ([f64; IDIM], [f64; IDIM]) {
            // ZIR of the cascade (zero input, live memories — clones).
            let mut lpc_mem = self.lpc_mem;
            let mut wfir = self.wfir;
            let mut wiir = self.wiir;
            let mut zir = [0.0; IDIM];
            for k in 0..IDIM {
                let mut y = 0.0;
                for j in 2..=(LPC + 1) {
                    y -= a[j - 1] * lpc_mem[j - 2];
                }
                for j in (1..LPC).rev() {
                    lpc_mem[j] = lpc_mem[j - 1];
                }
                lpc_mem[0] = y;
                let mut w = y;
                for j in 2..=(LPCW + 1) {
                    w += awz[j - 1] * wfir[j - 2];
                }
                for j in (1..LPCW).rev() {
                    wfir[j] = wfir[j - 1];
                }
                wfir[0] = y;
                for j in 2..=(LPCW + 1) {
                    w -= awp[j - 1] * wiir[j - 2];
                }
                for j in (1..LPCW).rev() {
                    wiir[j] = wiir[j - 1];
                }
                wiir[0] = w;
                zir[k] = w;
            }
            // Full response (ZIR + ZSR) with the real memories.
            let mut st = [0.0; IDIM];
            for k in 0..IDIM {
                let mut y = et[k];
                for j in 2..=(LPC + 1) {
                    y -= a[j - 1] * self.lpc_mem[j - 2];
                }
                for j in (1..LPC).rev() {
                    self.lpc_mem[j] = self.lpc_mem[j - 1];
                }
                self.lpc_mem[0] = y;
                st[k] = y;
                let mut w = y;
                for j in 2..=(LPCW + 1) {
                    w += awz[j - 1] * self.wfir[j - 2];
                }
                for j in (1..LPCW).rev() {
                    self.wfir[j] = self.wfir[j - 1];
                }
                self.wfir[0] = y;
                for j in 2..=(LPCW + 1) {
                    w -= awp[j - 1] * self.wiir[j - 2];
                }
                for j in (1..LPCW).rev() {
                    self.wiir[j] = self.wiir[j - 1];
                }
                self.wiir[0] = w;
            }
            (zir, st)
        }
    }

    #[test]
    fn zir_and_update_track_float_cascade() {
        // Drive the fixed blockzir + memory update against the float
        // direct-form cascade for many vectors with a stable synthesis
        // filter and realistic weighting coefficients; the quantized
        // speech ST and the ZIR must track within Q2-rounding noise.
        let lpc_taps = [-0.9f64, 0.2, -0.05, 0.01];
        let wz = [-0.8f64, 0.16, -0.03];
        let wp = [-0.5f64, 0.1, -0.02];
        let a = q14_lpc(&lpc_taps);
        let awz = q14_w(&wz);
        let awp = q14_w(&wp);
        let a_f: [f64; LPC + 1] = std::array::from_fn(|i| f64::from(a[i]) / 16384.0);
        let awz_f: [f64; LPCW + 1] = std::array::from_fn(|i| f64::from(awz[i]) / 16384.0);
        let awp_f: [f64; LPCW + 1] = std::array::from_fn(|i| f64::from(awp[i]) / 16384.0);

        let mut fixed = EncoderFiltersFixed::new();
        let mut float = FloatLoop::new();

        let mut seed = 0xdead_beef_1234_5678u64;
        for n in 0..40 {
            let zir_x = fixed.zero_input_response(&a, &awz, &awp);

            // A random excitation vector, |e| < 512, expressed as BFL
            // with nlset chosen so mantissas are 15-bit-normalized.
            let mut ef = [0.0f64; IDIM];
            for e in ef.iter_mut() {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                *e = f64::from(((seed >> 40) as i32 % 1024) - 512);
            }
            let peak = ef.iter().fold(0.0f64, |m, &v| m.max(v.abs())).max(1.0);
            let mut nlset = 0i32;
            let mut p = peak;
            while p < 8192.0 && nlset < 14 {
                p *= 2.0;
                nlset += 1;
            }
            let mut et = [0i16; IDIM];
            for k in 0..IDIM {
                et[k] = (ef[k] * (2.0f64).powi(nlset)).round() as i16;
                ef[k] = f64::from(et[k]) / (2.0f64).powi(nlset); // requantized
            }

            let (zir_f, st_f) = float.run(&ef, &a_f, &awz_f, &awp_f);
            let up = fixed.memory_update(&et, nlset, &a, &awz, &awp);

            if n >= 2 {
                for k in 0..IDIM {
                    let zx = f64::from(zir_x[k]) / 4.0; // Q2 → real
                    assert!(
                        (zx - zir_f[k]).abs() <= 2.0,
                        "n={n} ZIR[{k}] fixed {zx} vs float {}",
                        zir_f[k]
                    );
                    let sx = f64::from(up.st[k]) / (2.0f64).powi(up.nlsst);
                    assert!(
                        (sx - st_f[k]).abs() <= 2.0,
                        "n={n} ST[{k}] fixed {sx} vs float {}",
                        st_f[k]
                    );
                }
            }
        }
    }

    #[test]
    fn fresh_state_matches_table_g2() {
        let f = EncoderFiltersFixed::new();
        assert!(f.statelpc().iter().all(|&v| v == 0));
        assert!(f.nlsstate().iter().all(|&v| v == NLSSTATE_INIT));
    }
}
