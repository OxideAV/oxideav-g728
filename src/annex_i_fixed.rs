//! Annex I §I.5 — fixed-point frame-erasure building blocks.
//!
//! The fixed-point counterparts of [`crate::excitation_extrapolation`] /
//! [`crate::gain_growth_limiter`], transcribed from the fixed-point
//! pseudo-code of `T-REC-G.728-199905-AnnI.pdf`:
//!
//! * **Block 31SF** (§I.5.2) — erasure flags and the `FESCALE`
//!   attenuation schedule ([`EtPastFixed::on_erased_cycle`]);
//! * **Block 31FE** (§I.5.3) — voiced/unvoiced excitation extrapolation
//!   over the fixed-`Q2` `ETPAST()` array
//!   ([`EtPastFixed::extrapolate`]);
//! * **Block 31E** (§I.5.4) — the unconditional `ETPAST()` update from
//!   the 15-bit BFL `ET` ([`EtPastFixed::push_et`]);
//! * **Block 97FE** (§I.5.9) — the log-gain of the extrapolated
//!   excitation via the example `LIN2DB` polynomial
//!   ([`gstate1_of_extrapolated_et_q9`] / [`lin2db_gstate_q9`]);
//! * **Block 98AF** (§I.5.10) — the post-erasure log-gain growth
//!   limiter ([`limit_log_gain_after_fe_q9`]).
//!
//! §I.5.9 explicitly waives bit-exactness for the erasure path
//! ("Because bit exact compatibility is not required for frame
//! erasures, we have not included a bit exact description of the
//! LIN2DB function… Any fixed-point implementation of LIN2DB with
//! reasonable accuracy can be used here"); the polynomial below is the
//! annex's own example code.
//!
//! Clean-room: every constant and rule here is transcribed from the
//! §I.5 fixed-point pseudo-code, Table I.1/I.2 and Appendix I.I of the
//! staged Annex I PDF. No reference C / external implementation was
//! consulted.

use crate::annex_g_arith::{divide, rnd, shl_sat, vscale};
use crate::consts::{IDIM, KPMAX};
use crate::excitation_extrapolation::SlideSource;

/// `VTH = 7022` in Q14 (Table I.1/G.728) — the lowered frame-erasure
/// voicing threshold (`PPFTH / 1.4`), compared against the Q14 `PTAP`.
pub const VTH_Q14: i16 = 7022;

/// `FEGAINMAX = 1024` in Q9 (Table I.1/G.728) — the +2 dB/vector
/// post-erasure log-gain growth ceiling.
pub const FEGAINMAX_Q9: i32 = 1024;

/// `DIMINV = 13107` in Q16 (Table I.1/G.728) — `1/IDIM = 0.2`.
pub const DIMINV_Q16: i32 = 13_107;

/// `VOICEDFEGAIN(0:4)` in Q15 (§I.5.2 / Table I.I.1/G.728).
pub const VOICED_FE_GAIN_Q15: [i16; 5] = [26_214, 26_214, 19_661, 13_107, 6_554];

/// `UNVOICEDFEGAIN(0:5)` in Q15 (§I.5.2 / Table I.I.2/G.728).
pub const UNVOICED_FE_GAIN_Q15: [i16; 6] = [32_767, 32_767, 26_214, 19_661, 13_107, 6_554];

/// §I.5.10 block 98AF — the fixed-point log-gain limiter after a frame
/// erasure. Replaces block 98 (§G.3.16):
///
/// ```text
/// If LOGGAIN > 14336, set LOGGAIN = 14336
/// If LOGGAIN < −16384, set LOGGAIN = −16384
/// TMP = LOGGAIN − OGAINDB
/// If AFTERFE > 0 and TMP > FEGAINMAX, LOGGAIN = OGAINDB + FEGAINMAX
/// ```
///
/// `ogaindb_q9` is the §I.5.1 `OGAINDB` (Q9, offset-removed — the
/// block-98AF output on good vectors, the block-97FE output while
/// concealing). With `afterfe == 0` this is exactly block 98.
#[must_use]
pub fn limit_log_gain_after_fe_q9(loggain: i32, ogaindb_q9: i32, afterfe: usize) -> i16 {
    let mut lim = loggain.clamp(
        crate::annex_g_gain_adapter::LOGGAIN_MIN_Q9,
        crate::annex_g_gain_adapter::LOGGAIN_MAX_Q9,
    );
    if afterfe > 0 && lim - ogaindb_q9 > FEGAINMAX_Q9 {
        lim = ogaindb_q9 + FEGAINMAX_Q9;
    }
    lim as i16
}

// ---------------------------------------------------------------------
// Block 97FE (§I.5.9) — LIN2DB and the GSTATE(1) update.
// ---------------------------------------------------------------------

/// `LIN2DB` Taylor coefficients `c(0..=5)` in Q13 (§I.5.9): the log₂
/// expansion around 0.75.
const LIN2DB_C_Q13: [i32; 6] = [-3_400, 15_758, -10_505, 9_338, -9_338, 9_961];
/// `FAC1 = 24660` in Q13 (§I.5.9) — `10·log₁₀(2) = 3.0102999566`.
const LIN2DB_FAC1_Q13: i32 = 24_660;
/// `THRQTR = 24576` in Q15 (§I.5.9) — the 0.75 expansion point.
const LIN2DB_THRQTR_Q15: i32 = 24_576;

/// §I.5.9 example `LIN2DB`, with the block-42 gain-offset subtraction
/// folded in: converts a normalized linear energy `x · 2^(−nlsx)`
/// (mantissa `x ∈ [16384, 32767]`, i.e. a normalized value in
/// `[0.5, 1)` at Q15) to `10·log₁₀(·) − GOFF` in Q9 — the offset-removed
/// log-gain `GSTATE(1)` domain.
///
/// Note on the printed constant table: §I.5.9 lists
/// `GOFF = −16 384 = −32 in Q9` and the code line `AA1 = AA1 − GOFF`,
/// which taken literally *adds* 32 dB — inconsistent with the float
/// block 97FE (`GSTATE(1) = ETRMS − GOFF` with the base Recommendation's
/// `GOFF = +32 dB`) and with block 98AF's `[−16384, 14336]`
/// offset-removed clamp domain. The intent — subtract 32 dB — is
/// implemented here; §I.5.9 explicitly allows any reasonable `LIN2DB`.
#[must_use]
pub fn lin2db_gstate_q9(x: i16, nlsx: i32) -> i16 {
    // | AA0 = x                | normalized .5 <= x < 1, Q15
    // | AA0 = AA0 - THRQTR     | subtract .75 → [−.25, .25), Q15
    // | AA0 = AA0 << 1         | convert to Q16
    let aa0 = (i32::from(x) - LIN2DB_THRQTR_Q15) << 1;
    // | AA1 = 0
    // | For I = 5, 4, …, 1:
    // |   AA1 = AA1 + C(I)     | Q13
    // |   P = AA1 * AA0        | Q13 · Q16 = Q29
    // |   AA1 = P >> 16        | Q29 → Q13
    let mut aa1: i64 = 0;
    for i in (1..=5).rev() {
        aa1 += i64::from(LIN2DB_C_Q13[i]);
        let p = aa1 * i64::from(aa0);
        aa1 = p >> 16;
    }
    // | AA1 = AA1 + C(0)       | Q13
    // | AA1 = AA1 >> 3         | convert to Q10
    aa1 += i64::from(LIN2DB_C_Q13[0]);
    aa1 >>= 3;
    // | AA0 = 15 − NLSX        | exponent, Q0 → << 10 → Q10
    // | AA1 = AA0 + AA1        | log₂(x) in Q10
    aa1 += i64::from((15 - nlsx) << 10);
    // | AA1 = AA1 * FAC1       | Q10 · Q13 = Q23 → 10·log₁₀(x)
    // | AA1 = AA1 >> 14        | Q23 → Q9
    aa1 = (aa1 * i64::from(LIN2DB_FAC1_Q13)) >> 14;
    // Gain-offset subtraction (−32 dB in Q9; see the doc note above),
    // then the block-97 "lower limit" floor at −32 dB (§G.3.16 /
    // eq. G-9 — the pseudo-code has no upper clamp here; the maximum
    // Q2 excitation energy keeps the Q9 result well inside i16).
    aa1 -= 16_384;
    aa1.max(i64::from(crate::annex_g_gain_adapter::LOGGAIN_MIN_Q9)) as i16
}

/// §I.5.9 block 97FE — the offset-removed log-gain (Q9) of one
/// extrapolated excitation vector, from its fixed-`Q2` samples (the
/// `ETPAST(1..5)` scratch written by block 31FE):
///
/// ```text
/// AA0 = Σ ET(K)·ET(K)                  | ET is Q2 ⇒ AA0 is Q4
/// VSCALE(AA0, 30) ; ETRMS = AA0 >> 16 ; NLSETRMS = (4 + NLS) − 16
/// AA0 = ETRMS · DIMINV                 | DIMINV = 0.2 in Q16
/// VSCALE(AA0, 30) ; ETRMS = AA0 >> 16 ; NLSETRMS = NLSETRMS + NLS
/// If NLSETRMS > 14: ETRMS = 16384 ; NLSETRMS = 14   | floor at 1.0
/// GSTATE(1) = LIN2DB(ETRMS, NLSETRMS)
/// ```
#[must_use]
pub fn gstate1_of_extrapolated_et_q9(et_q2: &[i16; IDIM]) -> i16 {
    // | AA0 = Σ ET(K)² — Q2 samples ⇒ Q4 energy (fits 32 bits:
    // | 5 · 32768² < 2³¹ · 4).
    let mut aa0: i64 = 0;
    for &e in et_q2 {
        aa0 += i64::from(e) * i64::from(e);
    }
    let aa0 = aa0.clamp(0, i64::from(i32::MAX)) as i32;
    let (v, nls) = vscale(&[aa0], 1, 30);
    let etrms = v[0] >> 16;
    let mut nlsetrms = (4 + nls) - 16;

    let aa0 = etrms * DIMINV_Q16;
    let (v, nls) = vscale(&[aa0], 1, 30);
    let mut etrms = v[0] >> 16;
    nlsetrms += nls;

    if nlsetrms > 14 || etrms <= 0 {
        // ETRMS < 1 ⇒ clamp to exactly 1.0 (16384 at Q14) — the block-39
        // "check lower limit" floor.
        etrms = 16_384;
        nlsetrms = 14;
    }
    lin2db_gstate_q9(etrms as i16, nlsetrms)
}

// ---------------------------------------------------------------------
// Blocks 31SF / 31FE / 31E (§I.5.2 – §I.5.4) — fixed-Q2 ETPAST.
// ---------------------------------------------------------------------

/// `ETPAST(−139 ..= 5)` length: `KPMAX` past samples + the current
/// `IDIM`-sample scratch slot (Table I.2: array index range −139 to
/// IDIM). `ETPAST(k)` lives at `etpast[k + KPMAX − 1]`.
const ETPAST_LEN: usize = KPMAX + IDIM;
/// Offset of spec index 0 (the newest stored sample).
const ETPAST_ZERO: usize = KPMAX - 1;

/// Fixed-point frame-erasure excitation extrapolator: the `Q2`
/// `ETPAST()` array plus blocks 31SF / 31FE / 31E.
///
/// Mirrors [`crate::excitation_extrapolation::FrameErasureExcitation`]
/// but in the Annex G integer domain: `ETPAST` is fixed `Q2` ("there is
/// no point to use block floating-point for ETPAST", §I.5.3), the
/// extrapolated `ET` is returned in the 15-bit BFL form block 32
/// consumes, and the scales are the Q15 scalar-floating pairs of the
/// §I.5.2 pseudo-code.
#[derive(Debug, Clone)]
pub struct EtPastFixed {
    /// `ETPAST(−139 ..= 5)` in Q2 (initially zero, Table I.2).
    etpast: [i16; ETPAST_LEN],
    /// Last good `PTAP` (Q14, block 83) and `KP` (block 82).
    last_good_ptap_q14: i16,
    last_good_kp: usize,
    /// Erasure-onset latch.
    in_erasure: bool,
    /// `VOICED` — onset voicing decision.
    voiced: bool,
    /// `FEDELAY` — periodic slide-back (voiced) / last random slide.
    fedelay: usize,
    /// `AVMAG` scalar-floating pair (unvoiced onset magnitude).
    avmag: i16,
    nlsavmag: i32,
    /// `FESCALE` scalar-floating pair (current attenuation).
    fescale: i16,
    nlsfescale: i32,
}

impl Default for EtPastFixed {
    fn default() -> Self {
        Self::new()
    }
}

impl EtPastFixed {
    /// Fresh state: zeroed `ETPAST` (Table I.2 "initially zero").
    #[must_use]
    pub fn new() -> Self {
        Self {
            etpast: [0; ETPAST_LEN],
            last_good_ptap_q14: 0,
            last_good_kp: 0,
            in_erasure: false,
            voiced: false,
            fedelay: 0,
            avmag: 0,
            nlsavmag: 0,
            fescale: 0,
            nlsfescale: 15,
        }
    }

    /// Record the (Q14 `PTAP`, `KP`) pair of a good adaptation cycle and
    /// mark any erasure as ended (the next
    /// [`Self::on_erased_cycle`] re-latches the onset decision).
    pub fn observe_good_cycle(&mut self, ptap_q14: i16, kp: usize) {
        self.last_good_ptap_q14 = ptap_q14;
        self.last_good_kp = kp;
        self.in_erasure = false;
    }

    /// Block 31SF (§I.5.2) — set the erasure flags and the `FESCALE`
    /// scalar-floating attenuation. `n10msec = FECOUNT >> 2`.
    pub fn on_erased_cycle(&mut self, n10msec: usize) {
        if !self.in_erasure {
            self.in_erasure = true;
            // | If PTAP > VTH: FEDELAY = KP; VOICED = .TRUE.
            if self.last_good_ptap_q14 > VTH_Q14 {
                self.voiced = true;
                self.fedelay = self.last_good_kp;
            } else {
                // | VOICED = .FALSE.
                // | AA0 = Σ_{I=−39..0} |ETPAST(I)|   | 40 Q2 magnitudes
                // | VSCALE(AA0, 30); AVMAG = RND(AA0)
                // | NLSAVMAG = NLS − 16 + 3          | "+3" = "divided by 8"
                self.voiced = false;
                let mut aa0: i32 = 0;
                for k in (ETPAST_ZERO - 39)..=ETPAST_ZERO {
                    aa0 += i32::from(self.etpast[k]).abs();
                }
                let (v, nls) = vscale(&[aa0], 1, 30);
                self.avmag = rnd(v[0]);
                self.nlsavmag = nls - 16 + 3;
            }
        }

        if self.voiced {
            // | If N10MSEC < 5: FESCALE = VOICEDFEGAIN(N10MSEC) (Q15)
            self.fescale = VOICED_FE_GAIN_Q15.get(n10msec).copied().unwrap_or(0);
            self.nlsfescale = 15;
        } else if let Some(&g) = UNVOICED_FE_GAIN_Q15.get(n10msec) {
            // | FESCALE = UNVOICEDFEGAIN(N10MSEC)   | Q15
            // | AA0 = FESCALE * AVMAG ; VSCALE(AA0, 30)
            // | FESCALE = RND(AA0)
            // | NLSFESCALE = NLS + NLSAVMAG − 1     | −1 from (Q15 − 16)
            let aa0 = i32::from(g) * i32::from(self.avmag);
            let (v, nls) = vscale(&[aa0], 1, 30);
            self.fescale = rnd(v[0]);
            self.nlsfescale = nls + self.nlsavmag - 1;
        } else {
            self.fescale = 0;
            self.nlsfescale = 15;
        }
    }

    /// Block 31FE (§I.5.3) — extrapolate one excitation vector.
    ///
    /// Writes the extrapolated `Q2` samples into the `ETPAST(1..5)`
    /// scratch slot and returns them twice over: `(et, nlset, et_q2)` —
    /// the 15-bit BFL mantissas + shift for block 32, and the flat `Q2`
    /// copy consumed by block 97FE.
    pub fn extrapolate<S: SlideSource>(
        &mut self,
        slide: &mut S,
    ) -> ([i16; IDIM], i32, [i16; IDIM]) {
        let (temp, nlstemp);
        if self.voiced {
            // | TEMP = FESCALE ; NLSTEMP = 15       | FESCALE is Q15
            // | ETPAST(I) = ETPAST(I − FEDELAY)     | FEDELAY = pitch
            temp = self.fescale;
            nlstemp = 15;
            let d = self.fedelay.clamp(IDIM, KPMAX);
            for i in 1..=IDIM {
                self.etpast[ETPAST_ZERO + i] = self.etpast[ETPAST_ZERO + i - d];
            }
        } else {
            // | set FEDELAY = a random number between 5 and 140
            // | ETPAST(I) = ETPAST(I − FEDELAY) ; AA1 = Σ |ETPAST(1:5)|
            let d = slide.next_slide().clamp(IDIM, KPMAX);
            self.fedelay = d;
            let mut aa1: i32 = 0;
            for i in 1..=IDIM {
                let v = self.etpast[ETPAST_ZERO + i - d];
                self.etpast[ETPAST_ZERO + i] = v;
                aa1 += i32::from(v).abs();
            }
            if aa1 == 0 || self.fescale == 0 {
                // | TEMP = 0 ; NLSTEMP = 15
                temp = 0;
                nlstemp = 15;
            } else {
                // | VSCALE(AA1, 30) ; DEN = RND(AA1) ; NLSDEN = NLS − 16
                // | DIVIDE(FESCALE, NLSFESCALE, DEN, NLSDEN, TEMP, NLSTEMP)
                let (v, nls) = vscale(&[aa1], 1, 30);
                let den = rnd(v[0]);
                let nlsden = nls - 16;
                let (t, n) = divide(self.fescale, self.nlsfescale, den, nlsden);
                temp = t;
                nlstemp = n;
            }
        }

        // | For I = 1..IDIM:
        // |   AA0 = TEMP * ETPAST(I)     | Q(2 + NLSTEMP)
        // |   AA0 = AA0 >> NLSTEMP       | make Q2
        // |   clip to ±32767/−32768 ; ETPAST(I) = AA0
        let mut et_q2 = [0i16; IDIM];
        for i in 1..=IDIM {
            let p = i64::from(temp) * i64::from(self.etpast[ETPAST_ZERO + i]);
            let aa0 = if nlstemp >= 0 {
                p >> nlstemp
            } else {
                p << (-nlstemp).min(32)
            };
            let q2 = aa0.clamp(i64::from(i16::MIN), i64::from(i16::MAX)) as i16;
            self.etpast[ETPAST_ZERO + i] = q2;
            et_q2[i - 1] = q2;
        }

        // | VSCALE(ETPAST, IDIM, 13) → ET (15-bit BFL) ; NLSET = NLS + 2
        let vals: [i32; IDIM] = std::array::from_fn(|i| i32::from(et_q2[i]));
        let (v, nls) = vscale(&vals, IDIM, 13);
        let et: [i16; IDIM] = std::array::from_fn(|i| v[i] as i16);
        (et, nls + 2, et_q2)
    }

    /// Block 31E (§I.5.4) — shift `ETPAST()` by `IDIM` and append the
    /// current vector from its 15-bit BFL form. Must run on **every**
    /// vector, good or concealed.
    pub fn push_et(&mut self, et: &[i16; IDIM], nlset: i32) {
        // | For I = −KPMAX+1 .. −IDIM: ETPAST(I) = ETPAST(I + IDIM)
        self.etpast.copy_within(IDIM..=ETPAST_ZERO, 0);
        // | NRS = NLSET − 2 ; shift each ET mantissa to Q2, clip, store
        // | at ETPAST(I − IDIM) for I = 1..IDIM (i.e. −4..0).
        let nrs = nlset - 2;
        for i in 0..IDIM {
            let aa0 = i32::from(et[i]);
            let q2 = if nrs > 0 {
                aa0 >> nrs.min(31)
            } else {
                shl_sat(aa0, (-nrs) as u32)
            };
            self.etpast[ETPAST_ZERO - (IDIM - 1) + i] =
                q2.clamp(i32::from(i16::MIN), i32::from(i16::MAX)) as i16;
        }
    }

    /// Whether an erasure onset has been latched (for tests/audit).
    #[must_use]
    pub fn in_erasure(&self) -> bool {
        self.in_erasure
    }

    /// The onset `VOICED` decision (for tests/audit).
    #[must_use]
    pub fn voiced(&self) -> bool {
        self.voiced
    }

    /// The current `FESCALE` mantissa (0 ⇒ the §I.5.2 schedule has
    /// silenced the extrapolation) — for tests/audit.
    #[must_use]
    pub fn fescale(&self) -> i16 {
        self.fescale
    }

    /// Borrow the `Q2` `ETPAST` array (for tests/audit).
    #[must_use]
    pub fn etpast(&self) -> &[i16; ETPAST_LEN] {
        &self.etpast
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::{DIMINV, GOFF, VTH};

    #[test]
    fn table_i1_constants_match_float_values() {
        // Table I.1/G.728 lists both columns; cross-check the integer
        // renderings against the float constants.
        // 0.4285714 · 2¹⁴ = 7021.7 → the table prints 7022.
        assert_eq!(VTH_Q14, (VTH * 16384.0).round() as i16);
        assert_eq!(FEGAINMAX_Q9, (2.0f64 * 512.0) as i32);
        assert_eq!(DIMINV_Q16, (DIMINV * 65536.0).round() as i32);
        for (q, f) in VOICED_FE_GAIN_Q15
            .iter()
            .zip(crate::consts::VOICED_FE_GAIN.iter())
        {
            assert!((f64::from(*q) / 32768.0 - f).abs() < 2e-4);
        }
        for (q, f) in UNVOICED_FE_GAIN_Q15
            .iter()
            .zip(crate::consts::UNVOICED_FE_GAIN.iter())
        {
            assert!((f64::from(*q) / 32768.0 - f).abs() < 2e-4);
        }
    }

    #[test]
    fn block_98af_reduces_to_block_98_when_clamp_inactive() {
        assert_eq!(limit_log_gain_after_fe_q9(20_000, 0, 0), 14_336);
        assert_eq!(limit_log_gain_after_fe_q9(-20_000, 0, 0), -16_384);
        assert_eq!(limit_log_gain_after_fe_q9(1_234, -16_384, 0), 1_234);
    }

    #[test]
    fn block_98af_growth_clamp() {
        // OGAINDB = 0, LOGGAIN = 4096 (8 dB jump) with clamp active →
        // 0 + 1024 (= +2 dB).
        assert_eq!(limit_log_gain_after_fe_q9(4_096, 0, 3), 1_024);
        // Growth within budget passes through.
        assert_eq!(limit_log_gain_after_fe_q9(1_000, 0, 3), 1_000);
        // Decreases are never clamped.
        assert_eq!(limit_log_gain_after_fe_q9(-8_000, 0, 3), -8_000);
    }

    #[test]
    fn lin2db_tracks_float_log10() {
        // The example polynomial must track 10·log₁₀(x) − 32 within a
        // fraction of a dB across the mantissa range and several
        // exponents ("reasonable accuracy", §I.5.9).
        for nlsx in [8, 10, 14] {
            for xm in [16_384i16, 20_000, 24_576, 28_000, 32_000] {
                let x = f64::from(xm) / f64::powi(2.0, nlsx);
                let want = 10.0 * x.log10() - GOFF;
                let got = f64::from(lin2db_gstate_q9(xm, nlsx)) / 512.0;
                assert!(
                    (got - want).abs() < 0.05,
                    "x = {xm}·2^-{nlsx}: got {got} dB, want {want} dB"
                );
            }
        }
    }

    #[test]
    fn block_97fe_matches_float_gstate1() {
        // Fixed 97FE vs the float 97FE (ETRMS = 10·log₁₀(ΣET²/5),
        // floored at 1, minus 32 dB) on Q2 excitation vectors.
        for (name, q2) in [
            ("loud", [8_000i16, -6_000, 4_000, -2_000, 1_000]),
            ("medium", [400, -300, 200, -100, 50]),
            ("tiny", [4, -3, 2, -1, 1]),
            ("silence", [0, 0, 0, 0, 0]),
        ] {
            let e: f64 = q2
                .iter()
                .map(|&v| {
                    let s = f64::from(v) / 4.0;
                    s * s
                })
                .sum::<f64>()
                * DIMINV;
            let want = 10.0 * e.max(1.0).log10() - GOFF;
            let got = f64::from(gstate1_of_extrapolated_et_q9(&q2)) / 512.0;
            assert!(
                (got - want).abs() < 0.06,
                "{name}: got {got} dB, want {want} dB"
            );
        }
    }

    #[test]
    fn voiced_extrapolation_repeats_with_period_kp() {
        // Fill ETPAST with a ramp, latch a voiced onset (KP = 20) and
        // check the extrapolate→push loop repeats the last KP samples
        // scaled by FESCALE = 0.8 (Q15 26214) with period exactly KP.
        let mut fe = EtPastFixed::new();
        let mut val = 100i16;
        for _ in 0..30 {
            let mut et = [0i16; IDIM];
            for e in &mut et {
                *e = val;
                val += 4;
            }
            fe.push_et(&et, 2); // NLSET = 2 ⇒ mantissas ARE Q2
        }
        fe.observe_good_cycle(16_000, 20); // PTAP 0.977 ≫ VTH
        fe.on_erased_cycle(0);
        assert!(fe.voiced());
        assert_eq!(fe.fescale(), 26_214);

        let hist_tail: Vec<i16> = fe.etpast()[ETPAST_LEN - IDIM - 20..ETPAST_LEN - IDIM].to_vec();
        let mut slide = crate::excitation_extrapolation::LcgSlideSource::new(1);
        let mut stream = Vec::new();
        for _ in 0..8 {
            let (et, nlset, et_q2) = fe.extrapolate(&mut slide);
            fe.push_et(&et, nlset);
            stream.extend_from_slice(&et_q2);
        }
        for (n, &v) in stream.iter().enumerate() {
            let src = if n < 20 {
                f64::from(hist_tail[n])
            } else {
                f64::from(stream[n - 20])
            };
            let want = src * 0.8;
            assert!(
                (f64::from(v) - want).abs() <= 2.0,
                "sample {n}: {v} vs {want} (Q2 rounding budget)"
            );
        }
    }

    #[test]
    fn unvoiced_extrapolation_magnitude_matches_avmag() {
        // Unvoiced: each extrapolated 5-sample segment is scaled so its
        // magnitude sum ≈ FESCALE = AVMAG·UNVOICEDFEGAIN(0) — the
        // §I.4.1 magnitude match, here at UNVOICEDFEGAIN(0) = 1.0.
        struct Fixed(usize);
        impl SlideSource for Fixed {
            fn next_slide(&mut self) -> usize {
                self.0
            }
        }
        let mut fe = EtPastFixed::new();
        let mut seed = 0x1357_9bdfu32;
        for _ in 0..30 {
            let mut et = [0i16; IDIM];
            for e in &mut et {
                seed = seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                *e = ((seed >> 20) as i16 % 2_000) - 1_000;
            }
            fe.push_et(&et, 2);
        }
        fe.observe_good_cycle(0, 30); // PTAP = 0 → unvoiced
        fe.on_erased_cycle(0);
        assert!(!fe.voiced());

        // AVMAG (scalar float) ≈ Σ|last 40|/8.
        let manual: f64 = fe.etpast()[ETPAST_ZERO - 39..=ETPAST_ZERO]
            .iter()
            .map(|&v| f64::from(v).abs())
            .sum::<f64>()
            / 8.0;
        let mut slide = Fixed(60);
        let (_et, _nls, et_q2) = fe.extrapolate(&mut slide);
        let mag: f64 = et_q2.iter().map(|&v| f64::from(v).abs()).sum();
        assert!(
            (mag - manual).abs() / manual.max(1.0) < 0.01,
            "segment magnitude {mag} vs AVMAG {manual}"
        );
    }

    #[test]
    fn zero_history_extrapolates_silence() {
        let mut fe = EtPastFixed::new();
        fe.observe_good_cycle(0, 30);
        fe.on_erased_cycle(0);
        let mut slide = crate::excitation_extrapolation::LcgSlideSource::new(7);
        let (et, nlset, et_q2) = fe.extrapolate(&mut slide);
        assert_eq!(et_q2, [0; IDIM]);
        assert!(et.iter().all(|&v| v == 0));
        // Pushing the silent vector keeps the buffer finite/zero.
        fe.push_et(&et, nlset);
        assert!(fe.etpast().iter().all(|&v| v == 0));
    }
}
