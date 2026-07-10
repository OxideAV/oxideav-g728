//! Annex G §G.3.24 – §G.3.27 — the fixed-point long-term-postfilter
//! adaptation chain (blocks 81 / 82 / 83 / 84).
//!
//! These four decoder-side blocks derive the long-term (pitch)
//! postfilter parameters from the quantized speech, once per adaptation
//! cycle (at `ICOUNT = 3`):
//!
//! * **Block 81** (§G.3.24) — the 10th-order LPC inverse filter. Every
//!   vector it converts the BFL quantized-speech `ST` to the flat `Q2`
//!   postfilter form `SST` and FIR-filters it with the `Q13` predictor
//!   `APF` into the `Q1` prediction residual `D(IP+1 .. IP+IDIM)`.
//! * **Block 82** (§G.3.25) — pitch period extraction: 1 kHz low-pass +
//!   4:1 decimation of the newest `NFRSZ` residual samples, coarse
//!   correlation peak-picking over the decimated lag range, refinement
//!   `±3` in the undecimated domain, and the sub-multiple check against
//!   the previous frame's pitch `KP1 ± KPDELTA` (the fixed-point code
//!   compares `CMAX·SUM` against `CORMAX·TMP·ITAPTH` cross-multiplied,
//!   avoiding the two divisions of the floating-point version).
//! * **Block 83** (§G.3.26) — pitch predictor tap: the `Q14` ratio
//!   `PTAP = Σ SST(K)·SST(K−KP) / Σ SST(K−KP)²` over the last `NPWSZ`
//!   `Q0` quantized-speech samples, via the §G.1.3.4 `DIVIDE`.
//! * **Block 84** (§G.3.27) — long-term postfilter coefficients: the
//!   `PPFTH` on/off threshold, `B = PPFZCF·PTAP` (`Q16`),
//!   `GL = 1/(1+B)` (`Q14`, via `DIVIDE`) and the precomputed product
//!   `GLB = GL·B` (`Q16`) consumed by block 71
//!   ([`crate::annex_g_postfilter::PostfilterFixed`]).
//!
//! The §G.3.28 `LABEL` tail (re-normalizing the Levinson by-product
//! `APF` to `Q13` for block 81) lives here too ([`apf_to_q13`]), since
//! block 81 is its only consumer.
//!
//! Everything is transcribed from the §G.3.24 – §G.3.28 fixed-point
//! pseudo-code; buffer extents and initial values (`D(−139..100)` zero,
//! `DEC(−34..25)` zero, `IP = IPINIT = NPWSZ − NFRSZ + IDIM`,
//! `KP1 = 50`, `STLPCI` / `LPFFIR` / `LPFIIR` zero) come from Table
//! 2/G.728 and Table G.2/G.728.

use crate::annex_g_arith::{divide, findnls, rnd};
use crate::annex_g_postfilter::{LongTermCoeff, PF_ORDER};
use crate::consts::{IDIM, KPDELTA, KPMAX, KPMIN, NFRSZ, NPWSZ};

/// `D(k)` lives at `d[k + D_OFF]` for `k ∈ −(KPMAX−1) ..= NPWSZ`.
pub const D_OFF: usize = KPMAX - 1;
/// Length of the `Q1` LPC-prediction-residual buffer `D` (240).
pub const D_LEN: usize = D_OFF + NPWSZ + 1;
/// `DEC(n)` lives at `dec[n + DEC_OFF]` for `n ∈ −34 ..= NPWSZ/4`.
pub const DEC_OFF: usize = 34;
/// Length of the `Q1` decimated-residual buffer `DEC` (60).
pub const DEC_LEN: usize = DEC_OFF + NPWSZ / 4 + 1;

/// `IPINIT = NPWSZ − NFRSZ + IDIM` (Table 2/G.728 note b) — the initial
/// write pointer into `D`.
pub const IPINIT: usize = NPWSZ - NFRSZ + IDIM;

/// `KP1` initial value (Table 2/G.728): 50 samples.
pub const KP1_INIT: usize = 50;

/// §G.3.25 low-pass numerator `BL(0..=3)` in Q19 (the fixed-point
/// pseudo-code's inline constants; `BL(0) = BL(3)`, `BL(1) = BL(2)` —
/// the same symmetric 1 kHz FIR as the float Annex D `BL`).
pub const BL_Q19: [i32; 4] = [18721, -3668, -3668, 18721];
/// §G.3.25 low-pass denominator `AL(1..=3)` in Q13 (float Annex D `AL`
/// taps 2..4 × 2¹³).
pub const AL_Q13: [i32; 3] = [-19172, 16481, -5031];

/// `ITAPTH = 26214` in Q16 (Table G.1/G.728 `TAPTH = 0.4`) — the
/// fundamental-pitch replacement threshold used cross-multiplied in the
/// §G.3.25 sub-multiple decision.
pub const ITAPTH_Q16: i64 = 26214;

/// `PPFTH = 9830` in Q14 (Table G.1/G.728, tap threshold 0.6).
pub const PPFTH_Q14: i32 = 9830;
/// `PPFZCF = 9830` in Q16 (Table G.1/G.728, zero controlling factor
/// 0.15).
pub const PPFZCF_Q16: i32 = 9830;

/// §G.3.28 `LABEL` tail — re-normalize the order-10 Levinson by-product
/// `APF` (Q13/Q14/Q15 per `nlsatmp`) to the flat `Q13` block 81 consumes.
///
/// * `nlsatmp = 13`: already Q13 — copy through.
/// * `nlsatmp = 14`: `RND(APF(I) << 15)` (halve, round).
/// * `nlsatmp = 15`: `RND(APF(I) << 14)` (quarter, round).
///
/// `atil[0]` is the implicit unity head; the output's head is the exact
/// Q13 unity `8192`.
#[must_use]
pub fn apf_to_q13(atil: &[i16; PF_ORDER + 1], nlsatmp: i32) -> [i16; PF_ORDER + 1] {
    let mut out = [0i16; PF_ORDER + 1];
    out[0] = 1 << 13;
    let shift = match nlsatmp {
        13 => {
            out[1..].copy_from_slice(&atil[1..]);
            return out;
        }
        14 => 15,
        _ => 14, // nlsatmp = 15
    };
    for i in 1..=PF_ORDER {
        out[i] = rnd(i32::from(atil[i]) << shift);
    }
    out
}

/// The fixed-point long-term-postfilter adaptation chain (blocks 81 –
/// 84). Owns the `Q1` residual buffers, the low-pass state, the LPC
/// inverse-filter memory and the pitch bookkeeping (`KP` / `KP1`).
#[derive(Debug, Clone)]
pub struct PitchAdapterFixed {
    /// `D(−139 ..= 100)` — Q1 LPC prediction residual.
    d: [i16; D_LEN],
    /// `DEC(−34 ..= 25)` — Q1 decimated residual.
    dec: [i16; DEC_LEN],
    /// `LPFFIR(1..=3)` — Q1 low-pass FIR memory.
    lpffir: [i16; 3],
    /// `LPFIIR(1..=3)` — Q1 low-pass IIR memory.
    lpfiir: [i16; 3],
    /// `STLPCI(1..=10)` — Q2 LPC inverse-filter memory.
    stlpci: [i16; PF_ORDER],
    /// `IP` — the 1-based write pointer into `D`.
    ip: usize,
    /// `KP` — current pitch period (block 82 output).
    kp: usize,
    /// `KP1` — previous frame's pitch period.
    kp1: usize,
}

impl Default for PitchAdapterFixed {
    fn default() -> Self {
        Self::new()
    }
}

impl PitchAdapterFixed {
    /// Fresh state per Table 2/G.728: zero buffers, `IP = IPINIT`,
    /// `KP1 = 50`.
    #[must_use]
    pub fn new() -> Self {
        Self {
            d: [0; D_LEN],
            dec: [0; DEC_LEN],
            lpffir: [0; 3],
            lpfiir: [0; 3],
            stlpci: [0; PF_ORDER],
            ip: IPINIT,
            kp: KP1_INIT,
            kp1: KP1_INIT,
        }
    }

    /// The current pitch period `KP` (block 82's latest output).
    #[must_use]
    pub fn kp(&self) -> usize {
        self.kp
    }

    /// The block-81 residual write pointer `IP` (for tests/audit).
    ///
    /// Walks `IPINIT = 85 → 90 → 95 → 100`, then wraps to
    /// `NPWSZ − NFRSZ = 80` before the next write — the erratum-E3
    /// circular-buffer rewind (see [`Self::inverse_filter`]).
    #[must_use]
    pub fn ip(&self) -> usize {
        self.ip
    }

    /// Borrow the `Q1` residual buffer `D` (for tests/audit).
    #[must_use]
    pub fn d(&self) -> &[i16; D_LEN] {
        &self.d
    }

    /// §G.3.24 block 81 — convert the BFL quantized speech `ST`
    /// (mantissas + `nlsst`) to the flat `Q2` postfilter vector `SST`
    /// and run the 10th-order LPC inverse filter (`APF`, Q13) into the
    /// `Q1` residual `D`. Returns `SST(1..=IDIM)` (Q2) — the vector the
    /// caller hands to the §G.3.20 postfilter.
    pub fn inverse_filter(
        &mut self,
        st: &[i32; IDIM],
        nlsst: i32,
        apf_q13: &[i16; PF_ORDER + 1],
    ) -> [i32; IDIM] {
        // | NLS = 16 − NLSST + 2      | left shift for Q2
        // | For K: AA0 = ST(K) << NLS; SST(K) = RND(AA0)
        let nls = 16 - nlsst + 2;
        let mut sst = [0i32; IDIM];
        for k in 0..IDIM {
            let aa0 = shl_signed(i64::from(st[k]), nls);
            sst[k] = i32::from(rnd_i64(aa0));
        }

        // | If IP = NPWSZ, then set IP = NPWSZ − NFRSZ
        // (the float pseudo-code line; the D write pointer wraps inside
        // the newest-NFRSZ window 81..100 once the startup offset of
        // IPINIT = 85 has drained). The §G.3.24 *fixed-point* listing
        // mis-prints this reset target as `NFRSZ` alone; the paired
        // floating-point line (used here) is `NPWSZ − NFRSZ` and is what
        // the conformance vectors require — erratum E3 in
        // `docs/audio/g728/g728-errata.md`.
        if self.ip == NPWSZ {
            self.ip = NPWSZ - NFRSZ;
        }

        // | For K = 1..IDIM:
        // |   AA0 = SST(K); AA0 = AA0 << 13
        // |   For J = 10..2: AA0 += STLPCI(J)·APF(J+1); STLPCI(J) = STLPCI(J−1)
        // |   AA0 += STLPCI(1)·APF(2); STLPCI(1) = SST(K)
        // |   ITMP = IP + K; AA0 = AA0 << 2; D(ITMP) = RND(AA0)   | Q1
        for k in 0..IDIM {
            let mut aa0 = i64::from(sst[k]) << 13;
            for j in (2..=PF_ORDER).rev() {
                aa0 += i64::from(self.stlpci[j - 1]) * i64::from(apf_q13[j]);
                self.stlpci[j - 1] = self.stlpci[j - 2];
            }
            aa0 += i64::from(self.stlpci[0]) * i64::from(apf_q13[1]);
            self.stlpci[0] = clip_i16(sst[k]);
            let itmp = self.ip + k + 1; // ITMP = IP + K (1-based K)
            self.d[D_OFF + itmp] = rnd_i64(aa0 << 2);
        }

        // | IP = IP + IDIM
        self.ip += IDIM;

        sst
    }

    /// §G.3.25 block 82 — pitch period extraction. Call once per
    /// adaptation cycle (at `ICOUNT = 3`, after the cycle's residual
    /// has been written by block 81). Updates `KP` / `KP1` and shifts
    /// the `D` / `DEC` buffers.
    pub fn pitch_extract(&mut self) {
        // ---- 1 kHz low-pass + 4:1 decimation of the newest NFRSZ ----
        // | For K = NPWSZ − NFRSZ + 1 .. NPWSZ:
        // |   AA0 = D(K)·BL(0) + LPFFIR(1)·BL(1) + LPFFIR(2)·BL(2)
        // |       + LPFFIR(3)·BL(3)                | D Q1, BL Q19
        // |   shift LPFFIR; LPFFIR(1) = D(K)
        // |   AA0 = AA0 >> 6
        // |   AA0 −= LPFIIR(1)·AL(1) + LPFIIR(2)·AL(2) + LPFIIR(3)·AL(3)
        // |   shift LPFIIR; AA0 = AA0 << 3; LPFIIR(1) = RND(AA0)  | Q1
        // |   N = K >> 2; if K = N << 2, DEC(N) = LPFIIR(1)
        for k in (NPWSZ - NFRSZ + 1)..=NPWSZ {
            let dk = i64::from(self.d[D_OFF + k]);
            let mut aa0 = dk * i64::from(BL_Q19[0]);
            aa0 += i64::from(self.lpffir[0]) * i64::from(BL_Q19[1]);
            aa0 += i64::from(self.lpffir[1]) * i64::from(BL_Q19[2]);
            aa0 += i64::from(self.lpffir[2]) * i64::from(BL_Q19[3]);
            self.lpffir[2] = self.lpffir[1];
            self.lpffir[1] = self.lpffir[0];
            self.lpffir[0] = self.d[D_OFF + k];
            let mut aa0 = aa0 >> 6;
            aa0 -= i64::from(self.lpfiir[0]) * i64::from(AL_Q13[0]);
            aa0 -= i64::from(self.lpfiir[1]) * i64::from(AL_Q13[1]);
            aa0 -= i64::from(self.lpfiir[2]) * i64::from(AL_Q13[2]);
            self.lpfiir[2] = self.lpfiir[1];
            self.lpfiir[1] = self.lpfiir[0];
            self.lpfiir[0] = rnd_i64(aa0 << 3);
            if k % 4 == 0 {
                self.dec[DEC_OFF + k / 4] = self.lpfiir[0];
            }
        }

        // ---- Coarse peak-pick in the decimated domain ---------------
        // | M1 = KPMIN/4; M2 = KPMAX/4; AA1 = −2³¹
        let m1 = KPMIN / 4;
        let m2 = KPMAX / 4;
        let mut aa1 = i64::from(i32::MIN);
        let mut kmax = m1;
        for j in m1..=m2 {
            let mut aa0: i64 = 0;
            for n in 1..=(NPWSZ / 4) {
                aa0 += i64::from(self.dec[DEC_OFF + n]) * i64::from(self.dec[(DEC_OFF + n) - j]);
            }
            if aa0 > aa1 {
                aa1 = aa0;
                kmax = j;
            }
        }

        // | For N = −M2+1 .. (NPWSZ − NFRSZ)/4: DEC(N) = DEC(N + IDIM)
        for n in (1 - m2 as isize)..=((NPWSZ - NFRSZ) / 4) as isize {
            let dst = (DEC_OFF as isize + n) as usize;
            self.dec[dst] = self.dec[dst + IDIM];
        }

        // ---- Refine ±3 in the undecimated domain --------------------
        // | M1 = 4·KMAX − 3; M2 = 4·KMAX + 3; clamp to [KPMIN, KPMAX]
        let m1 = (4 * kmax - 3).max(KPMIN);
        let m2 = (4 * kmax + 3).min(KPMAX);
        let mut aa1 = i64::from(i32::MIN);
        let mut kp = m1;
        for j in m1..=m2 {
            // | AA0 = AA0 + D(K) · D(K − J)   | Correlation in undecimated domain
            //
            // The §G.3.25 *fixed-point* listing prints this inner MAC as
            // `AA0 = AA0 + D(K) * DEC(K – J)` — the second factor's array
            // name is a misprint: `DEC` is the *decimated* buffer, whose
            // index range (`−34 ..= NPWSZ/4`) cannot even hold `K − J` for
            // `K` up to `NPWSZ = 100`. The refined undecimated-domain
            // autocorrelation is `D(K)·D(K − J)`, exactly as the paired
            // floating-point listing on the previous page prints it
            // (`TMP = TMP + D(K) * D(K – J)`) and as the sibling
            // sub-multiple loop below prints it in the fixed-point text
            // itself. Conformance (outb4g) is bit-exact only with this
            // form — erratum E4 in `docs/audio/g728/g728-errata.md`.
            let mut aa0: i64 = 0;
            for k in 1..=NPWSZ {
                aa0 += i64::from(self.d[D_OFF + k]) * i64::from(self.d[(D_OFF + k) - j]);
            }
            if aa0 > aa1 {
                aa1 = aa0;
                kp = j;
            }
        }
        let cormax = aa1; // | Double precision save to CORMAX
        self.kp = kp;

        // ---- Sub-multiple check around the previous pitch -----------
        // | M1 = KP1 − KPDELTA; M2 = KP1 + KPDELTA
        // | If KP < M2 + 1, go to LABEL      | KP can't be a multiple
        // | clamp M1 to KPMIN, M2 to KPMAX (fixed-point-only clamp)
        let m2 = self.kp1 + KPDELTA;
        // | If KP < M2 + 1, go to LABEL  (i.e. only continue when KP > M2)
        if self.kp > m2 {
            let m1 = self.kp1.saturating_sub(KPDELTA).max(KPMIN);
            let m2 = m2.min(KPMAX);
            let mut aa1 = i64::from(i32::MIN);
            let mut kptmp = m1;
            for j in m1..=m2 {
                let mut aa0: i64 = 0;
                for k in 1..=NPWSZ {
                    aa0 += i64::from(self.d[D_OFF + k]) * i64::from(self.d[(D_OFF + k) - j]);
                }
                if aa0 > aa1 {
                    aa1 = aa0;
                    kptmp = j;
                }
            }
            let cmax = aa1;

            // | AA0 = Σ D(K−KP)²; AA1 = Σ D(K−KPTMP)²
            let mut aa0: i64 = 0;
            let mut aa1: i64 = 0;
            for k in 1..=NPWSZ {
                let a = i64::from(self.d[(D_OFF + k) - self.kp]);
                let b = i64::from(self.d[(D_OFF + k) - kptmp]);
                aa0 += a * a;
                aa1 += b * b;
            }

            // | Clip the correlations into [0, energy] (tap ∈ [0, 1]).
            let mut cormax = cormax;
            let mut cmax = cmax;
            if aa0 == 0 {
                cormax = 0;
            }
            if aa1 == 0 {
                cmax = 0;
            }
            cormax = cormax.clamp(0, aa0);
            cmax = cmax.clamp(0, aa1);

            // | Align both energies to 30 bits with one shared NLS.
            let (nls, aa0, aa1) = if aa0 > aa1 {
                let nls = findnls(&[clamp_i32(aa0)], 1, 30);
                (nls, shl_signed(aa0, nls), shl_signed(aa1, nls))
            } else {
                let nls = findnls(&[clamp_i32(aa1)], 1, 30);
                (nls, shl_signed(aa0, nls), shl_signed(aa1, nls))
            };

            // | SUM = AA0 >> 16; TMP = AA1 >> 16
            // | CORMAX = (CORMAX << NLS) >> 16; CMAX = (CMAX << NLS) >> 16
            let sum = aa0 >> 16;
            let tmp = aa1 >> 16;
            let cormax = shl_signed(cormax, nls) >> 16;
            let cmax = shl_signed(cmax, nls) >> 16;

            // | AA1 = ((CORMAX·TMP) >> 16) · ITAPTH
            // | AA0 = CMAX·SUM
            // | If AA0 > AA1, set KP = KPTMP     | TAP1 > TAPTH·TAP
            let aa1 = ((cormax * tmp) >> 16) * ITAPTH_Q16;
            let aa0 = cmax * sum;
            if aa0 > aa1 {
                self.kp = kptmp;
            }
        }

        // LABEL:
        // | KP1 = KP
        // | For K = −KPMAX+1 .. NPWSZ − NFRSZ: D(K) = D(K + NFRSZ)
        self.kp1 = self.kp;
        for k in (1 - KPMAX as isize)..=(NPWSZ - NFRSZ) as isize {
            let dst = (D_OFF as isize + k) as usize;
            self.d[dst] = self.d[dst + NFRSZ];
        }
    }

    /// §G.3.26 block 83 — the pitch predictor tap `PTAP` (Q14).
    ///
    /// `sst` is the postfilter's quantized-speech buffer (13-bit `Q0`
    /// deep past), laid out oldest-first with the previous vector's last
    /// sample `SST(0)` at `sst[sst_zero]`. The correlations run over
    /// `K = −NPWSZ+1 ..= 0`, reaching back a further `KP` samples.
    #[must_use]
    pub fn pitch_tap(&self, sst: &[i32], sst_zero: usize) -> i16 {
        let kp = self.kp;
        // | AA0 = Σ SST(K−KP)²; AA1 = Σ SST(K)·SST(K−KP)
        let mut aa0: i64 = 0;
        let mut aa1: i64 = 0;
        for k in (1 - NPWSZ as isize)..=0 {
            let idx = (sst_zero as isize + k) as usize;
            let cur = i64::from(sst[idx]);
            let lag = i64::from(sst[idx - kp]);
            aa0 += lag * lag;
            aa1 += cur * lag;
        }

        // | If AA0 = 0 or AA1 ≤ 0, PTAP = 0
        if aa0 == 0 || aa1 <= 0 {
            return 0;
        }
        // | If AA1 ≥ AA0, PTAP = 16384 (NLSPTAP = 14)
        if aa1 >= aa0 {
            return 16384;
        }
        // | VSCALE both to 30 bits, RND to 16, DIVIDE, re-normalize to Q14.
        let nlsden = findnls(&[clamp_i32(aa0)], 1, 30);
        let nlsnum = findnls(&[clamp_i32(aa1)], 1, 30);
        let den = rnd_i64(shl_signed(aa0, nlsden));
        let num = rnd_i64(shl_signed(aa1, nlsnum));
        let (ptap, nlsptap) = divide(num, nlsnum, den, nlsden);
        // | NRS = NLSPTAP − 14; PTAP = PTAP >> NRS
        //
        // `ptap` is a normalized 16-bit mantissa (`ptap ≥ 0` here — both
        // `num` and `den` are positive by the guards above), so on the
        // conformance path NRS stays in `[0, 15]` and the shift is exact.
        // A degenerate or corrupted codeword stream can push NLSPTAP far
        // enough that the raw `i16` shift amount reaches the type width
        // and panics ("shift with overflow"); guard the out-of-range
        // amounts (in-range arithmetic is left byte-identical):
        //   * NRS ≥ 16  → the mantissa shifts entirely below the LSB → 0;
        //   * −NRS ≥ 16 → the tap caps at the Q14 unit ceiling `16384`
        //     (= 1.0, matching the `AA1 ≥ AA0` early return; unreachable
        //     on valid input, where `AA1 < AA0` bounds the ratio below 1).
        let nrs = nlsptap - 14;
        if nrs >= 16 {
            0
        } else if nrs >= 0 {
            ptap >> nrs
        } else if -nrs < 16 {
            ptap << (-nrs)
        } else {
            16384
        }
    }

    /// §G.3.27 block 84 — the long-term postfilter coefficient
    /// calculator: threshold `PTAP` against `PPFTH`, form `B =
    /// PPFZCF·PTAP` (Q16), `GL = 1/(1+B)` (Q14 via `DIVIDE`) and the
    /// precomputed `GLB = GL·B` (Q16). Returns the [`LongTermCoeff`]
    /// (carrying the current `KP`) for blocks 71/72.
    #[must_use]
    pub fn long_term_coeff(&self, ptap: i16) -> LongTermCoeff {
        // | If PTAP < PPFTH, set PTAP = 0     | PPFTH = 9830 Q14
        let ptap = if i32::from(ptap) < PPFTH_Q14 {
            0
        } else {
            i32::from(ptap)
        };
        // | AA0 = PPFZCF·PTAP                 | Q16·Q14 = Q30
        let aa0 = PPFZCF_Q16 * ptap;
        // | B = AA0 >> 14                     | Q16
        let b = aa0 >> 14;
        // | AA0 = AA0 >> 16; AA0 += 16384; DEN = AA0    | 1 + B in Q14
        let den = (aa0 >> 16) + 16384;
        // | DIVIDE(16384, 14, DEN, 14, GL, NLS)
        let (gl, nls) = divide(16384, 14, den as i16, 14);
        // | AA0 = GL·B; GLB = AA0 >> NLS      | GLB Q16
        let glb = ((i64::from(gl) * i64::from(b)) >> nls) as i32;
        // | NRS = NLS − 14; if NRS > 0, GL = GL >> NRS  | GL Q14
        let nrs = nls - 14;
        let gl = if nrs > 0 {
            i32::from(gl) >> nrs
        } else {
            i32::from(gl)
        };

        LongTermCoeff {
            kp: self.kp,
            gl_q14: gl,
            glb_q16: glb,
        }
    }
}

// ---------------------------------------------------------------------
// Small helpers (§G.3 clip / shift idioms shared by the blocks above).
// ---------------------------------------------------------------------

/// `RND` on a 64-bit accumulator holding a ≤32-bit value (saturating
/// into the 32-bit range first, like the annex's DSP accumulator).
#[inline]
fn rnd_i64(v: i64) -> i16 {
    rnd(clamp_i32(v))
}

/// Saturate a 64-bit accumulator into the 32-bit range.
#[inline]
fn clamp_i32(v: i64) -> i32 {
    v.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32
}

/// Left shift by a possibly-negative count (negative ⇒ arithmetic right
/// shift) — the `VSCALE`-style alignment used on the 32-bit energies.
#[inline]
fn shl_signed(v: i64, k: i32) -> i64 {
    if k >= 0 {
        v << k
    } else {
        v >> (-k)
    }
}

/// Clip a 32-bit value to the signed 16-bit range.
#[inline]
fn clip_i16(v: i32) -> i16 {
    v.clamp(i32::from(i16::MIN), i32::from(i16::MAX)) as i16
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::{PPFTH, PPFZCF, TAPTH};

    #[test]
    fn constants_match_table_g1() {
        // Table G.1/G.728 integer columns vs the float Table 1 values.
        assert_eq!(PPFTH_Q14, (PPFTH * 16384.0).round() as i32);
        assert_eq!(PPFZCF_Q16, (PPFZCF * 65536.0).round() as i32);
        assert_eq!(ITAPTH_Q16, (TAPTH * 65536.0).round() as i64);
        // §G.3.25 low-pass constants vs the float Annex D tables.
        for i in 0..4 {
            let want = (crate::tables::BL[i] * (1 << 19) as f64).round() as i32;
            assert!(
                (BL_Q19[i] - want).abs() <= 1,
                "BL_Q19[{i}] = {} vs float {want}",
                BL_Q19[i]
            );
        }
        for i in 0..3 {
            let want = (crate::tables::AL[i + 1] * 8192.0).round() as i32;
            assert!(
                (AL_Q13[i] - want).abs() <= 1,
                "AL_Q13[{i}] = {} vs float {want}",
                AL_Q13[i]
            );
        }
    }

    #[test]
    fn apf_q13_conversion_all_three_formats() {
        // A predictor with |taps| < 1 expressed at each NLSATMP must
        // convert to (nearly) the same Q13 integers.
        let taps = [
            -0.5f64, 0.25, -0.125, 0.0625, 0.4, -0.3, 0.2, -0.1, 0.05, -0.025,
        ];
        for nls in [13, 14, 15] {
            let scale = (1i32 << nls) as f64;
            let mut atil = [0i16; PF_ORDER + 1];
            atil[0] = if nls == 15 { i16::MAX } else { 1 << nls };
            for (i, &t) in taps.iter().enumerate() {
                atil[i + 1] = (t * scale).round() as i16;
            }
            let q13 = apf_to_q13(&atil, nls);
            assert_eq!(q13[0], 8192);
            for i in 1..=PF_ORDER {
                let want = (taps[i - 1] * 8192.0).round() as i32;
                assert!(
                    (i32::from(q13[i]) - want).abs() <= 1,
                    "nls {nls} tap {i}: {} vs {want}",
                    q13[i]
                );
            }
        }
    }

    #[test]
    fn inverse_filter_tracks_float_residual() {
        // Block 81 vs the floating-point §4.6 LPC inverse filter on the
        // same Q2-quantized speech and Q13-quantized predictor.
        use crate::pitch_inverse_filter::PitchInverseFilter;

        let taps = [
            -0.7f64, 0.3, -0.15, 0.08, -0.04, 0.02, -0.01, 0.005, 0.0, 0.0,
        ];
        let mut apf = [0i16; PF_ORDER + 1];
        apf[0] = 8192;
        let mut a10 = [0.0f64; PF_ORDER + 1];
        a10[0] = 1.0;
        for (i, &t) in taps.iter().enumerate() {
            apf[i + 1] = (t * 8192.0).round() as i16;
            // The float filter sees the SAME requantized coefficients.
            a10[i + 1] = f64::from(apf[i + 1]) / 8192.0;
        }

        let mut fixed = PitchAdapterFixed::new();
        let mut float = PitchInverseFilter::new();
        float.set_from_synthesis_byproduct(&a10);

        let mut seed = 0x1122_3344_5566_7788u64;
        for _ in 0..12 {
            // A random Q2 speech vector (|s| < 4096) presented as 14-bit
            // BFL with nlsst = 2 (mantissa = Q2 value directly).
            let mut st = [0i32; IDIM];
            let mut sf = [0.0f64; IDIM];
            for k in 0..IDIM {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                let v = ((seed >> 40) as i32 % 8192) - 4096; // Q2 mantissa
                st[k] = v;
                sf[k] = f64::from(v) / 4.0;
            }
            let sst = fixed.inverse_filter(&st, 2, &apf);
            let dfloat = float.filter_vector(&sf);
            // SST must be the exact Q2 rendering of the input.
            for k in 0..IDIM {
                assert_eq!(sst[k], st[k], "SST Q2 passthrough");
            }
            // The newest 5 residual samples sit at D(IP−IDIM+1 .. IP).
            for k in 0..IDIM {
                let dq1 = f64::from(fixed.d[D_OFF + fixed.ip - IDIM + 1 + k]) / 2.0;
                assert!(
                    (dq1 - dfloat[k]).abs() <= 1.0,
                    "D[{k}] fixed {dq1} vs float {}",
                    dfloat[k]
                );
            }
        }
    }

    #[test]
    fn pitch_extract_finds_period_of_pulse_train() {
        // Feed block 81 a period-40 impulse train through an identity
        // predictor (APF = δ), then run block 82: the extracted pitch
        // must be exactly 40. Two cycles of warm-up fill D.
        let mut fixed = PitchAdapterFixed::new();
        let mut apf = [0i16; PF_ORDER + 1];
        apf[0] = 8192;

        let period = 40usize;
        let mut n = 0usize;
        for cycle in 0..6 {
            for _v in 0..4 {
                let mut st = [0i32; IDIM];
                for k in 0..IDIM {
                    if n % period == 0 {
                        st[k] = 4000; // Q2 pulse
                    }
                    n += 1;
                }
                let _ = fixed.inverse_filter(&st, 2, &apf);
            }
            fixed.pitch_extract();
            if cycle >= 2 {
                assert_eq!(fixed.kp(), period, "cycle {cycle}");
            }
        }
    }

    #[test]
    fn pitch_tap_matches_float_ratio() {
        // Block 83 vs the float Σ s(k)s(k−p) / Σ s(k−p)² on a decaying
        // sinusoid buffer with a known period.
        let kp = 50usize;
        let len = 400usize;
        let mut sst = vec![0i32; len];
        for (i, v) in sst.iter_mut().enumerate() {
            let t = i as f64;
            *v = ((t * std::f64::consts::TAU / kp as f64).sin() * 3000.0) as i32;
        }
        let sst_zero = len - 1;

        let mut a = PitchAdapterFixed::new();
        a.kp = kp;
        let ptap = a.pitch_tap(&sst, sst_zero);

        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for k in (1 - NPWSZ as isize)..=0 {
            let idx = (sst_zero as isize + k) as usize;
            num += f64::from(sst[idx]) * f64::from(sst[idx - kp]);
            den += f64::from(sst[idx - kp]) * f64::from(sst[idx - kp]);
        }
        let want = (num / den).clamp(0.0, 1.0);
        let got = f64::from(ptap) / 16384.0;
        assert!(
            (got - want).abs() < 2e-3,
            "PTAP fixed {got} vs float {want}"
        );
        // A perfectly periodic signal has PTAP ≈ 1.
        assert!(got > 0.9, "periodic signal must give a strong tap");
    }

    #[test]
    fn pitch_tap_zero_cases() {
        let a = PitchAdapterFixed::new();
        // All-zero history ⇒ AA0 = 0 ⇒ PTAP = 0.
        let sst = vec![0i32; 400];
        assert_eq!(a.pitch_tap(&sst, 399), 0);
        // Anti-correlated (negative AA1) ⇒ PTAP = 0.
        let mut sst = vec![0i32; 400];
        for (i, v) in sst.iter_mut().enumerate() {
            // Alternate sign every KP1_INIT = 50 samples ⇒ s(k)·s(k−50) < 0.
            *v = if (i / KP1_INIT) % 2 == 0 { 1000 } else { -1000 };
        }
        assert_eq!(a.pitch_tap(&sst, 399), 0);
    }

    #[test]
    fn long_term_coeff_matches_float_formula() {
        // Block 84 vs the float B = PPFZCF·β, GL = 1/(1+B) across the
        // whole tap range, including the PPFTH cut-off.
        let a = PitchAdapterFixed::new();
        for &beta in &[0.0f64, 0.3, 0.59, 0.6, 0.61, 0.75, 0.9, 1.0] {
            let ptap = (beta * 16384.0).round() as i16;
            let lt = a.long_term_coeff(ptap);
            // Same integer threshold the fixed point applies: PTAP (Q14)
            // against PPFTH_Q14 (9830 = the Q14 rendering of 0.6, which
            // round(0.6·2¹⁴) itself reaches — so β = 0.6 stays ON).
            let beta_eff = if i32::from(ptap) < PPFTH_Q14 {
                0.0
            } else {
                f64::from(ptap) / 16384.0
            };
            let b = PPFZCF * beta_eff;
            let gl = 1.0 / (1.0 + b);
            assert!(
                (f64::from(lt.gl_q14) / 16384.0 - gl).abs() < 2e-3,
                "β={beta}: GL fixed {} vs float {gl}",
                f64::from(lt.gl_q14) / 16384.0
            );
            assert!(
                (f64::from(lt.glb_q16) / 65536.0 - gl * b).abs() < 2e-3,
                "β={beta}: GLB fixed {} vs float {}",
                f64::from(lt.glb_q16) / 65536.0,
                gl * b
            );
            assert_eq!(lt.kp, a.kp());
        }
    }

    #[test]
    fn fresh_state_matches_table_2() {
        let a = PitchAdapterFixed::new();
        assert_eq!(a.ip, 85, "IPINIT = NPWSZ − NFRSZ + IDIM");
        assert_eq!(a.kp1, 50);
        assert!(a.d.iter().all(|&v| v == 0));
        assert!(a.dec.iter().all(|&v| v == 0));
        assert_eq!(D_LEN, 240);
        assert_eq!(DEC_LEN, 60);
    }

    #[test]
    fn erratum_e3_ip_wraps_to_npwsz_minus_nfrsz() {
        // Erratum E3 (docs/audio/g728/g728-errata.md): the §G.3.24
        // fixed-point listing prints the block-81 write-pointer reset as
        // `If IP = NPWSZ, then set IP = NFRSZ` — the reset target must be
        // `NPWSZ − NFRSZ` (= 80), as the paired floating-point listing
        // prints. Pin the whole IP trajectory: from IPINIT = 85 the
        // pointer walks 90 → 95 → 100, then rewinds to 80 before the
        // next write (so the post-write value is 85 again), and never
        // leaves the newest-NPWSZ window. With the as-printed `NFRSZ`
        // reset the post-write value after the wrap would be 25.
        let mut a = PitchAdapterFixed::new();
        let apf = apf_identity();
        let st = [0i32; IDIM];
        let expected = [90, 95, 100, 85, 90, 95, 100, 85, 90];
        for (call, &want) in expected.iter().enumerate() {
            let _ = a.inverse_filter(&st, 2, &apf);
            assert_eq!(a.ip(), want, "IP after call {call}");
            assert!(
                a.ip() > NPWSZ - NFRSZ && a.ip() <= NPWSZ,
                "IP must stay inside the newest-window (80, 100]"
            );
        }
    }

    /// Identity order-10 predictor in Q13 (`APF = δ`).
    fn apf_identity() -> [i16; PF_ORDER + 1] {
        let mut apf = [0i16; PF_ORDER + 1];
        apf[0] = 8192;
        apf
    }

    #[test]
    fn erratum_e4_refined_search_correlates_d_with_d() {
        // Erratum E4 (docs/audio/g728/g728-errata.md): the §G.3.25
        // fixed-point refined (undecimated) peak-pick prints its inner
        // MAC as `AA0 = AA0 + D(K) * DEC(K − J)`; the second factor's
        // array name is a misprint for `D(K − J)` (the paired
        // floating-point listing prints `TMP = TMP + D(K) * D(K − J)`
        // literally, and `DEC` cannot even be indexed at `K − J`).
        //
        // Pin the corrected behaviour with a two-pulse residual: pulses
        // at D(60) and D(100) (lag 40) plus a decimated-domain pulse at
        // DEC(15) so the coarse pick lands at KMAX = 10 (via the
        // lowpassed D(100) pulse entering DEC(25)). The refined search
        // over J = 37..43 must then pick the true undecimated lag 40 —
        // the argmax of Σ D(K)·D(K−J), recomputed independently here.
        let mut a = PitchAdapterFixed::new();
        a.d[D_OFF + 100] = 8000;
        a.d[D_OFF + 60] = 8000;
        a.dec[DEC_OFF + 15] = 8000;
        let d_before = a.d;

        a.pitch_extract();

        // Independent argmax of the refined-window autocorrelation on
        // the pre-extract residual snapshot.
        let mut best_j = 0usize;
        let mut best = i64::MIN;
        for j in 37..=43usize {
            let mut corr = 0i64;
            for k in 1..=NPWSZ {
                corr += i64::from(d_before[D_OFF + k]) * i64::from(d_before[(D_OFF + k) - j]);
            }
            if corr > best {
                best = corr;
                best_j = j;
            }
        }
        assert_eq!(best_j, 40, "test-vector self-check");
        assert_eq!(
            a.kp(),
            best_j,
            "block 82 must pick the argmax of the D·D refined correlation"
        );
    }
}
