//! Annex G §G.3.4 … §G.3.10 — fixed-point codebook search (the
//! analysis-by-synthesis coder of blocks 11 – 21).
//!
//! ITU-T G.728 Annex G (1994-11) §G.3 gives the bit-exact fixed-point
//! pseudo-code for the per-block modules of the floating-point coder.
//! This module transcribes the encoder-side search chain — the part of
//! Figure 2/G.728 between the perceptual-weighting target and the
//! emitted 10-bit channel index — one §G.3 sub-clause at a time:
//!
//! | §G.3 | Block(s) | Function |
//! |------|----------|----------|
//! | G.3.4  | 11    | [`vq_target`] — `TARGET = SW − ZIR` |
//! | G.3.5  | 12    | [`impulse_response`] — `H` from filter memory |
//! | G.3.6  | 13    | [`time_reversed_convolution`] — `PN` |
//! | G.3.7  | 14/15 | [`codevector_energy`] — `Y2[j]` |
//! | G.3.8  | 16    | [`normalize_target`] — `TARGET / GAIN` |
//! | G.3.9  | 17/18 | [`search`] — best `(IS, IG)`, `ICHAN` |
//! | G.3.10 | 19/21 | [`excitation`] — `ET = GAIN·GQ·Y` |
//!
//! The whole chain is driven end-to-end for one vector by
//! [`encode_vector`].
//!
//! All arithmetic is built on the §G.1.3 primitives in
//! [`crate::annex_g_arith`] ([`rnd`](crate::annex_g_arith::rnd),
//! [`vscale`](crate::annex_g_arith::vscale),
//! [`divide`](crate::annex_g_arith::divide)) and reads the fixed-point
//! ROM codebooks from [`crate::tables`] ([`Y_Q11`]). The Q-formats are
//! exactly those stated in the §G.3 prose:
//!
//! * `SW`, `ZIR`, `TARGET` enter in `Q2` (`NLSTARGET = 2`), the speech
//!   format;
//! * the predictor coefficients `A`, `AWZ`, `AWP` are `Q14`;
//! * the impulse response `H` is `Q13`;
//! * the shape codebook `Y` is `Q11`;
//! * `PN` is `Q7` (`NLSPN = 7`);
//! * the codevector energies `Y2` are `NLSY2 = 5`;
//! * the gain tables are `GB` (`NLS = 13`), `G2` (`NLS = 12`),
//!   `GSQ` (`NLS = 11`), `GQ` (`Q13`).
//!
//! The floating-point search already lives in
//! [`crate::codebook_search`]; this module is the §G.3 fixed-point
//! reformulation and is self-contained: it introduces no codec state.

use crate::annex_g_arith::{divide, rnd, vscale};
use crate::consts::{IDIM, NCWD, NG};
use crate::tables::Y_Q11;

/// `N1 = NG / 2` — the gain-codebook half size separating the positive
/// gain levels (indices `1 … N1`) from the sign-mirrored negatives
/// (indices `N1+1 … NG`) in the §G.3.9 search (1-based, as in the
/// pseudo-code).
pub const N1: usize = NG / 2;

/// `NLSPN` — the fixed left-shift count of the `PN` array produced by
/// the §G.3.6 time-reversed convolution (block 13). `PN` is held in
/// `Q7`, so `NLSPN = 7`.
pub const NLSPN: i32 = 7;

/// `NLSY2` — the fixed left-shift count of the codevector-energy array
/// `Y2` produced by §G.3.7 (block 14/15). With `NLSY = 11` and
/// `NLSH = 13`, the §G.3.7 formula
/// `NLSY2 = (NLSH + NLSY − 14)·2 − 15` gives `NLSY2 = 5`.
pub const NLSY2: i32 = 5;

// ---------------------------------------------------------------------
// §G.5 / §G.3.9 Q-format gain tables.
//
// The §G.3.9 search uses the gain-codebook quantities in fixed Q
// formats different from the floating-point [`crate::tables`] copies:
//   GB  in NLS = 13  (cell boundaries, eq. "NLS for P is 13 + 5 = 18")
//   G2  in NLS = 12  (= 2·GQ, eq. "NLSG2 = 12")
//   GSQ in NLS = 11  (= GQ², eq. "NLSGSQ = 11")
//   GQ  in Q13       (gain levels, §G.3.10)
// We re-derive each Q-format integer from the float [`crate::tables`]
// rationals (themselves derived from the spec footnote GQ(1)=33/64,
// GQ(i)=(7/4)·GQ(i−1)) by rounding to the stated Q quantum. The
// per-test cross-checks below prove the integer agrees with the float
// to one quantum.
// ---------------------------------------------------------------------

/// Gain levels `GQ[i]` in `Q13` (§G.3.10), `0 ≤ i < NG`.
pub const GQ_Q13: [i32; NG] = q13_gq();

/// Gain-quantizer cell boundaries `GB[k]` in `NLS = 13`
/// (§G.3.9), `0 ≤ k < NG − 1`. Only `k = 0,1,2` are read by the search
/// (the positive half-cell boundaries); the §G.3.9 fixed-point code
/// uses the absolute correlation, so the negative boundaries are not
/// needed.
pub const GB_Q13: [i32; NG - 1] = q13_gb();

/// `G2[i] = 2·GQ[i]` in `NLS = 12` (§G.3.9), `0 ≤ i < NG`.
pub const G2_Q12: [i32; NG] = q_scaled_gq(2.0, 12);

/// `GSQ[i] = GQ[i]²` in `NLS = 11` (§G.3.9), `0 ≤ i < NG`.
pub const GSQ_Q11: [i32; NG] = gsq_q11();

/// `NNGQ[i]` — `1 +` the number of left shifts needed to normalize the
/// `Q13` `GQ[i]` to a 16-bit word (§G.3.10): `3` for `i = 1,2,5,6`,
/// `2` for `i = 3,7`, `1` for `i = 4,8` (1-based). Stored 0-based.
pub const NNGQ: [u32; NG] = [3, 3, 2, 1, 3, 3, 2, 1];

/// The undivided floating gain levels, in the same order as
/// [`crate::tables::GQ`], used by the `const fn` table builders. Kept
/// local so the builders stay `const`.
const GQ_F: [f64; NG] = {
    let g0 = 33.0 / 64.0;
    let r = 7.0 / 4.0;
    let g1 = g0 * r;
    let g2 = g1 * r;
    let g3 = g2 * r;
    [g0, g1, g2, g3, -g0, -g1, -g2, -g3]
};

/// Build the `Q13` gain-level table.
const fn q13_gq() -> [i32; NG] {
    q_scaled_gq(1.0, 13)
}

/// Build `GB` in `Q13` from the mid-point definition
/// `GB[k] = (GQ[k] + GQ[k+1]) / 2`.
const fn q13_gb() -> [i32; NG - 1] {
    let scale = (1i64 << 13) as f64;
    // Mirror the float [`crate::tables::GB`] convention: d0,d1,d2 are
    // the positive half-cell mid-points; index 3 is the "any arbitrary
    // value (not used)" placeholder (held at 0); d4,d5,d6 are the
    // sign-mirrored negatives. The §G.3.9 fixed-point search reads only
    // the positive boundaries d0,d1,d2 (it uses |correlation|).
    [
        round_to_i32((GQ_F[0] + GQ_F[1]) / 2.0 * scale),
        round_to_i32((GQ_F[1] + GQ_F[2]) / 2.0 * scale),
        round_to_i32((GQ_F[2] + GQ_F[3]) / 2.0 * scale),
        0,
        round_to_i32((GQ_F[4] + GQ_F[5]) / 2.0 * scale),
        round_to_i32((GQ_F[5] + GQ_F[6]) / 2.0 * scale),
        round_to_i32((GQ_F[6] + GQ_F[7]) / 2.0 * scale),
    ]
}

/// Build `factor · GQ[i]` in the `q`-bit Q format.
const fn q_scaled_gq(factor: f64, q: u32) -> [i32; NG] {
    let mut out = [0i32; NG];
    let scale = (1i64 << q) as f64;
    let mut i = 0;
    while i < NG {
        out[i] = round_to_i32(GQ_F[i] * factor * scale);
        i += 1;
    }
    out
}

/// Build `GSQ[i] = GQ[i]²` in `Q11`.
const fn gsq_q11() -> [i32; NG] {
    let mut out = [0i32; NG];
    let scale = (1i64 << 11) as f64;
    let mut i = 0;
    while i < NG {
        out[i] = round_to_i32(GQ_F[i] * GQ_F[i] * scale);
        i += 1;
    }
    out
}

/// Round-half-away-from-zero of an `f64` to `i32`, usable in `const fn`.
const fn round_to_i32(x: f64) -> i32 {
    if x >= 0.0 {
        (x + 0.5) as i32
    } else {
        -((-x + 0.5) as i32)
    }
}

// ---------------------------------------------------------------------
// §G.3.4 — Block 11: VQ target vector computation.
// ---------------------------------------------------------------------

/// Block 11 (§G.3.4) — form the VQ target vector `TARGET = SW − ZIR`.
///
/// `SW` (weighted speech) and `ZIR` (zero-input response) both enter in
/// `Q2`, the speech format, so the output `TARGET` is also `Q2`
/// (`NLSTARGET = 2`). The subtraction is clipped to the 16-bit range
/// per the §G.3.4 saturation lines.
///
/// Returns the `IDIM`-length `TARGET` array; its `NLSTARGET` is the
/// constant `2`.
#[must_use]
pub fn vq_target(sw: &[i16; IDIM], zir: &[i16; IDIM]) -> [i16; IDIM] {
    let mut target = [0i16; IDIM];
    for k in 0..IDIM {
        // AA0 = SW(K) − ZIR(K), clipped to 16 bits.
        let aa0 = sw[k] as i32 - zir[k] as i32;
        target[k] = clip_i16(aa0);
    }
    target
}

// ---------------------------------------------------------------------
// §G.3.5 — Block 12: impulse response vector calculation.
// ---------------------------------------------------------------------

/// Block 12 (§G.3.5) — compute the impulse response `H(1..IDIM)` of the
/// cascade synthesis × weighting filter from the predictor coefficient
/// arrays.
///
/// `a` is the synthesis-filter coefficient array, `awz` the weighting
/// all-zero coefficients, `awp` the weighting all-pole coefficients,
/// **all in `Q14`** and in the crate's 0-based layout: `a[0]` is the
/// spec's `A(1)` (the implicit leading `1` = `16384`), `a[i - 1]` is
/// the spec's `A(i)` for `i = 2 … IDIM + 1` — exactly the head of the
/// live predictor / weighting arrays the §G.7 driver holds. Each slice
/// has length `IDIM + 1`.
///
/// The output `H` is `Q13` (`NLSH = 13`). The two `Q14` predictor
/// columns are merged into a single accumulator (the §G.3.5 note
/// "accumulators A1 and A2 … have been combined") and the `>> 14`
/// folds the `Q14` coefficient scale back out.
#[must_use]
pub fn impulse_response(
    a: &[i16; IDIM + 1],
    awz: &[i16; IDIM + 1],
    awp: &[i16; IDIM + 1],
) -> [i16; IDIM] {
    // TEMP = synthesis filter memory, WS = W(z) all-pole part memory,
    // both Q13 16-bit words seeded with 8192 = 1.0 in Q13.
    let mut temp = [0i16; IDIM + 1];
    let mut ws = [0i16; IDIM + 1];
    temp[1] = 8192;
    ws[1] = 8192;

    for k in 2..=IDIM {
        let mut aa0: i64 = 0; // synthesis branch
        let mut aa1: i64 = 0; // combined weighting branch
                              // For I = K, K−1, …, 2 (shift the memory then filter).
        let mut i = k;
        while i >= 2 {
            temp[i] = temp[i - 1];
            ws[i] = ws[i - 1];
            // A(I)·TEMP(I) etc., with the spec's 1-based A(I) at the
            // crate's 0-based slot `a[i − 1]` (`a[0]` = `A(1)` = 16384).
            aa0 -= a[i - 1] as i64 * temp[i] as i64;
            aa1 += awz[i - 1] as i64 * temp[i] as i64;
            aa1 -= awp[i - 1] as i64 * ws[i] as i64;
            i -= 1;
        }
        aa1 += aa0;
        // >> 14 because A, AWZ, AWP are Q14.
        let aa0 = aa0 >> 14;
        let aa1 = aa1 >> 14;
        temp[1] = aa0 as i16;
        ws[1] = aa1 as i16;
    }

    // Obtain h(n) by reversing the order of the W(z) all-pole memory.
    let mut h = [0i16; IDIM];
    let itmp = IDIM + 1;
    for k in 1..=IDIM {
        h[k - 1] = ws[itmp - k];
    }
    h
}

// ---------------------------------------------------------------------
// §G.3.6 — Block 13: time-reversed convolution.
// ---------------------------------------------------------------------

/// Block 13 (§G.3.6) — time-reversed convolution of `TARGET` with the
/// impulse response `H`, producing `PN` for the codebook search.
///
/// `target` is the (block-floating) VQ target with shift `nlstarget`
/// (its `NLSTARGET`, determined by block 16); `h` is the `Q13` impulse
/// response. The output `PN` is clipped to 16 bits and held at the
/// fixed `NLSPN = 7`.
#[must_use]
pub fn time_reversed_convolution(
    target: &[i16; IDIM],
    h: &[i16; IDIM],
    nlstarget: i32,
) -> [i16; IDIM] {
    let mut pn = [0i16; IDIM];
    // Right shift to make Q7: >> (13 + (NLSTARGET − 7)).
    let shift = 13 + (nlstarget - NLSPN);
    for k in 0..IDIM {
        let k1 = k; // K1 = K − 1 in 1-based → k in 0-based
        let mut aa0: i64 = 0;
        for j in k..IDIM {
            // P = TARGET(J) * H(J − K1)
            aa0 += target[j] as i64 * h[j - k1] as i64;
        }
        let aa0 = shr_signed(aa0, shift);
        pn[k] = clip_i16_64(aa0);
    }
    pn
}

// ---------------------------------------------------------------------
// §G.3.7 — Blocks 14 / 15: shape codevector convolution + energy.
// ---------------------------------------------------------------------

/// Blocks 14 / 15 (§G.3.7) — convolve each of the `NCWD` shape
/// codevectors with the impulse response `H` and compute its energy
/// `Y2[j]`.
///
/// `h` is the `Q13` impulse response; the shape codebook `Y` is read
/// from [`Y_Q11`] (`Q11`). The accumulator never overflows 32 bits
/// empirically (§G.3.7), so no overflow test is performed; the `>> 14`
/// puts the convolution `TEMP` into `Q10`, and the energy sum is
/// `>> 15` to give `Y2` at `NLSY2 = 5`.
///
/// Returns the `NCWD`-length `Y2` array (16-bit words).
#[must_use]
pub fn codevector_energy(h: &[i16; IDIM]) -> [i16; NCWD] {
    let mut y2 = [0i16; NCWD];
    for j in 0..NCWD {
        let mut temp = [0i16; IDIM];
        for k in 0..IDIM {
            // K1 = J1 + K + 1 in the 1-based spec; in 0-based the
            // codevector row is Y_Q11[j] and the convolution is
            // TEMP(k) = Σ_{i=1..k} H(i)·Y(k − i + 1).
            let mut aa0: i64 = 0;
            for i in 0..=k {
                aa0 += h[i] as i64 * Y_Q11[j][k - i] as i64;
            }
            let aa0 = aa0 >> 14;
            temp[k] = aa0 as i16; // lowest 16 bits only
        }
        let mut aa0: i64 = 0;
        for k in 0..IDIM {
            aa0 += temp[k] as i64 * temp[k] as i64;
        }
        let aa0 = aa0 >> 15;
        y2[j] = aa0 as i16; // lowest 16 bits only
    }
    y2
}

// ---------------------------------------------------------------------
// §G.3.8 — Block 16: VQ target vector normalization.
// ---------------------------------------------------------------------

/// Block 16 (§G.3.8) — normalize the VQ target by the excitation gain:
/// `TARGET ← TARGET / GAIN`, returning the new block-floating-point
/// shift `NLSTARGET`.
///
/// `target` enters at `NLSTARGET = 2`; `gain`/`nlsgain` is the scalar
/// floating-point gain (mantissa + NLS) from block 46. The division
/// reciprocal uses numerator `16384` at `NLS = 14`
/// (`DIVIDE(16384, 14, GAIN, NLSGAIN, …)`); the product is then
/// re-block-floated with [`vscale`].
///
/// Returns `(target, nlstarget)`.
#[must_use]
pub fn normalize_target(target: &[i16; IDIM], gain: i16, nlsgain: i32) -> ([i16; IDIM], i32) {
    // TMP = 1 / GAIN, numerator 16384 at NLS = 14.
    let (tmp, nlstmp) = divide(16384, 14, gain, nlsgain);

    let mut out = [0i32; IDIM];
    for k in 0..IDIM {
        // AA0 = TMP * TARGET(K); keep the lower 16 bits (>> 15).
        let aa0 = tmp as i32 * target[k] as i32;
        out[k] = aa0 >> 15;
    }

    // NLSTARGET = 2 + NLSTMP − 15, then re-block-float.
    let mut nlstarget = 2 + nlstmp - 15;
    let (scaled, nls) = vscale(&out, IDIM, 14);
    nlstarget += nls;

    let mut target_out = [0i16; IDIM];
    for k in 0..IDIM {
        target_out[k] = scaled[k] as i16;
    }
    (target_out, nlstarget)
}

// ---------------------------------------------------------------------
// §G.3.9 — Blocks 17 / 18: VQ search error calculator + best index.
// ---------------------------------------------------------------------

/// Result of the §G.3.9 codebook search: the 1-based shape index `is`
/// (`1 … NCWD`), the 1-based gain index `ig` (`1 … NG`), and the
/// packed 10-bit channel index `ichan = (IS − 1)·NG + (IG − 1)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SearchResult {
    /// 1-based best shape codebook index `IS` (`1 … NCWD`).
    pub is: usize,
    /// 1-based best gain codebook index `IG` (`1 … NG`).
    pub ig: usize,
    /// Packed 10-bit channel index `(IS − 1)·NG + (IG − 1)`.
    pub ichan: u16,
}

/// Blocks 17 / 18 (§G.3.9) — VQ search error calculator and best
/// codebook index selector.
///
/// For each of the `NCWD` shape codevectors this computes the inner
/// product `P_j = Σ_k PN(K)·Y(J1+K)`, picks the best gain level for its
/// magnitude via the `GB` cell boundaries, evaluates the distortion
/// `D = GSQ(IG)·Y2 − G2(IG)·|P_j|`, and tracks the minimum. After the
/// loop the sign of the winning inner product selects between the
/// positive (`IG`) and sign-mirrored negative (`IG + 4`) gain level.
///
/// The fixed-point code (per §G.3.9) uses the *absolute value* of the
/// correlation during the search to avoid a branch in the inner loop,
/// then re-derives the sign once at the end. `pn` is the `Q7` array
/// from block 13; `y2` is the `NLSY2 = 5` energy array from block
/// 14/15. Returns the [`SearchResult`].
#[must_use]
pub fn search(pn: &[i16; IDIM], y2: &[i16; NCWD]) -> SearchResult {
    // DISTM = largest representable (double precision).
    let mut distm: i64 = i64::from(i32::MAX);
    let mut ig: usize = 1; // 1-based positive gain index
    let mut is: usize = 1; // 1-based shape index

    for j in 0..NCWD {
        // AA0 = Σ PN(K)·Y(J1+K). NLS for AA0 is 7 + 11 = 18.
        let mut aa0: i64 = 0;
        for k in 0..IDIM {
            aa0 += pn[k] as i64 * Y_Q11[j][k] as i64;
        }
        // Take absolute value (gain search uses |correlation|).
        if aa0 < 0 {
            aa0 = -aa0;
        }

        // IDXG starts at 1 (1-based); bump once per crossed positive
        // boundary GB(1), GB(2), GB(3). NLS for P is 13 + 5 = 18,
        // matching AA0's NLS so the comparison is direct.
        let mut idxg: usize = 1;
        for kb in 0..(N1 - 1) {
            let p = GB_Q13[kb] as i64 * y2[j] as i64;
            if aa0 >= p {
                idxg += 1;
            }
        }

        // AA0 = AA0 >> 14 (NLS 18 → 4), then clip to 16 bits.
        let mut aa0s = aa0 >> 14;
        if aa0s > i16::MAX as i64 {
            aa0s = i16::MAX as i64;
        }

        // AA1 = GSQ(IDXG)·Y2(J) (NLS 11 + 5 = 16);
        // P   = G2(IDXG)·AA0   (NLS 12 + 4 = 16).
        let aa1 = GSQ_Q11[idxg - 1] as i64 * y2[j] as i64;
        let p = G2_Q12[idxg - 1] as i64 * aa0s;
        let dist = aa1 - p;

        if dist < distm {
            distm = dist;
            ig = idxg;
            is = j + 1; // store 1-based shape index
        }
    }

    // Re-determine the sign of the winning correlation: if it is ≤ 0,
    // select the sign-mirrored negative gain level IG + 4.
    let j1 = is - 1; // 0-based winning shape row
    let mut aa0: i64 = 0;
    for k in 0..IDIM {
        aa0 += pn[k] as i64 * Y_Q11[j1][k] as i64;
    }
    if aa0 <= 0 {
        ig += 4;
    }

    let ichan = ((is - 1) * NG + (ig - 1)) as u16;
    SearchResult { is, ig, ichan }
}

// ---------------------------------------------------------------------
// §G.3.10 — Block 19: excitation VQ codebook + block 21: gain scaling.
// ---------------------------------------------------------------------

/// Blocks 19 / 21 (§G.3.10) — reconstruct the gain-scaled excitation
/// `ET(K) = GAIN·GQ(IG)·Y(NN+K)` from the chosen indices.
///
/// `result` is the winning [`SearchResult`]; `gain`/`nlsgain` is the
/// scalar floating-point excitation gain (block 46). The product
/// `GQ(IG)·GAIN` is normalized to 32 bits (left-shift by `NNGQ(IG)`)
/// before rounding to a 16-bit `TMP`; the selected shape codevector is
/// block-floated to 16 bits; the per-sample product is rounded to a
/// 15-bit `ET`.
///
/// Returns `(et, nlset)`: the `IDIM`-length excitation and its NLS.
#[must_use]
pub fn excitation(result: &SearchResult, gain: i16, nlsgain: i32) -> ([i16; IDIM], i32) {
    let ig = result.ig; // 1-based gain index
    let is = result.is; // 1-based shape index

    // AA0 = GQ(IG)·GAIN (Q13·gain). The §G.3.10 note guarantees AA0 has
    // NNGQ(IG) leading zeros, so the product fits a 32-bit accumulator
    // and normalizing by << NNGQ(IG) keeps it in 32 bits. Round to the
    // upper 16 bits → TMP.
    let aa0 = (GQ_Q13[ig - 1] * gain as i32) << NNGQ[ig - 1];
    let tmp = rnd(aa0);

    // NLSAA0 = 13 + NLSGAIN; NLSTMP = NLSAA0 + NNGQ(IG) − 16.
    let nlsaa0 = 13 + nlsgain;
    let nlstmp = nlsaa0 + NNGQ[ig - 1] as i32 - 16;

    // Normalize the selected Q11 shape codevector to 16 bits (TEMP).
    let nn = is - 1; // 0-based shape row
    let row: [i32; IDIM] = {
        let mut r = [0i32; IDIM];
        for k in 0..IDIM {
            r[k] = Y_Q11[nn][k] as i32;
        }
        r
    };
    let (temp, nls) = vscale(&row, IDIM, 14);

    // ET(K) = round(TMP·TEMP(K)) — both are 16-bit normalized, so the
    // product has one leading zero; rounding to the high word gives a
    // 15-bit ET.
    let mut et = [0i16; IDIM];
    for k in 0..IDIM {
        let aa0 = tmp as i32 * temp[k];
        et[k] = rnd(aa0);
    }

    // NLSET = NLSTMP + 11 + NLS − 16.
    let nlset = nlstmp + 11 + nls - 16;
    (et, nlset)
}

// ---------------------------------------------------------------------
// End-to-end §G.3 fixed-point coder driver (blocks 11 → 21).
// ---------------------------------------------------------------------

/// The full per-vector output of the §G.3 fixed-point coder: the chosen
/// codebook indices / channel index ([`SearchResult`]) plus the
/// reconstructed gain-scaled excitation `ET` and its NLS.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CoderOutput {
    /// The codebook search result (`is`, `ig`, `ichan`).
    pub result: SearchResult,
    /// Gain-scaled excitation `ET(1..IDIM)` (block 19/21).
    pub et: [i16; IDIM],
    /// `NLSET` — the block-floating-point shift of `et`.
    pub nlset: i32,
}

/// Drive the entire §G.3.4 … §G.3.10 fixed-point coder for one
/// 5-sample vector, chaining the per-block modules of this file in the
/// Figure 2/G.728 order:
///
/// ```text
/// block 11 → TARGET = SW − ZIR
/// block 12 → H (impulse response of F(z)·W(z))
/// block 14 → Y2 (filtered-shape energies)
/// block 16 → TARGET / GAIN  (sets NLSTARGET)
/// block 13 → PN (time-reversed convolution)
/// block 17 → best (IS, IG), ICHAN
/// block 19/21 → ET (gain-scaled excitation)
/// ```
///
/// Inputs:
/// * `sw` / `zir` — the weighted-speech and zero-input-response vectors,
///   `Q2` (block 11 inputs).
/// * `a` / `awz` / `awp` — the `Q14` cascade predictor coefficients
///   indexed `[i] = coeff(i)` for `i = 2 … IDIM` (length `IDIM + 1`;
///   see [`impulse_response`]).
/// * `gain` / `nlsgain` — the scalar-floating-point excitation gain σ(n)
///   (block 46 output): mantissa + NLS.
///
/// Returns the [`CoderOutput`]. This is the fixed-point analogue of the
/// floating-point [`crate::codebook_search::CodebookSearch::search`]
/// path, with the impulse response / energy table recomputed per call
/// from the supplied coefficients (the per-cycle reuse and the backward
/// adapters that produce `a`/`awz`/`awp`/`gain` live in the
/// floating-point [`crate::encoder`] and the §G.2 fixed-point adapters).
#[must_use]
pub fn encode_vector(
    sw: &[i16; IDIM],
    zir: &[i16; IDIM],
    a: &[i16; IDIM + 1],
    awz: &[i16; IDIM + 1],
    awp: &[i16; IDIM + 1],
    gain: i16,
    nlsgain: i32,
) -> CoderOutput {
    // Block 11: VQ target.
    let target = vq_target(sw, zir);
    // Block 12: cascade impulse response.
    let h = impulse_response(a, awz, awp);
    // Blocks 14/15: filtered-shape energies.
    let y2 = codevector_energy(&h);
    // Block 16: normalize the target by the gain (sets NLSTARGET).
    let (target_n, nlstarget) = normalize_target(&target, gain, nlsgain);
    // Block 13: time-reversed convolution → PN (uses NLSTARGET).
    let pn = time_reversed_convolution(&target_n, &h, nlstarget);
    // Blocks 17/18: best codebook index.
    let result = search(&pn, &y2);
    // Blocks 19/21: gain-scaled excitation.
    let (et, nlset) = excitation(&result, gain, nlsgain);

    CoderOutput { result, et, nlset }
}

// ---------------------------------------------------------------------
// Small saturation / shift helpers (§G.3 "clip" lines).
// ---------------------------------------------------------------------

/// Clip a 32-bit value to the signed 16-bit range (the §G.3 "If AA0 >
/// 32767 … If AA0 < −32768" saturation lines).
#[inline]
fn clip_i16(value: i32) -> i16 {
    if value > i16::MAX as i32 {
        i16::MAX
    } else if value < i16::MIN as i32 {
        i16::MIN
    } else {
        value as i16
    }
}

/// Clip a 64-bit accumulator to the signed 16-bit range.
#[inline]
fn clip_i16_64(value: i64) -> i16 {
    if value > i16::MAX as i64 {
        i16::MAX
    } else if value < i16::MIN as i64 {
        i16::MIN
    } else {
        value as i16
    }
}

/// Arithmetic right shift of a 64-bit accumulator by a possibly
/// negative count (a negative `shift` means a left shift). The §G.3.6
/// `>> 13 + (NLSTARGET − 7)` exponent can go negative when
/// `NLSTARGET < 7`.
#[inline]
fn shr_signed(value: i64, shift: i32) -> i64 {
    if shift >= 0 {
        value >> shift
    } else {
        value << (-shift)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tables::{G2, GB, GQ, GSQ};

    // --- Q-format gain tables -------------------------------------

    #[test]
    fn gq_q13_matches_float() {
        for i in 0..NG {
            let want = (GQ[i] * 8192.0).round() as i32;
            assert_eq!(GQ_Q13[i], want, "GQ_Q13[{i}]");
        }
    }

    #[test]
    fn gb_q13_matches_float() {
        // Only the positive half-cell boundaries are read by the
        // search, but every derived entry must match the float to one
        // Q13 quantum.
        for k in 0..(NG - 1) {
            let want = (GB[k] * 8192.0).round() as i32;
            assert_eq!(GB_Q13[k], want, "GB_Q13[{k}]");
        }
    }

    #[test]
    fn g2_q12_matches_float() {
        for i in 0..NG {
            let want = (G2[i] * 4096.0).round() as i32;
            assert_eq!(G2_Q12[i], want, "G2_Q12[{i}]");
        }
    }

    #[test]
    fn gsq_q11_matches_float() {
        for i in 0..NG {
            let want = (GSQ[i] * 2048.0).round() as i32;
            assert_eq!(GSQ_Q11[i], want, "GSQ_Q11[{i}]");
        }
    }

    #[test]
    fn nngq_normalizes_q13_gq_to_16_bits() {
        // §G.3.10: NNGQ(i) = 1 + (left shifts to normalize the Q13
        // GQ(i)). So shifting the Q13 magnitude by NNGQ(i) − 1 must
        // land it in the normalized band [16384, 32767].
        for i in 0..NG {
            let shifted = (GQ_Q13[i].abs()) << (NNGQ[i] - 1);
            assert!(
                (16384..=32767).contains(&shifted),
                "NNGQ[{i}] = {} did not normalize {} → {}",
                NNGQ[i],
                GQ_Q13[i],
                shifted
            );
        }
    }

    // --- §G.3.4 Block 11 ------------------------------------------

    #[test]
    fn vq_target_subtracts_and_clips() {
        let sw = [100, -200, 30000, -30000, 0];
        let zir = [50, 50, -10000, 10000, 0];
        let t = vq_target(&sw, &zir);
        assert_eq!(t[0], 50);
        assert_eq!(t[1], -250);
        // 30000 − (−10000) = 40000 → clipped to 32767.
        assert_eq!(t[2], 32767);
        // −30000 − 10000 = −40000 → clipped to −32768.
        assert_eq!(t[3], -32768);
        assert_eq!(t[4], 0);
    }

    // --- §G.3.5 Block 12 ------------------------------------------

    #[test]
    fn impulse_response_first_tap_is_unity_q13() {
        // With all predictor coefficients zero, the cascade filter is
        // the identity: h(1) = 1.0 (= 8192 in Q13) and h(2..IDIM) = 0.
        let zero = [0i16; IDIM + 1];
        let h = impulse_response(&zero, &zero, &zero);
        assert_eq!(h[0], 8192, "h(1) should be unity in Q13");
        for k in 1..IDIM {
            assert_eq!(h[k], 0, "h({}) should be 0", k + 1);
        }
    }

    #[test]
    fn impulse_response_matches_reference_recursion() {
        // Cross-check the in-place block-12 recursion against an
        // independent re-implementation of the §G.3.5 pseudo-code on a
        // non-trivial coefficient set. This guards the in-place memory
        // shuffle and the combined-accumulator >> 14 against an
        // off-by-one without depending on a hand-computed magic value.
        // Crate layout: `a[0]` = A(1) = 16384, `a[i − 1]` = A(i).
        let a = [16384i16, 1000, -800, 600, -400, 0];
        let awz = [16384i16, 500, -300, 200, -100, 0];
        let awp = [16384i16, -700, 400, -250, 150, 0];
        let h = impulse_response(&a, &awz, &awp);

        // Reference: TEMP / WS as Q13 16-bit words with the spec's
        // exact 1-based loop (the truncation to i16 mirrors "Q13
        // 16-bit words"); `ref_x[i]` is the spec's X(i).
        let mut ref_a = [0i64; IDIM + 2];
        let mut ref_awz = [0i64; IDIM + 2];
        let mut ref_awp = [0i64; IDIM + 2];
        for i in 1..=(IDIM + 1) {
            ref_a[i] = a[i - 1] as i64;
            ref_awz[i] = awz[i - 1] as i64;
            ref_awp[i] = awp[i - 1] as i64;
        }
        let mut temp = [0i16; IDIM + 1];
        let mut ws = [0i16; IDIM + 1];
        temp[1] = 8192;
        ws[1] = 8192;
        for k in 2..=IDIM {
            let mut acc0: i64 = 0;
            let mut acc1: i64 = 0;
            let mut i = k;
            while i >= 2 {
                temp[i] = temp[i - 1];
                ws[i] = ws[i - 1];
                acc0 -= ref_a[i] * temp[i] as i64;
                acc1 += ref_awz[i] * temp[i] as i64;
                acc1 -= ref_awp[i] * ws[i] as i64;
                i -= 1;
            }
            acc1 += acc0;
            temp[1] = (acc0 >> 14) as i16;
            ws[1] = (acc1 >> 14) as i16;
        }
        for k in 1..=IDIM {
            assert_eq!(h[k - 1], ws[IDIM + 1 - k], "h({k})");
        }
    }

    // --- §G.3.6 Block 13 ------------------------------------------

    #[test]
    fn time_reversed_convolution_unit_impulse() {
        // With H = unit impulse in Q13 (h(1) = 8192, rest 0) and
        // NLSTARGET = 7, PN(K) = TARGET(K): the K-th tap reads
        // TARGET(K)·H(1) >> 13 = TARGET(K).
        let mut h = [0i16; IDIM];
        h[0] = 8192;
        let target = [10, -20, 30, -40, 50];
        let pn = time_reversed_convolution(&target, &h, 7);
        assert_eq!(pn, target);
    }

    // --- §G.3.7 Blocks 14/15 --------------------------------------

    #[test]
    fn codevector_energy_unit_impulse_matches_q11_self_energy() {
        // With H = unit impulse in Q13, TEMP(k) = (H(1)·Y(k)) >> 14 =
        // Y_Q11[j][k] >> 1, and Y2[j] = (Σ_k TEMP(k)²) >> 15.
        let mut h = [0i16; IDIM];
        h[0] = 8192;
        let y2 = codevector_energy(&h);
        for j in 0..NCWD {
            let mut acc: i64 = 0;
            for k in 0..IDIM {
                let t = (8192i64 * Y_Q11[j][k] as i64) >> 14;
                acc += t * t;
            }
            let want = (acc >> 15) as i16;
            assert_eq!(y2[j], want, "Y2[{j}]");
        }
    }

    // --- §G.3.8 Block 16 ------------------------------------------

    #[test]
    fn normalize_target_divides_by_unit_gain() {
        // GAIN = 1.0 represented as mantissa 16384 (Q14) with
        // NLSGAIN = 14. Then 1/GAIN = 1.0 and TARGET is unchanged in
        // value (modulo the block-float renormalization).
        let target = [4096i16, -2048, 1024, -512, 256];
        let (out, nls) = normalize_target(&target, 16384, 14);
        // Reconstruct the real value of each output and compare to the
        // Q2 input value (target / 4).
        for k in 0..IDIM {
            let got = out[k] as f64 / (1i64 << nls) as f64;
            let want = target[k] as f64 / 4.0;
            assert!(
                (got - want).abs() < 1e-2,
                "normalize_target[{k}]: got {got}, want {want}"
            );
        }
    }

    // --- §G.3.9 Blocks 17/18 --------------------------------------

    #[test]
    fn search_picks_codevector_aligned_with_target() {
        // Build a PN that is a positive multiple of one shape
        // codevector's Q11 row. The inner product with that row is the
        // largest positive correlation, so the search should select it
        // (IS = that 1-based index) with a positive gain (IG ≤ N1).
        let want_is0 = 17usize; // 0-based target row
        let mut pn = [0i16; IDIM];
        for k in 0..IDIM {
            // Scale the Q11 row into the PN range, staying within i16.
            pn[k] = Y_Q11[want_is0][k] / 4;
        }
        let h = {
            let mut hh = [0i16; IDIM];
            hh[0] = 8192;
            hh
        };
        let y2 = codevector_energy(&h);
        let r = search(&pn, &y2);
        assert_eq!(r.is, want_is0 + 1, "selected shape index");
        // Positive correlation ⇒ positive gain half (IG in 1..=N1).
        assert!(r.ig <= N1, "expected positive gain, got IG = {}", r.ig);
        assert_eq!(r.ichan, ((r.is - 1) * NG + (r.ig - 1)) as u16);
    }

    #[test]
    fn search_negated_target_flips_gain_sign() {
        // Negating PN flips the sign of every correlation, so the same
        // shape is selected but with the sign-mirrored negative gain
        // level (IG bumped by 4).
        let want_is0 = 42usize;
        let mut pn = [0i16; IDIM];
        for k in 0..IDIM {
            pn[k] = -(Y_Q11[want_is0][k] / 4);
        }
        let h = {
            let mut hh = [0i16; IDIM];
            hh[0] = 8192;
            hh
        };
        let y2 = codevector_energy(&h);
        let r = search(&pn, &y2);
        assert_eq!(r.is, want_is0 + 1, "selected shape index");
        assert!(r.ig > N1, "expected negative gain, got IG = {}", r.ig);
    }

    #[test]
    fn search_ichan_is_in_range() {
        // For an arbitrary PN the packed index must always fit the
        // 10-bit channel field (0 ..= 1023).
        let pn = [123i16, -456, 789, -321, 654];
        let h = {
            let mut hh = [0i16; IDIM];
            hh[0] = 8192;
            hh
        };
        let y2 = codevector_energy(&h);
        let r = search(&pn, &y2);
        assert!(r.ichan < 1024, "ICHAN out of 10-bit range: {}", r.ichan);
        assert!((1..=NCWD).contains(&r.is));
        assert!((1..=NG).contains(&r.ig));
    }

    // --- §G.3.10 Blocks 19/21 -------------------------------------

    #[test]
    fn excitation_reconstructs_gain_scaled_shape() {
        // With GAIN = 1.0 (mantissa 16384, NLSGAIN = 14) and a chosen
        // positive gain index, ET should be ≈ GQ(IG)·Y(IS) in real
        // value. Cross-check the reconstructed real ET against the
        // floating-point GQ·Y product to a small tolerance.
        let result = SearchResult {
            is: 10,
            ig: 2,
            ichan: ((10 - 1) * NG + (2 - 1)) as u16,
        };
        let (et, nlset) = excitation(&result, 16384, 14);
        let gq = crate::tables::GQ[result.ig - 1];
        for k in 0..IDIM {
            let got = et[k] as f64 / (1i64 << nlset) as f64;
            let y = Y_Q11[result.is - 1][k] as f64 / 2048.0;
            let want = gq * y;
            assert!((got - want).abs() < 5e-3, "ET[{k}]: got {got}, want {want}");
        }
    }

    // --- Cross-equivalence with the floating-point search --------

    #[test]
    fn fixed_point_search_agrees_with_float_on_shape_index() {
        // The §G.3.9 fixed-point reformulation must reproduce the
        // decisions of the existing floating-point coder
        // ([`crate::codebook_search`]). Drive both from the cold-start
        // identity cascade (H = unit impulse) so PN(k) = TARGET(k) and
        // Y2[j] = E_j, then feed a deterministic sweep of target
        // vectors and assert the selected *shape* index matches.
        //
        // The gain index can legitimately differ by one level for
        // targets that land within a Q-quantum of a gain-cell boundary
        // (the float path uses exact GB, the fixed path the Q13 GB), so
        // the strong invariant asserted here is shape-index agreement,
        // which is what dominates perceived quality.
        use crate::codebook_search::CodebookSearch;

        let float_search = CodebookSearch::new();
        let mut h = [0i16; IDIM];
        h[0] = 8192; // unit impulse, Q13
        let y2_fixed = codevector_energy(&h);

        // Deterministic LCG target sweep; values chosen well inside the
        // Q2 range so no clipping occurs and both paths see the same
        // real-valued target.
        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut mismatches = 0;
        let trials = 256;
        for _ in 0..trials {
            // PN / TARGET in Q7 for the fixed path (identity H ⇒
            // PN = TARGET), and the matching real value for the float.
            let mut pn = [0i16; IDIM];
            let mut target_f = [0f64; IDIM];
            for k in 0..IDIM {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let v = ((state >> 40) as i32 & 0x1fff) - 4096; // [-4096, 4095]
                pn[k] = v as i16;
                // PN is Q7, so the real target value is v / 128.
                target_f[k] = v as f64 / 128.0;
            }

            // Gain = 1.0 for both. (The float search divides the target
            // by GAIN before the inner product; with GAIN = 1 the
            // normalized target equals the raw target.)
            let fr = float_search.search(&target_f, 1.0);

            // Fixed path: PN already Q7, Y2 from block 14/15.
            let xr = search(&pn, &y2_fixed);

            // Float shape_index is 0-based; fixed `is` is 1-based.
            if (xr.is - 1) as u8 != fr.shape_index {
                mismatches += 1;
            }
        }
        // Allow a small number of near-boundary disagreements; the bulk
        // must agree, proving the reformulation tracks the float coder.
        assert!(
            mismatches <= trials / 20,
            "too many shape-index mismatches: {mismatches}/{trials}"
        );
    }

    // --- End-to-end §G.3 coder driver -----------------------------

    #[test]
    fn encode_vector_chain_agrees_with_float_search_identity_cascade() {
        // With the identity cascade (all predictor coefficients zero ⇒
        // H = unit impulse) and unit gain, encode_vector must reproduce
        // the floating-point coder's shape decision for a deterministic
        // target sweep, end-to-end through all of blocks 11 → 21.
        use crate::codebook_search::CodebookSearch;

        let float_search = CodebookSearch::new();
        let zero = [0i16; IDIM + 1];

        let mut state: u64 = 0x0f1e_2d3c_4b5a_6978;
        let mut mismatches = 0;
        let trials = 256;
        for _ in 0..trials {
            // SW − ZIR = TARGET in Q2; pick SW directly, ZIR = 0.
            let mut sw = [0i16; IDIM];
            let mut target_f = [0f64; IDIM];
            for k in 0..IDIM {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let v = ((state >> 40) as i32 & 0x3ff) - 512; // [-512, 511] Q2
                sw[k] = v as i16;
                target_f[k] = v as f64 / 4.0; // Q2 → real
            }
            let zir = [0i16; IDIM];

            // Gain = 1.0 (mantissa 16384, NLSGAIN = 14).
            let out = encode_vector(&sw, &zir, &zero, &zero, &zero, 16384, 14);
            let fr = float_search.search(&target_f, 1.0);

            if (out.result.is - 1) as u8 != fr.shape_index {
                mismatches += 1;
            }
            // ICHAN must always be a valid 10-bit field.
            assert!(out.result.ichan < 1024);
        }
        assert!(
            mismatches <= trials / 20,
            "too many end-to-end shape mismatches: {mismatches}/{trials}"
        );
    }

    #[test]
    fn encode_vector_excitation_matches_standalone_excitation() {
        // The driver's ET must equal calling excitation() directly with
        // the driver's own SearchResult and gain — i.e. the chaining is
        // faithful, not silently re-deriving a different index.
        let mut sw = [0i16; IDIM];
        for k in 0..IDIM {
            sw[k] = (k as i16 - 2) * 80;
        }
        let zir = [0i16; IDIM];
        let zero = [0i16; IDIM + 1];
        let out = encode_vector(&sw, &zir, &zero, &zero, &zero, 16384, 14);
        let (et, nlset) = excitation(&out.result, 16384, 14);
        assert_eq!(out.et, et);
        assert_eq!(out.nlset, nlset);
    }

    #[test]
    fn excitation_negative_gain_flips_sign() {
        // A negative gain index (IG in N1+1..=NG) yields a sign-flipped
        // excitation relative to the matching positive index.
        let pos = SearchResult {
            is: 5,
            ig: 3,
            ichan: 0,
        };
        let neg = SearchResult {
            is: 5,
            ig: 3 + 4,
            ichan: 0,
        };
        let (ep, np) = excitation(&pos, 16384, 14);
        let (en, nn) = excitation(&neg, 16384, 14);
        assert_eq!(np, nn, "matching positive/negative NLS");
        for k in 0..IDIM {
            assert_eq!(ep[k], -en[k], "ET sign flip at {k}");
        }
    }
}
