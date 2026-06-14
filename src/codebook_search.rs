//! Codebook search module — blocks 12 through 18 of Figure 2/G.728
//! (§3.9; pseudo-code §5.11).
//!
//! ## Scope
//!
//! This is the analysis-by-synthesis core of the LD-CELP encoder.
//! Given the VQ target vector `x(n)` (blocks 9 + 10 + 11, already in
//! `zero_input_response.rs`) and the predicted excitation gain `σ(n)`
//! (block 20, `gain_adapter.rs`), the module identifies the 10-bit
//! channel index whose excitation codevector minimises the
//! mean-squared error of eq. 3-16:
//!
//! ```text
//! D = ‖x(n) − x̃_ij‖² = σ²(n)·‖x̂(n) − g_i·H·y_j‖²        (3-16)
//! ```
//!
//! where `x̂(n) = x(n)/σ(n)` is the gain-normalised target, `y_j` is
//! the j-th shape codevector, `g_i` the i-th gain level and `H` the
//! lower-triangular impulse-response matrix of the cascaded filter
//! `F(z)·W(z)` (eq. 3-15). Per §3.9.1, expanding and dropping the
//! per-search constants reduces the criterion to
//!
//! ```text
//! D̂ = −b_i·P_j + c_i·E_j                                  (3-23)
//! ```
//!
//! with `b_i = 2·g_i` (eq. 3-21, table `G2`), `c_i = g_i²` (eq. 3-22,
//! table `GSQ`), `P_j = pᵀ(n)·y_j`, `p(n) = Hᵀ·x̂(n)` (eq. 3-19) and
//! `E_j = ‖H·y_j‖²` (eq. 3-20).
//!
//! ## Block split
//!
//! * **Block 12 — impulse response vector calculator.** Computes the
//!   first `IDIM = 5` samples of the impulse response of `F(z)·W(z)`
//!   by exciting the zero-memory cascade with `{1, 0, 0, 0, 0}`
//!   (§3.9.2). Executed once per adaptation cycle, when the synthesis
//!   and weighting coefficients are refreshed.
//! * **Blocks 14 + 15 — shape codevector convolution module + energy
//!   table calculator.** Convolves each of the 128 shape codevectors
//!   with `h(0..4)` (first five samples only) and stores the energy
//!   `E_j` of each result (eq. 3-20). Same cadence as block 12 —
//!   §3.9.2: "the computations in blocks 12, 14, and 15 are performed
//!   only once every four speech vectors".
//! * **Block 16 — VQ target vector normalization.** `x̂(n) = x(n)/σ(n)`,
//!   computed as a single reciprocal followed by five multiplies per
//!   the §5.11 pseudo-code (`TMP = 1/GAIN`).
//! * **Block 13 — time-reversed convolution module.** `p(n) = Hᵀ·x̂(n)`
//!   (eq. 3-19) — equivalent to reversing `x̂`, convolving with `h`,
//!   and reversing again (§3.9.2).
//! * **Blocks 17 + 18 — error calculator + best codebook index
//!   selector.** The division-free decision tree of §3.9.2 steps
//!   a)..n): for each shape index the best gain level is found by
//!   comparing `P_j` against `d_i·E_j` (the quantizer cell boundaries
//!   `GB` of Annex B, approach c) of §3.9.1), then `D̂` is evaluated
//!   once per shape and the running minimum tracked.
//!
//! Per the §5.11 pseudo-code header, blocks 12/14/15 run once per
//! adaptation cycle (`ICOUNT = 3`, after the new `A` / `AWZ` / `AWP`
//! coefficient sets are ready) while blocks 16/13/17/18 run on every
//! speech vector. [`CodebookSearch`] mirrors that split:
//! [`CodebookSearch::update_impulse_response`] is the per-cycle half,
//! [`CodebookSearch::search`] the per-vector half.
//!
//! ## Initial state
//!
//! Table 2/G.728 lists the initial `H` as `1, 0, 0, 0, 0` (the
//! impulse response of the cold-start all-pass cascade) and the
//! initial `Y2` as "Energy of `y_j`" — exactly the precomputed
//! [`Y_ENERGY`] table, since `H = identity` makes `‖H·y_j‖² = ‖y_j‖²`.

use crate::consts::{IDIM, LPC, NCWD, NG};
use crate::tables::{Y_ENERGY, Y_Q11};
use crate::weighting_filter_coeff::WeightingFilterCoeff;

/// Outcome of one §3.9 codebook search (blocks 17 + 18 outputs).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SearchResult {
    /// Best 3-bit gain codebook index, 0-based (`IG − 1` in the
    /// spec's 1-based pseudo-code; indexes [`crate::tables::GQ`]).
    pub gain_index: u8,
    /// Best 7-bit shape codebook index, 0-based (`IS − 1`; indexes
    /// [`crate::tables::Y_Q11`]).
    pub shape_index: u8,
    /// The concatenated 10-bit channel index
    /// `ICHAN = (IS − 1)·NG + (IG − 1)` (§5.11) — shape in the high
    /// seven bits, gain in the low three.
    pub channel_index: u16,
}

/// Blocks 12–18 — impulse response, filtered-shape energy table and
/// the per-vector best-codevector search.
///
/// Carries the two per-cycle arrays of Table 2/G.728: `H(1..IDIM)`
/// (impulse response of the cascade, initial `1,0,0,0,0`) and
/// `Y2(1..NCWD)` (energy of the convolved shape codevectors, initial
/// = energy of `y_j` since the cold-start cascade is the identity).
#[derive(Debug, Clone)]
pub struct CodebookSearch {
    /// `H(1..IDIM)` — first five samples of the impulse response of
    /// `F(z)·W(z)`. `h[0]` is the spec's `H(1) = h(0)`.
    h: [f64; IDIM],
    /// `Y2(1..NCWD)` — `E_j = ‖H·y_j‖²` per eq. 3-20. `y2[0]` is the
    /// spec's `Y2(1)`.
    y2: [f64; NCWD],
}

impl Default for CodebookSearch {
    fn default() -> Self {
        Self::new()
    }
}

impl CodebookSearch {
    /// Construct at the Table 2/G.728 initial state: `H = 1,0,0,0,0`
    /// and `Y2 = energy of y_j` (= [`Y_ENERGY`], because the identity
    /// impulse response makes the convolution a pass-through).
    pub fn new() -> Self {
        let mut h = [0.0f64; IDIM];
        h[0] = 1.0;
        Self { h, y2: Y_ENERGY }
    }

    /// Read-only view of the impulse response vector `H(1..IDIM)`.
    pub fn impulse_response(&self) -> &[f64; IDIM] {
        &self.h
    }

    /// Read-only view of the filtered-shape energy table
    /// `Y2(1..NCWD)` (eq. 3-20).
    pub fn energy_table(&self) -> &[f64; NCWD] {
        &self.y2
    }

    /// Blocks 12 + 14 + 15 — recompute the impulse response vector of
    /// the cascaded filter `F(z)·W(z)` and refresh the filtered-shape
    /// energy table. Per §5.11 this runs once per adaptation cycle,
    /// after the new `A` / `AWZ` / `AWP` coefficient sets are ready.
    ///
    /// `a` is the order-50 synthesis predictor in the crate's
    /// canonical `A` layout (`a[0] = A(1) = 1`, `A(J) → a[J − 1]`) —
    /// what [`crate::SynthesisAdapter::coefficients`] returns. `w`
    /// carries the block-38 weighting coefficients: `w.q_gamma1` is
    /// the spec's `AWZ` (numerator) and `w.q_gamma2` the spec's `AWP`
    /// (denominator), both with the implicit unity tap at index 0.
    pub fn update_impulse_response(&mut self, a: &[f64; LPC + 1], w: &WeightingFilterCoeff) {
        let awz = &w.q_gamma1; // spec AWZ(J) → awz[J - 1]
        let awp = &w.q_gamma2; // spec AWP(J) → awp[J - 1]

        // ===== Block 12 (§5.11): impulse response of F(z)·W(z) ======
        // Excite the zero-memory cascade with {1, 0, 0, 0, 0} and
        // collect the first IDIM output samples (§3.9.2). The
        // pseudo-code keeps two IDIM-deep histories with the most
        // recent sample at index 1: TEMP = synthesis-filter output
        // history, RC = cascade (weighting-filter) output history.
        //
        // | TEMP(1) = 1                     | first synthesis output
        // | RC(1) = 1                       | first cascade output
        // | For K = 2,3,..,IDIM:
        // |   A0 = 0; A1 = 0; A2 = 0
        // |   For I = K,K−1,..,3,2, do the next five lines
        // |     TEMP(I) = TEMP(I − 1)       | shift histories
        // |     RC(I) = RC(I − 1)
        // |     A0 = A0 − A(I)·TEMP(I)      | synthesis feedback
        // |     A1 = A1 + AWZ(I)·TEMP(I)    | weighting numerator
        // |     A2 = A2 − AWP(I)·RC(I)      | weighting denominator
        // |   TEMP(1) = A0
        // |   RC(1) = A0 + A1 + A2
        //
        // For K = 1 the impulse enters with no history: the synthesis
        // output is the impulse itself (A(1) = 1) and the cascade
        // output is likewise 1 (both unity leading taps) — the two
        // seed lines above. The zero input for K ≥ 2 makes A0 carry
        // only the feedback terms.
        let mut temp = [0.0f64; IDIM];
        let mut rc = [0.0f64; IDIM];
        temp[0] = 1.0;
        rc[0] = 1.0;
        for k in 2..=IDIM {
            let mut a0 = 0.0f64;
            let mut a1 = 0.0f64;
            let mut a2 = 0.0f64;
            for i in (2..=k).rev() {
                temp[i - 1] = temp[i - 2];
                rc[i - 1] = rc[i - 2];
                a0 -= a[i - 1] * temp[i - 1];
                a1 += awz[i - 1] * temp[i - 1];
                a2 -= awp[i - 1] * rc[i - 1];
            }
            temp[0] = a0;
            rc[0] = a0 + a1 + a2;
        }
        // | ITMP = IDIM + 1                 | h(n) is the cascade
        // | For K = 1,2,..,IDIM:            | output history read in
        // |   H(K) = RC(ITMP − K)           | reverse (oldest-first)
        for k in 1..=IDIM {
            self.h[k - 1] = rc[IDIM - k];
        }

        // ===== Blocks 14 + 15 (§5.11): convolution + energy ==========
        // | For J = 1,2,..,NCWD:
        // |   J1 = (J − 1)·IDIM
        // |   For K = 1,2,..,IDIM:
        // |     K1 = J1 + K + 1
        // |     TEMP(K) = 0
        // |     For I = 1,2,..,K:
        // |       TEMP(K) = TEMP(K) + H(I)·Y(K1 − I)    | convolution
        // |   Y2(J) = 0
        // |   For K = 1,2,..,IDIM:
        // |     Y2(J) = Y2(J) + TEMP(K)·TEMP(K)         | energy
        //
        // The flat Y(K1 − I) lookup resolves to component K + 1 − I
        // (1-based) of row J, i.e. y_row[k − i] in 0-based terms. The
        // shape rows are Annex B Q11 integers; divide by 2¹¹ = 2048
        // for the float value, the same convention every other
        // consumer of Y_Q11 uses.
        for (j, y2_slot) in self.y2.iter_mut().enumerate() {
            let row = &Y_Q11[j];
            let mut energy = 0.0f64;
            for k in 1..=IDIM {
                let mut conv = 0.0f64;
                for i in 1..=k {
                    conv += self.h[i - 1] * (row[k - i] as f64 / 2_048.0);
                }
                energy += conv * conv;
            }
            *y2_slot = energy;
        }
    }

    /// Blocks 16 + 13 + 17 + 18 — normalise the target, compute the
    /// correlation vector `p(n)` and run the §3.9.2 division-free
    /// search for the best (gain, shape) pair. Runs on every speech
    /// vector.
    ///
    /// `target` is the VQ target `x(n)` from block 11
    /// ([`crate::ZeroInputResponse::compute_target`]); `gain` is the
    /// predicted excitation gain `σ(n)` from block 20
    /// ([`crate::GainAdapter::predict_next`]). `σ(n)` is strictly
    /// positive by construction (block 48 is `10^(GAIN/20)` with the
    /// block-47 clamp to `[0, 60]` dB), so the block-16 reciprocal is
    /// always well-defined.
    pub fn search(&self, target: &[f64; IDIM], gain: f64) -> SearchResult {
        // Full codebook: shape indices 1..=NCWD (spec 1-based).
        self.search_range(target, gain, 1, NCWD)
    }

    /// §3.11 in-band-signalling half-codebook search (bit robbing).
    ///
    /// For an every-`N`-th "robbed" vector the encoder must search
    /// through **only half** of the shape codebook so that the leftmost
    /// (most-significant) bit of the 7-bit shape index is free to carry
    /// a synchronization / in-band-signalling bit (§3.11). Per the
    /// recommendation's convention:
    ///
    /// * desired bit `0` ⇒ search the **first** half (shape indices
    ///   `0..=63`), so the emitted index lies in `0..=63`;
    /// * desired bit `1` ⇒ search the **second** half (shape indices
    ///   `64..=127`), so the emitted index lies in `64..=127`.
    ///
    /// Because the seven shape bits precede the three sign-and-gain
    /// bits, that leftmost shape bit is the leftmost bit of the whole
    /// codeword (§3.11), which is the property in-band signalling needs.
    /// Both the encoder (here) and the decoder
    /// ([`crate::extract_sync_bit`]) must know which vectors are robbed,
    /// otherwise the two ends decode different excitation codevectors
    /// for those vectors (§3.11 "the encoder has to know which speech
    /// vectors will be robbed … Otherwise, the decoder will not have the
    /// same decoded excitation codevectors").
    ///
    /// The gain search is unaffected — the same blocks 17 + 18 decision
    /// tree runs; only the shape scan range narrows. The returned
    /// `channel_index` is therefore directly transmittable: its top
    /// shape bit equals `bit`.
    pub fn search_with_sync_bit(&self, target: &[f64; IDIM], gain: f64, bit: bool) -> SearchResult {
        // Spec 1-based shape range: first half = IS 1..=64 (0-based
        // 0..=63), second half = IS 65..=128 (0-based 64..=127). NCWD =
        // 128, so the split is at NCWD / 2 = 64.
        let half = NCWD / 2;
        let (start, end) = if bit { (half + 1, NCWD) } else { (1, half) };
        self.search_range(target, gain, start, end)
    }

    /// Core blocks 16 + 13 + 17 + 18 search over a 1-based shape index
    /// range `[start, end]` (inclusive). The full-codebook
    /// [`Self::search`] passes `1..=NCWD`; the §3.11 half-codebook
    /// [`Self::search_with_sync_bit`] passes one half.
    fn search_range(
        &self,
        target: &[f64; IDIM],
        gain: f64,
        start: usize,
        end: usize,
    ) -> SearchResult {
        // ===== Block 16 (§5.11): target normalization ================
        // | TMP = 1 / GAIN
        // | For K = 1,2,..,IDIM: TARGET(K) = TARGET(K) * TMP
        let tmp = 1.0 / gain;
        let mut xhat = [0.0f64; IDIM];
        for k in 0..IDIM {
            xhat[k] = target[k] * tmp;
        }

        // ===== Block 13 (§5.11): time-reversed convolution ===========
        // | For K = 1,2,..,IDIM:
        // |   K1 = K − 1
        // |   PN(K) = 0
        // |   For J = K,K+1,..,IDIM:
        // |     PN(K) = PN(K) + TARGET(J)·H(J − K1)
        //
        // H(J − K + 1) → h[j − k] in 0-based terms: p(n) = Hᵀ·x̂(n)
        // per eq. 3-19 (H is lower-triangular, eq. 3-15, so the
        // transpose pairs sample K with targets K..IDIM).
        let mut pn = [0.0f64; IDIM];
        for k in 1..=IDIM {
            for j in k..=IDIM {
                pn[k - 1] += xhat[j - 1] * self.h[j - k];
            }
        }

        // ===== Blocks 17 + 18 (§5.11): the search loop ===============
        // | Initialize DISTM to the largest representable number
        // | N1 = NG/2
        // | For J = 1,2,..,NCWD:
        // |   COR = Σ_K PN(K)·Y(J1 + K)                | P_j
        // |   If COR > 0: IDXG = N1; for K = 1..N1−1:
        // |     if COR < GB(K)·Y2(J) { IDXG = K; break }
        // |   If COR ≤ 0: IDXG = NG; for K = N1+1..NG−1:
        // |     if COR > GB(K)·Y2(J) { IDXG = K; break }
        // |   D = −G2(IDXG)·COR + GSQ(IDXG)·Y2(J)      | eq. 3-23
        // |   If D < DISTM { DISTM = D; IG = IDXG; IS = J }
        // | ICHAN = (IS − 1)·NG + (IG − 1)
        //
        // The two gain branches are §3.9.1 approach c): comparing
        // `P_j` against the precomputed cell boundaries `d_i·E_j`
        // avoids the per-shape division of the optimal-gain form
        // `ĝ = P_j/E_j`. The positive branch scans the three positive
        // boundaries d0..d2 (GB(1..3)); the negative branch scans the
        // mirrored d4..d6 (GB(5..7)); GB(4) is the spec's "any
        // arbitrary value (not used)" placeholder. Note the
        // pseudo-code routes the exact COR = 0 case to the negative
        // branch ("If COR ≤ 0") — the §3.9.2 prose step d) words the
        // boundary case the other way, but the two candidates have
        // mirrored gains and identical distortion at COR = 0, so the
        // choice only fixes which of the two equivalent indices is
        // emitted. We follow the pseudo-code.
        let n1 = NG / 2;
        let mut distm = f64::INFINITY;
        let mut ig = 1usize; // spec 1-based IG
        let mut is = start; // spec 1-based IS (range start)
        for j in start..=end {
            let row = &Y_Q11[j - 1];
            let mut cor = 0.0f64;
            for k in 0..IDIM {
                cor += pn[k] * (row[k] as f64 / 2_048.0);
            }

            let mut idxg = if cor > 0.0 { n1 } else { NG };
            if cor > 0.0 {
                for k in 1..=(n1 - 1) {
                    if cor < crate::tables::GB[k - 1] * self.y2[j - 1] {
                        idxg = k;
                        break;
                    }
                }
            } else {
                for k in (n1 + 1)..=(NG - 1) {
                    if cor > crate::tables::GB[k - 1] * self.y2[j - 1] {
                        idxg = k;
                        break;
                    }
                }
            }

            let d =
                -crate::tables::G2[idxg - 1] * cor + crate::tables::GSQ[idxg - 1] * self.y2[j - 1];
            if d < distm {
                distm = d;
                ig = idxg;
                is = j;
            }
        }

        let channel_index = ((is - 1) * NG + (ig - 1)) as u16;
        SearchResult {
            gain_index: (ig - 1) as u8,
            shape_index: (is - 1) as u8,
            channel_index,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::LPCW;
    use crate::decoder::pack_channel_index;
    use crate::tables::{G2, GB, GQ, GSQ};

    /// All-pass synthesis predictor (`A(1) = 1`, feedback taps zero).
    fn allpass_a() -> [f64; LPC + 1] {
        let mut a = [0.0f64; LPC + 1];
        a[0] = 1.0;
        a
    }

    /// Float view of one shape codebook row.
    fn y_row(j: usize) -> [f64; IDIM] {
        let mut out = [0.0f64; IDIM];
        for k in 0..IDIM {
            out[k] = Y_Q11[j][k] as f64 / 2_048.0;
        }
        out
    }

    #[test]
    fn fresh_state_matches_table_2_initials() {
        // Table 2/G.728: H initial = 1,0,0,0,0; Y2 initial = energy
        // of y_j (identity cascade).
        let cs = CodebookSearch::new();
        assert_eq!(cs.impulse_response(), &[1.0, 0.0, 0.0, 0.0, 0.0]);
        for j in 0..NCWD {
            assert_eq!(cs.energy_table()[j], Y_ENERGY[j], "j={j}");
        }
    }

    #[test]
    fn default_matches_new() {
        let a = CodebookSearch::default();
        let b = CodebookSearch::new();
        assert_eq!(a.impulse_response(), b.impulse_response());
        assert_eq!(a.energy_table(), b.energy_table());
    }

    #[test]
    fn allpass_cascade_keeps_identity_impulse_response() {
        // With F(z) all-pass and W(z) = 1, the cascade is the
        // identity: exciting it with {1,0,0,0,0} returns the impulse
        // unchanged, and the energy table stays at the bare shape
        // energies.
        let mut cs = CodebookSearch::new();
        cs.update_impulse_response(&allpass_a(), &WeightingFilterCoeff::disabled());
        assert_eq!(cs.impulse_response(), &[1.0, 0.0, 0.0, 0.0, 0.0]);
        for j in 0..NCWD {
            assert!((cs.energy_table()[j] - Y_ENERGY[j]).abs() < 1e-12, "j={j}");
        }
    }

    #[test]
    fn single_pole_synthesis_gives_geometric_impulse_response() {
        // F(z) = 1 / (1 − 0.5·z⁻¹) means A(2) = −0.5 in the spec's
        // `A` layout (output = input − Σ A(I)·past). Its impulse
        // response is 0.5^n; with W(z) = 1 the cascade inherits it.
        let mut a = allpass_a();
        a[1] = -0.5; // A(2)
        let mut cs = CodebookSearch::new();
        cs.update_impulse_response(&a, &WeightingFilterCoeff::disabled());
        let h = cs.impulse_response();
        for (n, &hv) in h.iter().enumerate() {
            let expect = 0.5f64.powi(n as i32);
            assert!(
                (hv - expect).abs() < 1e-15,
                "h[{n}] = {hv}, expected {expect}"
            );
        }
    }

    #[test]
    fn pure_fir_weighting_reproduces_numerator_taps() {
        // With F(z) all-pass and the weighting denominator forced to
        // unity, the cascade is the FIR numerator alone — its impulse
        // response IS the numerator tap vector. Build a
        // WeightingFilterCoeff by hand (fields are public) with a
        // known 4-tap numerator.
        let mut q_gamma1 = [0.0f64; LPCW + 1];
        q_gamma1[0] = 1.0;
        q_gamma1[1] = -0.9;
        q_gamma1[2] = 0.81;
        q_gamma1[3] = -0.729;
        let mut q_gamma2 = [0.0f64; LPCW + 1];
        q_gamma2[0] = 1.0;
        let w = WeightingFilterCoeff { q_gamma1, q_gamma2 };

        let mut cs = CodebookSearch::new();
        cs.update_impulse_response(&allpass_a(), &w);
        let h = cs.impulse_response();
        assert!((h[0] - 1.0).abs() < 1e-15);
        assert!((h[1] - (-0.9)).abs() < 1e-15);
        assert!((h[2] - 0.81).abs() < 1e-15);
        assert!((h[3] - (-0.729)).abs() < 1e-15);
        assert!((h[4] - 0.0).abs() < 1e-15);
    }

    #[test]
    fn energy_table_matches_direct_convolution() {
        // Cross-check blocks 14 + 15 against an independent
        // computation of ‖H·y_j‖² (eq. 3-20) done with the matrix
        // form of eq. 3-15 (lower-triangular convolution) on the
        // float view of the codebook.
        let mut a = allpass_a();
        a[1] = -0.6;
        a[2] = 0.2;
        let mut q = [0.0f64; LPCW + 1];
        q[0] = 1.0;
        q[1] = -0.4;
        q[2] = 0.1;
        let w = WeightingFilterCoeff::from_lpc(&q);

        let mut cs = CodebookSearch::new();
        cs.update_impulse_response(&a, &w);
        let h = cs.impulse_response();

        for j in 0..NCWD {
            let row = y_row(j);
            let mut energy = 0.0f64;
            for n in 0..IDIM {
                // (H·y)_n = Σ_{m ≤ n} h(n − m)·y(m)  (eq. 3-15)
                let mut s = 0.0f64;
                for m in 0..=n {
                    s += h[n - m] * row[m];
                }
                energy += s * s;
            }
            assert!(
                (cs.energy_table()[j] - energy).abs() < 1e-9 * energy.max(1.0),
                "j={j}: table {} direct {}",
                cs.energy_table()[j],
                energy
            );
        }
    }

    #[test]
    fn search_recovers_exact_codebook_entries_with_identity_cascade() {
        // With H = identity (fresh state) and σ(n) = 1, a target that
        // IS a scaled codevector g_i·y_j must come back as (i, j) with
        // zero distortion — eq. 3-16 has an exact zero at the true
        // pair, and no other (gain, shape) pair can beat zero.
        let cs = CodebookSearch::new();
        for &(ig, is) in &[(0usize, 0usize), (2, 17), (3, 127), (1, 64)] {
            let row = y_row(is);
            let mut target = [0.0f64; IDIM];
            for k in 0..IDIM {
                target[k] = GQ[ig] * row[k];
            }
            let res = cs.search(&target, 1.0);
            assert_eq!(res.shape_index as usize, is, "shape for ({ig},{is})");
            assert_eq!(res.gain_index as usize, ig, "gain for ({ig},{is})");
        }
    }

    #[test]
    fn search_recovers_negative_gain_levels() {
        // Indices 4..7 of GQ are the sign-mirrored negatives. A
        // target of −g·y_j with positive g must route through the
        // §5.11 negative branch (COR ≤ 0) and recover the mirrored
        // index.
        let cs = CodebookSearch::new();
        let row = y_row(9);
        for ig in 4usize..8 {
            let mut target = [0.0f64; IDIM];
            for k in 0..IDIM {
                target[k] = GQ[ig] * row[k];
            }
            let res = cs.search(&target, 1.0);
            assert_eq!(res.shape_index, 9, "shape for ig={ig}");
            assert_eq!(res.gain_index as usize, ig, "gain for ig={ig}");
        }
    }

    #[test]
    fn search_normalizes_by_gain_before_matching() {
        // Block 16: the target is divided by σ(n) before the search,
        // so a target of σ·g_i·y_j with σ ≠ 1 recovers the same pair
        // as the σ = 1 case.
        let cs = CodebookSearch::new();
        let row = y_row(42);
        let sigma = 137.5;
        let mut target = [0.0f64; IDIM];
        for k in 0..IDIM {
            target[k] = sigma * GQ[1] * row[k];
        }
        let res = cs.search(&target, sigma);
        assert_eq!(res.shape_index, 42);
        assert_eq!(res.gain_index, 1);
    }

    #[test]
    fn channel_index_packing_matches_decoder_convention() {
        // ICHAN = (IS − 1)·NG + (IG − 1) (§5.11) must agree with the
        // decoder-side pack_channel_index helper so the encoder's
        // output feeds Decoder::decode_vector unchanged.
        let cs = CodebookSearch::new();
        let row = y_row(100);
        let mut target = [0.0f64; IDIM];
        for k in 0..IDIM {
            target[k] = GQ[2] * row[k];
        }
        let res = cs.search(&target, 1.0);
        assert_eq!(
            res.channel_index,
            pack_channel_index(res.shape_index, res.gain_index)
        );
        assert!(res.channel_index < 1024);
    }

    #[test]
    fn decision_tree_matches_brute_force_over_all_1024_pairs() {
        // The §3.9.2 decision tree (boundary comparisons + eq. 3-23)
        // must select the same (i, j) as a brute-force minimisation
        // of the full distortion D̂ = −b_i·P_j + c_i·E_j over all
        // 1024 pairs — and that in turn must agree with the raw MSE
        // form ‖x̂ − g_i·H·y_j‖² of eq. 3-16 up to the constant
        // ‖x̂‖² term. Use a non-trivial cascade and a generic target
        // so no tie shows up.
        let mut a = allpass_a();
        a[1] = -0.55;
        a[2] = 0.18;
        a[3] = -0.05;
        let mut q = [0.0f64; LPCW + 1];
        q[0] = 1.0;
        q[1] = -0.35;
        q[2] = 0.12;
        let w = WeightingFilterCoeff::from_lpc(&q);
        let mut cs = CodebookSearch::new();
        cs.update_impulse_response(&a, &w);
        let h = *cs.impulse_response();

        let target = [13.7, -42.1, 8.9, 27.3, -19.6];
        let sigma = 3.25;
        let res = cs.search(&target, sigma);

        // Brute force on the raw eq. 3-16 distortion.
        let mut xhat = [0.0f64; IDIM];
        for k in 0..IDIM {
            xhat[k] = target[k] / sigma;
        }
        let mut best = f64::INFINITY;
        let mut best_pair = (0usize, 0usize);
        for j in 0..NCWD {
            let row = y_row(j);
            // filtered = H·y_j (lower-triangular, eq. 3-15)
            let mut filtered = [0.0f64; IDIM];
            for n in 0..IDIM {
                for m in 0..=n {
                    filtered[n] += h[n - m] * row[m];
                }
            }
            for i in 0..NG {
                let mut d = 0.0f64;
                for k in 0..IDIM {
                    let e = xhat[k] - GQ[i] * filtered[k];
                    d += e * e;
                }
                if d < best {
                    best = d;
                    best_pair = (i, j);
                }
            }
        }
        assert_eq!(
            (res.gain_index as usize, res.shape_index as usize),
            best_pair,
            "decision tree disagrees with brute-force eq. 3-16 minimum"
        );
    }

    #[test]
    fn gain_boundary_tables_satisfy_annex_b_relations() {
        // The decision tree leans on GB being the mid-points and
        // G2/GSQ the doubled/squared gain levels (eqs. 3-21/3-22).
        // Pin the relations the search assumes.
        for i in 0..3 {
            assert!((GB[i] - (GQ[i] + GQ[i + 1]) / 2.0).abs() < 1e-15, "GB[{i}]");
            assert!(
                (GB[i + 4] - (GQ[i + 4] + GQ[i + 5]) / 2.0).abs() < 1e-15,
                "GB[{}]",
                i + 4
            );
        }
        for i in 0..NG {
            assert!((G2[i] - 2.0 * GQ[i]).abs() < 1e-15, "G2[{i}]");
            assert!((GSQ[i] - GQ[i] * GQ[i]).abs() < 1e-15, "GSQ[{i}]");
        }
    }

    #[test]
    fn zero_target_routes_to_negative_branch_per_pseudocode() {
        // §5.11 sends COR = 0 to the "If COR ≤ 0" branch. With a zero
        // target every COR is 0, every distortion is GSQ(IDXG)·Y2(J),
        // and the negative branch scans GB(5..7) (all negative), so
        // 0 > GB(5)·Y2 fires immediately → IDXG = 5 (1-based) =
        // gain_index 4. The shape minimising GSQ(5)·Y2(J) is the
        // minimum-energy filtered codevector.
        let cs = CodebookSearch::new();
        let res = cs.search(&[0.0; IDIM], 1.0);
        assert_eq!(res.gain_index, 4);
        let mut min_j = 0usize;
        for j in 1..NCWD {
            if cs.energy_table()[j] < cs.energy_table()[min_j] {
                min_j = j;
            }
        }
        assert_eq!(res.shape_index as usize, min_j);
    }

    #[test]
    fn half_codebook_search_confines_shape_index_to_the_requested_half() {
        // §3.11: bit 0 ⇒ shape ∈ 0..=63, bit 1 ⇒ shape ∈ 64..=127.
        // Drive a generic non-trivial cascade + target so the unrobbed
        // optimum could fall in either half, then confirm each robbed
        // search stays in its half and the top shape bit equals the bit.
        let mut a = allpass_a();
        a[1] = -0.55;
        a[2] = 0.18;
        a[3] = -0.05;
        let mut q = [0.0f64; LPCW + 1];
        q[0] = 1.0;
        q[1] = -0.35;
        q[2] = 0.12;
        let w = WeightingFilterCoeff::from_lpc(&q);
        let mut cs = CodebookSearch::new();
        cs.update_impulse_response(&a, &w);

        let target = [13.7, -42.1, 8.9, 27.3, -19.6];
        let sigma = 3.25;

        let r0 = cs.search_with_sync_bit(&target, sigma, false);
        assert!(
            r0.shape_index < 64,
            "bit0 shape {} not in 0..64",
            r0.shape_index
        );
        assert!(!crate::extract_sync_bit(r0.channel_index));

        let r1 = cs.search_with_sync_bit(&target, sigma, true);
        assert!(
            r1.shape_index >= 64,
            "bit1 shape {} not in 64..128",
            r1.shape_index
        );
        assert!(crate::extract_sync_bit(r1.channel_index));
    }

    #[test]
    fn half_codebook_search_finds_the_best_in_half_by_brute_force() {
        // The half-codebook search must select the same (i, j) as a
        // brute-force minimisation of eq. 3-16 restricted to the chosen
        // half — only the shape scan range narrows; the gain decision
        // tree is unchanged.
        let mut a = allpass_a();
        a[1] = -0.4;
        a[2] = 0.2;
        let w = WeightingFilterCoeff::from_lpc(&{
            let mut q = [0.0f64; LPCW + 1];
            q[0] = 1.0;
            q[1] = -0.3;
            q
        });
        let mut cs = CodebookSearch::new();
        cs.update_impulse_response(&a, &w);
        let h = *cs.impulse_response();
        let target = [5.1, -8.3, 2.7, 11.0, -4.4];
        let sigma = 2.0;

        for bit in [false, true] {
            let res = cs.search_with_sync_bit(&target, sigma, bit);
            let (lo, hi) = if bit { (64usize, 128usize) } else { (0, 64) };
            let mut xhat = [0.0f64; IDIM];
            for k in 0..IDIM {
                xhat[k] = target[k] / sigma;
            }
            let mut best = f64::INFINITY;
            let mut best_pair = (0usize, 0usize);
            for j in lo..hi {
                let row = y_row(j);
                let mut filtered = [0.0f64; IDIM];
                for n in 0..IDIM {
                    for m in 0..=n {
                        filtered[n] += h[n - m] * row[m];
                    }
                }
                for i in 0..NG {
                    let mut d = 0.0f64;
                    for k in 0..IDIM {
                        let e = xhat[k] - GQ[i] * filtered[k];
                        d += e * e;
                    }
                    if d < best {
                        best = d;
                        best_pair = (i, j);
                    }
                }
            }
            assert_eq!(
                (res.gain_index as usize, res.shape_index as usize),
                best_pair,
                "half-search bit={bit} disagrees with brute force"
            );
        }
    }
}
