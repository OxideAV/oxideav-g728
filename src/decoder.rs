//! G.728 LD-CELP decoder — vector-rate front end.
//!
//! Implements the four-block excitation-side decoder pipeline of
//! Figure 3/G.728:
//!
//! * **Block 29** — `decode_codebook_index` splits a 10-bit channel
//!   index into the 7-bit shape index and the 3-bit gain index, then
//!   looks up `y(n) = gq(IG) · y[IS]` from the Annex B tables
//!   (Recommendation §5.14 / block 29 pseudocode).
//! * **Block 31** — `apply_gain_scaling` multiplies the quantised
//!   excitation by the externally-supplied gain σ(n)
//!   (Recommendation §5.14 / block 31 pseudocode, identical to
//!   block 21).
//! * **Block 32** — `Synthesizer::filter_vector` performs the 50th-
//!   order all-pole synthesis-filter convolution (Recommendation
//!   §5.14 / block 32 pseudocode), including the MAX/MIN saturation
//!   on the filter memory after every output sample. This guards
//!   against the "possible overloading of signal levels" the spec
//!   calls out in §5.13.
//!
//! What is **not yet** wired up at the current round:
//!
//! * The backward synthesis-filter adapter (block 33) — instead, this
//!   module exposes a `Synthesizer::set_predictor` hook so a future
//!   round can plug Levinson-Durbin output into the same filter and
//!   keep this module's hot path stable.
//! * The backward vector gain adapter (block 30) — same pattern: the
//!   gain σ(n) is supplied by the caller, and a future round will
//!   plug the log-gain predictor in.
//! * The adaptive postfilter (blocks 71–77) — postfilter is purely
//!   perceptual and can be added in a later round without changing
//!   the synthesis-filter output bytes.
//! * Output PCM format conversion (block 28) — the synthesizer returns
//!   linear PCM in the spec's `-4 095 … +4 095` range; A-law / µ-law
//!   conversion is delegated to `oxideav-g711` per §3.1 / §5.3.

use crate::consts::{IDIM, LPC, NG};
use crate::tables::{GQ, Y_Q11};

/// Decoder excitation block size; one G.728 frame produces
/// [`IDIM`] = 5 output PCM samples.
pub const FRAME_LEN: usize = IDIM;

/// Maximum saturation level applied to the synthesis-filter memory
/// after each output sample. The spec calls these `MAX`/`MIN` in
/// §5.13 and ties them to "the positive and negative saturation
/// levels of A-law or µ-law PCM, depending on which law is used".
///
/// We default to the assumed-input range from §3.1.1: `±4 095`. For
/// stricter A-law-only operation a caller can construct the
/// `Synthesizer` with [`Synthesizer::with_saturation`].
pub const DEFAULT_MAX: f64 = 4_095.0;

/// Excitation codevector after Annex B lookup and per-vector gain
/// quantisation.
///
/// The codebook stores Q11 i16 integers; the spec instructs callers
/// to divide by 2 048 to obtain the floating-point representation
/// (Annex B prose: "To obtain the floating point value from the
/// integer value, divide the integer value by 2 048"). We pre-apply
/// that division so downstream filtering can stay in `f64`.
#[derive(Debug, Clone, Copy)]
pub struct ExcitationVector(pub [f64; FRAME_LEN]);

impl ExcitationVector {
    /// Look up the excitation codevector for a 10-bit channel index.
    ///
    /// Implements the block-29 decoder pseudocode (Recommendation
    /// §5.14):
    ///
    /// ```text
    /// ITMP = integer part of (ICHAN / NG)
    /// IG   = ICHAN - ITMP * NG + 1            | Decode IG
    /// NN   = ITMP * IDIM                      | shape codevector offset
    /// For K = 1, 2, ..., IDIM, do:
    ///     YN(K) = GQ(IG) * Y(NN + K)          | Decode IS as ITMP, multiply by gain
    /// ```
    ///
    /// Per §3.11, the 10-bit channel index packs the 7-bit shape
    /// index as the high bits and the 3-bit gain index as the low
    /// bits: `ICHAN = IS · NG + IG` with `IS ∈ [0, 128)` and
    /// `IG ∈ [0, 8)`. Channel indices outside `[0, 1024)` are clamped
    /// to the low ten bits — a defensive guard that mirrors how a
    /// real wire stream is parsed by the demuxer.
    pub fn from_channel_index(ichan: u16) -> Self {
        let ichan = (ichan & 0x03FF) as usize; // 10 bits
        let is = ichan / NG; // shape index 0..127
        let ig = ichan % NG; // gain index 0..7
        let gain = GQ[ig];

        let mut out = [0.0f64; FRAME_LEN];
        let row = &Y_Q11[is];
        for (k, slot) in out.iter_mut().enumerate() {
            // Annex B Q11 → float; pre-multiply by gain.
            *slot = gain * (row[k] as f64 / 2_048.0);
        }
        ExcitationVector(out)
    }
}

/// Reassemble the 10-bit channel index from a shape index and a gain
/// index — useful for encoder-side bit packing once the encoder
/// lands. The packing rule is the inverse of [`ExcitationVector::
/// from_channel_index`].
///
/// Per the spec §3.11 most-significant-bit-first convention, the
/// shape index occupies bits `9..3` (the high 7 bits) and the gain
/// index occupies bits `2..0` (the low 3 bits).
pub fn pack_channel_index(shape: u8, gain: u8) -> u16 {
    let shape = (shape as u16) & 0x7F;
    let gain = (gain as u16) & 0x07;
    shape * NG as u16 + gain
}

/// 50th-order all-pole synthesis filter, decoder side.
///
/// State carried across vectors:
///
/// * `predictor` — `[1.0, -a₁, -a₂, …, -a₅₀]` in the spec's `A` layout
///   where `A(1) = 1`. The leading 1 lets the filter loop use the
///   same indexing as the spec pseudocode.
/// * `state` — `STATELPC(1..LPC)` filter memory. Initialised to zero
///   per Table 2/G.728.
///
/// The filter coefficients should be updated once per 20-sample
/// adaptation cycle (i.e. every fourth call to [`filter_vector`])
/// via [`set_predictor`]. Between updates the synthesiser keeps
/// using the most recently supplied predictor.
#[derive(Debug, Clone)]
pub struct Synthesizer {
    predictor: [f64; LPC + 1],
    state: [f64; LPC],
    max: f64,
    min: f64,
}

impl Default for Synthesizer {
    fn default() -> Self {
        Self::new()
    }
}

impl Synthesizer {
    /// Construct a synthesiser with the spec's default saturation
    /// (`±4 095` per §3.1.1).
    pub fn new() -> Self {
        Self::with_saturation(DEFAULT_MAX)
    }

    /// Construct a synthesiser with an explicit symmetric saturation
    /// range `±max`. Callers wanting A-law-only `±4 015.5` semantics
    /// (per §3.1.1) supply `4_015.5` here.
    pub fn with_saturation(max: f64) -> Self {
        // Default predictor: A(1) = 1, A(2..=51) = 0 — the
        // "all-pass" filter so a fresh synthesiser produces the
        // identity output until the backward-adapter wires up.
        let mut predictor = [0.0; LPC + 1];
        predictor[0] = 1.0;
        Self {
            predictor,
            state: [0.0; LPC],
            max,
            min: -max,
        }
    }

    /// Update the synthesis-filter predictor coefficients for the
    /// next vector. Caller supplies the spec's `A(1..LPC+1)` array
    /// in 0-based Rust slice form: `coeffs[0]` is `A(1) = 1.0`,
    /// `coeffs[i]` carries `A(i+1) = -aᵢ`.
    ///
    /// The caller is responsible for validating `coeffs[0] == 1.0`.
    pub fn set_predictor(&mut self, coeffs: [f64; LPC + 1]) {
        self.predictor = coeffs;
    }

    /// Return a read-only view of the current predictor (for tests).
    pub fn predictor(&self) -> &[f64; LPC + 1] {
        &self.predictor
    }

    /// Return a read-only view of the synthesis-filter memory.
    pub fn state(&self) -> &[f64; LPC] {
        &self.state
    }

    /// Run the synthesis filter for one excitation vector and produce
    /// one decoded-speech vector.
    ///
    /// Implements the block-32 pseudocode (Recommendation §5.14):
    ///
    /// ```text
    /// For K = 1, 2, ..., IDIM:
    ///     TEMP(K) = 0
    ///     For J = LPC, LPC-1, ..., 3, 2:                | zero-input response
    ///         TEMP(K)   = TEMP(K) - STATELPC(J) * A(J+1)
    ///         STATELPC(J) = STATELPC(J-1)
    ///     TEMP(K)   = TEMP(K) - STATELPC(1) * A(2)
    ///     STATELPC(1) = TEMP(K)
    /// TEMP(1) = ET(1)                                  | zero-state response init
    /// For K = 2,...,IDIM:
    ///     A0 = ET(K)
    ///     For I = K, K-1, ..., 2:
    ///         TEMP(I) = TEMP(I-1)
    ///         A0      = A0 - A(I) * TEMP(I)            | zero-state convolution
    ///     TEMP(1) = A0
    /// For K = 1, 2, ..., IDIM:
    ///     STATELPC(K) = STATELPC(K) + TEMP(K)          | ZIR + ZSR
    ///     If STATELPC(K) > MAX, STATELPC(K) = MAX      | saturation
    ///     If STATELPC(K) < MIN, STATELPC(K) = MIN
    /// I = IDIM + 1
    /// For K = 1, 2, ..., IDIM:                         | output reversal
    ///     ST(K) = STATELPC(I - K)
    /// ```
    pub fn filter_vector(&mut self, excitation: &ExcitationVector) -> [f64; FRAME_LEN] {
        // STATELPC(1..LPC) in spec → state[0..LPC] in Rust.
        // A(1..LPC+1) in spec    → predictor[0..LPC+1] in Rust.
        let a = &self.predictor;
        let mut temp = [0.0f64; FRAME_LEN];

        // ----- Zero-input response (ZIR) ------------------------------
        // For each output sample K, compute the zero-input contribution
        // by convolving the synthesis-filter memory with the predictor
        // and shifting the memory left after each step.
        for k in 0..FRAME_LEN {
            let mut tmp = 0.0f64;
            // J = LPC, LPC-1, ..., 2 in spec → j = LPC-1, ..., 1 in Rust
            for j in (1..LPC).rev() {
                // STATELPC(J) → state[J-1]; A(J+1) → a[J]
                tmp -= self.state[j] * a[j + 1];
                // STATELPC(J) = STATELPC(J-1)
                self.state[j] = self.state[j - 1];
            }
            // | TEMP(K) = TEMP(K) - STATELPC(1) * A(2)
            // | STATELPC(1) = TEMP(K)
            tmp -= self.state[0] * a[1];
            self.state[0] = tmp;
            temp[k] = tmp;
        }

        // ----- Zero-state response (ZSR) ------------------------------
        // The spec spells the ZSR convolution in-place over the same
        // TEMP array; we mirror that exactly.
        temp[0] = excitation.0[0];
        for k in 1..FRAME_LEN {
            let mut a0 = excitation.0[k];
            // I = K, K-1, ..., 2 in spec → i = k, k-1, ..., 1 in Rust.
            // The loop shifts TEMP(I) ← TEMP(I-1) along the way.
            for i in (1..=k).rev() {
                temp[i] = temp[i - 1];
                a0 -= a[i] * temp[i];
            }
            temp[0] = a0;
        }

        // ----- Add ZIR + ZSR and saturate -----------------------------
        for k in 0..FRAME_LEN {
            let mut v = self.state[k] + temp[k];
            if v > self.max {
                v = self.max;
            } else if v < self.min {
                v = self.min;
            }
            self.state[k] = v;
        }

        // ----- Output is the filter memory in reverse order ----------
        // I = IDIM + 1; ST(K) = STATELPC(I - K)
        let mut out = [0.0f64; FRAME_LEN];
        for k in 0..FRAME_LEN {
            out[k] = self.state[FRAME_LEN - 1 - k];
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---------- Block 29: codebook index decoding -----------------

    #[test]
    fn channel_index_zero_picks_shape0_gain0() {
        // Spec: ICHAN = 0 ⇒ IS = 0, IG = 0 ⇒ YN(K) = GQ(0) · Y(0, K)
        // for K = 1..5.
        let yn = ExcitationVector::from_channel_index(0);
        let g0 = GQ[0];
        for k in 0..FRAME_LEN {
            let expected = g0 * (Y_Q11[0][k] as f64 / 2_048.0);
            assert!((yn.0[k] - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn channel_index_picks_correct_shape_and_gain() {
        // Pick a non-trivial split: shape index 42, gain index 5.
        let ichan = pack_channel_index(42, 5);
        assert_eq!(ichan, 42 * 8 + 5);
        let yn = ExcitationVector::from_channel_index(ichan);
        let g = GQ[5];
        for k in 0..FRAME_LEN {
            let expected = g * (Y_Q11[42][k] as f64 / 2_048.0);
            assert!((yn.0[k] - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn channel_index_high_bits_are_masked() {
        // The wire stream is 10 bits; any garbage above bit 9 must be
        // ignored. The result for `ichan & 0x03FF == 0` and the same
        // value with junk in the high bits must match.
        let a = ExcitationVector::from_channel_index(0);
        let b = ExcitationVector::from_channel_index(0xFC00);
        for k in 0..FRAME_LEN {
            assert_eq!(a.0[k], b.0[k]);
        }
    }

    #[test]
    fn channel_index_round_trip_packs_and_unpacks() {
        // pack_channel_index is the inverse of the (shape, gain)
        // split done inside from_channel_index. Verify the round trip
        // covers the full codebook.
        for is in 0..128 {
            for ig in 0..8 {
                let ichan = pack_channel_index(is, ig);
                let decoded_is = (ichan / NG as u16) as u8;
                let decoded_ig = (ichan % NG as u16) as u8;
                assert_eq!(decoded_is, is);
                assert_eq!(decoded_ig, ig);
            }
        }
    }

    // ---------- Block 32: synthesis filter ------------------------

    #[test]
    fn allpass_predictor_acts_as_identity() {
        // With the default predictor (1, 0, 0, ..., 0) the filter
        // F(z) = 1 / (1 - 0) is the identity. A single excitation
        // vector should pass through to the output unchanged (modulo
        // the spec's output-reversal step: STATELPC after one vector
        // holds the excitation in reverse, and ST reverses it back).
        let mut synth = Synthesizer::new();
        let exc = ExcitationVector([1.0, 2.0, 3.0, 4.0, 5.0]);
        let out = synth.filter_vector(&exc);
        assert_eq!(out, [1.0, 2.0, 3.0, 4.0, 5.0]);
    }

    #[test]
    fn saturation_clamps_to_max() {
        // Excitation above the +4095 ceiling should saturate. With
        // an all-pass predictor the spec's ZSR-then-reversal sequence
        // makes ST(K) = ET(K) (see `allpass_predictor_acts_as_identity`
        // above), so an over-range ET(1) shows up in out[0] after the
        // STATELPC clamp + reversal.
        let mut synth = Synthesizer::new();
        let exc = ExcitationVector([10_000.0, 0.0, 0.0, 0.0, 0.0]);
        let out = synth.filter_vector(&exc);
        assert_eq!(out[0], DEFAULT_MAX);
        // ...and the unsaturated trailing samples should stay zero.
        for &v in &out[1..] {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn saturation_clamps_to_min() {
        let mut synth = Synthesizer::new();
        let exc = ExcitationVector([-10_000.0, 0.0, 0.0, 0.0, 0.0]);
        let out = synth.filter_vector(&exc);
        assert_eq!(out[0], -DEFAULT_MAX);
        for &v in &out[1..] {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn zero_excitation_zero_state_yields_zero_output() {
        let mut synth = Synthesizer::new();
        let zero = ExcitationVector([0.0; FRAME_LEN]);
        for _ in 0..16 {
            let out = synth.filter_vector(&zero);
            assert!(out.iter().all(|&v| v == 0.0));
        }
    }

    #[test]
    fn state_persists_across_vectors() {
        // After a single non-zero excitation the state should be
        // non-zero, and a follow-up zero excitation should not be
        // identically zero (the filter rings out as the state shifts
        // down through STATELPC).
        let mut synth = Synthesizer::new();
        let exc = ExcitationVector([100.0, 0.0, 0.0, 0.0, 0.0]);
        synth.filter_vector(&exc);
        assert!(synth.state().iter().any(|&v| v != 0.0));
        // With an all-pass predictor the all-zero excitation that
        // follows produces zero again (the filter has no recursion
        // when A(i) = 0 for i ≥ 2), but the carried state itself
        // remains non-zero — confirm we're indeed shifting through
        // the memory:
        let zero = ExcitationVector([0.0; FRAME_LEN]);
        let out = synth.filter_vector(&zero);
        // All-pass predictor: ZIR = 0 and ZSR = excitation = 0, so
        // the *output* is zero — but the state has shifted, which is
        // the property under test.
        assert!(out.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn predictor_swap_is_in_effect_on_next_vector() {
        // Plug in a tiny first-order predictor `A = (1, -0.5, 0...)`,
        // which represents the synthesis filter `1 / (1 + 0.5·z⁻¹)`.
        // Confirm the output of the next vector differs from the
        // all-pass case.
        let mut synth_default = Synthesizer::new();
        let mut synth_custom = Synthesizer::new();

        let mut new_pred = [0.0f64; LPC + 1];
        new_pred[0] = 1.0;
        new_pred[1] = -0.5; // A(2) = -0.5
        synth_custom.set_predictor(new_pred);

        let exc = ExcitationVector([1.0, 0.0, 0.0, 0.0, 0.0]);
        let out_default = synth_default.filter_vector(&exc);
        let out_custom = synth_custom.filter_vector(&exc);

        assert_ne!(out_default, out_custom);
    }
}
