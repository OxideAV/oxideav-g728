//! ITU-T G.728 LD-CELP encoder — typed front-end scaffold.
//!
//! Round 235 lands the **typed encoder surface** ahead of the
//! analysis-by-synthesis pipeline that §3.9 of the Recommendation
//! describes. No encoder logic is implemented yet; this module exposes
//! a stable [`Encoder`] type and a [`make_encoder`] factory consistent
//! with the workspace's dual-API convention so consumers can wire the
//! encoder symbol now and the pipeline can fill in behind it.
//!
//! ## Scope
//!
//! The G.728 encoder reads 16-bit linear PCM (or A/µ-law per §3.1), runs
//! the per-vector analysis-by-synthesis loop of §3.9, and emits one
//! 10-bit channel index per `IDIM = 5`-sample vector. The shared
//! backward adapters (synthesis-filter adapter §3.7, gain adapter §3.8,
//! both already transcribed in the decoder modules) are reused
//! unchanged — see §4.4 / §4.5 of the Recommendation.
//!
//! What §3.9 splits into precomputable constants:
//!
//! - `Y_ENERGY[j] = Σ_k y_j(k)²` (the per-shape codevector energy) —
//!   landed this round in `tables.rs` as a const-derived
//!   `[f64; NCWD]` array (eq. 3-23). The per-test in `tables.rs`
//!   asserts it equals the direct dot product on `y_f64()` to machine
//!   precision.
//! - `G2[i] = 2·gq[i]` and `GSQ[i] = gq[i]²` — already in `tables.rs`
//!   from r189, used by the `b_i = 2·g_i`, `c_i = g_i²` substitution
//!   in eq. 3-21 / 3-22 of the gain-quantiser decision tree.
//! - `GB[i]` — the mid-point boundaries already in `tables.rs`, used
//!   by the division-free gain decision in §3.9.2.
//!
//! All four tables are constants of the codebook; the encoder will
//! consume them by reference without recomputation.
//!
//! ## What is NOT implemented
//!
//! The full §3.9 search loop (perceptual weighting filter + impulse
//! response + zero-state response correlation + shape-codebook scan +
//! gain-quantiser decision tree) is intentionally absent. Every
//! `encode_*` call returns [`Error::NotImplemented`]. The shared
//! backward adapters ([`crate::SynthesisAdapter`], [`crate::GainAdapter`])
//! are accessible through the encoder's typed surface so future rounds
//! can wire blocks 1..28 + 67..70 against them without restructuring
//! this module.

use crate::consts::NCWD;
use crate::decoder::FRAME_LEN;
use crate::gain_adapter::GainAdapter;
use crate::synthesis_adapter::SynthesisAdapter;
use crate::tables::Y_ENERGY;
use crate::weighting_filter_coeff::WeightingFilterCoeff;
use crate::Error;

/// Number of bits in one G.728 channel index (§2.1). Exposed as a
/// constant so encoder consumers can size their bitstream buffers
/// without re-deriving the rate-bit-count relationship.
pub const CHANNEL_INDEX_BITS: u32 = 10;

/// LD-CELP encoder — typed scaffold only as of round 235.
///
/// Carries the two backward adapters (synthesis-filter adapter §3.7,
/// log-gain adapter §3.8) plus the previously gain-scaled excitation
/// `ET(1..IDIM)` slot — the same per-vector state the decoder uses,
/// because §4.4 / §4.5 of the Recommendation require the encoder and
/// decoder backward adapters to be bit-for-bit identical so that the
/// quantised speech at both ends evolves in lockstep.
///
/// The [`Encoder::shape_energy`] accessor exposes the precomputed
/// `E_j` table the analysis-by-synthesis search will consume per §3.9
/// equation 3-23.
#[derive(Debug)]
pub struct Encoder {
    /// Backward synthesis-filter adapter, shared block-33 logic.
    synth_adapter: SynthesisAdapter,
    /// Backward vector gain adapter, shared block-30 logic.
    gain_adapter: GainAdapter,
    /// Perceptual weighting filter coefficients (block 38 output).
    /// Driven from the order-10 LPC predictor produced by block 37 of
    /// the perceptual weighting filter adapter (block 3) per
    /// equations 3-4b / 3-4c. Initialised to the §3.4 / §3.4.1
    /// "filter disabled" all-pass `W(z) = 1` per Table 2/G.728 —
    /// before the first adaptation cycle has run, no `q_i` is yet
    /// available, and the spec's start-up convention is to leave the
    /// filter as a pass-through.
    weighting_filter_coeff: WeightingFilterCoeff,
    /// Previous gain-scaled excitation vector — block-67 1-vector
    /// delay. Initialised to zero per Table 2/G.728.
    et_prev: [f64; FRAME_LEN],
}

impl Default for Encoder {
    fn default() -> Self {
        Self::new()
    }
}

impl Encoder {
    /// Construct a fresh encoder with all internal state at the
    /// Table 2/G.728 initial values. Identical structural shape to
    /// [`crate::Decoder::new`] because §4.4 / §4.5 require the two
    /// adapters to share their start-up state.
    pub fn new() -> Self {
        Self {
            synth_adapter: SynthesisAdapter::new(),
            gain_adapter: GainAdapter::new(),
            // §3.4 / §3.4.1: before the first adaptation cycle the
            // weighting filter is the all-pass `W(z) = 1`. The
            // disabled-mode constructor produces exactly this state.
            weighting_filter_coeff: WeightingFilterCoeff::disabled(),
            et_prev: [0.0; FRAME_LEN],
        }
    }

    /// Borrow the precomputed shape-codevector energy table
    /// `E_j = Σ_k y_j(k)²` (§3.9, eq. 3-23). The encoder's
    /// analysis-by-synthesis search reads this once per shape-codebook
    /// scan; exposing it through the encoder is the convenience hook
    /// future encoder rounds will use.
    pub fn shape_energy(&self) -> &'static [f64; NCWD] {
        &Y_ENERGY
    }

    /// Borrow the encoder's synthesis-filter adapter (§3.7). The same
    /// type is used by the decoder; §4.5 of the Recommendation
    /// mandates that the two backward adapters be bit-for-bit
    /// identical.
    pub fn synthesis_adapter(&self) -> &SynthesisAdapter {
        &self.synth_adapter
    }

    /// Borrow the encoder's log-gain adapter (§3.8). Same constraint
    /// as [`Self::synthesis_adapter`].
    pub fn gain_adapter(&self) -> &GainAdapter {
        &self.gain_adapter
    }

    /// Borrow the previous gain-scaled excitation `ET(1..IDIM)`. Block
    /// 67 (1-vector delay) reads this on every encoded vector.
    pub fn previous_excitation(&self) -> &[f64; FRAME_LEN] {
        &self.et_prev
    }

    /// Borrow the current perceptual-weighting filter coefficients
    /// `(Q̃(z/γ₁), Q̃(z/γ₂))` — output of block 38 (§3.3 eqs. 3-4b,
    /// 3-4c). Before the first adaptation cycle has run (round-248
    /// scaffold state), this is the all-pass `W(z) = 1`.
    pub fn weighting_filter_coeff(&self) -> &WeightingFilterCoeff {
        &self.weighting_filter_coeff
    }

    /// Set the perceptual-weighting filter coefficients from the
    /// order-10 LPC predictor `q_i` that block 37 (the Levinson-Durbin
    /// recursion of the weighting filter adapter, §3.3) produces. This
    /// is block 38's spec contract: the substitutions `z ← z/γ₁`,
    /// `z ← z/γ₂` of eq. 3-4b / 3-4c.
    ///
    /// `q[0]` must equal `1.0` per the spec's eq. 3-3a leading-tap
    /// convention; `q[1..=LPCW]` carries the 10 broadened predictor
    /// coefficients.
    pub fn set_weighting_filter_coeff_from_lpc(&mut self, q: &[f64; crate::consts::LPCW + 1]) {
        self.weighting_filter_coeff = WeightingFilterCoeff::from_lpc(q);
    }

    /// §3.4.1 non-speech-mode entry. Force the perceptual weighting
    /// filter to the all-pass `W(z) = 1` state regardless of the
    /// last adaptation cycle's `q_i`. Exposed so encoder callers
    /// integrating G.728 into a modem path can flip the spec's
    /// "disable weighting filter" switch without re-running the
    /// transform.
    pub fn disable_weighting_filter(&mut self) {
        self.weighting_filter_coeff = WeightingFilterCoeff::disabled();
    }

    /// Encode one `FRAME_LEN`-sample input vector into one 10-bit
    /// channel index.
    ///
    /// The full §3.9 analysis-by-synthesis search (blocks 1..28 +
    /// 67..70) is not yet implemented — the call returns
    /// [`Error::NotImplemented`] in round 235. The signature is final
    /// and matches the symmetry of [`crate::Decoder::decode_vector`]:
    /// one input vector in, one channel index out.
    pub fn encode_vector(&mut self, _input: &[f64; FRAME_LEN]) -> Result<u16, Error> {
        Err(Error::NotImplemented)
    }
}

/// Direct factory entry for the typed encoder scaffold, mirroring
/// [`crate::make_decoder`]. Returns a fresh [`Encoder`] with all
/// internal state at Table 2/G.728 initial values.
///
/// The decoder/encoder symmetry of this dual-API convention lets
/// downstream consumers wire the encoder type now; round-235's
/// `encode_*` calls still surface [`Error::NotImplemented`] because
/// the §3.9 analysis-by-synthesis pipeline is not yet implemented.
pub fn make_encoder() -> Encoder {
    Encoder::new()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::{IDIM, NCWD};

    #[test]
    fn fresh_encoder_state_matches_table_2_initials() {
        // Table 2/G.728 lists ET (1..IDIM) initial as all zero. The
        // backward adapters' own initial-state tests live in their
        // modules; here we only confirm the encoder's own field
        // (et_prev) follows the same Table-2 invariant.
        let enc = Encoder::new();
        for &v in enc.previous_excitation() {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn default_and_new_are_equivalent_initial_state() {
        // Both construction paths must produce the same et_prev block.
        let a = Encoder::default();
        let b = Encoder::new();
        assert_eq!(a.previous_excitation(), b.previous_excitation());
    }

    #[test]
    fn shape_energy_accessor_returns_table_constant() {
        // The accessor is a thin borrow over the same const-derived
        // table the encoder will consume in eq. 3-23 of §3.9. Compare
        // value-for-value: Rust does not guarantee a single
        // monomorphic backing storage for a const reference (the
        // compiler may rematerialise the const at the use site), so
        // the contract is element-wise equality, not pointer identity.
        let enc = Encoder::new();
        let table = enc.shape_energy();
        assert_eq!(table.len(), NCWD);
        for j in 0..NCWD {
            assert_eq!(table[j], Y_ENERGY[j]);
        }
    }

    #[test]
    fn encode_vector_is_not_implemented_in_r235() {
        // Round 235 lands only the typed surface. The full §3.9
        // analysis-by-synthesis search is not yet implemented; every
        // encode_vector call must surface Error::NotImplemented so
        // downstream callers can detect the pre-encoder-pipeline state.
        let mut enc = Encoder::new();
        let zero_in = [0.0f64; IDIM];
        let err = enc.encode_vector(&zero_in);
        assert_eq!(err, Err(Error::NotImplemented));
    }

    #[test]
    fn channel_index_bits_matches_spec_value() {
        // §2.1: G.728 emits a 10-bit channel index per vector. Sanity
        // guard against a future typo in the constant.
        assert_eq!(CHANNEL_INDEX_BITS, 10);
    }

    #[test]
    fn make_encoder_returns_fresh_state() {
        // The dual-API factory entry must produce the same all-zero
        // initial state as Encoder::new().
        let enc = make_encoder();
        for &v in enc.previous_excitation() {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn encoder_exposes_synthesis_adapter_at_fresh_state() {
        // The encoder must expose the same synthesis-filter adapter
        // type as the decoder (§4.4 / §4.5 bit-for-bit-identical rule).
        // The leading tap of the bandwidth-expanded predictor must be
        // 1.0 — the standing `A(1) = 1` invariant of the synthesis
        // adapter, valid at fresh-state construction.
        let enc = Encoder::new();
        let a = enc.synthesis_adapter().coefficients();
        assert_eq!(a[0], 1.0);
    }

    #[test]
    fn encoder_exposes_gain_adapter_at_fresh_state() {
        // The encoder's gain adapter must start at the Table 2/G.728
        // initial state (GSTATE = -32 dB filler, ICOUNT = 1) — same as
        // the decoder's. The gain-adapter module tests cover the field
        // values; here we only confirm the encoder's accessor surface
        // returns one.
        let enc = Encoder::new();
        assert_eq!(enc.gain_adapter().icount(), 1);
    }

    #[test]
    fn fresh_encoder_weighting_filter_is_all_pass() {
        // §3.4 / §3.4.1: before the first adaptation cycle the
        // perceptual weighting filter is the all-pass W(z) = 1, i.e.
        // both q_gamma1 and q_gamma2 are [1, 0, 0, ..., 0]. The
        // module-level tests cover the spec invariants in detail;
        // here we just confirm the encoder's fresh state matches
        // WeightingFilterCoeff::disabled().
        let enc = Encoder::new();
        let coeff = enc.weighting_filter_coeff();
        assert_eq!(*coeff, WeightingFilterCoeff::disabled());
    }

    #[test]
    fn set_weighting_filter_coeff_drives_block_38() {
        // Round-trip: feeding the encoder a known q_i vector through
        // set_weighting_filter_coeff_from_lpc must yield the same
        // (q_gamma1, q_gamma2) as a direct WeightingFilterCoeff::from_lpc
        // call. This is the spec's block 38 dataflow — the encoder
        // doesn't transform the q_i further on its way through.
        let q = [
            1.0f64,
            0.5,
            -0.25,
            0.125,
            -0.0625,
            0.03125,
            -0.015_625,
            0.007_812_5,
            -0.003_906_25,
            0.001_953_125,
            -0.000_976_562_5,
        ];

        let mut enc = Encoder::new();
        enc.set_weighting_filter_coeff_from_lpc(&q);

        let direct = WeightingFilterCoeff::from_lpc(&q);
        assert_eq!(*enc.weighting_filter_coeff(), direct);
    }

    #[test]
    fn disable_weighting_filter_returns_to_all_pass() {
        // §3.4.1 manual disable path: after the calculator has been
        // populated by a non-trivial q_i, the encoder must still be
        // able to flip back to W(z) = 1 on demand.
        let q = [
            1.0f64,
            0.5,
            -0.25,
            0.125,
            -0.0625,
            0.03125,
            -0.015_625,
            0.007_812_5,
            -0.003_906_25,
            0.001_953_125,
            -0.000_976_562_5,
        ];

        let mut enc = Encoder::new();
        enc.set_weighting_filter_coeff_from_lpc(&q);
        assert_ne!(
            *enc.weighting_filter_coeff(),
            WeightingFilterCoeff::disabled()
        );

        enc.disable_weighting_filter();
        assert_eq!(
            *enc.weighting_filter_coeff(),
            WeightingFilterCoeff::disabled()
        );
    }
}
