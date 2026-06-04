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
}
