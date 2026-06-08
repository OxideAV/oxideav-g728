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

use crate::consts::{NCWD, NFRSZ};
use crate::decoder::FRAME_LEN;
use crate::gain_adapter::GainAdapter;
use crate::levinson::LevinsonError;
use crate::synthesis_adapter::SynthesisAdapter;
use crate::tables::Y_ENERGY;
use crate::weighting_filter::PerceptualWeightingFilter;
use crate::weighting_filter_adapter::WeightingFilterAdapter;
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
    /// Block 4 — the perceptual weighting filter `W(z)` applied to the
    /// current input speech vector `s(n)` to produce the weighted
    /// speech vector `v(n)` (§3.4). Per §3.4 the per-sample memory
    /// must not be reset to zero outside initialisation; the filter
    /// is constructed once here at the all-pass state and thereafter
    /// only receives new coefficients via
    /// [`set_weighting_filter_coeff_from_lpc`] (which mirrors them
    /// into [`weighting_filter`]) — never a state reset.
    weighting_filter: PerceptualWeightingFilter,
    /// Upstream weighting-filter adapter (blocks 36 + 37 of §3.3).
    /// Owns the hybrid-window state (`SBW`, `REXPW` per Table 2) and
    /// the most recently produced order-10 predictor `q_i`. Driven
    /// per adaptation cycle by [`Encoder::adapt_weighting_filter`];
    /// at the third vector of each cycle the encoder reads the
    /// predictor back out and pushes it through block 38 via
    /// [`Encoder::commit_weighting_filter_coefficients`].
    weighting_filter_adapter: WeightingFilterAdapter,
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
            // §3.4 spec note: filter memory is zeroed only at
            // initialisation. PerceptualWeightingFilter::new()
            // produces that state — coefficients = disabled() (so
            // both q_gamma1 and q_gamma2 collapse to [1, 0, …, 0]),
            // both delay lines = [0; LPCW]. The encoder then never
            // calls a reset on this field for the lifetime of the
            // Encoder struct.
            weighting_filter: PerceptualWeightingFilter::new(),
            // §3.3 cold start: the upstream block-36/37 adapter
            // owns its own Table-2/G.728 initial state (SBW = 0,
            // REXPW = 0, predictor = all-pass). No coefficient swap
            // happens until the first cycle completes and the
            // encoder explicitly calls
            // commit_weighting_filter_coefficients.
            weighting_filter_adapter: WeightingFilterAdapter::new(),
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
    ///
    /// The new coefficient set is also mirrored into the live
    /// [`PerceptualWeightingFilter`] (block 4) — per §3.3 / §3.4 the
    /// freeze-and-swap convention only swaps the taps, never the
    /// per-sample memory.
    pub fn set_weighting_filter_coeff_from_lpc(&mut self, q: &[f64; crate::consts::LPCW + 1]) {
        self.weighting_filter_coeff = WeightingFilterCoeff::from_lpc(q);
        self.weighting_filter
            .set_coefficients(self.weighting_filter_coeff);
    }

    /// §3.4.1 non-speech-mode entry. Force the perceptual weighting
    /// filter to the all-pass `W(z) = 1` state regardless of the
    /// last adaptation cycle's `q_i`. Exposed so encoder callers
    /// integrating G.728 into a modem path can flip the spec's
    /// "disable weighting filter" switch without re-running the
    /// transform.
    ///
    /// Per §3.4 the filter memory is not touched — only the
    /// coefficient set is swapped to the all-pass state.
    pub fn disable_weighting_filter(&mut self) {
        self.weighting_filter_coeff = WeightingFilterCoeff::disabled();
        self.weighting_filter
            .set_coefficients(self.weighting_filter_coeff);
    }

    /// Borrow the live block-4 perceptual weighting filter state
    /// (current coefficients + per-sample delay-line memory).
    pub fn weighting_filter(&self) -> &PerceptualWeightingFilter {
        &self.weighting_filter
    }

    /// Borrow the upstream weighting-filter adapter (blocks 36 + 37
    /// of §3.3). The adapter is the encoder's source of truth for the
    /// order-10 predictor `q_i` that drives block 38. Read-only so
    /// callers see what the next coefficient swap will commit, but
    /// cannot bypass [`Self::commit_weighting_filter_coefficients`]
    /// (which is what propagates the new `q_i` to the live block-4
    /// filter).
    pub fn weighting_filter_adapter(&self) -> &WeightingFilterAdapter {
        &self.weighting_filter_adapter
    }

    /// Run blocks 36 + 37 of the §3.3 weighting-filter adapter on
    /// one adaptation cycle of input speech (`NFRSZ = 20` samples =
    /// `NUPDATE = 4` vectors of `IDIM = 5` samples each). The new
    /// order-10 predictor lands in the adapter's cache; the encoder
    /// does **not** yet propagate it into the live block-38 /
    /// block-4 path — that step is gated on the §3.3 "third vector
    /// of each cycle" timing rule and is exposed separately as
    /// [`Self::commit_weighting_filter_coefficients`].
    ///
    /// On Levinson-Durbin failure the cached predictor is left
    /// untouched and the error is propagated so the caller can log
    /// or trace it; the block-33 mirror policy is documented on
    /// [`WeightingFilterAdapter::adapt`].
    pub fn adapt_weighting_filter(&mut self, speech: &[f64; NFRSZ]) -> Result<(), LevinsonError> {
        self.weighting_filter_adapter.adapt(speech)?;
        Ok(())
    }

    /// Push the upstream adapter's cached predictor through block 38
    /// (the §3.3 eq. 3-4b / 3-4c substitutions) and into the live
    /// block-4 perceptual weighting filter. Per §3.3 the encoder
    /// only does this at the **third vector of each adaptation
    /// cycle**; the timing gate is the caller's responsibility, the
    /// same way the synthesis-filter adapter's commit timing is the
    /// caller's responsibility for block 33.
    ///
    /// Per §3.4 spec rule the per-sample memory of block 4 is
    /// preserved across the swap — only the coefficients change.
    pub fn commit_weighting_filter_coefficients(&mut self) {
        let q = *self.weighting_filter_adapter.predictor();
        self.set_weighting_filter_coeff_from_lpc(&q);
    }

    /// Apply the block-4 perceptual weighting filter to one
    /// `IDIM`-sample input speech vector `s(n)` and return the
    /// corresponding weighted vector `v(n)` (§3.4). The filter's
    /// per-sample memory advances as a side effect.
    ///
    /// The encoder's analysis-by-synthesis search of §3.9 consumes
    /// `v(n)` (eventually combined with the zero-input response from
    /// block 10 to form the VQ target `x(n)`, §3.5 / §3.6). Block 10
    /// is intentionally still NOT wired up — its cascade with the
    /// synthesis filter requires the §3.10 pre-/post-save memory
    /// handling and lands in a future round.
    pub fn apply_weighting_filter(&mut self, s: &[f64; FRAME_LEN]) -> [f64; FRAME_LEN] {
        self.weighting_filter.filter_vector(s)
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
    fn fresh_encoder_weighting_filter_passes_input_through() {
        // §3.4 / §3.4.1 cold start: at the all-pass W(z) = 1 state
        // with zeroed memory, block 4 emits v(n) = s(n) bit-for-bit.
        let mut enc = Encoder::new();
        let s = [42.0, -17.0, 0.0, 100.5, -250.25];
        let v = enc.apply_weighting_filter(&s);
        for k in 0..IDIM {
            assert_eq!(v[k], s[k], "k={k}");
        }
    }

    #[test]
    fn set_weighting_filter_coeff_propagates_to_block_4_taps() {
        // The setter must mirror its newly-computed coefficients into
        // the live block-4 filter so that the next apply_weighting_filter
        // sees them. Without this propagation the encoder would carry
        // a "pending" coefficient set that block 4 ignores — a silent
        // spec violation.
        let q = [
            1.0_f64,
            -0.5,
            0.25,
            -0.125,
            0.0625,
            -0.03125,
            0.015_625,
            -0.007_812_5,
            0.003_906_25,
            -0.001_953_125,
            0.000_976_562_5,
        ];
        let mut enc = Encoder::new();
        enc.set_weighting_filter_coeff_from_lpc(&q);

        let live = enc.weighting_filter().coefficients();
        let expected = WeightingFilterCoeff::from_lpc(&q);
        assert_eq!(*live, expected);
    }

    #[test]
    fn weighting_filter_memory_survives_coefficient_swap() {
        // §3.3 / §3.4 invariant: when block 38 commits a new
        // coefficient set, block 4's per-sample memory must carry
        // continuously across the swap.
        let mut enc = Encoder::new();
        let s = [10.0, 20.0, 30.0, 40.0, 50.0];
        let _ = enc.apply_weighting_filter(&s);
        let mem_before = *enc.weighting_filter().input_memory();

        let mut q = [0.0f64; crate::consts::LPCW + 1];
        q[0] = 1.0;
        q[1] = -0.4;
        enc.set_weighting_filter_coeff_from_lpc(&q);

        assert_eq!(*enc.weighting_filter().input_memory(), mem_before);
    }

    #[test]
    fn disable_weighting_filter_keeps_per_sample_memory() {
        // §3.4 spec rule: even when the encoder flips the §3.4.1
        // "non-speech mode" disable switch, the per-sample memory of
        // the filter is left alone. The coefficients return to the
        // all-pass state but block 4's delay lines carry whatever
        // they carried before the flip.
        let mut enc = Encoder::new();
        let s = [10.0, 20.0, 30.0, 40.0, 50.0];
        let _ = enc.apply_weighting_filter(&s);
        let mem_before = *enc.weighting_filter().input_memory();

        enc.disable_weighting_filter();

        assert_eq!(*enc.weighting_filter().input_memory(), mem_before);
        assert_eq!(
            *enc.weighting_filter().coefficients(),
            WeightingFilterCoeff::disabled()
        );
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

    #[test]
    fn fresh_encoder_weighting_filter_adapter_is_allpass() {
        // §3.3 cold start: the upstream block-36/37 adapter starts
        // with the all-pass q_i = (1, 0, ..., 0) so block 38 on the
        // same vector collapses to W(z) = 1. The encoder must expose
        // exactly that state at construction time so the live block-4
        // filter and the cached upstream predictor agree on "no
        // weighting yet".
        let enc = Encoder::new();
        let q = enc.weighting_filter_adapter().predictor();
        assert_eq!(q[0], 1.0);
        assert!(q[1..].iter().all(|&v| v == 0.0));
    }

    #[test]
    fn adapt_weighting_filter_returns_levinson_error_on_zero_input() {
        // The encoder must surface a Levinson failure verbatim — its
        // job is to relay the upstream adapter's contract, not to
        // mask failures behind a generic Error::NotImplemented or a
        // silent success. Zero input forces R(1) ≤ 0; we expect
        // either ZeroSignal or TrailingZero per the adapter docs.
        let mut enc = Encoder::new();
        let zero = [0.0f64; NFRSZ];
        let result = enc.adapt_weighting_filter(&zero);
        assert!(
            matches!(
                result,
                Err(LevinsonError::ZeroSignal | LevinsonError::TrailingZero)
            ),
            "expected Levinson failure on zero input, got {:?}",
            result
        );
    }

    #[test]
    fn adapt_weighting_filter_on_nontrivial_input_succeeds() {
        // A speech-like input must drive the adapter cleanly through
        // hybrid window + Levinson. Run 6 cycles (the same budget
        // the synthesis adapter's tests use) and assert at least one
        // adapt() succeeds — exactly one would be enough but the
        // multi-cycle smoke gives the hybrid window's recursive part
        // a chance to populate.
        let mut enc = Encoder::new();
        let mut input = [0.0f64; NFRSZ];
        for k in 0..NFRSZ {
            let t = k as f64;
            input[k] = 100.0 * (0.4 * t).sin() + 50.0 * (0.9 * t).cos() + (t * 13.0).fract();
        }
        let mut ok_count = 0;
        for _ in 0..6 {
            if enc.adapt_weighting_filter(&input).is_ok() {
                ok_count += 1;
            }
        }
        assert!(
            ok_count >= 1,
            "expected at least one successful weighting-filter adaptation"
        );
    }

    #[test]
    fn commit_weighting_filter_coefficients_propagates_to_block_38_and_block_4() {
        // End-to-end §3.3 wiring: adapt on a speech-like cycle,
        // commit the resulting q_i, then assert (a) block 38's
        // (q_gamma1, q_gamma2) reflect the new predictor (not the
        // disabled all-pass), and (b) the live block-4 filter sees
        // the same coefficients. This is the spec's full block 36 →
        // block 37 → block 38 → block 4 chain.
        let mut enc = Encoder::new();
        let mut input = [0.0f64; NFRSZ];
        for k in 0..NFRSZ {
            let t = k as f64;
            input[k] = 200.0 * (0.3 * t).sin() + 50.0 * (1.1 * t).cos() + 7.0;
        }

        // Drive several cycles to land Levinson on a non-degenerate
        // RTMP, then commit.
        for _ in 0..6 {
            let _ = enc.adapt_weighting_filter(&input);
        }
        enc.commit_weighting_filter_coefficients();

        let committed = *enc.weighting_filter_coeff();
        // After a successful adaptation + commit the encoder must
        // not be at the disabled all-pass shape.
        assert_ne!(committed, WeightingFilterCoeff::disabled());
        // And block 4 must have picked up the same coefficients —
        // the setter's invariant, exercised here through the
        // commit path.
        assert_eq!(*enc.weighting_filter().coefficients(), committed);
    }

    #[test]
    fn commit_weighting_filter_preserves_block_4_memory() {
        // §3.4 spec rule: coefficient swap MUST NOT touch block 4's
        // per-sample delay-line memory. The encoder's commit path
        // goes through set_weighting_filter_coeff_from_lpc, which
        // already honours that rule for the manual setter path; this
        // test pins the same rule for the upstream-driven commit
        // path so a future refactor cannot regress it.
        let mut enc = Encoder::new();
        // Push some samples to populate block 4's memory.
        let s = [10.0, 20.0, 30.0, 40.0, 50.0];
        let _ = enc.apply_weighting_filter(&s);
        let mem_before = *enc.weighting_filter().input_memory();

        // Adapt and commit a non-trivial coefficient set.
        let mut input = [0.0f64; NFRSZ];
        for k in 0..NFRSZ {
            input[k] = 100.0 * (0.5 * k as f64).sin() + 30.0;
        }
        for _ in 0..6 {
            let _ = enc.adapt_weighting_filter(&input);
        }
        enc.commit_weighting_filter_coefficients();

        assert_eq!(*enc.weighting_filter().input_memory(), mem_before);
    }

    #[test]
    fn weighting_filter_adapter_accessor_is_read_only_view() {
        // The accessor must return the encoder's own adapter (same
        // structural state as it was set up at construction). At
        // fresh-state the predictor is all-pass and exactly equal
        // to a standalone WeightingFilterAdapter's initial state.
        let enc = Encoder::new();
        let q_enc = enc.weighting_filter_adapter().predictor();
        let standalone = WeightingFilterAdapter::new();
        let q_solo = standalone.predictor();
        assert_eq!(q_enc, q_solo);
    }
}
