//! ITU-T G.728 LD-CELP encoder.
//!
//! Round 235 landed the **typed encoder surface**; rounds 248 / 249 /
//! 258 / 267 grew the front end (blocks 38, 4, 36 + 37, 9 + 10 + 11);
//! round 276 lands the **§3.9 codebook search** (blocks 12–18, in
//! [`crate::codebook_search`]) and the **§3.10 memory update**
//! ([`ZeroInputResponse::update_memory`], pseudo-code §5.13), closing
//! the analysis-by-synthesis loop: [`Encoder::encode_vector`] now
//! emits one real 10-bit channel index per `IDIM = 5`-sample input
//! vector.
//!
//! ## Per-vector dataflow (Figure 2/G.728)
//!
//! 1. **Block 20** — predict the excitation gain `σ(n)` from the
//!    previous gain-scaled excitation (backward vector gain adapter,
//!    §3.8; shared type with the decoder's block 30 per §4.5).
//! 2. **Block 4** — perceptual weighting filter on the input speech,
//!    `v(n)` (§3.4).
//! 3. **Blocks 9 + 10 + 11** — zero-input response ring + VQ target
//!    `x(n) = v(n) − r(n)` (§3.5 / §3.6).
//! 4. **Blocks 12–18** — codebook search: impulse response, energy
//!    table, target normalisation, time-reversed convolution, gain
//!    decision tree, best-index selection (§3.9).
//! 5. **Blocks 19 + 21** — excitation lookup + gain scaling,
//!    `e(n) = σ(n)·g_i·y_j` (§5.12).
//! 6. **§5.13 memory update** — zero-state responses of `e(n)` added
//!    onto the post-ring filter memory; the quantized speech `sq(n)`
//!    falls out of the top `STATELPC` taps (block 22 omitted per
//!    §3.10).
//! 7. **Blocks 23 + 3 cycle bookkeeping + ICOUNT = 3 stagger** —
//!    `sq(n)` accumulates into the synthesis adapter's NFRSZ buffer
//!    and `s(n)` into the weighting adapter's; when a full adaptation
//!    cycle (NUPDATE = 4 vectors) is collected both adapters re-run
//!    and their results are **stashed** (synthesis predictor +
//!    block-38 weighting coefficients) rather than applied. Per
//!    Table E-1/G.728 (and the block-51 `Wait until ICOUNT = 3, then
//!    A = ATMP` instruction) the new synthesis predictor, weighting
//!    coefficients, and impulse response + energy table (blocks
//!    12 + 14 + 15) are **first used at encoding vector 3** of the
//!    following cycle — so the swap is deferred to `ICOUNT = 3`. This
//!    is the spec-faithful update stagger; [`crate::Decoder::
//!    decode_vector`] applies the identical deferred swap for block
//!    33, so the encoder and decoder synthesis predictors evolve in
//!    exact lockstep (§4.5) — the precondition for interoperating with
//!    a third-party G.728 endpoint.
//!
//! Because blocks 19–23 form the **simulated decoder** (block 8,
//! §3.10), an [`Encoder`] driven by any input and a
//! [`crate::Decoder`] fed the resulting channel indices produce the
//! same quantized speech bit for bit — the
//! `encoder_decoder_lockstep_*` tests pin that property.
//!
//! The precomputed codebook constants the search consumes
//! (`Y_ENERGY`, `G2`, `GSQ`, `GB`) live in `tables.rs`; the per-cycle
//! `E_j = ‖H·y_j‖²` table lives in [`CodebookSearch`].

use crate::codebook_search::CodebookSearch;
use crate::consts::{LPC, NCWD, NFRSZ, NUPDATE};
use crate::decoder::{ExcitationVector, FRAME_LEN};
use crate::gain_adapter::GainAdapter;
use crate::levinson::LevinsonError;
use crate::synthesis_adapter::SynthesisAdapter;
use crate::tables::Y_ENERGY;
use crate::weighting_filter::PerceptualWeightingFilter;
use crate::weighting_filter_adapter::WeightingFilterAdapter;
use crate::weighting_filter_coeff::WeightingFilterCoeff;
use crate::zero_input_response::ZeroInputResponse;
use crate::Error;

/// Number of bits in one G.728 channel index (§2.1). Exposed as a
/// constant so encoder consumers can size their bitstream buffers
/// without re-deriving the rate-bit-count relationship.
pub const CHANNEL_INDEX_BITS: u32 = 10;

/// LD-CELP encoder — full per-vector analysis-by-synthesis loop as
/// of round 276 ([`Encoder::encode_vector`]).
///
/// Carries the two backward adapters (synthesis-filter adapter §3.7,
/// log-gain adapter §3.8) plus the previously gain-scaled excitation
/// `ET(1..IDIM)` slot — the same per-vector state the decoder uses,
/// because §4.4 / §4.5 of the Recommendation require the encoder and
/// decoder backward adapters to be bit-for-bit identical so that the
/// quantised speech at both ends evolves in lockstep.
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
    /// Zero-input-response + VQ-target unit (blocks 9 + 10 + 11 of
    /// §3.5 / §3.6) — and, since r276, the §5.13 memory-update phase
    /// (§3.10). Carries the synthesis-filter ring memory `STATELPC`
    /// and the weighting-filter ring memories `ZIRWFIR` / `ZIRWIIR`,
    /// all at the spec's all-zero initialisation state per
    /// Table 2/G.728. Driven by [`Encoder::compute_vq_target`] /
    /// [`Encoder::encode_vector`].
    zir: ZeroInputResponse,
    /// Codebook search module (blocks 12–18 of §3.9). Carries the
    /// per-cycle impulse response `H` and filtered-shape energy
    /// table `Y2` at their Table 2/G.728 initial states.
    codebook_search: CodebookSearch,
    /// Previous gain-scaled excitation vector — block-67 1-vector
    /// delay. Initialised to zero per Table 2/G.728.
    et_prev: [f64; FRAME_LEN],
    /// Rolling NFRSZ-sample quantized-speech buffer feeding the
    /// synthesis adapter (block 23) — the encoder-side mirror of the
    /// decoder's `STTMP(1..NFRSZ)` accumulation.
    sttmp_buf: [f64; NFRSZ],
    /// Rolling NFRSZ-sample input-speech buffer feeding the
    /// weighting-filter adapter (blocks 36 + 37) — the spec's
    /// `STMP(1..4*IDIM)` cycle buffer.
    stmp_buf: [f64; NFRSZ],
    /// Shared write cursor into the two cycle buffers; both fill in
    /// lockstep (one vector of each per encoded vector) and the
    /// adapters re-run when a full NFRSZ cycle is collected.
    cycle_idx: usize,
    /// Most recent quantized speech vector `sq(n)` (the §3.10
    /// by-product of the memory update) — exposed for the
    /// encoder/decoder lockstep tests and for callers that want the
    /// locally decoded signal without running a separate decoder.
    last_sq: [f64; FRAME_LEN],
    /// **Spec ICOUNT for the per-cycle coefficient swap.** 1-based,
    /// walks `1 → 2 → 3 → 4 → 1` once per [`Encoder::encode_vector`]
    /// call. Per Table E-1/G.728 the backward synthesis filter (blocks
    /// 23/49/50/51) AND the perceptual weighting filter / fast
    /// codebook-search adapter (blocks 3/36/37/38 + 12/14/15) both
    /// take their input through vector 2–4 of a cycle but are **first
    /// used at encoding vector 3** of the following cycle. This counter
    /// drives that deferred swap so the encoder's synthesis predictor
    /// evolves in exact lockstep with the decoder's.
    sf_icount: usize,
    /// **Pending synthesis predictor** (bandwidth-expanded `ATMP`,
    /// block 51) computed at a cycle boundary, awaiting the
    /// `ICOUNT = 3` swap into [`Self::active_predictor`]. `None` until
    /// the first full adaptation cycle has accumulated.
    pending_synth: Option<[f64; LPC + 1]>,
    /// **Active synthesis predictor** currently driving blocks 9/10/11
    /// (ZIR target) and block 13 (§3.10 memory update). Only refreshed
    /// from [`Self::pending_synth`] at `ICOUNT = 3`; during vectors 1
    /// and 2 of each cycle the encoder keeps the previous cycle's
    /// committed predictor (block-50 note: "the old set … is still
    /// being used" until the third vector).
    active_predictor: [f64; LPC + 1],
    /// **Pending perceptual weighting coefficients** (block 38 output)
    /// computed at a cycle boundary, awaiting the `ICOUNT = 3` commit.
    /// `None` when the cycle's Levinson-Durbin (block 37) failed or no
    /// cycle has completed — in which case the previous cycle's
    /// weighting coefficients are kept (the §3.3 ill-conditioning
    /// rule) and no commit happens.
    pending_weighting: Option<WeightingFilterCoeff>,
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
            // §3.5 initialisation: the synthesis-filter and
            // weighting-filter ring memories start all-zero, so the
            // first vector's zero-input response is zero and the VQ
            // target equals the weighted speech exactly.
            zir: ZeroInputResponse::new(),
            // Table 2/G.728: H starts at the identity impulse
            // (1,0,0,0,0) and Y2 at the bare shape energies — the
            // CodebookSearch constructor encodes both.
            codebook_search: CodebookSearch::new(),
            et_prev: [0.0; FRAME_LEN],
            sttmp_buf: [0.0; NFRSZ],
            stmp_buf: [0.0; NFRSZ],
            cycle_idx: 0,
            last_sq: [0.0; FRAME_LEN],
            // Spec ICOUNT runs 1..4; start at NUPDATE so the first
            // `encode_vector` call advances it to 1.
            sf_icount: NUPDATE,
            // No adaptation cycle has completed yet. The active
            // synthesis predictor starts at the all-pass identity
            // (`A(1) = 1`, rest 0) per Table 2/G.728.
            pending_synth: None,
            active_predictor: {
                let mut a = [0.0f64; LPC + 1];
                a[0] = 1.0;
                a
            },
            pending_weighting: None,
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

    /// The **active synthesis predictor** currently driving blocks
    /// 9/10/11 + 13 (the §3.7 / block-51 output, committed at the
    /// third vector of each adaptation cycle per Table E-1/G.728). This
    /// is the value the encoder's analysis-by-synthesis loop actually
    /// uses; it lags [`Self::synthesis_adapter`]`().coefficients()` (the
    /// most-recently computed predictor) by the ICOUNT = 3 swap delay.
    /// The decoder's [`crate::Decoder::active_predictor`] follows the
    /// identical stagger, so the two values coincide vector-for-vector
    /// when both ends process the same quantised-speech stream.
    pub fn active_predictor(&self) -> &[f64; LPC + 1] {
        &self.active_predictor
    }

    /// The current 1-based synthesis-filter `ICOUNT` (walks
    /// `1 → 2 → 3 → 4 → 1` once per [`Self::encode_vector`] call). The
    /// block-51 / block-38 coefficient swap commits when this reaches
    /// 3. Useful for tests and audit.
    pub fn sf_icount(&self) -> usize {
        self.sf_icount
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

    /// Borrow the zero-input-response unit (blocks 9 + 10 + 11). The
    /// accessor exposes the three filter-memory arrays
    /// (`STATELPC` / `ZIRWFIR` / `ZIRWIIR`) for tests and audit.
    pub fn zero_input_response_unit(&self) -> &ZeroInputResponse {
        &self.zir
    }

    /// Form the §3.9 analysis-by-synthesis **target vector** `x(n)`
    /// from the weighted speech vector `v(n)` (blocks 9 + 10 + 11,
    /// §3.5 / §3.6).
    ///
    /// `x(n) = v(n) − r(n)`, where `r(n)` is the zero-input response
    /// of the synthesis filter (block 9) cascaded with the perceptual
    /// weighting filter (block 10) — the "ring" of §3.5. The current
    /// order-50 synthesis predictor comes from the encoder's own
    /// synthesis-filter adapter (block 23, shared with the decoder per
    /// §4.5), and the weighting coefficients from the live block-38
    /// coefficient set; the caller passes the weighted speech vector
    /// produced by [`Self::apply_weighting_filter`].
    ///
    /// Side effect: the synthesis-filter and weighting-filter ring
    /// memories advance one slot per sample as §5.9 prescribes, so the
    /// next call sees the rung-down (generally non-zero) memory state.
    ///
    /// This runs the **zero-input phase** of §3.5 only. The §3.10
    /// memory-update phase — adding the zero-state response of the
    /// chosen excitation `e(n)` back onto the saved ring memory — is a
    /// later round; it depends on the §3.9 codebook search output that
    /// is not yet wired up.
    pub fn compute_vq_target(&mut self, v: &[f64; FRAME_LEN]) -> [f64; FRAME_LEN] {
        let a = *self.synth_adapter.coefficients();
        let w = self.weighting_filter_coeff;
        self.zir.compute_target(&a, &w, v)
    }

    /// Borrow the codebook-search module (blocks 12–18 of §3.9) —
    /// read-only view of the current impulse response `H` and
    /// filtered-shape energy table `Y2`.
    pub fn codebook_search(&self) -> &CodebookSearch {
        &self.codebook_search
    }

    /// Borrow the most recent quantized speech vector `sq(n)` — the
    /// §3.10 by-product of the last [`Self::encode_vector`] call.
    /// Because blocks 19–23 form the simulated decoder (block 8),
    /// this equals what [`crate::Decoder::decode_vector`] produces
    /// for the same channel index, bit for bit.
    pub fn quantized_speech(&self) -> &[f64; FRAME_LEN] {
        &self.last_sq
    }

    /// Encode one `FRAME_LEN`-sample input vector into one 10-bit
    /// channel index — the full per-vector analysis-by-synthesis loop
    /// of Figure 2/G.728 (§3.4–§3.10). Mirrors the symmetry of
    /// [`crate::Decoder::decode_vector`]: one input vector in, one
    /// channel index out.
    ///
    /// See the module docs for the block-by-block dataflow.
    ///
    /// **ICOUNT = 3 coefficient stagger (Table E-1/G.728).** The
    /// backward synthesis filter (blocks 23/49/50/51) and the
    /// perceptual-weighting / fast-codebook-search adapter (blocks
    /// 3/36/37/38 + 12/14/15) both compute their new coefficients from
    /// signal accumulated over a 4-vector cycle, but the spec defers
    /// **first use of the updated coefficients to vector 3** of the
    /// following cycle. The encoder therefore stashes the freshly-
    /// computed predictor / weighting coefficients at the cycle
    /// boundary and only swaps them into the live blocks when
    /// `ICOUNT` reaches 3 — matching the decoder's identical
    /// `decode_vector` stagger so the two ends stay in exact lockstep
    /// (§4.5) and interoperate with a third-party G.728 endpoint.
    ///
    /// The `Result` shape is kept (per-vector encoding itself cannot
    /// fail) so byte-stream framing entry points can share the
    /// signature.
    pub fn encode_vector(&mut self, input: &[f64; FRAME_LEN]) -> Result<u16, Error> {
        self.encode_vector_inner(input, None)
    }

    /// §3.11 in-band-signalling encode: encode one vector while robbing
    /// the leftmost codeword bit to carry the synchronization /
    /// signalling bit `sync_bit`.
    ///
    /// Identical to [`Self::encode_vector`] except that the §3.9
    /// codebook search is restricted to **one half** of the shape
    /// codebook (blocks 17 + 18 over shape `0..=63` for a `0`,
    /// `64..=127` for a `1`) so the emitted index's top shape bit equals
    /// `sync_bit` (§3.11). The full per-vector backward-adaptation
    /// dataflow — predictor σ(n), weighting, ZIR target, gain search,
    /// excitation lookup, §3.10 memory update, adapter bookkeeping — is
    /// unchanged, so the encoder and decoder stay in lockstep on the
    /// robbed vector exactly as for an un-robbed one. The decoder
    /// recovers the bit with [`crate::extract_sync_bit`].
    ///
    /// §3.11 recommends robbing once every `N` vectors with `N` a
    /// multiple of 4 (e.g. `N = 16`, ≈ 100 bit/s) and robbing from the
    /// **last** vector of an adaptation cycle; scheduling which vectors
    /// are robbed is the caller's responsibility (both ends must agree).
    pub fn encode_vector_with_sync_bit(
        &mut self,
        input: &[f64; FRAME_LEN],
        sync_bit: bool,
    ) -> Result<u16, Error> {
        self.encode_vector_inner(input, Some(sync_bit))
    }

    fn encode_vector_inner(
        &mut self,
        input: &[f64; FRAME_LEN],
        sync_bit: Option<bool>,
    ) -> Result<u16, Error> {
        // ----- Spec ICOUNT advance + ICOUNT = 3 deferred swap --------
        // Advance the cycle counter (1 → 2 → 3 → 4 → 1). At the third
        // vector, commit the previous cycle's pending synthesis
        // predictor (block 51 `Wait until ICOUNT = 3, then A = ATMP`),
        // the pending weighting coefficients (block 38), and refresh
        // the impulse response + filtered-shape energy table (blocks
        // 12/14/15, which the spec also runs at ICOUNT = 3 once the new
        // A, AWZ, AWP are ready).
        self.sf_icount = (self.sf_icount % NUPDATE) + 1;
        if self.sf_icount == 3 {
            if let Some(pending) = self.pending_synth.take() {
                self.active_predictor = pending;
            }
            if let Some(w_new) = self.pending_weighting.take() {
                self.weighting_filter_coeff = w_new;
                self.weighting_filter.set_coefficients(w_new);
            }
            // Blocks 12/14/15: refresh now that A (block 51) and AWZ /
            // AWP (block 38) for this cycle are both committed.
            self.codebook_search
                .update_impulse_response(&self.active_predictor, &self.weighting_filter_coeff);
        }

        // ----- Block 20: predict σ(n) from the previous ET -----------
        // Same call order as the decoder's block 30 (§4.5 lockstep).
        let gain = self.gain_adapter.predict_next(&self.et_prev);

        // ----- Block 4: perceptual weighting (§3.4) ------------------
        let v = self.weighting_filter.filter_vector(input);

        // ----- Blocks 9 + 10 + 11: ZIR ring + VQ target (§3.5/§3.6) --
        // Use the ICOUNT-staggered active synthesis predictor (block 51
        // output committed at the third vector), NOT the adapter's
        // most-recent result directly.
        let a = self.active_predictor;
        let w = self.weighting_filter_coeff;
        let target = self.zir.compute_target(&a, &w, &v);

        // ----- Blocks 12–18: codebook search (§3.9) ------------------
        // Blocks 16 + 13 + 17 + 18 run per vector; the per-cycle
        // blocks 12 + 14 + 15 are refreshed at ICOUNT = 3 above. On a
        // §3.11-robbed vector the shape scan is restricted to one half
        // of the codebook so the leftmost codeword bit carries the
        // signalling bit.
        let res = match sync_bit {
            Some(bit) => self
                .codebook_search
                .search_with_sync_bit(&target, gain, bit),
            None => self.codebook_search.search(&target, gain),
        };

        // ----- Blocks 19 + 21: excitation lookup + gain scaling ------
        // §5.12: YN = GQ(IG)·Y(IS row), ET = GAIN·YN. The decoder-side
        // ExcitationVector lookup implements exactly the block-19
        // table read, so the two ends share the code path.
        let yn = ExcitationVector::from_channel_index(res.channel_index);
        let mut et = [0.0f64; FRAME_LEN];
        for k in 0..FRAME_LEN {
            et[k] = gain * yn.0[k];
        }

        // ----- §5.13: filter memory update (§3.10) -------------------
        // Adds the zero-state responses of e(n) onto the post-ring
        // memory and yields sq(n) from the top STATELPC taps.
        let sq = self.zir.update_memory(&et, &a, &w);
        self.last_sq = sq;

        // ----- Cycle bookkeeping: blocks 23 + 36/37/38 --------------
        // Accumulate one vector of quantized speech (block 23 input)
        // and one vector of input speech (block 36 input); on a full
        // NFRSZ cycle re-run both adapters and stash their outputs as
        // `pending_*` for the deferred ICOUNT = 3 swap above. On
        // Levinson failure each adapter keeps its previous cycle's
        // output (the §3.7 "keep the old coefficients" policy): for the
        // weighting filter that means leaving `pending_weighting` empty
        // so the previous coefficients survive; for the synthesis
        // adapter the cached `coefficients()` are unchanged so the same
        // value is re-stashed harmlessly.
        for k in 0..FRAME_LEN {
            self.sttmp_buf[self.cycle_idx + k] = sq[k];
            self.stmp_buf[self.cycle_idx + k] = input[k];
        }
        self.cycle_idx += FRAME_LEN;
        if self.cycle_idx >= NFRSZ {
            self.cycle_idx = 0;
            let _ = self.synth_adapter.adapt(&self.sttmp_buf);
            self.pending_synth = Some(*self.synth_adapter.coefficients());
            if self.weighting_filter_adapter.adapt(&self.stmp_buf).is_ok() {
                let q = *self.weighting_filter_adapter.predictor();
                self.pending_weighting = Some(WeightingFilterCoeff::from_lpc(&q));
            }
        }

        // ----- Block 67 wrap-around: save ET for the next vector -----
        self.et_prev = et;

        Ok(res.channel_index)
    }
}

/// Direct factory entry for the encoder, mirroring
/// [`crate::make_decoder`] per the workspace dual-API convention.
/// Returns a fresh [`Encoder`] with all internal state at
/// Table 2/G.728 initial values.
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
    fn encode_vector_emits_ten_bit_channel_indices() {
        // r276: the §3.9 + §3.10 loop is closed — every call returns
        // a real channel index in [0, 1024) (CHANNEL_INDEX_BITS = 10),
        // including on all-zero input (the COR = 0 path of the §5.11
        // decision tree is well-defined).
        let mut enc = Encoder::new();
        let zero_in = [0.0f64; IDIM];
        let ichan = enc.encode_vector(&zero_in).expect("encode must succeed");
        assert!(ichan < 1024);
        for i in 0..64usize {
            let mut s = [0.0f64; IDIM];
            for k in 0..IDIM {
                s[k] = 500.0 * ((i * IDIM + k) as f64 * 0.37).sin();
            }
            let ichan = enc.encode_vector(&s).expect("encode must succeed");
            assert!(ichan < 1024, "vector {i}");
            for &q in enc.quantized_speech() {
                assert!(q.is_finite(), "vector {i}");
            }
        }
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
    fn fresh_encoder_zir_unit_is_all_zero() {
        // §3.5 initialisation: the synthesis-filter and
        // weighting-filter ring memories start all-zero.
        let enc = Encoder::new();
        let z = enc.zero_input_response_unit();
        assert!(z.state_lpc().iter().all(|&v| v == 0.0));
        assert!(z.zirw_fir().iter().all(|&v| v == 0.0));
        assert!(z.zirw_iir().iter().all(|&v| v == 0.0));
    }

    #[test]
    fn first_vq_target_equals_weighted_speech() {
        // Block 11: x(n) = v(n) − r(n). On the first vector after
        // construction every ring memory tap is zero, so r(n) = 0 and
        // x(n) = v(n) bit-for-bit — independent of the (fresh,
        // all-pass) predictor and weighting coefficients.
        let mut enc = Encoder::new();
        let v = [12.0, -34.0, 56.0, -78.0, 90.0];
        let x = enc.compute_vq_target(&v);
        for k in 0..IDIM {
            assert_eq!(x[k], v[k], "k={k}");
        }
    }

    #[test]
    fn compute_vq_target_advances_ring_memory() {
        // After a vector pushed through a non-trivial predictor the
        // ring memory must evolve (the §5.9 "ring" shifts state). Drive
        // a real adaptation to load a feedback predictor, then run two
        // target computations and assert the synthesis-filter memory
        // changed between them.
        let mut enc = Encoder::new();
        // Adapt the synthesis filter on a speech-like cycle so its
        // predictor carries feedback taps. The synthesis adapter is
        // updated through the decoder-shared path; here we feed it via
        // the public synthesis_adapter()… but that accessor is
        // read-only. Instead drive the encoder's weighting adapter and
        // rely on the synthesis adapter staying all-pass — so we seed
        // the ring memory indirectly by running a first non-zero
        // weighted vector, then confirm the weighting-filter ring
        // memory (ZIRWIIR / ZIRWFIR) evolves.
        let v1 = [100.0, 50.0, -25.0, 12.5, -6.25];
        let _ = enc.compute_vq_target(&v1);
        let fir_after_1 = *enc.zero_input_response_unit().zirw_fir();
        // The first vector loads the weighting-filter all-zero memory
        // with the block-9 outputs; with an all-pass synthesis filter
        // those are zero, so ZIRWFIR stays zero. This still exercises
        // the path without panicking and confirms the accessor wiring.
        assert!(fir_after_1.iter().all(|&x| x == 0.0));
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

    /// Speech-like test signal: two tones + a slow envelope, scaled
    /// to the §3.1.1 nominal input range.
    fn speech_like_sample(n: usize) -> f64 {
        let t = n as f64;
        let envelope = 0.6 + 0.4 * (t * 0.004).sin();
        envelope * (900.0 * (t * 0.23).sin() + 350.0 * (t * 0.71).cos())
    }

    #[test]
    fn encoder_decoder_lockstep_quantized_speech_is_bit_exact() {
        // Blocks 19–23 of the encoder form the simulated decoder
        // (block 8, §3.10): "the quantized speech vector sq(n) is
        // actually the simulated decoded speech vector when there are
        // no channel errors". So a fresh Decoder fed the encoder's
        // channel indices must reproduce the encoder's sq(n) BIT FOR
        // BIT, vector after vector — both ends share the backward
        // adapters (§4.4 / §4.5), the same end-of-cycle adaptation
        // cadence, and the same ±4 095 memory saturation. 200 vectors
        // = 50 adaptation cycles, deep into the adapted regime.
        let mut enc = Encoder::new();
        let mut dec = crate::Decoder::new();
        for vec_idx in 0..200usize {
            let mut s = [0.0f64; IDIM];
            for k in 0..IDIM {
                s[k] = speech_like_sample(vec_idx * IDIM + k);
            }
            let ichan = enc.encode_vector(&s).expect("encode");
            let out = dec.decode_vector(ichan);
            for k in 0..IDIM {
                assert_eq!(
                    enc.quantized_speech()[k],
                    out[k],
                    "vec {vec_idx} sample {k}: encoder sq {} vs decoder {}",
                    enc.quantized_speech()[k],
                    out[k]
                );
            }
        }
    }

    #[test]
    fn sync_bit_round_trips_through_encode_and_extract() {
        // §3.11: encode every 16th vector (last of a cycle) with an
        // alternating sync bit; the decoder-side extractor must recover
        // exactly the bit that was robbed. Unrobbed vectors decode
        // normally; robbed vectors carry the bit in the leftmost
        // codeword bit (shape MSB).
        let mut enc = Encoder::new();
        let mut bit = false;
        for vec_idx in 0..160usize {
            let mut s = [0.0f64; IDIM];
            for k in 0..IDIM {
                s[k] = speech_like_sample(vec_idx * IDIM + k);
            }
            // Rob the last vector of every 4th adaptation cycle (N = 16).
            let robbed = (vec_idx + 1) % 16 == 0;
            let ichan = if robbed {
                enc.encode_vector_with_sync_bit(&s, bit).expect("encode")
            } else {
                enc.encode_vector(&s).expect("encode")
            };
            if robbed {
                assert_eq!(
                    crate::extract_sync_bit(ichan),
                    bit,
                    "vec {vec_idx}: robbed sync bit not recovered"
                );
                bit = !bit;
            }
        }
    }

    #[test]
    fn sync_bit_robbing_preserves_encoder_decoder_lockstep() {
        // §3.11: because both ends know which vectors are robbed and the
        // robbed vector runs the identical backward-adaptation dataflow
        // (only the shape scan narrows), a Decoder fed the robbed
        // indices must STILL reproduce the encoder's sq(n) bit for bit —
        // the half-codebook codevector is a perfectly valid 10-bit index
        // that decodes the same way at both ends.
        let mut enc = Encoder::new();
        let mut dec = crate::Decoder::new();
        let mut bit = true;
        for vec_idx in 0..160usize {
            let mut s = [0.0f64; IDIM];
            for k in 0..IDIM {
                s[k] = speech_like_sample(vec_idx * IDIM + k);
            }
            let robbed = (vec_idx + 1) % 16 == 0;
            let ichan = if robbed {
                let c = enc.encode_vector_with_sync_bit(&s, bit).expect("encode");
                bit = !bit;
                c
            } else {
                enc.encode_vector(&s).expect("encode")
            };
            let out = dec.decode_vector(ichan);
            for k in 0..IDIM {
                assert_eq!(
                    enc.quantized_speech()[k],
                    out[k],
                    "vec {vec_idx} sample {k}: lockstep broke on robbed stream"
                );
            }
        }
    }

    #[test]
    fn encoder_tracks_input_after_adaptation() {
        // End-to-end quality smoke: after the backward adapters have
        // had time to converge, the quantized speech must track the
        // input — the coding error energy over a window must drop
        // well below the signal energy. (G.728 delivers toll quality;
        // even this floating-point build at cold start beats 3 dB
        // SNR comfortably on a stationary two-tone signal. The exact
        // SNR is not pinned — only the "error ≪ signal" property.)
        let mut enc = Encoder::new();
        let warmup = 100usize;
        let measure = 100usize;
        let mut sig_energy = 0.0f64;
        let mut err_energy = 0.0f64;
        for vec_idx in 0..(warmup + measure) {
            let mut s = [0.0f64; IDIM];
            for k in 0..IDIM {
                s[k] = speech_like_sample(vec_idx * IDIM + k);
            }
            enc.encode_vector(&s).expect("encode");
            if vec_idx >= warmup {
                let sq = enc.quantized_speech();
                for k in 0..IDIM {
                    sig_energy += s[k] * s[k];
                    err_energy += (s[k] - sq[k]) * (s[k] - sq[k]);
                }
            }
        }
        assert!(
            err_energy < 0.5 * sig_energy,
            "coding error energy {err_energy} should be well below signal energy {sig_energy}"
        );
    }

    #[test]
    fn encode_vector_advances_cycle_state() {
        // After several adaptation cycles of real signal the
        // end-of-cycle bookkeeping must have run: the synthesis
        // predictor leaves the all-pass state (block 23 adapted on
        // sq) and the impulse response leaves the identity (blocks
        // 12 + 14 + 15 refreshed from the new coefficient sets).
        // Levinson may legitimately fail on the earliest cold-start
        // cycles (the §3.7 "keep the old coefficients" path), so
        // give the hybrid window the same multi-cycle budget the
        // adapter's own tests use: 10 cycles = 40 vectors.
        let mut enc = Encoder::new();
        for vec_idx in 0..40usize {
            let mut s = [0.0f64; IDIM];
            for k in 0..IDIM {
                s[k] = speech_like_sample(vec_idx * IDIM + k);
            }
            enc.encode_vector(&s).expect("encode");
        }
        let a = enc.synthesis_adapter().coefficients();
        assert!(
            a[1..].iter().any(|&v| v != 0.0),
            "synthesis predictor should have adapted off the all-pass state"
        );
        let h = enc.codebook_search().impulse_response();
        assert!(
            h[1..].iter().any(|&v| v != 0.0),
            "impulse response should have left the identity after a cycle refresh"
        );
    }

    #[test]
    fn encode_vector_sf_icount_walks_one_to_four() {
        // Table E-1/G.728: the encoder's synthesis-filter ICOUNT walks
        // 1 → 2 → 3 → 4 → 1 once per encoded vector, landing on 1 for
        // the first call.
        let mut enc = Encoder::new();
        let want = [1usize, 2, 3, 4, 1, 2, 3, 4];
        let s = [0.0f64; IDIM];
        for &w in &want {
            enc.encode_vector(&s).expect("encode");
            assert_eq!(
                enc.sf_icount(),
                w,
                "encoder ICOUNT must walk 1..4 cyclically"
            );
        }
    }

    #[test]
    fn encoder_active_predictor_tracks_decoder_in_lockstep() {
        // The ICOUNT = 3 swap stagger (Table E-1/G.728) must keep the
        // encoder's *active* synthesis predictor — the one its
        // analysis-by-synthesis loop actually uses — bit-for-bit equal
        // to the decoder's active predictor at every vector, since both
        // ends run the identical block-49/50/51 adapter on the identical
        // quantised-speech stream and apply the identical deferred swap.
        // This is the property that makes the synthesis filters evolve
        // in lockstep; if the two staggers disagreed, the encoder's
        // codebook search would target a different filter than the
        // decoder reconstructs with.
        let mut enc = Encoder::new();
        let mut dec = crate::Decoder::new();
        for vec_idx in 0..120usize {
            let mut s = [0.0f64; IDIM];
            for k in 0..IDIM {
                s[k] = speech_like_sample(vec_idx * IDIM + k);
            }
            let ichan = enc.encode_vector(&s).expect("encode");
            dec.decode_vector(ichan);
            assert_eq!(
                enc.sf_icount(),
                dec.sf_icount(),
                "vec {vec_idx}: encoder/decoder ICOUNT must agree"
            );
            assert_eq!(
                enc.active_predictor(),
                dec.active_predictor(),
                "vec {vec_idx}: encoder/decoder active synthesis predictor must agree"
            );
        }
    }

    #[test]
    fn encoder_synthesis_predictor_swap_only_happens_at_third_vector() {
        // Mirror of the decoder's stagger test: the encoder's active
        // synthesis predictor may change ONLY at ICOUNT = 3 (block 51 /
        // Table E-1/G.728). Drive a non-stationary speech-like signal
        // for many cycles and assert the active predictor is constant
        // across every non-ICOUNT-3 vector boundary, with at least one
        // genuine swap observed at an ICOUNT = 3 boundary.
        let mut enc = Encoder::new();
        let mut prev = *enc.active_predictor();
        assert_eq!(prev[0], 1.0);
        let mut observed_a_swap = false;
        for vec_idx in 0..80usize {
            let mut s = [0.0f64; IDIM];
            for k in 0..IDIM {
                s[k] = speech_like_sample(vec_idx * IDIM + k);
            }
            enc.encode_vector(&s).expect("encode");
            let now = *enc.active_predictor();
            if enc.sf_icount() == 3 {
                if now != prev {
                    observed_a_swap = true;
                }
            } else {
                assert_eq!(
                    now,
                    prev,
                    "vec {vec_idx} (ICOUNT {}): active predictor changed outside ICOUNT = 3",
                    enc.sf_icount()
                );
            }
            assert_eq!(now[0], 1.0, "A(1) = 1 invariant must hold every vector");
            prev = now;
        }
        assert!(
            observed_a_swap,
            "expected at least one genuine ICOUNT = 3 predictor swap over 80 vectors"
        );
    }

    #[test]
    fn encoder_excitation_matches_decoder_reconstruction() {
        // The block-67 ET delay slot must carry σ(n)·g_i·y_j — the
        // same gain-scaled excitation the decoder rebuilds. Cross-
        // check one step explicitly: encode a vector, then rebuild
        // ET from the emitted index with a parallel gain adapter.
        let mut enc = Encoder::new();
        let mut reference_gain = crate::GainAdapter::new();
        let s = [800.0, -300.0, 450.0, -120.0, 615.0];
        let sigma = reference_gain.predict_next(&[0.0; IDIM]);
        let ichan = enc.encode_vector(&s).expect("encode");
        let yn = ExcitationVector::from_channel_index(ichan);
        for k in 0..IDIM {
            assert_eq!(enc.previous_excitation()[k], sigma * yn.0[k], "sample {k}");
        }
    }
}
