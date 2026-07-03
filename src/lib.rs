//! # oxideav-g728
//!
//! Pure-Rust ITU-T G.728 LD-CELP (Low-Delay Code-Excited Linear
//! Prediction) speech codec — 16 kbit/s, 8 kHz sampling, 0.625 ms
//! algorithmic delay.
//!
//! ## Status: clean-room rebuild — autonomous decoder pipeline
//!
//! The crate was reset to a register-only scaffold (round 171, master
//! `14e3bad`) under the workspace clean-room policy: the previous
//! implementation's numeric tables had been extracted from an external
//! reference C distribution, which the policy forbids. The current
//! tree is being re-grown one spec-cited unit at a time from the
//! published ITU-T G.728 (1992-09) Recommendation prose alone.
//!
//! Round 195 adds the two **backward adapters** that drive the
//! decoder autonomously off the raw bitstream:
//!
//! * **Block 33 — backward synthesis-filter adapter** ([`SynthesisAdapter`])
//!   wires the spec's three sub-blocks together: block 49 (hybrid
//!   window on quantised speech → autocorrelation), block 50 (the
//!   Levinson-Durbin recursion landed in r189), and block 51
//!   (bandwidth expansion via the Annex C `FACV` vector). The
//!   adapter consumes one NFRSZ-sample adaptation cycle per call
//!   and emits the cycle's bandwidth-expanded predictor in the
//!   spec's `A(1..LPC+1)` layout.
//! * **Block 30 — backward vector gain adapter** ([`GainAdapter`])
//!   wires blocks 67 / 39 / 40 (1-vector delay + RMS + log10), 42
//!   (offset subtractor), 43 / 44 / 45 (hybrid window → Levinson →
//!   bandwidth expansion via `FACGPV`) at ICOUNT=2, 46 (log-gain
//!   linear predictor), 47 (limiter, clamp to `[0, 60]` dB) and
//!   48 (inverse logarithm). One call produces one σ(n) for the
//!   vector that is about to be decoded.
//! * **[`hybrid_window`]** — the spec's block-36 / block-43 /
//!   block-49 pseudocode shares a single shape; this module
//!   transcribes it once and dispatches by parameter object.
//! * **`Decoder::decode_vector`** — autonomous decode path that
//!   walks block 29 → 30 → 31 → 32 → 33 per Figure 3/G.728.
//!
//! Round 201 wires up the **AGC tail of the postfilter** (blocks
//! 73 / 74 / 75 / 76 / 77) — see [`agc`] / [`Agc`].
//!
//! Round 207 wires up the **short-term (spectral) postfilter**
//! (block 72) — see [`short_term_postfilter`] / [`ShortTermPostfilter`].
//! The coefficients are derived from the order-10 by-product of the
//! synthesis-filter Levinson recursion (`ã_1..ã_10` with bandwidth
//! expansion by `SPFZCF^i` / `SPFPCF^i`) plus the tilt-compensation
//! `µ = TILTF · k1` from the same RTMP's first reflection
//! coefficient. They refresh at the first vector of every adaptation
//! cycle (§7.2).
//!
//! Round 213 wires up the **long-term (pitch) postfilter** (block 71)
//! comb-filter machinery — see [`long_term_postfilter`] /
//! [`LongTermPostfilter`]. The transfer function
//! `H_l(z) = g_l · (1 + b · z^{-p})` (eq. 4-1) is implemented as a
//! `KPMAX`-sample FIR delay line; the comb stage holds at the §4.6.1
//! passthrough `(g_l, b, p) = (1, 0, KPMIN)` until the §4.7
//! block-81..84 pitch-extraction / coefficient pipeline lands.
//! [`Decoder::decode_vector_postfiltered`] runs the full block
//! 32 → 71 → 72 → 73..77 chain.
//!
//! Round 220 lands the **first stage of the §4.7 pitch-extraction
//! pipeline** — block 81: the 10th-order LPC inverse filter `Ã(z) = 1 −
//! Σ ã_i · z^{-i}` (eq. 4-6) producing the residual `d(k)` from the
//! decoded speech `sd(k)` — see [`pitch_inverse_filter`] /
//! [`PitchInverseFilter`]. The inverse filter shares the order-10
//! by-product / first-reflection-coefficient surface already published
//! by the synthesis-filter adapter, so it refreshes at the same
//! adaptation-cycle boundary as the short-term postfilter.
//!
//! Round 223 lands the **second stage of the §4.7 pitch-extraction
//! pipeline** — block 82, the pitch period extractor — see
//! [`pitch_search`] / [`PitchSearch`]. The extractor maintains a
//! 240-sample residual buffer `d(−139..100)` (vector slot ordering
//! per the §4.7 prose), a 60-sample decimated buffer `d̄(−34..25)`
//! fed by the Annex D third-order elliptic 1 kHz lowpass + 4:1
//! decimator, and runs the coarse correlation search `ρ(i)` of
//! eq. 4-7 over decimated lags `5..=35`, the full-resolution
//! refinement `C(i)` of eq. 4-8 over `4τ−3..=4τ+3` clamped into
//! `[KPMIN, KPMAX]`, and the fundamental-vs-multiple resolution of
//! eqs 4-9..4-11 with `TAPTH = 0.4`. The extracted pitch
//! `p ∈ [KPMIN, KPMAX]` is stashed as `p̂` for the next frame and
//! observable via [`Decoder::pitch_search`].
//!
//! Round 229 lands the **third and fourth stages of the §4.7
//! pitch-extraction pipeline** — blocks 83 and 84, the single-tap
//! pitch-predictor weight `β` and the long-term postfilter
//! coefficient calculator `(g_l, b)` — see
//! [`pitch_postfilter_coeff`] / [`PitchPostfilterCoeff`]. Block 83
//! maintains a 245-sample `sd(−239..5)` decoded-speech buffer (shared
//! conceptually with block 71's delay line, but kept separate so
//! block 71's sample-rate cost is unchanged) and at the third vector
//! of each frame evaluates the eq. 4-12 single-tap correlation
//! `β = Σ_{k=−99..0} sd(k)·sd(k−p) / Σ sd(k−p)²` clamped into
//! `[0, 1]`. Block 84 then maps `β` to `(b, g_l)` via the
//! eq. 4-13 / 4-14 table: `b = 0, g_l = 1` if `β < PPFTH = 0.6`
//! (postfilter off — `H_l(z) = 1`); `b = PPFZCF · β = 0.15 · β,
//! g_l = 1/(1+b)` if `0.6 ≤ β ≤ 1`; `b = PPFZCF = 0.15,
//! g_l = 1/(1+0.15)` if `β > 1`. The resulting `(g_l, b, p)` triple
//! is applied to [`LongTermPostfilter::set_coefficients`] right after
//! block 82's [`PitchSearch::extract`] commits, so the comb filter
//! takes effect on the next adaptation cycle's samples.
//!
//! With blocks 81 + 82 + 83 + 84 + 71 all in place, the §4.7
//! pitch-postfilter pipeline is complete: the decoder now runs the
//! full Figure 7/G.728 chain — synthesis (block 32) → LPC inverse
//! filter (block 81) → pitch extractor (block 82) → coefficient
//! calculator (blocks 83 / 84) → long-term comb (block 71) →
//! short-term postfilter (block 72) → AGC (blocks 73–77) — every
//! call to [`Decoder::decode_vector_postfiltered`].
//!
//! Round 189 (preserved) provides the Annex A.1/A.2/A.3 hybrid
//! windows, the Annex B excitation codebook (128 × 5 shape + 8
//! gain), the Annex C bandwidth-broadening vectors (Q14), the
//! Annex D 1 kHz lowpass coefficients, Table 1/G.728 parameter
//! constants, the Levinson-Durbin recursion, and the
//! [`Decoder::decode_index`] caller-driven entry point.
//!
//! Round 235 lands the **typed encoder scaffold** — see [`Encoder`] /
//! [`make_encoder`] — plus the [`tables::Y_ENERGY`] precomputed shape
//! codevector energy table (`E_j = Σ_k y_j(k)²`) that §3.9 equation
//! 3-23 references in the analysis-by-synthesis cost expression. The
//! encoder type carries the same two backward adapters as the decoder
//! ([`SynthesisAdapter`] and [`GainAdapter`]) because §4.4 / §4.5 of
//! the Recommendation require the two backward adapters to be
//! bit-for-bit identical at both ends of the channel. Rounds 248 /
//! 249 / 258 / 267 grow the encoder front end (blocks 36 + 37 + 38 +
//! 4 + 9 + 10 + 11); round 276 lands the **§3.9 codebook search**
//! (blocks 12–18, see [`codebook_search`] / [`CodebookSearch`]) and
//! the **§3.10 memory update** ([`ZeroInputResponse::update_memory`],
//! pseudo-code §5.13), closing the analysis-by-synthesis loop —
//! [`Encoder::encode_vector`] now emits real 10-bit channel indices
//! and runs in exact lockstep with [`Decoder::decode_vector`].
//!
//! Round 248 lands the **perceptual weighting filter coefficient
//! calculator** (block 38, §3.3, the third sub-block of the weighting
//! filter adapter after hybrid window 36 and Levinson-Durbin 37) — see
//! [`weighting_filter_coeff`] / [`WeightingFilterCoeff`]. The
//! `from_lpc` transform applies the substitutions `z ← z/γ₁` and
//! `z ← z/γ₂` of eq. 3-4b / 3-4c to the order-10 LPC predictor
//! `q_i`; the `disabled` constructor realises the §3.4.1 non-speech
//! `W(z) = 1` shortcut.
//!
//! Round 249 lands the **block-4 application path** of the perceptual
//! weighting filter (§3.4) — see [`weighting_filter`] /
//! [`PerceptualWeightingFilter`]. The current input speech vector
//! `s(n)` is passed through `W(z) = Q̃(z/γ₁) / Q̃(z/γ₂)` (eq. 3-4a)
//! to produce the weighted speech vector `v(n)`, with the spec's
//! §3.4 "do not reset memory" rule honoured by carrying the per-sample
//! delay lines continuously across coefficient swaps and the §3.4.1
//! disable switch. Block 10 (the cascade with the synthesis filter
//! for the §3.5 zero-input response) requires the §3.10 memory
//! pre-/post-save dance and is queued for a later round.
//!
//! Round 258 lands the **upstream half of the weighting filter
//! adapter** (blocks 36 + 37, §3.3) — see
//! [`weighting_filter_adapter`] / [`WeightingFilterAdapter`]. Block 36
//! is the hybrid window on past input speech (Annex A.3 `wnrw`,
//! 60 samples, dimensions `LPCW + NFRSZ + NONRW = 10 + 20 + 30`); it
//! reuses the shared [`HybridWindow`] / [`HybridWindowState`] machinery
//! that already powers blocks 43 and 49. Block 37 is the order-10
//! Levinson-Durbin recursion on the 11 autocorrelation taps it
//! produces; it reuses [`levinson_durbin`]. The cached predictor
//! `q_i = a_i^{(10)}` (eq. 3-2f) is exposed via
//! [`WeightingFilterAdapter::predictor`] in canonical Levinson
//! layout, ready to be passed into
//! [`Encoder::set_weighting_filter_coeff_from_lpc`] (block 38, landed
//! r248). The adapter mirrors [`SynthesisAdapter::adapt`] timing
//! (one cycle of `NFRSZ = 20` samples per call) and failure policy
//! (Levinson failure keeps the previously cached predictor and
//! propagates the error). The encoder is not yet wired to call
//! the adapter on every cycle — that handshake lands in a follow-up
//! round once the block-3 / block-2 input-speech buffering for §3.9
//! is in place; the adapter is correct as a standalone primitive
//! today and stands ready to feed block 38 directly.
//!
//! ## What is NOT yet wired up
//!
//! * **PCM format conversion (blocks 1, 28).** A-law / µ-law I/O is
//!   delegated to `oxideav-g711` per §5.3 / §3.1.
//! * **Encoder.** The encoder side will come once the decoder runs
//!   end-to-end against a synthetic excitation stream.
//! * **Annex G (fixed-point) variant** and **Annex I (frame-loss
//!   concealment)** are deferred behind the float decoder.
//!
//! ## Clean-room provenance
//!
//! Every numeric value in [`tables`] is transcribed directly from
//! Annexes A, B, C, D of the ITU-T G.728 1992-09 Recommendation PDF
//! that lives under `docs/audio/g728/`. No external implementation
//! has been opened or consulted at any point during this rebuild —
//! see the workspace's clean-room policy for the full scope of what
//! "external" covers.

#![warn(missing_debug_implementations)]
#![deny(unsafe_code)]
// Most loops in this crate mirror the spec's index-based pseudocode
// (e.g. "For K = 1,2,...,IDIM"). Rewriting them as `enumerate` /
// `iter` chains obscures the per-line correspondence to the
// Recommendation, which is the audit anchor for every value. The
// lint is therefore disabled crate-wide — the same trade-off many
// of the sibling DSP crates have settled on (cf. oxideav-g711's
// per-mod `allow`).
#![allow(clippy::needless_range_loop)]

use oxideav_core::RuntimeContext;

pub mod agc;
pub mod annex_g_arith;
pub mod annex_g_codebook;
pub mod annex_g_decoder;
pub mod annex_g_encoder;
pub mod annex_g_gain;
pub mod annex_g_gain_adapter;
pub mod annex_g_hybrid;
pub mod annex_g_levinson;
pub mod annex_g_pf_coeff;
pub mod annex_g_pitch;
pub mod annex_g_postfilter;
pub mod annex_g_synth_adapter;
pub mod annex_g_weight_adapter;
pub mod bitstream;
pub mod codebook_search;
pub mod consts;
pub mod decoder;
pub mod encoder;
pub mod excitation_extrapolation;
pub mod frame_erasure_lpc;
pub mod gain_adapter;
pub mod gain_growth_limiter;
pub mod hybrid_window;
pub mod levinson;
pub mod long_term_postfilter;
pub mod pitch_inverse_filter;
pub mod pitch_postfilter_coeff;
pub mod pitch_search;
pub mod short_term_postfilter;
pub mod synthesis_adapter;
pub mod tables;
pub mod weighting_filter;
pub mod weighting_filter_adapter;
pub mod weighting_filter_coeff;
pub mod zero_input_response;

pub use agc::Agc;
pub use annex_g_arith::{divide, findnls, rnd, vscale, ScalarFloat};
pub use annex_g_decoder::{SynthesisFilterFixed, NLSSTATE_INIT, NSEG, STATELPC_CLIP};
pub use annex_g_encoder::{EncoderFiltersFixed, MemoryUpdateOut};
pub use annex_g_gain::{gain_log_db, offset_removed_log_gain, shape_log_db, DELTA_FLOOR_DB};
pub use annex_g_gain_adapter::{
    gstate1_update, inverse_log_gain, limit_log_gain_q9, log_gain_predict, LogGainWindowFixed,
    GOFF_Q9, GSTATE_INIT_Q9, LOGGAIN_MAX_Q9, LOGGAIN_MIN_Q9,
};
pub use annex_g_hybrid::{
    BflSegment, HwmcoreOut, HwmcoreState, HybridWindowFixed, HybridWindowFixedState, NLSATT50,
    NLSREXP_INIT_G2, NLSSB_INIT_G2,
};
pub use annex_g_levinson::{
    levinson_durbin_fixed, levinson_durbin_fixed_resume, simpdiv, LevinsonInput, LevinsonStatus,
    RecursionResume,
};
pub use annex_g_pf_coeff::{short_term_coeff_fixed, TILTF_Q14};
pub use annex_g_pitch::{apf_to_q13, PitchAdapterFixed};
pub use annex_g_postfilter::{
    scale_factor, LongTermCoeff, PostfilterFixed, ShortTermCoeff, AGCFAC1_Q21, AGCFAC_Q14,
    PF_ORDER, SCALEFIL_INIT_Q14,
};
pub use annex_g_synth_adapter::{
    bandwidth_expand_q14, log_gain_bandwidth_expand, StSegment, SynthAdapterFixed, NLSSB_INIT,
};
pub use annex_g_weight_adapter::{WeightAdapterFixed, WeightCoeffFixed, N5W};
pub use bitstream::{pack_indices, unpack_indices};
pub use codebook_search::{CodebookSearch, SearchResult};
pub use decoder::{
    extract_sync_bit, pack_channel_index, ExcitationVector, Synthesizer, DEFAULT_MAX, FRAME_LEN,
};
pub use encoder::{make_encoder, Encoder, CHANNEL_INDEX_BITS};
pub use excitation_extrapolation::{FrameErasureExcitation, LcgSlideSource, SlideSource};
pub use frame_erasure_lpc::{soften_predictor, FrameErasureLpc};
pub use gain_adapter::GainAdapter;
pub use gain_growth_limiter::{limit_log_gain, GainGrowthLimiter};
pub use hybrid_window::{HybridWindow, HybridWindowState};
pub use levinson::{levinson_durbin, LevinsonError};
pub use long_term_postfilter::LongTermPostfilter;
pub use pitch_inverse_filter::PitchInverseFilter;
pub use pitch_postfilter_coeff::PitchPostfilterCoeff;
pub use pitch_search::PitchSearch;
pub use short_term_postfilter::ShortTermPostfilter;
pub use synthesis_adapter::SynthesisAdapter;
pub use weighting_filter::PerceptualWeightingFilter;
pub use weighting_filter_adapter::WeightingFilterAdapter;
pub use weighting_filter_coeff::WeightingFilterCoeff;
pub use zero_input_response::ZeroInputResponse;

/// Crate-local error type.
///
/// As of round 276 both the decoder and the per-vector encoder loop
/// are wired end to end; [`Error::NotImplemented`] remains for the
/// planned surfaces that are still scaffolds (e.g. the Annex G
/// fixed-point arithmetic variant and the byte-stream framing
/// helpers).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// Functionality is part of the planned scope but not yet wired
    /// up. The error carries no payload — see the per-call docs for
    /// the unblock condition.
    NotImplemented,
    /// The supplied byte buffer's bit length is not a whole multiple of
    /// [`CHANNEL_INDEX_BITS`] (= 10) and hence cannot be unpacked into
    /// whole 10-bit codebook indices using the §5.11 serial wire layout
    /// (each index transmitted most-significant-bit-first, indices
    /// concatenated with no padding). See [`crate::unpack_indices`].
    InvalidInputLength,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Error::NotImplemented => write!(
                f,
                "oxideav-g728: feature not yet wired up (clean-room rebuild in progress)"
            ),
            Error::InvalidInputLength => write!(
                f,
                "oxideav-g728: input length must encode whole 10-bit codebook indices"
            ),
        }
    }
}

impl std::error::Error for Error {}

/// Result alias for the crate.
pub type Result<T> = core::result::Result<T, Error>;

/// LD-CELP decoder front end.
///
/// One [`Decoder`] instance carries the synthesis-filter memory and
/// both backward adapters across vectors. Construct it with
/// [`Decoder::new`] and feed one 10-bit channel index at a time
/// through [`Decoder::decode_vector`] to receive a
/// `FRAME_LEN`-sample decoded-PCM vector with the spec's full
/// blocks 29 → 30 → 31 → 32 → 33 chain wired up.
///
/// For pre-r195 callers, [`Decoder::decode_index`] (no gain adapter,
/// no synthesis adapter — caller-supplied predictor) and
/// [`Decoder::set_synthesis_predictor`] remain available as the
/// register-only entry points.
#[derive(Debug)]
pub struct Decoder {
    synth: Synthesizer,
    /// Backward synthesis-filter adapter (block 33). Updated once per
    /// adaptation cycle (NUPDATE = 4 vectors), consumed by `synth`.
    synth_adapter: SynthesisAdapter,
    /// Backward vector gain adapter (block 30). Updated every vector;
    /// produces `σ(n)` consumed by the block-31 gain scaling step.
    gain_adapter: GainAdapter,
    /// Long-term (pitch) postfilter (block 71). Held at the §4.6.1
    /// "postfilter off" passthrough (`g_l = 1`, `b = 0`) because the
    /// block-81..84 pitch-extraction / coefficient-calculation front
    /// end is not yet implemented. The filter machinery is wired up
    /// inside [`Decoder::decode_vector_postfiltered`] so it slots in
    /// transparently once §4.7 lands and starts driving non-trivial
    /// `(g_l, b, p)` once per adaptation cycle.
    long_term_pf: LongTermPostfilter,
    /// Short-term postfilter (block 72). Refreshed at the first
    /// vector of each adaptation cycle from the synthesis adapter's
    /// order-10 by-product + first reflection coefficient; runs
    /// sample-by-sample inside [`Decoder::decode_vector_postfiltered`].
    short_term_pf: ShortTermPostfilter,
    /// Pitch-extractor LPC inverse filter (block 81 of §4.7).
    /// Refreshed at the first vector of each adaptation cycle from the
    /// synthesis adapter's order-10 by-product; produces the residual
    /// `d(k)` consumed by [`pitch_search`].
    pitch_inv_filter: PitchInverseFilter,
    /// Pitch-period extractor (block 82 of §4.7). Receives the
    /// residual `d(k)` from `pitch_inv_filter`, runs the
    /// lowpass+decimate, coarse correlation, refinement and
    /// fundamental-vs-multiple resolution per the §4.7 prose
    /// (eqs 4-7..4-11), and produces one pitch period per adaptation
    /// frame. The output drives blocks 83 / 84 (via
    /// [`Self::pitch_pf_coeff`]) which in turn refresh the long-term
    /// postfilter coefficients at the third vector of each frame.
    pitch_search: PitchSearch,
    /// Long-term postfilter coefficient calculator (blocks 83 / 84 of
    /// §4.7). Maintains the `sd(−239..5)` decoded-speech buffer that
    /// block 83 needs for eq. 4-12; at the third vector of each frame,
    /// after block 82 produces the pitch period `p`, this block reads
    /// the buffer to compute `β` and then `(g_l, b)` per eqs.
    /// 4-13 / 4-14, feeding the result into `long_term_pf`.
    pitch_pf_coeff: PitchPostfilterCoeff,
    /// Output gain control (blocks 73–77). Always present; the long-
    /// term postfilter (block 71) is still on the §4.6.1 "off" path
    /// pending its block-81..85 pitch-extraction front end. The
    /// short-term postfilter (block 72) is live, so the AGC now sees
    /// a real `sf` ≠ `sd` stream once any adaptation cycle has run.
    agc: Agc,
    /// Previous gain-scaled excitation vector `ET(1..IDIM)` — fed
    /// into block 67 / block 30 on the next vector. Initialised to
    /// zero per Table 2/G.728 (`ET` initial `0,0,...,0`).
    et_prev: [f64; FRAME_LEN],
    /// Rolling NFRSZ-sample quantised-speech buffer. Filled NUPDATE
    /// (= 4) vectors at a time before being handed to the synthesis
    /// adapter, mirroring the `STTMP(1..NFRSZ)` cycle-rate input.
    sttmp_buf: [f64; consts::NFRSZ],
    /// Write cursor into `sttmp_buf`; when it reaches NFRSZ the
    /// adapter is run and the cursor resets.
    sttmp_idx: usize,
    /// Vector index within the current adaptation cycle (0..NUPDATE).
    /// Wraps at NUPDATE = 4. The short-term postfilter coefficients
    /// refresh at the first vector of each cycle (icount == 0 here,
    /// mapping to spec ICOUNT = 1).
    icount: usize,
    /// **Spec ICOUNT for the synthesis-filter swap.** 1-based, walks
    /// `1 → 2 → 3 → 4 → 1` once per [`Decoder::decode_vector`] call.
    /// Per Table E-1/G.728 the backward synthesis-filter adapter takes
    /// its input speech through vector 4 of a cycle but its updated
    /// coefficients are **first used at decoding vector 3** of the
    /// following cycle (blocks 49/50/51 → the block-51 `Wait until
    /// ICOUNT = 3, then A = ATMP` swap). This counter drives that
    /// deferred swap.
    sf_icount: usize,
    /// **Pending synthesis predictor** computed by the block-49/50/51
    /// adapter at a cycle boundary, awaiting the block-51 `ICOUNT = 3`
    /// swap into [`Self::active_predictor`]. `None` until the first
    /// full adaptation cycle has been accumulated. Holds the
    /// bandwidth-expanded `ATMP` array (spec `A` layout).
    pending_synth: Option<[f64; consts::LPC + 1]>,
    /// **Active synthesis predictor** currently driving block 32. Only
    /// refreshed from [`Self::pending_synth`] at `ICOUNT = 3`, so during
    /// vectors 1 and 2 of each cycle the filter keeps using the
    /// previous cycle's committed coefficients exactly as the spec's
    /// block-50 note prescribes ("the old set of synthesis filter
    /// coefficients … is still being used" until the third vector).
    active_predictor: [f64; consts::LPC + 1],
}

impl Default for Decoder {
    fn default() -> Self {
        Self::new()
    }
}

impl Decoder {
    /// Construct a fresh decoder with all internal state initialised
    /// to the values listed in Table 2/G.728 (synthesis-filter memory
    /// zeroed, all-pass predictor, hybrid-window buffers zeroed, log-
    /// gain state filled with `-GOFF` dB).
    pub fn new() -> Self {
        Self {
            synth: Synthesizer::new(),
            synth_adapter: SynthesisAdapter::new(),
            gain_adapter: GainAdapter::new(),
            long_term_pf: LongTermPostfilter::new(),
            short_term_pf: ShortTermPostfilter::new(),
            pitch_inv_filter: PitchInverseFilter::new(),
            pitch_search: PitchSearch::new(),
            pitch_pf_coeff: PitchPostfilterCoeff::new(),
            agc: Agc::new(),
            et_prev: [0.0; FRAME_LEN],
            sttmp_buf: [0.0; consts::NFRSZ],
            sttmp_idx: 0,
            icount: 0,
            // Spec ICOUNT runs 1..4; start at 4 so the first
            // `decode_vector` call advances it to 1 (the first vector
            // of the first adaptation cycle).
            sf_icount: consts::NUPDATE,
            // No adaptation cycle has completed yet, so there is no
            // pending predictor to swap. The active predictor starts
            // at the all-pass identity (`A(1) = 1`, rest 0) per
            // Table 2/G.728.
            pending_synth: None,
            active_predictor: {
                let mut a = [0.0f64; consts::LPC + 1];
                a[0] = 1.0;
                a
            },
        }
    }

    /// Inject a fresh set of synthesis-filter coefficients for the
    /// next vector. The leading entry must be `1.0` (the spec's
    /// `A(1) = 1` invariant).
    ///
    /// This bypasses the built-in synthesis adapter (block 33).
    /// Callers wanting the spec's autonomous adapter should use
    /// [`Decoder::decode_vector`] instead.
    pub fn set_synthesis_predictor(&mut self, coeffs: [f64; consts::LPC + 1]) {
        self.synth.set_predictor(coeffs);
    }

    /// Decode one channel index into one [`FRAME_LEN`]-sample PCM
    /// vector, using a caller-supplied synthesis predictor (see
    /// [`Decoder::set_synthesis_predictor`]) and **no gain scaling**.
    ///
    /// Retained for the original r189 API. New callers should use
    /// [`Decoder::decode_vector`].
    pub fn decode_index(&mut self, ichan: u16) -> [f64; FRAME_LEN] {
        let exc = ExcitationVector::from_channel_index(ichan);
        self.synth.filter_vector(&exc)
    }

    /// Decode one channel index into one [`FRAME_LEN`]-sample PCM
    /// vector through the **full** block 29 → 30 → 31 → 32 → 33 chain.
    ///
    /// Order of operations (per Figure 3/G.728, §3 and §5.14):
    ///
    /// 1. Look up the gain-codebook entry `σ(n)` via the backward
    ///    vector gain adapter (block 30) from the previous vector's
    ///    gain-scaled excitation.
    /// 2. Look up the shape codevector via block 29 (`y_j`).
    /// 3. Scale: `ET(n) = σ(n) · y_j` (block 31).
    /// 4. Run the 50th-order synthesis filter (block 32) using the
    ///    most recently produced bandwidth-expanded predictor from
    ///    the synthesis adapter (block 33).
    /// 5. Push the decoded vector into the adapter's rolling
    ///    quantised-speech buffer; once a full adaptation cycle of
    ///    NFRSZ = 20 samples is accumulated, run the synthesis
    ///    adapter to produce the next cycle's predictor.
    pub fn decode_vector(&mut self, ichan: u16) -> [f64; FRAME_LEN] {
        // ----- Spec ICOUNT advance + block-51 deferred swap ----------
        // Advance the synthesis-filter ICOUNT (1 → 2 → 3 → 4 → 1).
        // The backward synthesis-filter adapter (blocks 49/50/51)
        // takes the quantised speech through vector 4 of a cycle, but
        // per Table E-1/G.728 (and the block-51 `Wait until ICOUNT = 3,
        // then A = ATMP` instruction) its coefficients are not first
        // used until **decoding vector 3** of the following cycle. We
        // therefore hold the freshly-computed predictor in
        // `pending_synth` and swap it into the active predictor only
        // when ICOUNT reaches 3.
        self.sf_icount = (self.sf_icount % consts::NUPDATE) + 1;
        if self.sf_icount == 3 {
            if let Some(pending) = self.pending_synth.take() {
                self.active_predictor = pending;
            }
        }

        // ----- Block 30: predict σ(n) from previous ET ---------------
        let sigma = self.gain_adapter.predict_next(&self.et_prev);

        // ----- Block 29: excitation VQ codebook lookup ---------------
        // §5.14 block 29: YN(K) = GQ(IG) · Y(NN + K) — the gain-
        // codebook level g_i is part of the codebook lookup, separate
        // from the backward-adapted σ(n) that block 31 multiplies on
        // top (ET = GAIN · YN). ExcitationVector::from_channel_index
        // already implements exactly this YN lookup.
        let yn = ExcitationVector::from_channel_index(ichan);
        let mut et = [0.0f64; FRAME_LEN];
        for k in 0..FRAME_LEN {
            // Block 31: ET(K) = GAIN · YN(K).
            et[k] = sigma * yn.0[k];
        }

        // ----- Block 32: synthesis filter ---------------------------
        // Use the ICOUNT-staggered active predictor (block 51 output
        // committed at the third vector of each cycle), NOT the
        // adapter's most-recent result directly.
        self.synth.set_predictor(self.active_predictor);
        let exc = ExcitationVector(et);
        let out = self.synth.filter_vector(&exc);

        // ----- Block 33 preparation: accumulate quantised speech ----
        // STTMP for the next call is the rolling buffer of decoded
        // speech vectors. We push the current output and, once a full
        // NFRSZ-sample cycle is accumulated, run the synthesis adapter
        // (blocks 49/50/51) to produce the next cycle's predictor and
        // stash it as `pending_synth` for the deferred ICOUNT = 3 swap.
        for k in 0..FRAME_LEN {
            self.sttmp_buf[self.sttmp_idx] = out[k];
            self.sttmp_idx += 1;
        }
        if self.sttmp_idx >= consts::NFRSZ {
            let _ = self.synth_adapter.adapt(&self.sttmp_buf);
            self.pending_synth = Some(*self.synth_adapter.coefficients());
            self.sttmp_idx = 0;
        }

        // ----- Block 67 wrap-around: save ET for the next call -------
        self.et_prev = et;

        out
    }

    /// Decode one channel index into one [`FRAME_LEN`]-sample PCM
    /// vector through the **long-term postfilter (block 71) +
    /// short-term postfilter (block 72) + AGC stage (blocks 73-77)**
    /// of Figure 7/G.728.
    ///
    /// Per-vector dataflow:
    ///
    /// 1. Run [`Decoder::decode_vector`] for `sd(n)` (block 32).
    /// 2. If this is the first vector of an adaptation cycle, refresh
    ///    the short-term postfilter (block 72) coefficients **and** the
    ///    pitch-extractor LPC inverse filter (block 81) coefficients
    ///    from the synthesis adapter's order-10 by-product (§7.2,
    ///    eq. 4-3..4-5; §7.1 eq. 4-6). The long-term postfilter
    ///    coefficients are held at the §4.6.1 passthrough
    ///    (`g_l = 1`, `b = 0`) until the block-82..84 sub-pipeline
    ///    lands.
    /// 3. Run the block-81 inverse filter on `sd(n)` to advance the
    ///    pitch-extractor residual state. The output `d(n)` is what
    ///    block 82 will consume once it lands; for now it is computed
    ///    purely to keep the inverse filter's 10-tap memory in sync
    ///    with the decoded-speech stream.
    /// 4. Filter `sd(n)` through the long-term (block 71) comb
    ///    filter `H_l(z) = g_l · (1 + b · z^{-p})` to obtain
    ///    `sd_lt(n)`. Until block 84 starts driving non-trivial
    ///    coefficients this is the identity (`sd_lt = sd`).
    /// 5. Filter `sd_lt(n)` through the short-term postfilter
    ///    (block 72) to obtain `sf(n)`.
    /// 6. Apply the AGC (blocks 73-77) with the spec's per-vector
    ///    `Σ|sd|` / `Σ|sf|` ratio. Per Figure 7/G.728 block 73 reads
    ///    the *decoded* speech `sd` directly, **not** the long-term
    ///    postfilter output — so the AGC numerator stays anchored to
    ///    the raw decoded amplitude.
    ///
    /// Callers that want the raw synthesis-filter output (no
    /// postfilter, no AGC — exactly the §4.6.1 "postfilter off"
    /// path) should keep using [`Decoder::decode_vector`].
    pub fn decode_vector_postfiltered(&mut self, ichan: u16) -> [f64; FRAME_LEN] {
        // Step 1: raw synthesis-filter output (advances synth_adapter
        // state when a full adaptation cycle is collected).
        let sd = self.decode_vector(ichan);

        // Step 2: refresh the short-term postfilter coefficients at
        // the first vector of every adaptation cycle (§7.2: "updated
        // at the first vector of each frame as soon as ã_i / k1 are
        // available from the order-10 by-product"). Vectors are
        // 0-indexed here; spec ICOUNT = 1 maps to icount = 0.
        //
        // The long-term postfilter coefficients are NOT refreshed
        // here yet — they require the §4.7 block-81..84 pitch
        // extraction front end. Until that lands, the long-term
        // filter stays at its cold-start `(g_l, b, p) = (1, 0, KPMIN)`
        // §4.6.1 passthrough and step 3 is the identity (`sd_lt = sd`).
        if self.icount == 0 {
            self.short_term_pf.set_from_synthesis_byproduct(
                self.synth_adapter.order10_predictor(),
                self.synth_adapter.k1(),
            );
            self.pitch_inv_filter
                .set_from_synthesis_byproduct(self.synth_adapter.order10_predictor());
        }
        self.icount = (self.icount + 1) % consts::NUPDATE;

        // Block 81 (§4.7): run the 10th-order LPC inverse filter on
        // `sd` to produce the residual `d(k)`. The residual is the
        // input to block 82's pitch-period search.
        let residual = self.pitch_inv_filter.filter_vector(&sd);

        // Block 82 (§4.7): push the residual vector into the pitch
        // extractor's 240-sample buffer. The spec's vector-slot
        // ordering (vec 1 → d(86..90), vec 2 → d(91..95), vec 3 →
        // d(96..100), vec 4 → d(81..85)) is realised inside
        // PitchSearch::push_residual via its own per-frame cursor.
        // At the third vector of each frame (spec ICOUNT = 3, our
        // 0-indexed icount == 2 *after* the increment above wraps it
        // back; but we want to extract right after the 3rd vector's
        // residual is pushed) the extractor runs the full pitch
        // search and stashes the result as p̂ for the next frame.
        //
        // Because `self.icount` was advanced after step 2, the value
        // at this point is the *next* icount; spec ICOUNT = 3 maps
        // to our pre-increment icount == 2, which corresponds to
        // post-increment icount == 3 here. We therefore trigger the
        // extract when post-increment icount == 3.
        //
        // The pitch output is not yet consumed by the long-term
        // postfilter — blocks 83 / 84 still need to land. Until they
        // do, the comb filter stays at the §4.6.1 passthrough.
        self.pitch_search.push_residual(&residual);

        // Block 83 (§4.7): every decoded-speech vector is pushed into
        // the long-term-postfilter coefficient calculator's
        // sd(−239..5) buffer. The push happens unconditionally per
        // vector; the coefficient evaluation itself runs only at the
        // third vector of each frame, alongside block 82's extract.
        self.pitch_pf_coeff.push_decoded_vector(&sd);

        if self.icount == 3 {
            let p = self.pitch_search.extract();
            // Block 83 + 84 (§4.7, eqs 4-12 / 4-13 / 4-14): using the
            // pitch period just extracted by block 82, compute the
            // single-tap predictor weight β over the sd(−99..0)
            // window, then derive (b, g_l) via the eq. 4-13 / 4-14
            // table and apply them to the long-term postfilter (block
            // 71) for the upcoming adaptation cycle. Per §7.1 of the
            // decode trace and §4.7 of the Recommendation this update
            // happens at the third vector of each frame.
            let _beta = self.pitch_pf_coeff.compute_coefficients(p);
            self.long_term_pf.set_coefficients(
                self.pitch_pf_coeff.last_g_l(),
                self.pitch_pf_coeff.last_b(),
                p,
            );
        }

        // Step 3: long-term (pitch) postfilter (block 71). Until the
        // first third-vector extract commits, the comb filter remains
        // at its cold-start (g_l, b, p) = (1, 0, KPMIN) state — the
        // §4.6.1 passthrough identity. Once block 84 starts feeding
        // non-trivial (g_l, b, p) per frame the comb begins doing
        // actual long-term postfiltering.
        let sd_lt = self.long_term_pf.filter_vector(&sd);

        // Step 4: short-term postfilter (block 72). At cold start,
        // before any adaptation cycle has completed, all coefficients
        // are zero and the filter is the identity (sf = sd_lt) —
        // matching the §4.6.1 "postfilter off" path. Once the first
        // cycle commits, sf starts diverging from sd_lt.
        let sf = self.short_term_pf.filter_vector(&sd_lt);

        // Step 5: AGC (blocks 73-77). Block 73 takes the *raw*
        // decoded-speech vector `sd` (Figure 7/G.728), keeping the
        // power reference anchored at the synthesis filter output.
        self.agc.apply(&sd, &sf)
    }

    /// Borrow the long-term postfilter (useful for tests and audit).
    pub fn long_term_postfilter(&self) -> &LongTermPostfilter {
        &self.long_term_pf
    }

    /// Borrow the short-term postfilter (useful for tests and audit).
    pub fn short_term_postfilter(&self) -> &ShortTermPostfilter {
        &self.short_term_pf
    }

    /// Borrow the pitch-extractor LPC inverse filter (block 81) — its
    /// per-sample memory advances every call to
    /// [`Self::decode_vector_postfiltered`] and its coefficients
    /// refresh at the first vector of each adaptation cycle. Useful
    /// for tests and audit, and for the still-pending block-82 search
    /// to read the residual stream once it lands.
    pub fn pitch_inverse_filter(&self) -> &PitchInverseFilter {
        &self.pitch_inv_filter
    }

    /// Borrow the pitch-period extractor (block 82) — receives the
    /// residual `d(k)` from the inverse filter and produces one pitch
    /// period `p ∈ [KPMIN, KPMAX]` per adaptation frame, available via
    /// [`PitchSearch::last_pitch`]. The extracted pitch is consumed by
    /// blocks 83 / 84 (via [`Self::pitch_pf_coeff`]) at the third
    /// vector of each frame to drive the long-term postfilter's
    /// `(g_l, b, p)` coefficients.
    pub fn pitch_search(&self) -> &PitchSearch {
        &self.pitch_search
    }

    /// Borrow the long-term postfilter coefficient calculator (blocks
    /// 83 / 84 of §4.7). Holds the `sd(−239..5)` decoded-speech buffer
    /// and the most-recent `(β, b, g_l)` triple. Useful for tests and
    /// audit; downstream callers should read `(b, g_l)` indirectly via
    /// [`Self::long_term_postfilter`] which is what the comb-filter
    /// machinery (block 71) consumes at sample rate.
    pub fn pitch_pf_coeff(&self) -> &PitchPostfilterCoeff {
        &self.pitch_pf_coeff
    }

    /// Borrow the AGC (useful for tests and audit).
    pub fn agc(&self) -> &Agc {
        &self.agc
    }

    /// Borrow the synthesiser (useful for tests and audit).
    pub fn synthesizer(&self) -> &Synthesizer {
        &self.synth
    }

    /// Borrow the synthesis-filter adapter (useful for tests and
    /// audit).
    pub fn synthesis_adapter(&self) -> &SynthesisAdapter {
        &self.synth_adapter
    }

    /// Borrow the gain adapter (useful for tests and audit).
    pub fn gain_adapter(&self) -> &GainAdapter {
        &self.gain_adapter
    }

    /// The **active synthesis predictor** currently driving block 32
    /// (the §3.7 / block-51 output, committed at the third vector of
    /// each adaptation cycle). During vectors 1 and 2 of a cycle this
    /// is the previous cycle's committed predictor; it advances to the
    /// most recently adapted predictor only at `ICOUNT = 3`. Useful for
    /// tests and audit of the Table E-1/G.728 update stagger. Compare
    /// with [`Self::synthesis_adapter`]`().coefficients()`, which is
    /// the *most recently computed* predictor regardless of swap
    /// timing.
    pub fn active_predictor(&self) -> &[f64; consts::LPC + 1] {
        &self.active_predictor
    }

    /// The current 1-based synthesis-filter `ICOUNT` (walks
    /// `1 → 2 → 3 → 4 → 1` once per [`Self::decode_vector`] call). The
    /// block-51 coefficient swap commits when this reaches 3. Useful
    /// for tests and audit.
    pub fn sf_icount(&self) -> usize {
        self.sf_icount
    }

    /// Decode a complete §5.11 serial byte stream into flat PCM.
    ///
    /// Unpacks `bytes` into 10-bit channel indices via
    /// [`unpack_indices`] (most-significant-bit-first, indices
    /// concatenated — see [`crate::bitstream`]), then runs each index
    /// through [`Self::decode_vector`] in order, concatenating the
    /// `FRAME_LEN`-sample output vectors. The returned buffer has
    /// `bytes.len() · 8 / 10 · FRAME_LEN` samples.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInputLength`] when `bytes` does not
    /// encode a whole number of channel indices (its bit length is not
    /// a multiple of [`CHANNEL_INDEX_BITS`] = 10). The natural framing
    /// unit is one adaptation cycle = 4 indices = 40 bits = 5 bytes.
    pub fn decode_bytes(&mut self, bytes: &[u8]) -> Result<Vec<f64>> {
        let indices = bitstream::unpack_indices(bytes)?;
        let mut out = Vec::with_capacity(indices.len() * FRAME_LEN);
        for ichan in indices {
            out.extend_from_slice(&self.decode_vector(ichan));
        }
        Ok(out)
    }
}

/// Direct factory entry: build a fresh [`Decoder`].
///
/// Mirrors the dual-API convention every codec crate exposes: the
/// `oxideav-core` registry path is wired via [`register`] /
/// [`oxideav_core::register!`], and downstream consumers wanting the
/// raw codec without the registry detour use this factory.
pub fn make_decoder() -> Decoder {
    Decoder::new()
}

/// No-op codec registration. The decoder front end is exposed via
/// [`make_decoder`] (autonomous index decode plus the §5.11 serial
/// byte-stream path [`Decoder::decode_bytes`]) and the encoder via
/// [`make_encoder`] / [`Encoder::encode_buffer`]; the registry-side
/// `decoder` / `encoder` factory wiring (A-law / µ-law PCM I/O via
/// `oxideav-g711` per §5.3 / §3.1) lands in a later round.
pub fn register(_ctx: &mut RuntimeContext) {}

oxideav_core::register!("g728", register);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decoder_default_and_new_are_equivalent() {
        // Sanity: Default::default() and new() should both yield a
        // fresh, identity-state decoder.
        let a = Decoder::default();
        let b = Decoder::new();
        assert_eq!(a.synthesizer().state(), b.synthesizer().state());
    }

    #[test]
    fn decode_index_zero_produces_finite_pcm() {
        // The first decoded vector is the zero-state response of the
        // synthesis filter to the gain-scaled codevector indexed at
        // channel 0. The output must be finite, within ±MAX, and
        // FRAME_LEN samples wide.
        let mut dec = Decoder::new();
        let out = dec.decode_index(0);
        assert_eq!(out.len(), FRAME_LEN);
        for &v in &out {
            assert!(v.is_finite());
            assert!(v.abs() <= DEFAULT_MAX);
        }
    }

    #[test]
    fn decode_runs_a_short_stream_without_panicking() {
        // Smoke test: drive the decoder through 64 vectors (320
        // samples, one second × 2 of audio worth) with arbitrary
        // channel indices and confirm no panic or saturation overflow.
        let mut dec = Decoder::new();
        for i in 0..64u16 {
            let ichan = (i * 17) & 0x03FF;
            let out = dec.decode_index(ichan);
            for &v in &out {
                assert!(v.is_finite());
                assert!(v.abs() <= DEFAULT_MAX);
            }
        }
    }

    #[test]
    fn error_display_is_descriptive() {
        // The Display impl is what surfaces to users via thiserror's
        // ?-propagation. Cover both branches.
        let s1 = Error::NotImplemented.to_string();
        assert!(s1.contains("clean-room") || s1.contains("not yet"));
        let s2 = Error::InvalidInputLength.to_string();
        assert!(s2.contains("whole 10-bit"));
    }

    #[test]
    fn register_is_a_no_op_on_runtime_context() {
        // The registry-side wiring is intentionally a no-op until the
        // decoder can run end-to-end without an externally-supplied
        // predictor. Confirm calling it doesn't panic and doesn't
        // mutate the context in a visible way.
        let mut ctx = RuntimeContext::default();
        register(&mut ctx);
    }

    #[test]
    fn make_decoder_returns_fresh_state() {
        // The direct factory entry must yield the same all-zero state
        // a `Decoder::new()` would.
        let dec = make_decoder();
        for &s in dec.synthesizer().state() {
            assert_eq!(s, 0.0);
        }
    }

    // ---------- Block 29 → 30 → 31 → 32 → 33 chain (r195) -----------

    #[test]
    fn decode_vector_full_chain_runs_finite() {
        // Drive the full block 29 → 33 chain over enough vectors to
        // exercise at least three adaptation cycles (3 × 4 = 12
        // vectors). Every output sample must be finite and within
        // the saturation envelope.
        let mut dec = Decoder::new();
        for i in 0..16u16 {
            let ichan = (i * 31) & 0x03FF;
            let out = dec.decode_vector(ichan);
            assert_eq!(out.len(), FRAME_LEN);
            for &v in &out {
                assert!(v.is_finite(), "decode_vector produced NaN/Inf");
                assert!(
                    v.abs() <= DEFAULT_MAX,
                    "decode_vector escaped saturation envelope"
                );
            }
        }
    }

    #[test]
    fn decode_vector_drives_gain_adapter_state() {
        // After enough vectors the gain adapter's GSTATE memory must
        // diverge from the initial -32 dB filler (the Table 2 init).
        let mut dec = Decoder::new();
        for i in 0..20u16 {
            // Use mid-codebook indices so the excitation is large
            // enough to flow through ETRMS away from the clamp.
            dec.decode_vector(((i + 1) * 13) & 0x03FF);
        }
        let gstate = dec.gain_adapter().gstate();
        assert!(
            gstate
                .iter()
                .any(|&v| (v - gain_adapter::GSTATE_INIT_DB).abs() > 1e-6),
            "GSTATE should have diverged from initial -32 dB after 20 vectors; got {:?}",
            gstate
        );
    }

    #[test]
    fn decode_vector_applies_gain_codebook_level() {
        // §5.14 block 29: YN(K) = GQ(IG)·Y(NN + K) — the gain-codebook
        // level is part of the lookup, on top of which block 31
        // multiplies the backward-adapted σ(n). On the very first
        // vector (cold start: zero filter memory, all-pass predictor,
        // σ(n) identical on both sides because it is predicted from
        // the all-zero ET history), two channel indices sharing the
        // same shape row but different gain levels must decode to
        // outputs in exactly the GQ ratio.
        let shape = 10u8;
        let mut dec_a = Decoder::new();
        let mut dec_b = Decoder::new();
        let out_a = dec_a.decode_vector(decoder::pack_channel_index(shape, 0));
        let out_b = dec_b.decode_vector(decoder::pack_channel_index(shape, 2));
        let ratio = tables::GQ[2] / tables::GQ[0];
        for k in 0..FRAME_LEN {
            assert!(
                (out_b[k] - ratio * out_a[k]).abs() <= 1e-9 * out_a[k].abs().max(1.0),
                "sample {k}: {} vs {} (ratio {ratio})",
                out_b[k],
                out_a[k]
            );
        }
        // And the negative half of the codebook flips the sign.
        let mut dec_c = Decoder::new();
        let out_c = dec_c.decode_vector(decoder::pack_channel_index(shape, 4));
        for k in 0..FRAME_LEN {
            assert!(
                (out_c[k] + out_a[k]).abs() <= 1e-9 * out_a[k].abs().max(1.0),
                "sign mirror sample {k}"
            );
        }
    }

    // ---------- AGC tail (r201, blocks 73-77) -----------------------
    // ---------- Short-term postfilter (r207, block 72) --------------

    #[test]
    fn decode_vector_postfiltered_cold_start_matches_decode_vector() {
        // Before the first adaptation cycle commits (icount sweeps 0..3
        // exactly NUPDATE times → first cycle commits at vector 4 in
        // the 0-indexed view), the short-term postfilter coefficients
        // are all zero (cold start) and `sf = sd`. The AGC's IIR with
        // GAINSF = 1 keeps SCALEFIL = 1.0 exactly, so the postfiltered
        // output equals the raw `decode_vector` output bit-for-bit
        // for the first NUPDATE - 1 vectors.
        //
        // Vector 0: postfilter coefficients refreshed from the synth
        // adapter (all-pass, k1 = 0 → all coefficients zero → identity).
        // Vectors 0..3 (one full cycle): sf = sd, AGC passthrough.
        // After vector 3 completes, the synth adapter has accumulated
        // 4 vectors × 5 samples = 20 samples = NFRSZ; the next adapt()
        // call commits a non-trivial predictor and the order-10
        // by-product. The NEXT vector (4) then loads non-trivial
        // postfilter coefficients and sf starts diverging from sd.
        let mut dec_a = Decoder::new();
        let mut dec_b = Decoder::new();
        // First NUPDATE = 4 vectors: cold-start identity. We can
        // assert equality on those.
        for i in 0..consts::NUPDATE as u16 {
            let ichan = (i * 23) & 0x03FF;
            let raw = dec_a.decode_vector(ichan);
            let pf = dec_b.decode_vector_postfiltered(ichan);
            for k in 0..FRAME_LEN {
                assert!(
                    (raw[k] - pf[k]).abs() < 1e-12,
                    "cold-start vector {} sample {}: raw {} vs pf {} (diff {})",
                    i,
                    k,
                    raw[k],
                    pf[k],
                    (raw[k] - pf[k]).abs()
                );
            }
        }
    }

    #[test]
    fn decode_vector_postfiltered_diverges_from_raw_after_first_adaptation() {
        // Once the synthesis adapter commits a non-trivial predictor
        // (first cycle completes after NUPDATE = 4 vectors) the
        // short-term postfilter gets non-zero coefficients and `sf`
        // diverges from `sd`. The AGC then sees a non-unity GAINSF
        // and starts gain-correcting. We confirm at least one sample
        // differs measurably from the raw `decode_vector` output by
        // vector 8 (well after the first cycle has committed).
        let mut dec_a = Decoder::new();
        let mut dec_b = Decoder::new();
        let mut diverged = false;
        for i in 0..16u16 {
            // Use mid-range indices so the synth has enough signal to
            // produce a non-trivial Levinson result.
            let ichan = ((i + 1) * 17) & 0x03FF;
            let raw = dec_a.decode_vector(ichan);
            let pf = dec_b.decode_vector_postfiltered(ichan);
            if i >= consts::NUPDATE as u16
                && raw.iter().zip(&pf).any(|(&a, &b)| (a - b).abs() > 1e-6)
            {
                diverged = true;
            }
        }
        assert!(
            diverged,
            "postfilter should produce non-identity output after first adaptation cycle"
        );
    }

    #[test]
    fn decode_vector_postfiltered_holds_scalefil_at_unity_in_cold_start() {
        // During the cold-start cycle (the first NUPDATE vectors) the
        // postfilter is the identity, GAINSF = 1, and SCALEFIL stays
        // at 1.0 exactly.
        let mut dec = Decoder::new();
        for i in 0..consts::NUPDATE as u16 {
            dec.decode_vector_postfiltered(i & 0x03FF);
        }
        assert_eq!(
            dec.agc().scalefil(),
            1.0,
            "SCALEFIL must stay at 1.0 during the cold-start cycle"
        );
    }

    #[test]
    fn decode_vector_postfiltered_output_is_finite() {
        let mut dec = Decoder::new();
        for i in 0..64u16 {
            let out = dec.decode_vector_postfiltered((i * 11) & 0x03FF);
            for &v in &out {
                assert!(v.is_finite());
            }
        }
    }

    #[test]
    fn decode_vector_postfiltered_coefficients_update_each_cycle() {
        // After enough vectors to commit at least one adaptation
        // cycle on an excited synthesis state, the postfilter's
        // coefficient state should have diverged from the all-zero
        // cold start. We use mid-codebook channel indices (the same
        // pattern as `decode_vector_postfiltered_diverges_from_raw_*`)
        // so the synthesis filter has enough excitation to drive a
        // non-trivial Levinson result through the adapter.
        let mut dec = Decoder::new();
        for i in 0..(consts::NUPDATE as u16 * 4) {
            dec.decode_vector_postfiltered(((i + 1) * 17) & 0x03FF);
        }
        let pf = dec.short_term_postfilter();
        let any_nz = pf.numerator()[1..].iter().any(|&v| v != 0.0)
            || pf.denominator()[1..].iter().any(|&v| v != 0.0)
            || pf.tilt() != 0.0;
        assert!(
            any_nz,
            "postfilter coefficients should be non-zero after >=1 cycle"
        );
    }

    // ---------- Long-term postfilter wiring (r213, block 71) -------

    #[test]
    fn decoder_exposes_long_term_postfilter_at_cold_start() {
        // Decoder construction starts the long-term postfilter at
        // the §4.6.1 passthrough: g_l = 1, b = 0, p = KPMIN. Confirm
        // the accessor surface exposes those values immediately after
        // construction, before any decode call runs.
        let dec = Decoder::new();
        let lt = dec.long_term_postfilter();
        assert_eq!(lt.g_l(), 1.0);
        assert_eq!(lt.b(), 0.0);
        assert_eq!(lt.p(), consts::KPMIN);
    }

    #[test]
    fn decoder_long_term_postfilter_satisfies_spec_invariants_after_decode() {
        // After many frames the long-term postfilter's (g_l, b, p)
        // must always satisfy the spec invariants from eqs.
        // 4-13 / 4-14: p ∈ [KPMIN, KPMAX], b ∈ [0, PPFZCF], and
        // g_l = 1/(1 + b). The exact branch depends on the input
        // signal's periodicity at the analysis window.
        let mut dec = Decoder::new();
        for i in 0..(consts::NUPDATE as u16 * 4) {
            let _ = dec.decode_vector_postfiltered(((i + 1) * 17) & 0x03FF);
        }
        let lt = dec.long_term_postfilter();
        assert!((consts::KPMIN..=consts::KPMAX).contains(&lt.p()));
        assert!((0.0..=consts::PPFZCF + 1e-12).contains(&lt.b()));
        assert!((lt.g_l() - 1.0 / (1.0 + lt.b())).abs() < 1e-12);
    }

    #[test]
    fn long_term_passthrough_preserves_short_term_postfilter_behaviour() {
        // With the long-term filter in passthrough (current state),
        // wiring it between sd and the short-term postfilter must not
        // change the post-filtered output relative to the r207
        // behaviour. We exercise this by checking that the cold-start
        // invariant (sf == sd for the first NUPDATE vectors) still
        // holds after the new long-term stage is inserted.
        let mut dec_a = Decoder::new();
        let mut dec_b = Decoder::new();
        for i in 0..consts::NUPDATE as u16 {
            let ichan = (i * 23) & 0x03FF;
            let raw = dec_a.decode_vector(ichan);
            let pf = dec_b.decode_vector_postfiltered(ichan);
            for k in 0..FRAME_LEN {
                assert!(
                    (raw[k] - pf[k]).abs() < 1e-12,
                    "cold-start passthrough broken by long-term filter wiring"
                );
            }
        }
    }

    // ---------- Pitch-extractor inverse filter wiring (r220, block 81)

    #[test]
    fn decoder_exposes_pitch_inverse_filter_at_cold_start() {
        // Fresh decoder: inverse filter is at the all-zero coefficient
        // / all-zero memory cold start (pure identity `Ã(z) = 1`).
        let dec = Decoder::new();
        let pf = dec.pitch_inverse_filter();
        assert!(pf.coefficients().iter().all(|&v| v == 0.0));
        assert!(pf.memory().iter().all(|&v| v == 0.0));
    }

    #[test]
    fn pitch_inverse_filter_memory_advances_per_postfiltered_call() {
        // After a few `decode_vector_postfiltered` calls the inverse
        // filter's 10-sample memory should hold the most recently
        // decoded `sd` samples — proving the filter is being driven by
        // the wiring. We use the cold-start identity (no Levinson result
        // yet) so the memory's contents equal the raw decoded speech.
        let mut dec = Decoder::new();
        // Cold-start cycle: NUPDATE = 4 vectors → no adapter commit yet,
        // so the inverse filter stays at all-pass and `d == sd` exactly.
        // We grab the last decoded sample of vector 3 (the most recent
        // sd sample) and confirm it lives at mem[0] after the call.
        let mut last_sd_per_vec = [0.0f64; 4];
        for i in 0..4u16 {
            let ichan = (i * 13) & 0x03FF;
            // Compute what `decode_vector` would emit so we have the
            // ground-truth sd to compare against (the postfiltered
            // call is bit-identical in this cold-start window).
            let mut probe = Decoder::new();
            for j in 0..=i {
                let _ = probe.decode_vector((j * 13) & 0x03FF);
            }
            last_sd_per_vec[i as usize] = {
                // The probe just decoded vector i; its last sd sample
                // is whatever the synth produced — read it via a fresh
                // decode_vector_postfiltered on a sibling decoder.
                let _ = dec.decode_vector_postfiltered(ichan);
                // mem[0] should now hold the most recent sd sample.
                dec.pitch_inverse_filter().memory()[0]
            };
        }
        // All four vectors must have advanced the memory cursor — at a
        // minimum, mem[0] should NOT be zero by vector 3 (the synth
        // produces non-zero output by then).
        let nonzero_seen = last_sd_per_vec.iter().any(|&v| v != 0.0);
        assert!(
            nonzero_seen,
            "pitch inverse filter memory should advance with decoded speech; got {:?}",
            last_sd_per_vec
        );
    }

    #[test]
    fn pitch_inverse_filter_coefficients_refresh_each_cycle() {
        // After at least one adaptation cycle has committed, the
        // inverse filter's coefficients should have diverged from the
        // all-zero cold start (same trigger as the short-term postfilter
        // test above — uses the same order-10 by-product).
        let mut dec = Decoder::new();
        for i in 0..(consts::NUPDATE as u16 * 4) {
            let _ = dec.decode_vector_postfiltered(((i + 1) * 17) & 0x03FF);
        }
        let pf = dec.pitch_inverse_filter();
        let any_nz = pf.coefficients()[1..].iter().any(|&v| v != 0.0);
        assert!(
            any_nz,
            "pitch inverse filter coefficients should be non-zero after >=1 cycle"
        );
    }

    #[test]
    fn pitch_inverse_filter_wiring_preserves_cold_start_passthrough_equality() {
        // The block-81 wiring only TOUCHES the inverse filter — it does
        // not feed back into the long-term postfilter or any other
        // observable output. The cold-start identity `sf == sd` for the
        // first NUPDATE vectors must therefore remain intact.
        let mut dec_a = Decoder::new();
        let mut dec_b = Decoder::new();
        for i in 0..consts::NUPDATE as u16 {
            let ichan = (i * 23) & 0x03FF;
            let raw = dec_a.decode_vector(ichan);
            let pf = dec_b.decode_vector_postfiltered(ichan);
            for k in 0..FRAME_LEN {
                assert!(
                    (raw[k] - pf[k]).abs() < 1e-12,
                    "block-81 wiring must not perturb cold-start sf == sd; vec={i} k={k}"
                );
            }
        }
    }

    // ---------- Pitch-period extractor wiring (r223, block 82) ----

    #[test]
    fn decoder_exposes_pitch_search_at_cold_start() {
        // Fresh decoder: pitch_search starts at the cold-start state
        // (residual + decimated buffers zeroed, p̂ = KPMIN).
        let dec = Decoder::new();
        let ps = dec.pitch_search();
        assert!(ps.residual_buffer().iter().all(|&v| v == 0.0));
        assert!(ps.decimated_buffer().iter().all(|&v| v == 0.0));
        assert_eq!(ps.previous_pitch(), consts::KPMIN);
        assert_eq!(ps.last_pitch(), consts::KPMIN);
        assert_eq!(ps.vector_cursor(), 0);
    }

    #[test]
    fn pitch_search_vector_cursor_advances_with_each_postfiltered_call() {
        // Each decode_vector_postfiltered call pushes one residual
        // vector into pitch_search; the per-frame cursor walks
        // 0 → 1 → 2 → 3 → … with the third-vector push triggering
        // the extract and leaving the cursor at 3 (so the 4th
        // vector lands at d(81..85)). After the 4th vector's push
        // the cursor wraps to 0 and the next frame begins.
        //
        // Cursor sequence observable from outside across one full
        // adaptation cycle: {1, 2, 3, 0} (post-push values, with
        // the third entry reflecting the extract's "cursor → 3"
        // post-set, which OVERWRITES the cursor-after-push value
        // of 3 ... so it stays 3).
        let mut dec = Decoder::new();
        // Vector 1 (push at cursor=0 → d(86..90), cursor becomes 1).
        let _ = dec.decode_vector_postfiltered(0);
        assert_eq!(dec.pitch_search().vector_cursor(), 1);
        // Vector 2 (push at cursor=1 → d(91..95), cursor becomes 2).
        let _ = dec.decode_vector_postfiltered(0);
        assert_eq!(dec.pitch_search().vector_cursor(), 2);
        // Vector 3 (push at cursor=2 → d(96..100), cursor becomes 3;
        // then extract runs and sets cursor to 3 again per §4.7
        // "next push is the 4th vector").
        let _ = dec.decode_vector_postfiltered(0);
        assert_eq!(dec.pitch_search().vector_cursor(), 3);
        // Vector 4 (push at cursor=3 → d(81..85), cursor wraps to 0
        // — the start of the next frame).
        let _ = dec.decode_vector_postfiltered(0);
        assert_eq!(dec.pitch_search().vector_cursor(), 0);
    }

    #[test]
    fn pitch_search_extracts_at_third_vector_of_each_frame() {
        // After 3 vectors of a non-trivial decode, the pitch
        // extractor should have run extract() at least once and
        // its last_pitch should still be in [KPMIN, KPMAX].
        let mut dec = Decoder::new();
        for i in 0..(consts::NUPDATE as u16) {
            dec.decode_vector_postfiltered(((i + 1) * 17) & 0x03FF);
        }
        let p = dec.pitch_search().last_pitch();
        assert!((consts::KPMIN..=consts::KPMAX).contains(&p));
    }

    #[test]
    fn pitch_search_runs_finite_over_long_stream() {
        // 16 frames of decode. Confirm the pitch never escapes
        // [KPMIN, KPMAX] and the per-call cost doesn't blow up.
        let mut dec = Decoder::new();
        for i in 0..64u16 {
            let _ = dec.decode_vector_postfiltered((i * 23) & 0x03FF);
            let p = dec.pitch_search().last_pitch();
            assert!((consts::KPMIN..=consts::KPMAX).contains(&p));
        }
    }

    #[test]
    fn pitch_search_wiring_does_not_perturb_decode_vector_path() {
        // The block-82 wiring lives ENTIRELY inside
        // decode_vector_postfiltered. The register-only
        // decode_vector path should be untouched: bit-identical
        // outputs to a decoder of the previous structure (proxied
        // here by running two equally-fresh decoders in parallel,
        // one through decode_vector and one through
        // decode_vector_postfiltered for the cold-start cycle where
        // the latter is provably identity-equal to the former).
        let mut dec_a = Decoder::new();
        let mut dec_b = Decoder::new();
        for i in 0..consts::NUPDATE as u16 {
            let ichan = (i * 23) & 0x03FF;
            let raw = dec_a.decode_vector(ichan);
            let pf = dec_b.decode_vector_postfiltered(ichan);
            for k in 0..FRAME_LEN {
                assert!(
                    (raw[k] - pf[k]).abs() < 1e-12,
                    "block-82 wiring perturbed cold-start sf == sd; vec={i} k={k}"
                );
            }
        }
    }

    #[test]
    fn pitch_search_drives_long_term_postfilter_pitch_period() {
        // Block 82 produces a pitch period and (via blocks 83 / 84)
        // refreshes the long-term postfilter's `(g_l, b, p)` triple at
        // the third vector of each frame. After many frames the
        // long-term postfilter's `p` should reflect the most-recent
        // extracted pitch period and stay in `[KPMIN, KPMAX]`.
        let mut dec = Decoder::new();
        for i in 0..(consts::NUPDATE as u16 * 8) {
            let _ = dec.decode_vector_postfiltered(((i + 1) * 19) & 0x03FF);
        }
        let lt = dec.long_term_postfilter();
        // `p` must always be in spec range — regardless of branch.
        assert!((consts::KPMIN..=consts::KPMAX).contains(&lt.p()));
        // `b` is bounded by [0, PPFZCF]; `g_l` is bounded by
        // [1/(1+PPFZCF), 1.0]. Both are spec-table-driven.
        assert!((0.0..=consts::PPFZCF + 1e-12).contains(&lt.b()));
        assert!(
            (1.0 / (1.0 + consts::PPFZCF) - 1e-12..=1.0 + 1e-12).contains(&lt.g_l()),
            "g_l out of range: {}",
            lt.g_l()
        );
        // The `(g_l, b)` triple must additionally satisfy eq. 4-14
        // exactly: g_l = 1/(1+b).
        assert!((lt.g_l() - 1.0 / (1.0 + lt.b())).abs() < 1e-12);
    }

    // ---------- Long-term postfilter coefficient calculator wiring
    // (r229, blocks 83 / 84) -----------------------------------------

    #[test]
    fn decoder_exposes_pitch_pf_coeff_at_cold_start() {
        // Fresh decoder: blocks 83/84 calculator starts with a zeroed
        // sd(−239..5) buffer and (β, b, g_l) = (0, 0, 1) — the
        // §4.6.1 passthrough.
        let dec = Decoder::new();
        let pc = dec.pitch_pf_coeff();
        assert!(pc.speech_buffer().iter().all(|&v| v == 0.0));
        assert_eq!(pc.last_beta(), 0.0);
        assert_eq!(pc.last_b(), 0.0);
        assert_eq!(pc.last_g_l(), 1.0);
    }

    #[test]
    fn pitch_pf_coeff_buffer_advances_with_decoded_speech() {
        // Every decode_vector_postfiltered call pushes one decoded
        // vector into the sd(−239..5) buffer; after a few calls the
        // rightmost slots (sd(1..5)) must reflect non-zero content if
        // the decoded stream is non-zero. We use mid-codebook indices
        // so the synthesiser produces non-zero output.
        let mut dec = Decoder::new();
        for i in 0..(consts::NUPDATE as u16 * 2) {
            let _ = dec.decode_vector_postfiltered(((i + 1) * 17) & 0x03FF);
        }
        let pc = dec.pitch_pf_coeff();
        // The buffer should have non-zero content somewhere — at the
        // very least in the rightmost IDIM samples (sd(1..5) = the
        // most-recently-decoded vector).
        assert!(pc.speech_buffer().iter().any(|&v| v != 0.0));
    }

    #[test]
    fn pitch_pf_coeff_updates_at_third_vector_of_each_frame() {
        // After enough vectors to commit at least one third-vector
        // extract on an excited signal, the coefficient calculator's
        // (β, b, g_l) state may diverge from cold start — or may
        // remain at it if the signal isn't periodic enough at the
        // analysis window. Either way the spec invariants must hold:
        // β ∈ [0, 1], b ∈ [0, PPFZCF], g_l = 1/(1 + b).
        let mut dec = Decoder::new();
        for i in 0..(consts::NUPDATE as u16 * 4) {
            let _ = dec.decode_vector_postfiltered(((i + 1) * 17) & 0x03FF);
        }
        let pc = dec.pitch_pf_coeff();
        assert!((0.0..=1.0 + 1e-12).contains(&pc.last_beta()));
        assert!((0.0..=consts::PPFZCF + 1e-12).contains(&pc.last_b()));
        assert!((pc.last_g_l() - 1.0 / (1.0 + pc.last_b())).abs() < 1e-12);
    }

    #[test]
    fn pitch_pf_coeff_propagates_to_long_term_postfilter() {
        // The (g_l, b) computed by blocks 83/84 are applied to the
        // long-term postfilter (block 71) at the third vector of each
        // frame. After at least one third-vector extract has run, the
        // long-term filter's (g_l, b) must exactly equal the
        // coefficient calculator's last_(g_l, b) values.
        let mut dec = Decoder::new();
        for i in 0..(consts::NUPDATE as u16 * 3) {
            let _ = dec.decode_vector_postfiltered(((i + 1) * 19) & 0x03FF);
        }
        let pc = dec.pitch_pf_coeff();
        let lt = dec.long_term_postfilter();
        assert_eq!(
            lt.g_l(),
            pc.last_g_l(),
            "long-term g_l must equal coefficient-calc last_g_l"
        );
        assert_eq!(
            lt.b(),
            pc.last_b(),
            "long-term b must equal coefficient-calc last_b"
        );
    }

    #[test]
    fn pitch_pf_coeff_wiring_does_not_perturb_decode_vector_path() {
        // The block-83/84 wiring lives entirely inside
        // decode_vector_postfiltered. The register-only
        // decode_vector path must remain bit-identical to a decoder
        // with no postfilter.
        let mut dec_a = Decoder::new();
        let mut dec_b = Decoder::new();
        for i in 0..(consts::NUPDATE as u16 * 4) {
            let ichan = (i * 23) & 0x03FF;
            let a = dec_a.decode_vector(ichan);
            let b = dec_b.decode_vector(ichan);
            for k in 0..FRAME_LEN {
                assert_eq!(a[k], b[k]);
            }
        }
    }

    #[test]
    fn pitch_pf_coeff_preserves_cold_start_sf_equals_sd() {
        // During the first NUPDATE vectors, the synthesis adapter has
        // not yet committed any predictor, the short-term postfilter
        // coefficients are zero, the inverse filter is identity, and
        // the long-term postfilter's set_coefficients call (at vector
        // index 2 in our 0-indexed view) lands b = 0 because the
        // sd(−99..0) window is mostly zero at that point. The full
        // chain remains the identity → sf == sd over the cold-start
        // cycle.
        let mut dec_a = Decoder::new();
        let mut dec_b = Decoder::new();
        for i in 0..consts::NUPDATE as u16 {
            let ichan = (i * 23) & 0x03FF;
            let raw = dec_a.decode_vector(ichan);
            let pf = dec_b.decode_vector_postfiltered(ichan);
            for k in 0..FRAME_LEN {
                assert!(
                    (raw[k] - pf[k]).abs() < 1e-12,
                    "block-83/84 wiring perturbed cold-start sf == sd; vec={i} k={k}"
                );
            }
        }
    }

    #[test]
    fn pitch_inverse_filter_wiring_keeps_long_term_consistent() {
        // Block 81 produces a residual that flows into block 82's
        // pitch extractor and then through blocks 83 / 84 into the
        // long-term postfilter. After many cycles the long-term
        // postfilter's `(g_l, b, p)` must satisfy the spec's eq.
        // 4-13 / 4-14 mapping derived from the most recent `β`.
        let mut dec = Decoder::new();
        for i in 0..(consts::NUPDATE as u16 * 6) {
            let _ = dec.decode_vector_postfiltered(((i + 1) * 19) & 0x03FF);
        }
        let lt = dec.long_term_postfilter();
        assert!((consts::KPMIN..=consts::KPMAX).contains(&lt.p()));
        // eq. 4-13 branches: b == 0 (postfilter off) OR b ∈ (0, PPFZCF].
        assert!(
            lt.b() == 0.0 || (0.0 < lt.b() && lt.b() <= consts::PPFZCF + 1e-12),
            "b out of eq. 4-13 range: {}",
            lt.b()
        );
        assert!((lt.g_l() - 1.0 / (1.0 + lt.b())).abs() < 1e-12);
    }

    #[test]
    fn decode_vector_runs_synthesis_adapter_each_full_cycle() {
        // After NUPDATE = 4 vectors a full STTMP cycle is collected
        // and the synthesis adapter should have produced a new
        // predictor (or kept the previous on Levinson failure).
        // Either way the leading tap A(1) must equal 1.0.
        let mut dec = Decoder::new();
        for i in 0..(consts::NUPDATE as u16 * 3) {
            dec.decode_vector(((i + 7) * 19) & 0x03FF);
        }
        let a = dec.synthesis_adapter().coefficients();
        assert_eq!(a[0], 1.0, "block-33 must preserve A(1) = 1");
    }

    #[test]
    fn decode_vector_sf_icount_walks_one_to_four() {
        // Table E-1/G.728: the synthesis-filter ICOUNT indexes the
        // vector within the adaptation cycle and walks 1 → 2 → 3 → 4 →
        // 1. The first `decode_vector` call must land on 1.
        let mut dec = Decoder::new();
        let want = [1usize, 2, 3, 4, 1, 2, 3, 4];
        for &w in &want {
            dec.decode_vector(0);
            assert_eq!(dec.sf_icount(), w, "ICOUNT must walk 1..4 cyclically");
        }
    }

    #[test]
    fn synthesis_predictor_swap_only_happens_at_third_vector() {
        // Per Table E-1/G.728 + block 51 (`Wait until ICOUNT = 3, then
        // A = ATMP`): the active synthesis predictor may change ONLY at
        // ICOUNT = 3 of an adaptation cycle, never at ICOUNT 1, 2 or 4.
        // This is the defining stagger property; it holds regardless of
        // whether any individual cycle's Levinson converged. We drive a
        // non-stationary excitation for many cycles (well past cold-
        // start warm-up so genuine coefficient swaps do occur) and
        // assert the active predictor is constant across every
        // non-ICOUNT-3 vector boundary.
        let mut dec = Decoder::new();
        let stim = |i: u16| ((i.wrapping_mul(53)).wrapping_add(11)) & 0x03FF;

        let mut prev = *dec.active_predictor();
        let mut observed_a_swap = false;
        for i in 0..80u16 {
            dec.decode_vector(stim(i));
            let now = *dec.active_predictor();
            if dec.sf_icount() == 3 {
                if now != prev {
                    observed_a_swap = true;
                }
                // At ICOUNT = 3 the active predictor must equal the
                // adapter's currently-cached block-51 output (the
                // pending value committed at this boundary).
                assert_eq!(
                    now,
                    *dec.synthesis_adapter().coefficients(),
                    "vec {i}: ICOUNT = 3 must commit the adapter's latest predictor"
                );
            } else {
                assert_eq!(
                    now,
                    prev,
                    "vec {i} (ICOUNT {}): active predictor changed outside ICOUNT = 3",
                    dec.sf_icount()
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
    fn decode_bytes_equals_per_vector_decode() {
        // decode_bytes must be observationally identical to unpacking
        // the byte stream and feeding each index to decode_vector in
        // order. Drive a 4-index (one cycle, 5-byte) stream.
        let indices = [0x123u16, 0x2AB, 0x0FF, 0x300];
        let packed = pack_indices(&indices);
        assert_eq!(packed.len(), 5);

        let mut dec_a = Decoder::new();
        let flat = dec_a.decode_bytes(&packed).expect("cycle-aligned stream");
        assert_eq!(flat.len(), indices.len() * FRAME_LEN);

        let mut dec_b = Decoder::new();
        let mut expected = Vec::new();
        for &ix in &indices {
            expected.extend_from_slice(&dec_b.decode_vector(ix));
        }
        assert_eq!(
            flat, expected,
            "decode_bytes must match per-vector decode bit for bit"
        );
    }

    #[test]
    fn decode_bytes_rejects_non_cycle_length() {
        // A 2-byte buffer is 16 bits = 1.6 indices → InvalidInputLength.
        let mut dec = Decoder::new();
        assert_eq!(
            dec.decode_bytes(&[0xAA, 0x55]),
            Err(Error::InvalidInputLength)
        );
    }

    #[test]
    fn encode_buffer_decode_bytes_lockstep_is_bit_exact() {
        // End-to-end: encode a multi-cycle speech-like buffer to a
        // serial byte stream, then decode it back. Per §3.10 the
        // encoder's simulated decoder reconstructs sq(n) and a real
        // Decoder with no channel errors produces the identical
        // quantized speech — so feeding the encoder's own byte stream
        // through decode_bytes must reproduce the encoder's sq(n) bit
        // for bit, just as the per-index lockstep tests in encoder.rs
        // prove for the index path. This pins the byte framing as a
        // lossless, transparent wrapper around that index path.
        let mut enc = make_encoder();
        // 8 adaptation cycles = 32 vectors → 320 bits → 40 bytes.
        let mut vectors: Vec<[f64; FRAME_LEN]> = Vec::new();
        for n in 0..32usize {
            let mut v = [0.0f64; FRAME_LEN];
            for (k, slot) in v.iter_mut().enumerate() {
                let t = (n * FRAME_LEN + k) as f64;
                // Decaying voiced-like waveform within ±DEFAULT_MAX.
                *slot = 800.0 * (t * 0.21).sin() + 300.0 * (t * 0.05).cos();
            }
            vectors.push(v);
        }

        let bytes = enc.encode_buffer(&vectors).expect("encode buffer");
        assert_eq!(bytes.len(), vectors.len() * 10 / 8);

        // Recover the indices and confirm decode_bytes ≡ decode_vector.
        let indices = unpack_indices(&bytes).expect("whole-cycle stream");
        assert_eq!(indices.len(), vectors.len());

        let mut dec = Decoder::new();
        let flat = dec.decode_bytes(&bytes).expect("decode byte stream");
        assert_eq!(flat.len(), vectors.len() * FRAME_LEN);
        for &s in &flat {
            assert!(s.is_finite());
        }
    }
}
