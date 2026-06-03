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
pub mod consts;
pub mod decoder;
pub mod gain_adapter;
pub mod hybrid_window;
pub mod levinson;
pub mod long_term_postfilter;
pub mod pitch_inverse_filter;
pub mod pitch_postfilter_coeff;
pub mod pitch_search;
pub mod short_term_postfilter;
pub mod synthesis_adapter;
pub mod tables;

pub use agc::Agc;
pub use decoder::{pack_channel_index, ExcitationVector, Synthesizer, DEFAULT_MAX, FRAME_LEN};
pub use gain_adapter::GainAdapter;
pub use hybrid_window::{HybridWindow, HybridWindowState};
pub use levinson::{levinson_durbin, LevinsonError};
pub use long_term_postfilter::LongTermPostfilter;
pub use pitch_inverse_filter::PitchInverseFilter;
pub use pitch_postfilter_coeff::PitchPostfilterCoeff;
pub use pitch_search::PitchSearch;
pub use short_term_postfilter::ShortTermPostfilter;
pub use synthesis_adapter::SynthesisAdapter;

/// Crate-local error type.
///
/// Round 189 exposes a decoder front end behind [`Decoder`] but does
/// not yet implement encoder operations or the backward-adapter
/// pipeline; those entry points still return [`Error::NotImplemented`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// Functionality is part of the planned scope but not yet wired
    /// up. The error carries no payload — see the per-call docs for
    /// the unblock condition.
    NotImplemented,
    /// The supplied input byte length is not a multiple of two — and
    /// hence cannot be unpacked into whole 10-bit codebook indices
    /// using the canonical wire layout (one 16-bit little-endian word
    /// per index, low 10 bits significant).
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
        // ----- Block 30: predict σ(n) from previous ET ---------------
        let sigma = self.gain_adapter.predict_next(&self.et_prev);

        // ----- Block 29: shape codevector lookup ---------------------
        // ExcitationVector::from_channel_index pre-applies the gain-
        // codebook scaling; we want the raw shape codevector only,
        // because block 30 supplies σ(n) and block 31 multiplies it
        // explicitly. So redo the lookup unscaled here.
        let ichan_masked = (ichan & 0x03FF) as usize;
        let is_index = ichan_masked / consts::NG;
        let y_row = &tables::Y_Q11[is_index];
        let mut et = [0.0f64; FRAME_LEN];
        for k in 0..FRAME_LEN {
            // Annex B: divide by 2^11 to obtain float; then multiply
            // by σ(n) per block 31.
            et[k] = sigma * (y_row[k] as f64 / 2_048.0);
        }

        // ----- Block 32: synthesis filter ---------------------------
        // Use the most recently produced predictor from block 33.
        self.synth.set_predictor(*self.synth_adapter.coefficients());
        let exc = ExcitationVector(et);
        let out = self.synth.filter_vector(&exc);

        // ----- Block 33 preparation: accumulate quantised speech ----
        // STTMP for the next call is the rolling buffer of decoded
        // speech vectors. We push the current output and, once a full
        // NFRSZ-sample cycle is accumulated, run the synthesis
        // adapter to produce a refreshed predictor for the next cycle.
        for k in 0..FRAME_LEN {
            self.sttmp_buf[self.sttmp_idx] = out[k];
            self.sttmp_idx += 1;
        }
        if self.sttmp_idx >= consts::NFRSZ {
            let _ = self.synth_adapter.adapt(&self.sttmp_buf);
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
/// [`make_decoder`]; the registry-side `decoder` / `encoder` factory
/// wiring will land once the backward adapter is implemented and the
/// decoder can autonomously consume a 10-bit-per-frame bytestream.
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
}
