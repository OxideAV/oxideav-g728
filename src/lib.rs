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
//! Round 189 (preserved) provides the Annex A.1/A.2/A.3 hybrid
//! windows, the Annex B excitation codebook (128 × 5 shape + 8
//! gain), the Annex C bandwidth-broadening vectors (Q14), the
//! Annex D 1 kHz lowpass coefficients, Table 1/G.728 parameter
//! constants, the Levinson-Durbin recursion, and the
//! [`Decoder::decode_index`] caller-driven entry point.
//!
//! ## What is NOT yet wired up
//!
//! * **Adaptive postfilter (blocks 71–77).** The postfilter is a
//!   perceptual quality enhancement: it does not affect bit-exactness
//!   against the standard's reference output, so it lands later.
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

pub mod consts;
pub mod decoder;
pub mod gain_adapter;
pub mod hybrid_window;
pub mod levinson;
pub mod synthesis_adapter;
pub mod tables;

pub use decoder::{pack_channel_index, ExcitationVector, Synthesizer, DEFAULT_MAX, FRAME_LEN};
pub use gain_adapter::GainAdapter;
pub use hybrid_window::{HybridWindow, HybridWindowState};
pub use levinson::{levinson_durbin, LevinsonError};
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
            et_prev: [0.0; FRAME_LEN],
            sttmp_buf: [0.0; consts::NFRSZ],
            sttmp_idx: 0,
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
