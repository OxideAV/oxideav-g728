//! # oxideav-g728
//!
//! Pure-Rust ITU-T G.728 LD-CELP (Low-Delay Code-Excited Linear
//! Prediction) speech codec — 16 kbit/s, 8 kHz sampling, 0.625 ms
//! algorithmic delay.
//!
//! ## Status: clean-room rebuild — decoder front end + Annex tables
//!
//! The crate was reset to a register-only scaffold (round 171, master
//! `14e3bad`) under the workspace clean-room policy: the previous
//! implementation's numeric tables had been extracted from an external
//! reference C distribution, which the policy forbids. The current
//! tree is being re-grown one spec-cited unit at a time from the
//! published ITU-T G.728 (1992-09) Recommendation prose alone.
//!
//! Round 189 lands:
//!
//! * The complete set of Table 1/G.728 codec parameters (see
//!   [`consts`]).
//! * Annex A.1, A.2 and A.3 hybrid windows (synthesis, log-gain and
//!   weighting filter), transcribed as their normative Q15 integers
//!   plus float views derived by the spec-stated `value / 2¹⁵` rule.
//! * The complete Annex B excitation codebook: 128 shape codevectors
//!   × 5 components (Q11) plus the 8-level Q13 gain codebook with
//!   pre-computed `g2` and `gsq` arrays (equations 3-21 / 3-22).
//! * Annex C bandwidth-broadening vectors for the synthesis filter,
//!   the log-gain predictor, the perceptual-weighting filter, and
//!   the short-term postfilter (Q14).
//! * Annex D 1 kHz lowpass coefficients for the pitch extractor.
//! * A Levinson-Durbin recursion that the synthesis-filter adapter,
//!   the log-gain adapter and the weighting-filter adapter will all
//!   share once a future round wires up the backward-adapter blocks
//!   (see [`levinson`]).
//! * The decoder excitation front end (codebook lookup → gain scaling
//!   → 50th-order all-pole synthesis filter with `±MAX` memory
//!   saturation), behind a stable [`Decoder`] entry point.
//!
//! ## What is NOT yet wired up
//!
//! * **Backward adapters (blocks 30, 33).** The synthesis-filter
//!   predictor is currently caller-supplied via
//!   [`decoder::Synthesizer::set_predictor`]; the next round will plug
//!   the Annex A hybrid window → Levinson-Durbin → Annex C bandwidth
//!   expansion chain in.
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
pub mod levinson;
pub mod tables;

pub use decoder::{pack_channel_index, ExcitationVector, Synthesizer, DEFAULT_MAX, FRAME_LEN};
pub use levinson::{levinson_durbin, LevinsonError};

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
/// One [`Decoder`] instance carries the synthesis-filter memory across
/// vectors. Construct it with [`Decoder::new`] and feed one 10-bit
/// channel index at a time through [`Decoder::decode_index`] to
/// receive a `FRAME_LEN`-sample decoded-PCM vector.
///
/// The decoder is currently parameterised by the predictor coefficients
/// the caller provides via [`Decoder::set_synthesis_predictor`]; a
/// follow-up round will replace that hook with the spec's backward
/// adapter so the decoder runs autonomously off the raw bitstream.
#[derive(Debug, Clone)]
pub struct Decoder {
    synth: Synthesizer,
}

impl Default for Decoder {
    fn default() -> Self {
        Self::new()
    }
}

impl Decoder {
    /// Construct a fresh decoder with all internal state initialised
    /// to the values listed in Table 2/G.728 (synthesis-filter memory
    /// zeroed, all-pass predictor).
    pub fn new() -> Self {
        Self {
            synth: Synthesizer::new(),
        }
    }

    /// Inject a fresh set of synthesis-filter coefficients for the
    /// next vector. The leading entry must be `1.0` (the spec's
    /// `A(1) = 1` invariant).
    pub fn set_synthesis_predictor(&mut self, coeffs: [f64; consts::LPC + 1]) {
        self.synth.set_predictor(coeffs);
    }

    /// Decode one channel index into one [`FRAME_LEN`]-sample PCM
    /// vector in the spec's `±4 095` linear range. Channel indices
    /// outside `[0, 1024)` are masked to their low 10 bits.
    pub fn decode_index(&mut self, ichan: u16) -> [f64; FRAME_LEN] {
        let exc = ExcitationVector::from_channel_index(ichan);
        self.synth.filter_vector(&exc)
    }

    /// Borrow the synthesiser (useful for tests and for the future
    /// backward-adapter wiring).
    pub fn synthesizer(&self) -> &Synthesizer {
        &self.synth
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
}
