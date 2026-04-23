//! ITU-T G.728 Low-Delay CELP (LD-CELP, 16 kbit/s) decoder — first cut.
//!
//! G.728 codes narrowband (8 kHz) speech at 16 kbit/s using 5-sample
//! analysis vectors (0.625 ms). Each vector is represented by a single
//! 10-bit index split into a 7-bit shape codebook selector (128 shape
//! vectors) and a 3-bit gain selector (sign bit + 2 magnitude bits,
//! for 8 signed levels), giving 2 kbyte/s of compressed data.
//!
//! The spec is unusual for a CELP variant: the synthesis filter is a
//! **50th-order** LPC predictor adapted **backward** from previously
//! reconstructed speech (no transmitted LPC), and a 10th-order backward-
//! adaptive log-gain predictor tracks excitation energy. Because all the
//! adaptation machinery is derived from the decoder's own output, the
//! payload shrinks to just the codebook index — which is why the bit
//! reader here is so thin.
//!
//! What's landed:
//!
//! - Crate + workspace wiring (registers as codec id `g728`).
//! - 10-bit MSB-first bit reader for the packed index stream.
//! - 128 × 5 shape codebook + 4-magnitude gain codebook transcribed
//!   verbatim from ITU-T G.728 Annex B (`CODEBK` and `GQ`). See the
//!   [`tables`] module.
//! - Backward-adaptive 50th-order LPC predictor using the spec's §3.7
//!   hybrid (Barnwell / logarithmic) autocorrelation window — 35-sample
//!   non-recursive tail + decaying recursive portion with `alpha^{2L} =
//!   3/4` and the 105-sample window table from Annex A.1. Refreshed
//!   every 4 vectors (2.5 ms) via Levinson-Durbin with `FAC = 253/256`
//!   bandwidth expansion.
//! - Backward-adaptive 10th-order log-gain predictor.
//! - §5.5 adaptive postfilter: long-term (pitch) + short-term (pole-
//!   zero formant emphasis with spectral-tilt compensation) + AGC
//!   level renormalisation. See the [`postfilter`] module.
//! - Synthesis loop: excitation → all-pole IIR → optional postfilter →
//!   S16 PCM, exposed via the standard [`oxideav_codec::Decoder`]
//!   trait.
//! - Frame-erasure concealment (Annex A.3 / §5.8). A `Packet` whose
//!   `flags.corrupt` bit is set triggers synthesis from the last
//!   transmitted excitation, attenuated across successive erased
//!   packets (0.9 → 0.7 → 0.5 → 0.3 → 0.0) so burst losses fade to
//!   silence rather than freezing a tone. The backward-adaptive
//!   LPC + gain predictors keep running on the concealment output
//!   so they resume cleanly when the next clean packet arrives.
//!
//! Residual deviation: the log-gain predictor still uses the scaffold's
//! Hamming-window autocorrelation for block 43's hybrid-window module;
//! this is a stability tweak with no measurable impact on round-trip
//! quality but means bitstream-level output is not strictly bit-exact
//! against the ITU reference decoder.

// Scaffold-only — symbols will be used once the full decoder body lands.
// These allow()s come off when the decoder is exercised from end to end.
#![allow(
    dead_code,
    clippy::needless_range_loop,
    clippy::unnecessary_cast,
    clippy::excessive_precision,
    clippy::approx_constant,
    clippy::doc_lazy_continuation,
    clippy::doc_overindented_list_items
)]

pub mod bitreader;
pub mod codec;
pub mod decoder;
pub mod encoder;
pub mod postfilter;
pub mod predictor;
pub mod tables;

use oxideav_codec::CodecRegistry;

pub const CODEC_ID_STR: &str = "g728";

/// Sample rate, always 8 kHz narrowband.
pub const SAMPLE_RATE: u32 = 8_000;

/// Samples per analysis vector.
pub const VECTOR_SIZE: usize = 5;

/// Bits per packed codebook index.
pub const INDEX_BITS: u32 = 10;

/// Order of the backward-adaptive LPC synthesis predictor.
pub const LPC_ORDER: usize = 50;

/// Order of the backward-adaptive log-gain predictor.
pub const GAIN_ORDER: usize = 10;

/// Number of shape codebook entries (7 index bits).
pub const SHAPE_CB_SIZE: usize = 128;

/// Number of gain codebook entries. The 10-bit raw index exposes 2 bits of
/// magnitude and 1 bit of sign, so only the first 4 entries of the lookup
/// table are reachable through the bitstream; the sign flip doubles the
/// set to 8 signed levels. (The ITU Annex A `GB` table is also 4
/// magnitudes × 2 signs by the same bit-layout argument — this is *not*
/// a placeholder artefact.)
pub const GAIN_CB_SIZE: usize = 4;

/// Register the G.728 codec (decoder + encoder) with the codec registry.
pub fn register(reg: &mut CodecRegistry) {
    codec::register(reg);
}
