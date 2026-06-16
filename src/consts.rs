//! Basic coder parameters for ITU-T G.728 LD-CELP.
//!
//! Transcribed from Table 1/G.728 (Recommendation 09/92, p. 25). The
//! third column ("Value") of that table is the source for every named
//! constant defined here.

/// AGC adaptation speed controlling factor (postfilter; Table 1/G.728).
pub const AGCFAC: f64 = 0.99;

/// Bandwidth expansion factor of the 50th-order synthesis filter
/// (`λ = 253/256`, equation 3-7 / Table 1/G.728).
pub const FAC: f64 = 253.0 / 256.0;

/// Bandwidth expansion factor of the 10th-order log-gain predictor
/// (`λ_g = 29/32`, equation 3-12 / Table 1/G.728).
pub const FACGP: f64 = 29.0 / 32.0;

/// Reciprocal of vector dimension; 1/IDIM (Table 1/G.728).
pub const DIMINV: f64 = 0.2;

/// Vector dimension — five samples per excitation block (Table 1/G.728).
pub const IDIM: usize = 5;

/// Log-gain offset value, in dB (Table 1/G.728; §3.8).
pub const GOFF: f64 = 32.0;

/// Allowed deviation from previous pitch period (Table 1/G.728).
pub const KPDELTA: usize = 6;

/// Minimum pitch period, in samples (Table 1/G.728).
pub const KPMIN: usize = 20;

/// Maximum pitch period, in samples (Table 1/G.728).
pub const KPMAX: usize = 140;

/// 50th-order synthesis-filter LPC order (Table 1/G.728).
pub const LPC: usize = 50;

/// 10th-order log-gain predictor order (Table 1/G.728).
pub const LPCLG: usize = 10;

/// 10th-order perceptual-weighting-filter order (Table 1/G.728).
pub const LPCW: usize = 10;

/// Shape codebook size; 128 codevectors (Table 1/G.728).
pub const NCWD: usize = 128;

/// Frame (adaptation-cycle) size; 20 samples = 4 × IDIM (Table 1/G.728).
pub const NFRSZ: usize = 20;

/// Gain codebook size; 8 levels (Table 1/G.728).
pub const NG: usize = 8;

/// Non-recursive window length for the synthesis-filter hybrid window
/// (Table 1/G.728).
pub const NONR: usize = 35;

/// Non-recursive window length for the log-gain hybrid window
/// (Table 1/G.728).
pub const NONRLG: usize = 20;

/// Non-recursive window length for the weighting-filter hybrid window
/// (Table 1/G.728).
pub const NONRW: usize = 30;

/// Pitch analysis window size, in samples (Table 1/G.728).
pub const NPWSZ: usize = 100;

/// Predictor update period, expressed in vectors per adaptation cycle
/// (Table 1/G.728).
pub const NUPDATE: usize = 4;

/// Tap threshold below which the long-term postfilter is disabled
/// (Table 1/G.728).
pub const PPFTH: f64 = 0.6;

/// Pitch-postfilter zero controlling factor (Table 1/G.728).
pub const PPFZCF: f64 = 0.15;

/// Short-term postfilter pole controlling factor (Table 1/G.728).
pub const SPFPCF: f64 = 0.75;

/// Short-term postfilter zero controlling factor (Table 1/G.728).
pub const SPFZCF: f64 = 0.65;

/// Tap threshold for fundamental pitch replacement (Table 1/G.728).
pub const TAPTH: f64 = 0.4;

/// Spectral tilt-compensation controlling factor (Table 1/G.728).
pub const TILTF: f64 = 0.15;

/// White-noise correction factor; `257/256` ≈ 24 dB below average
/// speech power (equation 3-1i / Table 1/G.728).
pub const WNCF: f64 = 257.0 / 256.0;

/// Perceptual-weighting filter pole controlling factor `γ₂`
/// (Table 1/G.728).
pub const WPCF: f64 = 0.6;

/// Perceptual-weighting filter zero controlling factor `γ₁`
/// (Table 1/G.728).
pub const WZCF: f64 = 0.9;

/// Frame-erasure LPC-softening bandwidth-expansion factor (Annex I,
/// §I.4.2). During erased frames the last good 50th-order LPC predictor
/// is "softened" by bandwidth expansion using this factor instead of
/// the normal `FAC = 253/256 ≈ 0.9883` of block 51; the spec spells it
/// out as "the bandwidth expansion factor FAC (5.1/G.728) is 0.97
/// rather than 253/256".
pub const FEFAC: f64 = 0.97;

/// Frame-erasure LPC re-softening cadence, in adaptation cycles per
/// 10 ms (Annex I, §I.4.2). The softened predictor is further expanded
/// by another factor of `FEFAC` every additional 10 ms of erasure; one
/// adaptation cycle is `NFRSZ = 20` samples = 2.5 ms, so 10 ms is four
/// adaptation cycles ("the LPC coefficients are updated again at the
/// third vector in the 5th adaptation cycle").
pub const FE_LPC_CYCLES_PER_STEP: usize = 4;

/// Frame-erasure gain-growth limit, in dB per 5-sample vector (Annex I,
/// §I.4.5 / Table I.1/G.728, `FEGAINMAX`). After a frame erasure the
/// backward-adapted log-gain is not allowed to grow by more than this
/// many dB per vector relative to the last log-gain, for the first few
/// good frames; the spec spells it out as "we limit the gain growth to
/// no more than `FEGAINMAX = +2 dB` per 5-sample vector".
pub const FEGAINMAX: f64 = 2.0;

/// Maximum gain-clamping duration after a frame erasure, in 2.5 ms
/// adaptation cycles (Annex I, §I.4.5 / Table I.1/G.728, `AFTERFEMAX`).
/// The gain-growth clamp lasts for the length of the erasure but never
/// longer than `AFTERFEMAX` cycles (= 40 ms): "If the frame erasure is
/// longer than 40 ms long, the gain clamping is limited to the first
/// 40 ms after the erasure."
pub const AFTERFEMAX: usize = 16;

/// Initial / reset value of the "last predicted gain in dB" tracker
/// `OGAINDB` (Annex I, §I.5.1 — `OGAINDB = -32 | Last predicted gain in
/// dB`). The gain-growth limiter clamps the new log-gain to at most
/// `OGAINDB + FEGAINMAX`; before any gain has been predicted the
/// reference value is the log-gain floor (`GAIN = 0` maps to
/// `0 - GOFF = -32` dB after the §3.8 offset subtraction).
pub const OGAINDB_INIT: f64 = -32.0;

/// Frame-erasure voicing threshold for excitation extrapolation (Annex
/// I, §I.4.1). At the first erased frame the decoder reuses the last
/// good single-tap pitch-predictor weight `PTAP` (block 83) to decide
/// voiced vs unvoiced extrapolation, but against a *lowered* threshold
/// `VTH = PPFTH / 1.4` (rather than the normal `PPFTH = 0.6` used to
/// gate the long-term postfilter): "we use a voicing threshold,
/// `VTH = PPFTH/1.4`, that is lower than `PPFTH`".
pub const VTH: f64 = PPFTH / 1.4;

/// Frame-erasure voiced excitation-extrapolation scaling factors
/// (Annex I, §I.5.2 block 31SF, `VOICEDFEGAIN`). Indexed by
/// `N10MSEC = FECOUNT >> 2` (the number of completed 10 ms intervals of
/// the erasure). The floating-point pseudo-code spells out
/// `VOICEDFEGAIN(0:4) = .8, .8, .6, .4, .2`; beyond index 4 the scale
/// `FESCALE` is forced to 0 (no extrapolation past 50 ms) to avoid
/// "long, unnatural, periodic artifacts".
pub const VOICED_FE_GAIN: [f64; 5] = [0.8, 0.8, 0.6, 0.4, 0.2];

/// Frame-erasure unvoiced excitation-extrapolation attenuation factors
/// (Annex I, §I.5.2 block 31SF, `UNVOICEDFEGAIN`). Indexed by
/// `N10MSEC = FECOUNT >> 2`. The floating-point pseudo-code spells out
/// `UNVOICEDFEGAIN(0:5) = 1., 1., .8, .6, .4, .2`; beyond index 5 the
/// scale `FESCALE` is forced to 0 (no extrapolation past 60 ms). These
/// multiply the magnitude-match scale `AVMAG` so the extrapolated
/// excitation level decays as the erasure lengthens.
pub const UNVOICED_FE_GAIN: [f64; 6] = [1.0, 1.0, 0.8, 0.6, 0.4, 0.2];

/// Length of the `ETPAST()` look-back window for the unvoiced random
/// slide-back (Annex I, §I.4.1 / block 31FE). The unvoiced extrapolator
/// randomly slides the 5-sample window back "between 5 and 140 samples"
/// into `ETPAST()`; `KPMAX = 140` is the maximum slide-back and equals
/// the maximum pitch period, so the look-back buffer must hold at least
/// `KPMAX` past samples (plus the current `IDIM` vector).
pub const FE_UNVOICED_MAX_SLIDE: usize = KPMAX;

/// Minimum slide-back for the unvoiced random extrapolator (Annex I,
/// §I.4.1 / block 31FE): "a random integer between 5 and 140". Five is
/// the vector dimension `IDIM`, so the random segment never overlaps
/// the vector currently being written.
pub const FE_UNVOICED_MIN_SLIDE: usize = IDIM;

/// Number of past `ET()` samples used to compute the unvoiced
/// magnitude-match reference `AVMAG` (Annex I, §I.5.2 block 31SF). The
/// sum is taken over "the last 40 samples of `ETPAST()`" (the last
/// 5 ms before the erasure) and then divided by 8 — the spec keeps the
/// `/8` rather than `/40` because only the *ratio* of the 5-sample
/// segment magnitude to this reference matters in block 31FE.
pub const FE_AVMAG_SAMPLES: usize = 40;

/// Total number of vectors per encoded second at 16 kbit/s.
///
/// G.728 emits one 10-bit codebook index per IDIM-sample vector; 8 kHz
/// sampling × IDIM = 5 yields 1600 vectors per second = 16 000 bit/s.
pub const VECTORS_PER_SECOND: usize = 8_000 / IDIM;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diminv_matches_idim() {
        // Cross-check: Table 1 lists both IDIM = 5 and DIMINV = 0.2
        // independently; they must agree (sanity guard against a
        // future typo in either constant).
        assert!((DIMINV - 1.0 / IDIM as f64).abs() < 1e-12);
    }

    #[test]
    fn bitrate_matches_16kbits() {
        // Sanity: G.728 is a 16 kbit/s codec. Each vector is IDIM = 5
        // samples at 8 kHz with a 10-bit codebook index, giving exactly
        // 16 000 bit/s when both sides agree on IDIM.
        let bits_per_second = VECTORS_PER_SECOND * 10;
        assert_eq!(bits_per_second, 16_000);
    }

    #[test]
    fn nfrsz_is_four_vectors() {
        // Table 1 lists NFRSZ = 20 and NUPDATE = 4 independently; they
        // must satisfy NFRSZ = NUPDATE × IDIM.
        assert_eq!(NFRSZ, NUPDATE * IDIM);
    }

    #[test]
    fn fac_matches_paragraph_value() {
        // Spec §3.7 spells out `λ = 253/256 = 0.98828125`; round-trip
        // the rational value through f64 to guard the literal.
        assert!((FAC - 0.98828125).abs() < 1e-15);
    }

    #[test]
    fn facgp_matches_paragraph_value() {
        // Spec §3.8 spells out `(29/32)¹ = 0.90625` for i = 1; round
        // trip the rational the same way.
        assert!((FACGP - 0.90625).abs() < 1e-15);
    }
}
