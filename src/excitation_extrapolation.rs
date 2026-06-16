//! Annex I §I.4.1 — frame-erasure excitation extrapolation.
//!
//! During an erased frame the G.728 decoder has no channel index, so it
//! cannot look up an excitation codevector (block 29) or backward-adapt
//! the gain (block 20). Annex I instead **extrapolates the gain-scaled
//! excitation** `ET()` — the output of the gain-scaling unit (block 31),
//! *not* the raw codebook output — directly from a history buffer
//! `ETPAST()` (§I.4.1):
//!
//! > We do not extrapolate the output of the excitation VQ codebook
//! > (Block 29). Instead, we extrapolate the output of the gain scaling
//! > unit (Block 31), or the gain-scaled excitation, which is denoted by
//! > the `ET()` array. […] we need to maintain an additional array
//! > `ETPAST()` which stores the present and past `ET()` vectors.
//!
//! The extrapolation has two regimes, chosen once at erasure onset from
//! the last good single-tap pitch weight `PTAP` (block 83) against the
//! lowered voicing threshold [`VTH`] `= PPFTH/1.4` (§I.4.1):
//!
//! - **Voiced** (`PTAP > VTH`): periodically repeat a scaled-down copy
//!   of the last `KP` (= `FEDELAY`) samples of `ETPAST()`. Each output
//!   sample is `ET(i) = FESCALE · ETPAST(i − FEDELAY)`. The scale
//!   `FESCALE` starts at `0.8` (first 20 ms) and steps down per
//!   [`VOICED_FE_GAIN`], reaching `0` past 50 ms.
//! - **Unvoiced** (`PTAP ≤ VTH`): randomly slide a 5-sample window back
//!   `5..=140` samples into `ETPAST()` and magnitude-match the segment
//!   to the average excitation magnitude `AVMAG` of the last 40 samples
//!   before the erasure, attenuated per [`UNVOICED_FE_GAIN`]: each
//!   vector is rescaled by `FESCALE / MAG` where `MAG` is the segment's
//!   own magnitude sum.
//!
//! These three operations correspond to the §I.5 pseudo-code blocks:
//!
//! - **31SF** ([`FrameErasureExcitation::on_erased_cycle`]) — set the
//!   voiced/unvoiced flag, save `FEDELAY = KP` (voiced) or compute
//!   `AVMAG` (unvoiced) at onset, and recompute `FESCALE` every 10 ms.
//! - **31FE** ([`FrameErasureExcitation::extrapolate`]) — produce one
//!   extrapolated `ET()` vector.
//! - **31E** (the buffer shift inside [`FrameErasureExcitation::push`])
//!   — slide `ETPAST()` by `IDIM` and append the new vector. Block 31E
//!   runs in **every** frame (good or erased) so the history is ready
//!   when an erasure begins; in good frames the *real* decoded `ET()` is
//!   pushed, in erased frames the extrapolated one.
//!
//! Only the **floating-point** pseudo-code of §I.5 is realised here; the
//! fixed-point Q2 `ETPAST` representation is part of the deferred Annex G
//! fixed-point build.
//!
//! The "random" slide-back of the unvoiced path is supplied by the
//! caller as a [`SlideSource`]. The spec explicitly licenses any
//! aperiodic integer sequence in `5..=140` here — "interoperability with
//! the G.728 encoder is not an issue during frame erasures […] The
//! numbers also do not have to be truly random" — so the source is a
//! trait, letting tests drive it deterministically while production uses
//! a small in-crate LCG.
//!
//! Clean-room: every value and rule here is transcribed from §I.4.1 and
//! the §I.5.2/§I.5.3/§I.5.4 floating-point pseudo-code of
//! `T-REC-G.728-199905-AnnI.pdf`; no reference C / external
//! implementation was consulted.

use crate::consts::{
    FE_AVMAG_SAMPLES, FE_UNVOICED_MAX_SLIDE, FE_UNVOICED_MIN_SLIDE, IDIM, KPMAX, UNVOICED_FE_GAIN,
    VOICED_FE_GAIN, VTH,
};

/// Source of the unvoiced "random" slide-back distance (block 31FE).
///
/// Each call returns how far back in `ETPAST()` to slide the 5-sample
/// window, as a `FEDELAY` in `FE_UNVOICED_MIN_SLIDE..=FE_UNVOICED_MAX_SLIDE`
/// (`5..=140`). Implementations need only stay inside that closed range
/// and avoid being periodic; the spec requires nothing more.
pub trait SlideSource {
    /// Next slide-back distance in samples, clamped to `5..=140`.
    fn next_slide(&mut self) -> usize;
}

/// Default aperiodic slide-back source: a 32-bit linear congruential
/// generator mapped onto the `5..=140` range.
///
/// The spec only requires an aperiodic integer sequence in `5..=140`, so
/// the exact generator is unconstrained; this is a self-contained LCG
/// (Numerical-Recipes constants) with no external dependency. It is
/// **not** part of the bit-exact decode path — during an erasure the
/// output is concealment, not a reproducible reconstruction.
#[derive(Debug, Clone)]
pub struct LcgSlideSource {
    state: u32,
}

impl LcgSlideSource {
    /// New generator seeded with `seed` (any value is valid).
    #[must_use]
    pub fn new(seed: u32) -> Self {
        Self { state: seed }
    }
}

impl Default for LcgSlideSource {
    fn default() -> Self {
        // Arbitrary non-zero seed; the sequence's only requirement is
        // aperiodicity over the erasure, which any LCG seed satisfies.
        Self::new(0x1234_5678)
    }
}

impl SlideSource for LcgSlideSource {
    fn next_slide(&mut self) -> usize {
        // Numerical Recipes LCG: state = state·1664525 + 1013904223.
        self.state = self
            .state
            .wrapping_mul(1_664_525)
            .wrapping_add(1_013_904_223);
        // Map the top bits onto 5..=140 (span = 136 values).
        let span = (FE_UNVOICED_MAX_SLIDE - FE_UNVOICED_MIN_SLIDE + 1) as u32;
        FE_UNVOICED_MIN_SLIDE + ((self.state >> 16) % span) as usize
    }
}

/// One extrapolated gain-scaled excitation vector (`ET()`), `IDIM`
/// samples, ready to be fed to the (softened, §I.4.2) synthesis filter.
pub type ExcitationVector = [f64; IDIM];

/// Frame-erasure excitation extrapolator — the §I.4.1 mechanism wrapping
/// the `ETPAST()` history buffer plus blocks 31SF / 31FE / 31E.
///
/// Usage mirrors the decoder loop:
///
/// - In **good** frames, push the real decoded `ET()` vector via
///   [`Self::push`] (block 31E), and call [`Self::observe_good_cycle`]
///   once per adaptation cycle with the freshly computed `PTAP`/`KP`
///   from the postfilter's pitch extractor so the extrapolator knows the
///   last-good voicing decision when an erasure begins.
/// - At each **erased** adaptation cycle boundary, call
///   [`Self::on_erased_cycle`] (block 31SF) to refresh the scale.
/// - For each **erased** vector, call [`Self::extrapolate`] (block 31FE)
///   to synthesise `ET()`, then [`Self::push`] it (block 31E) so the
///   buffer keeps growing with the concealed excitation.
#[derive(Debug, Clone)]
pub struct FrameErasureExcitation {
    /// `ETPAST()` ring of past gain-scaled excitation samples, newest at
    /// the high index. Length `KPMAX + IDIM` so a slide-back of up to
    /// `KPMAX` from the most-recent vector still lands in-buffer.
    etpast: Vec<f64>,
    /// Last good single-tap pitch weight `PTAP` (block 83 of the last
    /// good adaptation cycle).
    last_good_ptap: f64,
    /// Last good pitch period `KP` (block 82 of the last good cycle).
    last_good_kp: usize,
    /// Voiced flag fixed at erasure onset (block 31SF, `VOICED`).
    voiced: bool,
    /// Periodic slide-back `FEDELAY = KP` saved at onset (voiced path).
    fedelay: usize,
    /// Magnitude reference `AVMAG` saved at onset (unvoiced path) —
    /// `(1/8)·Σ|ETPAST|` over the last 40 samples.
    avmag: f64,
    /// Current excitation scale `FESCALE` (block 31SF output).
    fescale: f64,
    /// Whether we are currently inside an erasure (onset not yet seen
    /// when `false`).
    in_erasure: bool,
}

/// Buffer capacity: `KPMAX` past samples plus the current `IDIM` vector.
const ETPAST_LEN: usize = KPMAX + IDIM;

impl Default for FrameErasureExcitation {
    fn default() -> Self {
        Self::new()
    }
}

impl FrameErasureExcitation {
    /// Fresh extrapolator with a zeroed `ETPAST()` history, not in an
    /// erasure.
    #[must_use]
    pub fn new() -> Self {
        Self {
            etpast: vec![0.0; ETPAST_LEN],
            last_good_ptap: 0.0,
            last_good_kp: 0,
            voiced: false,
            fedelay: 0,
            avmag: 0.0,
            fescale: 0.0,
            in_erasure: false,
        }
    }

    /// Record the pitch analysis (`PTAP`, `KP`) of a **good** adaptation
    /// cycle (block 83 / block 82). At the start of the next erasure
    /// block 31SF reads these last-good values to decide voiced vs
    /// unvoiced and to set `FEDELAY = KP`.
    ///
    /// Calling this also marks the erasure as ended: the next
    /// [`Self::on_erased_cycle`] is treated as a fresh onset (`N10MSEC`
    /// restarts at 0).
    pub fn observe_good_cycle(&mut self, ptap: f64, kp: usize) {
        self.last_good_ptap = ptap;
        self.last_good_kp = kp;
        self.in_erasure = false;
    }

    /// Block 31E — append a freshly produced `ET()` vector to
    /// `ETPAST()`, sliding the buffer back by `IDIM` samples.
    ///
    /// In good frames `et` is the real decoded gain-scaled excitation;
    /// in erased frames it is the [`Self::extrapolate`] output. The block
    /// runs unconditionally so the history is always current.
    pub fn push(&mut self, et: &ExcitationVector) {
        // "For I = -KPMAX+1 .. -IDIM: ETPAST(I) = ETPAST(I+IDIM)" then
        // "For I = 1 .. IDIM: ETPAST(I-IDIM) = ET(I)": shift the whole
        // window down by IDIM and append the new vector at the top.
        self.etpast.copy_within(IDIM.., 0);
        let tail = ETPAST_LEN - IDIM;
        self.etpast[tail..].copy_from_slice(et);
    }

    /// Block 31SF — set the frame-erasure flags and the excitation
    /// scaling factor at an erased adaptation-cycle boundary.
    ///
    /// `n10msec` is `FECOUNT >> 2` — the number of completed 10 ms
    /// intervals (4 adaptation cycles) of the current erasure. At onset
    /// (`n10msec == 0`) the voiced/unvoiced flag, `FEDELAY` and `AVMAG`
    /// are fixed from the last-good pitch analysis; thereafter only the
    /// scale `FESCALE` is re-derived from `n10msec`.
    pub fn on_erased_cycle(&mut self, n10msec: usize) {
        if !self.in_erasure {
            // First erased cycle: latch the onset decision (block 31SF,
            // "If N10MSEC == 0 …"). We treat the first on_erased_cycle of
            // an erasure as onset regardless of the caller's count, so a
            // stale in_erasure flag can never leak a prior decision.
            self.in_erasure = true;
            if self.last_good_ptap > VTH {
                self.voiced = true;
                self.fedelay = self.last_good_kp;
            } else {
                self.voiced = false;
                self.avmag = self.compute_avmag();
            }
        }

        // FESCALE schedule (block 31SF). VOICED uses VOICEDFEGAIN
        // directly; UNVOICED multiplies AVMAG by UNVOICEDFEGAIN.
        self.fescale = if self.voiced {
            VOICED_FE_GAIN.get(n10msec).copied().unwrap_or(0.0)
        } else {
            match UNVOICED_FE_GAIN.get(n10msec) {
                Some(&g) => self.avmag * g,
                None => 0.0,
            }
        };
    }

    /// Average-magnitude reference `AVMAG` (block 31SF unvoiced branch):
    /// `(1/8)·Σ|ETPAST(I)|` over the last `FE_AVMAG_SAMPLES` (= 40)
    /// samples before the erasure. The `/8` (rather than `/40`) is kept
    /// verbatim from the spec — only the ratio against the 5-sample
    /// segment magnitude in block 31FE matters.
    fn compute_avmag(&self) -> f64 {
        let start = ETPAST_LEN - FE_AVMAG_SAMPLES;
        let sum: f64 = self.etpast[start..].iter().map(|s| s.abs()).sum();
        0.125 * sum
    }

    /// Block 31FE — extrapolate one gain-scaled excitation vector.
    ///
    /// Must be called after [`Self::on_erased_cycle`] has set the scale
    /// for the current cycle. For the unvoiced path the `slide` source
    /// supplies the random slide-back distance.
    ///
    /// Returns the `IDIM`-sample `ET()` vector; the caller is expected to
    /// [`Self::push`] it immediately (block 31E).
    pub fn extrapolate<S: SlideSource>(&self, slide: &mut S) -> ExcitationVector {
        let mut et = [0.0f64; IDIM];
        // The vector currently being written occupies ETPAST positions
        // [tail .. tail+IDIM); a slide-back of `d` reads from
        // [tail-d .. tail-d+IDIM).
        let tail = ETPAST_LEN - IDIM;

        if self.voiced {
            // Voiced: ET(I) = FESCALE * ETPAST(I - FEDELAY).
            let d = self.fedelay.max(1).min(tail);
            for (i, out) in et.iter_mut().enumerate() {
                *out = self.fescale * self.etpast[tail - d + i];
            }
        } else {
            // Unvoiced: pull a random 5-sample segment, magnitude-match
            // it to FESCALE (= AVMAG·UNVOICEDFEGAIN) via RATIO = FESCALE
            // / MAG where MAG is the segment's own magnitude sum.
            let d = slide
                .next_slide()
                .clamp(FE_UNVOICED_MIN_SLIDE, FE_UNVOICED_MAX_SLIDE)
                .min(tail);
            let mut mag = 0.0;
            for (i, out) in et.iter_mut().enumerate() {
                let v = self.etpast[tail - d + i];
                *out = v;
                mag += v.abs();
            }
            if mag == 0.0 {
                mag = 1.0; // block 31FE: "If MAG = 0 do MAG = 1."
            }
            let ratio = self.fescale / mag;
            for out in &mut et {
                *out *= ratio;
            }
        }
        et
    }

    /// Whether the extrapolator is currently inside an erasure (onset
    /// latched, no good cycle observed since). Exposed for tests / audit.
    #[must_use]
    pub fn in_erasure(&self) -> bool {
        self.in_erasure
    }

    /// Onset voicing decision (block 31SF `VOICED`). Meaningful only
    /// while [`Self::in_erasure`] is `true`.
    #[must_use]
    pub fn voiced(&self) -> bool {
        self.voiced
    }

    /// Current excitation scale `FESCALE` (block 31SF output).
    #[must_use]
    pub fn fescale(&self) -> f64 {
        self.fescale
    }

    /// Periodic slide-back `FEDELAY` latched at onset (voiced path).
    #[must_use]
    pub fn fedelay(&self) -> usize {
        self.fedelay
    }

    /// Unvoiced magnitude reference `AVMAG` latched at onset.
    #[must_use]
    pub fn avmag(&self) -> f64 {
        self.avmag
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic slide source that yields a fixed schedule of
    /// distances, so block 31FE's unvoiced path is reproducible in tests.
    struct FixedSlide {
        seq: Vec<usize>,
        idx: usize,
    }
    impl FixedSlide {
        fn new(seq: &[usize]) -> Self {
            Self {
                seq: seq.to_vec(),
                idx: 0,
            }
        }
    }
    impl SlideSource for FixedSlide {
        fn next_slide(&mut self) -> usize {
            let d = self.seq[self.idx % self.seq.len()];
            self.idx += 1;
            d
        }
    }

    /// Fill ETPAST with a recognisable ramp so slide-backs are checkable:
    /// the most-recent sample (top of buffer) is the largest.
    fn ramped() -> FrameErasureExcitation {
        let mut fe = FrameErasureExcitation::new();
        // Push 30 distinct vectors (150 samples) of a rising ramp.
        let mut next = 1.0;
        for _ in 0..30 {
            let mut v = [0.0; IDIM];
            for s in &mut v {
                *s = next;
                next += 1.0;
            }
            fe.push(&v);
        }
        fe
    }

    #[test]
    fn vth_is_ppfth_over_1_4() {
        // §I.4.1: VTH = PPFTH/1.4, lower than the postfilter PPFTH.
        assert!((VTH - crate::consts::PPFTH / 1.4).abs() < 1e-15);
        const { assert!(VTH < crate::consts::PPFTH) };
    }

    #[test]
    fn voiced_gain_schedule_matches_spec() {
        // VOICEDFEGAIN(0:4) = .8, .8, .6, .4, .2 ; 0 beyond.
        assert_eq!(VOICED_FE_GAIN, [0.8, 0.8, 0.6, 0.4, 0.2]);
    }

    #[test]
    fn unvoiced_gain_schedule_matches_spec() {
        // UNVOICEDFEGAIN(0:5) = 1., 1., .8, .6, .4, .2 ; 0 beyond.
        assert_eq!(UNVOICED_FE_GAIN, [1.0, 1.0, 0.8, 0.6, 0.4, 0.2]);
    }

    #[test]
    fn push_shifts_and_appends() {
        // Block 31E: pushing a vector slides the window so the previous
        // top samples move down by IDIM and the new vector lands on top.
        let mut fe = FrameErasureExcitation::new();
        let a = [1.0, 2.0, 3.0, 4.0, 5.0];
        let b = [6.0, 7.0, 8.0, 9.0, 10.0];
        fe.push(&a);
        fe.push(&b);
        // Top IDIM samples are `b`; the IDIM below them are `a`.
        let top = ETPAST_LEN - IDIM;
        assert_eq!(&fe.etpast[top..], &b);
        assert_eq!(&fe.etpast[top - IDIM..top], &a);
    }

    #[test]
    fn voiced_onset_latches_kp_and_voiced_flag() {
        // PTAP > VTH at the last good cycle → voiced, FEDELAY = KP.
        let mut fe = ramped();
        fe.observe_good_cycle(0.9, 37); // strongly voiced
        fe.on_erased_cycle(0);
        assert!(fe.voiced());
        assert!(fe.in_erasure());
        assert_eq!(fe.fedelay(), 37);
        // First-cycle FESCALE = VOICEDFEGAIN(0) = 0.8.
        assert!((fe.fescale() - 0.8).abs() < 1e-15);
    }

    #[test]
    fn voiced_extrapolation_is_periodic_scaled_copy() {
        // ET(i) = FESCALE * ETPAST(i - FEDELAY): a scaled copy of the
        // segment FEDELAY samples back from the current vector position.
        let mut fe = ramped();
        let kp = 20;
        fe.observe_good_cycle(0.9, kp);
        fe.on_erased_cycle(0);
        let mut slide = FixedSlide::new(&[5]); // unused on voiced path
        let et = fe.extrapolate(&mut slide);
        let tail = ETPAST_LEN - IDIM;
        for i in 0..IDIM {
            let expected = 0.8 * fe.etpast[tail - kp + i];
            assert!((et[i] - expected).abs() < 1e-12, "sample {i}");
        }
    }

    #[test]
    fn voiced_fescale_decays_with_erasure_length() {
        // N10MSEC schedule: 0.8, 0.8, 0.6, 0.4, 0.2, then 0.
        let mut fe = ramped();
        fe.observe_good_cycle(0.9, 30);
        let want = [0.8, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0];
        for (n, &w) in want.iter().enumerate() {
            fe.on_erased_cycle(n);
            assert!((fe.fescale() - w).abs() < 1e-15, "n10msec={n}");
        }
    }

    #[test]
    fn unvoiced_onset_latches_avmag() {
        // PTAP <= VTH → unvoiced, AVMAG = (1/8)·Σ|last 40 ETPAST|.
        let mut fe = ramped();
        fe.observe_good_cycle(0.1, 30); // unvoiced (0.1 < VTH)
        fe.on_erased_cycle(0);
        assert!(!fe.voiced());
        // Recompute the reference independently: last 40 samples of the
        // ramp are the 40 largest values pushed.
        let start = ETPAST_LEN - FE_AVMAG_SAMPLES;
        let manual: f64 = fe.etpast[start..].iter().map(|s| s.abs()).sum::<f64>() * 0.125;
        assert!((fe.avmag() - manual).abs() < 1e-9);
        assert!(fe.avmag() > 0.0);
    }

    #[test]
    fn unvoiced_extrapolation_matches_target_magnitude() {
        // After magnitude-matching, Σ|ET| of the output equals FESCALE
        // (= AVMAG·UNVOICEDFEGAIN(n)) exactly, because RATIO = FESCALE/MAG
        // and the output is the segment times RATIO.
        let mut fe = ramped();
        fe.observe_good_cycle(0.1, 30);
        fe.on_erased_cycle(0); // n10msec=0 → UNVOICEDFEGAIN(0)=1.0
        let mut slide = FixedSlide::new(&[17, 42, 88]);
        let et = fe.extrapolate(&mut slide);
        let out_mag: f64 = et.iter().map(|s| s.abs()).sum();
        // FESCALE here = AVMAG * 1.0 = AVMAG.
        assert!((out_mag - fe.fescale()).abs() < 1e-9);
        assert!((fe.fescale() - fe.avmag()).abs() < 1e-15);
    }

    #[test]
    fn unvoiced_zero_segment_does_not_divide_by_zero() {
        // Block 31FE: "If MAG = 0 do MAG = 1." A slide into the zeroed
        // pre-history yields a zero segment; output must be finite (all
        // zero), not NaN/Inf.
        let mut fe = FrameErasureExcitation::new(); // all-zero ETPAST
        fe.observe_good_cycle(0.0, 30); // unvoiced, AVMAG = 0
        fe.on_erased_cycle(0);
        let mut slide = FixedSlide::new(&[100]);
        let et = fe.extrapolate(&mut slide);
        assert!(et.iter().all(|s| s.is_finite()));
        assert_eq!(et, [0.0; IDIM]);
    }

    #[test]
    fn good_cycle_ends_erasure_and_reonsets() {
        // A good cycle clears in_erasure; the next erased cycle is a
        // fresh onset that re-reads the (new) last-good voicing.
        let mut fe = ramped();
        fe.observe_good_cycle(0.9, 25); // voiced
        fe.on_erased_cycle(0);
        assert!(fe.voiced() && fe.in_erasure());
        // Recovery: a good cycle that is now unvoiced.
        fe.observe_good_cycle(0.05, 25);
        assert!(!fe.in_erasure());
        fe.on_erased_cycle(0); // fresh onset
        assert!(!fe.voiced(), "re-onset must use the new last-good PTAP");
    }

    #[test]
    fn lcg_slide_source_stays_in_range_and_is_not_constant() {
        // The production slide source must always land in 5..=140 and
        // must not be trivially periodic/constant over a short run.
        let mut s = LcgSlideSource::new(1);
        let mut seen = std::collections::BTreeSet::new();
        for _ in 0..200 {
            let d = s.next_slide();
            assert!((FE_UNVOICED_MIN_SLIDE..=FE_UNVOICED_MAX_SLIDE).contains(&d));
            seen.insert(d);
        }
        assert!(seen.len() > 50, "slide source is too repetitive");
    }

    #[test]
    fn voiced_then_extended_erasure_silences_after_50ms() {
        // Past index 4 (i.e. > 50 ms) the voiced FESCALE is 0, so the
        // extrapolated vector is all-zero regardless of history.
        let mut fe = ramped();
        fe.observe_good_cycle(0.9, 20);
        fe.on_erased_cycle(5); // beyond VOICEDFEGAIN length
        let mut slide = FixedSlide::new(&[10]);
        let et = fe.extrapolate(&mut slide);
        assert_eq!(et, [0.0; IDIM]);
    }
}
