//! Pitch period extraction — block 82 of Figure 7/G.728 (§4.7).
//!
//! Reads the LPC-residual stream produced by block 81 and produces a
//! single integer pitch period `p ∈ [KPMIN, KPMAX] = [20, 140]` per
//! adaptation frame. The output drives the long-term postfilter
//! (block 71) coefficient calculator (blocks 83 / 84, eq. 4-12..4-14)
//! in later rounds.
//!
//! ## Spec dataflow (§4.7, prose around eqs 4-7..4-11)
//!
//! 1. **Residual buffer (240 samples).** Block 82 holds the most recent
//!    240 LPC-residual samples in a buffer indexed
//!    `d(−139), d(−138), …, d(100)`. The spec's update rule places the
//!    four 5-sample vectors of each adaptation frame into
//!    `d(81..85)` (4th vector of the *previous* frame),
//!    `d(86..90)` (1st vector of the *current* frame),
//!    `d(91..95)` (2nd vector of the current frame) and
//!    `d(96..100)` (3rd vector of the current frame). Older samples in
//!    `d(−139..80)` are shifted-in from previous frames.
//! 2. **1 kHz lowpass + 4:1 decimation (Annex D).** The last 20 samples
//!    `d(81..100)` are filtered through the third-order elliptic lowpass
//!    `H(z) = (b_0 + b_1 z⁻¹ + b_2 z⁻² + b_3 z⁻³) / (1 + a_1 z⁻¹ + a_2 z⁻²
//!    + a_3 z⁻³)` (coefficients in [`crate::tables::BL`] /
//!    [`crate::tables::AL`]) and 4:1 down-sampled into the five new
//!    decimated samples `d̄(21..25)`. The decimated buffer holds 60
//!    samples `d̄(−34..25)`.
//! 3. **Coarse correlation search (eq. 4-7).** Compute
//!    `ρ(i) = Σ_{n=1..25} d̄(n)·d̄(n−i)` for `i ∈ {5,…,35}`. The 31
//!    lags span the decimated equivalents of pitch periods
//!    `4·5 = 20 ≤ p ≤ 4·35 = 140`. Pick the lag `τ` with the largest
//!    `ρ(τ)`.
//! 4. **Full-resolution refinement (eq. 4-8).** Compute
//!    `C(i) = Σ_{k=1..100} d(k)·d(k−i)` over the seven lags
//!    `i ∈ {4τ−3, …, 4τ+3}` and pick `p0` as the lag of largest `C(i)`.
//! 5. **Fundamental-vs-multiple resolution (eqs 4-9..4-11).** If the
//!    previous frame's pitch `p̂` falls outside the neighbourhood
//!    `[p0 − KPDELTA, p0 + KPDELTA]` (`KPDELTA = 6`), evaluate `C(i)`
//!    over the 13 lags `p̂−6..p̂+6` and pick `p1` as that
//!    neighbourhood's argmax. Compute the single-tap predictor weights
//!    `β0 = C(p0) / Σ d(k−p0)²` and `β1 = C(p1) / Σ d(k−p1)²`,
//!    clamped into `[0, 1]`. Then `p = p0` if `β1 ≤ TAPTH·β0` else
//!    `p = p1` (`TAPTH = 0.4`).
//! 6. **Carry-over.** The chosen `p` becomes the next frame's `p̂`.
//!
//! ## Cold start
//!
//! Before the decoder has decoded enough vectors to fill the
//! residual buffer the lowpass-decimate pass and the correlation sums
//! still operate on the (zero-initialised) buffer. The all-zero
//! correlations make `ρ(i) = C(i) = 0` for every lag, the argmax
//! defaults to the lowest index, and the resulting `p` is `KPMIN = 20`.
//! `p̂` is initialised to `KPMIN` so the fundamental-vs-multiple
//! check is also stable from the first frame.
//!
//! ## Output domain
//!
//! `p` is always clamped into `[KPMIN, KPMAX]`. The clamp is redundant
//! with the search range (coarse 4·5 = 20 ≤ p ≤ 4·35 = 140) but is
//! applied defensively, matching the decode-trace `KPMIN ≤ p ≤ KPMAX`
//! invariant called out at §7.1.

use crate::consts::{IDIM, KPDELTA, KPMAX, KPMIN, NFRSZ, NPWSZ, NUPDATE, TAPTH};
use crate::tables::{AL, BL};

/// Length of the spec's `d(−139..100)` residual buffer: 240 samples.
///
/// The spec indices run `−139 ≤ k ≤ 100`, i.e. 240 entries. We store
/// them contiguously with the in-buffer index `i = k + 139` so that
/// `d(−139)` lives at slot `0` and `d(100)` at slot `239`.
pub const RESIDUAL_BUF_LEN: usize = (KPMAX - 1) + NPWSZ + 1; // 139 + 100 + 1 = 240 — see eq above

/// Length of the spec's decimated `d̄(−34..25)` buffer: 60 samples.
///
/// The decimated indices run `−34 ≤ n ≤ 25`. We store them
/// contiguously with `n = −34` at slot `0` and `n = 25` at slot `59`.
/// The decimation factor is 4, so the decimated buffer is one-fourth
/// the size of the residual buffer.
pub const DECIMATED_BUF_LEN: usize = RESIDUAL_BUF_LEN / 4; // 60

/// Number of new residual samples per adaptation frame: `NFRSZ = 20`.
pub const FRAME_RESIDUAL: usize = NFRSZ;

/// Number of new decimated samples per adaptation frame: `5`
/// (= `NFRSZ / 4`).
pub const FRAME_DECIMATED: usize = FRAME_RESIDUAL / 4;

/// Coarse-search lag range (inclusive) over the decimated buffer.
const COARSE_LAG_MIN: usize = 5;
const COARSE_LAG_MAX: usize = 35;

/// Block 82 — pitch period extractor.
///
/// One instance per decoder. Carries:
///
/// * The 240-sample LPC-residual buffer covering `d(−139..100)`.
/// * The 60-sample decimated-residual buffer covering `d̄(−34..25)`.
/// * The 3-sample lowpass-filter input and output memories (Annex D
///   third-order elliptic, eq. of `H(z)` in module docs).
/// * The previous frame's chosen pitch period `p̂`.
/// * A cursor counting how many vector-residuals have been ingested
///   in the current adaptation frame (0..NUPDATE).
#[derive(Debug, Clone)]
pub struct PitchSearch {
    /// `d(−139..100)` LPC-residual buffer. `residual[k + 139] = d(k)`.
    residual: [f64; RESIDUAL_BUF_LEN],
    /// `d̄(−34..25)` decimated-residual buffer. `decimated[n + 34] =
    /// d̄(n)`.
    decimated: [f64; DECIMATED_BUF_LEN],
    /// Annex D lowpass-filter input memory.
    /// `lpf_x[0] = x(n−1)`, `lpf_x[1] = x(n−2)`, `lpf_x[2] = x(n−3)`.
    lpf_x: [f64; 3],
    /// Annex D lowpass-filter output memory.
    /// `lpf_y[0] = y(n−1)`, `lpf_y[1] = y(n−2)`, `lpf_y[2] = y(n−3)`.
    lpf_y: [f64; 3],
    /// Previous frame's chosen pitch period `p̂`. Initialised to
    /// `KPMIN` so the fundamental-vs-multiple resolution is stable
    /// from the first frame onward.
    p_prev: usize,
    /// How many vectors of LPC residual have been ingested in the
    /// current adaptation frame (`0..NUPDATE = 0..4`). Used by the
    /// caller to decide when to call [`Self::extract`].
    vec_idx: usize,
    /// Last pitch period the extractor reported (or `KPMIN` before
    /// the first `extract` call). Exposed via [`Self::last_pitch`]
    /// for tests and audit.
    last_p: usize,
}

impl Default for PitchSearch {
    fn default() -> Self {
        Self::new()
    }
}

impl PitchSearch {
    /// Construct a fresh pitch-search state with all buffers and
    /// filter memories zero-initialised, `p̂ = KPMIN`, and the
    /// per-frame vector cursor at 0.
    pub fn new() -> Self {
        Self {
            residual: [0.0; RESIDUAL_BUF_LEN],
            decimated: [0.0; DECIMATED_BUF_LEN],
            lpf_x: [0.0; 3],
            lpf_y: [0.0; 3],
            p_prev: KPMIN,
            vec_idx: 0,
            last_p: KPMIN,
        }
    }

    /// The previous frame's pitch period `p̂` — the value that the
    /// fundamental-vs-multiple resolution will compare against on the
    /// next [`Self::extract`] call. Exposed for tests and audit.
    pub fn previous_pitch(&self) -> usize {
        self.p_prev
    }

    /// The most recent pitch period the extractor reported. Equal to
    /// [`Self::previous_pitch`] outside of the per-extract update
    /// window; exposed separately so callers can observe the value
    /// without ambiguity.
    pub fn last_pitch(&self) -> usize {
        self.last_p
    }

    /// Per-frame vector cursor (`0..NUPDATE`). Resets to 0 each time
    /// [`Self::extract`] runs. Exposed for tests.
    pub fn vector_cursor(&self) -> usize {
        self.vec_idx
    }

    /// Read-only access to the residual buffer for tests and audit.
    /// `residual()[k + 139]` is `d(k)`.
    pub fn residual_buffer(&self) -> &[f64; RESIDUAL_BUF_LEN] {
        &self.residual
    }

    /// Read-only access to the decimated buffer for tests and audit.
    /// `decimated_buffer()[n + 34]` is `d̄(n)`.
    pub fn decimated_buffer(&self) -> &[f64; DECIMATED_BUF_LEN] {
        &self.decimated
    }

    /// Reset every internal piece of state to the cold-start defaults.
    /// Useful for tests.
    pub fn reset(&mut self) {
        self.residual = [0.0; RESIDUAL_BUF_LEN];
        self.decimated = [0.0; DECIMATED_BUF_LEN];
        self.lpf_x = [0.0; 3];
        self.lpf_y = [0.0; 3];
        self.p_prev = KPMIN;
        self.vec_idx = 0;
        self.last_p = KPMIN;
    }

    /// Push one IDIM-sample residual vector into the spec-position of
    /// the residual buffer.
    ///
    /// The spec stores the four 5-sample vectors of each frame into
    /// the slots `d(81..85)` (4th vector of the *previous* frame),
    /// `d(86..90)`, `d(91..95)` and `d(96..100)` in left-to-right
    /// order. We honour the same ordering via the per-frame counter
    /// `vec_idx` — caller responsibility is just "feed one residual
    /// vector per `decode_vector_postfiltered` call, then call
    /// [`Self::extract`] at the third vector of each frame."
    ///
    /// The frame's vector cursor advances `0 → 1 → 2 → 3` and wraps
    /// back to `0` at the start of each new frame inside
    /// [`Self::extract`]. The spec's special vector ordering is
    /// realised by writing into the slots:
    /// * cursor `0` (current frame's 1st vector) → `d(86..90)`,
    /// * cursor `1` (2nd)                          → `d(91..95)`,
    /// * cursor `2` (3rd)                          → `d(96..100)`,
    /// * cursor `3` (4th)                          → `d(81..85)`.
    ///
    /// The 4th vector lands into `d(81..85)` because by the time the
    /// 4th vector of the *current* frame is decoded, the next frame
    /// has not started yet — `d(81..85)` is therefore semantically
    /// "4th vector of the most-recently-completed frame," exactly the
    /// spec's intent.
    pub fn push_residual(&mut self, residual_vec: &[f64; IDIM]) {
        // Map the per-frame vector cursor to the d() target index.
        // Slots 81..85 / 86..90 / 91..95 / 96..100 — spec absolute
        // indices, with buffer offset +139.
        let target_d_start = match self.vec_idx {
            0 => 86, // current frame's 1st vector
            1 => 91, // current frame's 2nd vector
            2 => 96, // current frame's 3rd vector
            3 => 81, // current frame's 4th vector (= previous frame
            // 's "last vector" in the spec's view at the *next*
            // frame's extract time)
            _ => unreachable!(),
        };
        let buf_start = target_d_start + (KPMAX - 1); // +139 offset
        self.residual[buf_start..buf_start + IDIM].copy_from_slice(residual_vec);
        self.vec_idx = (self.vec_idx + 1) % NUPDATE;
    }

    /// Extract the pitch period for the current frame.
    ///
    /// Should be called once per adaptation frame, immediately after
    /// the third vector's residual has been pushed in
    /// [`Self::push_residual`] (spec: "at the third vector of each
    /// frame"). The caller is expected to keep feeding the 4th
    /// vector's residual *after* this call so the extractor's buffer
    /// reflects the spec's d(81..85) "4th vector of last frame"
    /// placement at the next extract.
    ///
    /// Per the §4.7 dataflow:
    ///
    /// 1. Lowpass-filter and 4:1 decimate the last 20 residual
    ///    samples `d(81..100)` into the five new decimated samples
    ///    `d̄(21..25)`. Slide the decimated buffer left by 5 (so
    ///    `d̄(−34..20)` is what `d̄(−29..25)` used to be).
    /// 2. Coarse correlation search over decimated lags
    ///    `5..=35` (eq. 4-7) → `τ`.
    /// 3. Full-resolution refinement over `4τ−3..=4τ+3` (eq. 4-8)
    ///    → `p0`. Pick `p0` as the argmax of `C(i)`.
    /// 4. If `p̂` is not in `[p0 − KPDELTA, p0 + KPDELTA]`, run a
    ///    second eq. 4-8 search over `p̂−6..=p̂+6` → `p1`. Compute
    ///    `β0 = C(p0) / Σ d(k − p0)²` (clamped to `[0, 1]`) and
    ///    `β1 = C(p1) / Σ d(k − p1)²` (also clamped) per eqs.
    ///    4-9 / 4-10. Select `p = p0` if `β1 ≤ TAPTH·β0`, else `p1`
    ///    (eq. 4-11).
    /// 5. Slide the residual buffer left by `NFRSZ = 20` so the next
    ///    frame's pushes land in the right slots (this realises the
    ///    spec's "d(−139..80) are previous LPC residuals shifted to
    ///    the correct time order" sentence). Reset the per-frame
    ///    vector cursor to 0 and stash `p` as `p̂` for the next call.
    ///
    /// Returns the chosen pitch period, clamped into `[KPMIN, KPMAX]`.
    pub fn extract(&mut self) -> usize {
        // ----- Step 1: lowpass-filter + 4:1 decimate ----------------
        // The Annex D filter is causal IIR with denominator order 3
        // and numerator order 3. The spec writes
        //   y(n) = Σ_{i=0..3} b_i·x(n−i) − Σ_{i=1..3} a_i·y(n−i)
        // We run this over the 20 input samples d(81..100), keep
        // only every 4th output (the 4:1 decimation), and use the
        // lpf_x / lpf_y memories so the filter state carries across
        // frame boundaries.
        //
        // The spec stores the filter coefficients in tables/AL, BL.
        // Per `consts.rs` / `tables.rs` the layout is
        // `AL = [1.0, a_1, a_2, a_3]` (with the implicit leading 1)
        // and `BL = [b_0, b_1, b_2, b_3]`.
        //
        // Implementation note: the buffer slot d(81) is at index
        // (81 + 139) = 220; d(100) is at 239. After the loop the
        // decimated samples are slotted into d̄(21..25) which after
        // the post-extract shift become d̄(-34..-30) of the new frame
        // and so on. We compute the 5 outputs into local variables
        // then push them.

        // Slide decimated buffer left by 5 first so the new samples
        // land in positions [DECIMATED_BUF_LEN - 5 ..]. The pre-slide
        // buffer holds d̄(−34..25); after this slide d̄(−39..20) sits
        // in slots [0..55] (the leftmost 5 entries d̄(−39..−35) drop
        // off the left, the next slide-in puts d̄(21..25) at the
        // rightmost 5 slots).
        for i in 0..(DECIMATED_BUF_LEN - FRAME_DECIMATED) {
            self.decimated[i] = self.decimated[i + FRAME_DECIMATED];
        }

        // Now run the lowpass over d(81..100) and decimate.
        let d81_idx = 81 + (KPMAX - 1);
        let mut sample_idx = 0usize;
        for n in 0..FRAME_RESIDUAL {
            let x = self.residual[d81_idx + n];
            // y(n) = b_0·x + b_1·x(n-1) + b_2·x(n-2) + b_3·x(n-3)
            //        − a_1·y(n-1) − a_2·y(n-2) − a_3·y(n-3)
            let y =
                BL[0] * x + BL[1] * self.lpf_x[0] + BL[2] * self.lpf_x[1] + BL[3] * self.lpf_x[2]
                    - AL[1] * self.lpf_y[0]
                    - AL[2] * self.lpf_y[1]
                    - AL[3] * self.lpf_y[2];
            // Shift LPF input memory: lpf_x[0] becomes most-recent
            // past input x(n-1); lpf_x[2] is x(n-3).
            self.lpf_x[2] = self.lpf_x[1];
            self.lpf_x[1] = self.lpf_x[0];
            self.lpf_x[0] = x;
            // Shift LPF output memory likewise.
            self.lpf_y[2] = self.lpf_y[1];
            self.lpf_y[1] = self.lpf_y[0];
            self.lpf_y[0] = y;
            // 4:1 decimation — keep every 4th output starting at n=3
            // so that the decimated rate is one-fourth the input rate.
            // The choice of phase (3 vs 0) shifts each decimated
            // sample by one input sample; we use phase 3 so the last
            // output sample corresponds to the most recent input
            // (d̄(25) lines up with d(100)). The spec doesn't pin
            // the phase explicitly; phase 3 is the only choice that
            // produces exactly 5 outputs from 20 inputs while keeping
            // the rightmost decimated sample co-aligned with the
            // rightmost input sample (the natural reading of "the
            // last 20 samples … resulting in five lowpass filtered
            // and decimated LPC residual samples, denoted
            // d̄(21), d̄(22), …, d̄(25), which are stored as the last
            // five samples in a decimated LPC residual buffer").
            if n % 4 == 3 {
                self.decimated[DECIMATED_BUF_LEN - FRAME_DECIMATED + sample_idx] = y;
                sample_idx += 1;
            }
        }
        debug_assert_eq!(sample_idx, FRAME_DECIMATED);

        // ----- Step 2: coarse correlation search (eq. 4-7) ----------
        //
        // ρ(i) = Σ_{n=1..25} d̄(n)·d̄(n−i)  for i = 5..=35
        //
        // In our buffer-offset terms d̄(n) lives at slot (n + 34) so
        // d̄(1) is slot 35 and d̄(25) is slot 59. The sum window
        // `n = 1..25` becomes slot range `[35, 60)`. For lag i,
        // d̄(n−i) lives at slot `35 + (n−1) − i = 34 + n − i`.
        let dbar_n_start = (KPMAX - 1) / 4 + 1; // 35
        let dbar_n_end = DECIMATED_BUF_LEN; // 60 (exclusive)
        let mut best_rho = f64::NEG_INFINITY;
        let mut tau = COARSE_LAG_MIN;
        for i in COARSE_LAG_MIN..=COARSE_LAG_MAX {
            let mut rho = 0.0f64;
            for n_slot in dbar_n_start..dbar_n_end {
                rho += self.decimated[n_slot] * self.decimated[n_slot - i];
            }
            if rho > best_rho {
                best_rho = rho;
                tau = i;
            }
        }

        // ----- Step 3: full-resolution refinement (eq. 4-8) ---------
        //
        // C(i) = Σ_{k=1..100} d(k)·d(k−i)  for i = 4τ−3..=4τ+3
        //
        // d(k) lives at slot k + 139; the sum window k = 1..100
        // becomes slot range [140, 240). For lag i, d(k−i) is at
        // slot 139 + k − i = 140 + (k−1) − i.
        let d_k_start = KPMAX - 1 + 1; // 140
        let d_k_end = RESIDUAL_BUF_LEN; // 240 (exclusive)

        // Clamp the refinement window into [KPMIN, KPMAX] — the
        // spec's pitch range. For τ = 5 (the lowest coarse lag) the
        // raw window 17..23 reaches below KPMIN = 20; for τ = 35
        // (the highest) the raw window 137..143 reaches above
        // KPMAX = 140 (and a lag of 143 would index past the
        // residual buffer's start, since the sum window k = 1..100
        // would otherwise hit d(-43) for k = 100, which is fine for
        // the buffer but k=1, d(1-143) = d(-142) IS past d(-139)).
        // The clamp protects both the buffer access and the
        // semantic "p must be in [KPMIN, KPMAX]" invariant the long-
        // term postfilter requires.
        let p0_min = (tau * 4).saturating_sub(3).max(KPMIN);
        let p0_max = (tau * 4 + 3).min(KPMAX);
        let (p0, _c_p0) = self.argmax_correlation(d_k_start, d_k_end, p0_min, p0_max);

        // ----- Step 4: fundamental-vs-multiple resolution ----------
        //
        // If the previous frame's p̂ is INSIDE the neighbourhood
        // [p0 − KPDELTA, p0 + KPDELTA] then there is no separate
        // candidate to test: spec says "if the time lag p0 obtained
        // above is not in the neighbourhood of p̂". We interpret
        // "not in the neighbourhood" symmetrically — the test is
        // whether the *previous* pitch falls inside the p0-centred
        // window. If it does, no second search is run; if it does
        // not, search a 13-lag window around p̂.
        let p_hat = self.p_prev;
        let p_in_neighbourhood = p_hat + KPDELTA >= p0 && p_hat <= p0 + KPDELTA;

        let p = if p_in_neighbourhood {
            p0
        } else {
            let p1_min = p_hat.saturating_sub(KPDELTA).max(KPMIN);
            let p1_max = (p_hat + KPDELTA).min(KPMAX);
            let (p1, _c_p1) = self.argmax_correlation(d_k_start, d_k_end, p1_min, p1_max);

            // ----- β0, β1 (eqs 4-9, 4-10) -----------------------------
            // β = C(p) / Σ d(k − p)²  for k = 1..100, clamped to
            // [0, 1].
            let beta0 = self.beta(d_k_start, d_k_end, p0);
            let beta1 = self.beta(d_k_start, d_k_end, p1);

            // Eq. 4-11: p = p0 if β1 ≤ TAPTH·β0 else p1.
            if beta1 <= TAPTH * beta0 {
                p0
            } else {
                p1
            }
        };

        // ----- Step 5: post-extract buffer slide, stash p̂ -----------
        //
        // After the extract, prepare for the next frame:
        //   - Slide d(−139..100) left by NFRSZ=20 so the next push of
        //     vec_idx=0 (current frame's 1st vector) lands at d(86..90)
        //     of the SLIDED buffer, which spec-wise is now d(86..90)
        //     of the NEW frame. The 4th vector of the current frame
        //     (still to be decoded after the extract) will land at
        //     d(81..85) of the new frame (vec_idx=3 path in
        //     push_residual). Per the spec d(81..85) at extract time
        //     is exactly "4th vector of the previous frame," so the
        //     placement is consistent.
        for i in 0..(RESIDUAL_BUF_LEN - FRAME_RESIDUAL) {
            self.residual[i] = self.residual[i + FRAME_RESIDUAL];
        }
        for i in (RESIDUAL_BUF_LEN - FRAME_RESIDUAL)..RESIDUAL_BUF_LEN {
            self.residual[i] = 0.0;
        }

        // Clamp into [KPMIN, KPMAX] defensively (the search range
        // already guarantees this; the clamp is a guard against future
        // refactors).
        let p = p.clamp(KPMIN, KPMAX);
        self.p_prev = p;
        self.last_p = p;

        // Per §4.7 prose, after the third-vector extract the
        // current frame's *4th* vector is the next residual to be
        // pushed; it must land at d(81..85) of the SLID buffer. In
        // our cursor table that mapping is vec_idx == 3. Setting
        // the cursor to 3 here aligns the next push with the spec's
        // "4th vector of the just-extracted frame" slot. After that
        // push the cursor wraps to 0 — which then maps the new
        // frame's 1st vector to d(86..90), matching the spec.
        self.vec_idx = 3;

        p
    }

    /// Compute the argmax of `C(i) = Σ_{k=k_start..k_end} d(k)·d(k−i)`
    /// over the closed integer range `i ∈ [i_min, i_max]` and return
    /// `(i_max_arg, C(i_max_arg))`. Helper for [`Self::extract`].
    fn argmax_correlation(
        &self,
        k_start: usize,
        k_end: usize,
        i_min: usize,
        i_max: usize,
    ) -> (usize, f64) {
        let mut best_c = f64::NEG_INFINITY;
        let mut best_i = i_min;
        for i in i_min..=i_max {
            let mut c = 0.0f64;
            for k_slot in k_start..k_end {
                // k_slot - i is safe: i ≤ KPMAX = 140, k_slot ≥ 140.
                c += self.residual[k_slot] * self.residual[k_slot - i];
            }
            if c > best_c {
                best_c = c;
                best_i = i;
            }
        }
        (best_i, best_c)
    }

    /// Compute `β = C(p) / Σ d(k − p)²` for `k = k_start..k_end`,
    /// clamped into `[0, 1]` per the §4.7 prose after eq. 4-9 / 4-10
    /// ("The value of β0 is then clamped between 0 and 1"). On a
    /// zero denominator return `0.0` (the energy at lag `p` is zero,
    /// so the predictor tap is undefined; the spec doesn't spell out
    /// this corner explicitly, but `0` matches the "no useful
    /// prediction here" intent and falls through to the `p0` branch
    /// of eq. 4-11 since `0 ≤ TAPTH·β0` is always true).
    fn beta(&self, k_start: usize, k_end: usize, p: usize) -> f64 {
        let mut num = 0.0f64;
        let mut denom = 0.0f64;
        for k_slot in k_start..k_end {
            let d_k = self.residual[k_slot];
            let d_kp = self.residual[k_slot - p];
            num += d_k * d_kp;
            denom += d_kp * d_kp;
        }
        if denom == 0.0 {
            0.0
        } else {
            (num / denom).clamp(0.0, 1.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---------- Construction + state invariants ----------------------

    #[test]
    fn fresh_state_is_zeroed_and_p_hat_is_kpmin() {
        let ps = PitchSearch::new();
        assert!(ps.residual_buffer().iter().all(|&v| v == 0.0));
        assert!(ps.decimated_buffer().iter().all(|&v| v == 0.0));
        assert_eq!(ps.previous_pitch(), KPMIN);
        assert_eq!(ps.last_pitch(), KPMIN);
        assert_eq!(ps.vector_cursor(), 0);
    }

    #[test]
    fn buf_dimensions_match_spec() {
        // Decode-trace §7.1 / §4.7 prose: 240 LPC-residual samples
        // and 60 decimated samples.
        assert_eq!(RESIDUAL_BUF_LEN, 240);
        assert_eq!(DECIMATED_BUF_LEN, 60);
        assert_eq!(FRAME_RESIDUAL, 20);
        assert_eq!(FRAME_DECIMATED, 5);
    }

    #[test]
    fn reset_returns_to_cold_start() {
        let mut ps = PitchSearch::new();
        ps.push_residual(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        ps.push_residual(&[6.0, 7.0, 8.0, 9.0, 10.0]);
        ps.reset();
        assert!(ps.residual_buffer().iter().all(|&v| v == 0.0));
        assert_eq!(ps.previous_pitch(), KPMIN);
        assert_eq!(ps.vector_cursor(), 0);
    }

    // ---------- Vector cursor + push placement ------------------------

    #[test]
    fn vector_cursor_walks_zero_through_three_then_wraps() {
        let mut ps = PitchSearch::new();
        assert_eq!(ps.vector_cursor(), 0);
        ps.push_residual(&[0.0; IDIM]);
        assert_eq!(ps.vector_cursor(), 1);
        ps.push_residual(&[0.0; IDIM]);
        assert_eq!(ps.vector_cursor(), 2);
        ps.push_residual(&[0.0; IDIM]);
        assert_eq!(ps.vector_cursor(), 3);
        ps.push_residual(&[0.0; IDIM]);
        assert_eq!(ps.vector_cursor(), 0);
    }

    #[test]
    fn push_residual_places_vectors_in_spec_d_slots() {
        // Per §4.7 prose: 1st current-frame vec → d(86..90),
        // 2nd → d(91..95), 3rd → d(96..100), 4th → d(81..85).
        let mut ps = PitchSearch::new();
        let v1 = [11.0, 12.0, 13.0, 14.0, 15.0];
        let v2 = [21.0, 22.0, 23.0, 24.0, 25.0];
        let v3 = [31.0, 32.0, 33.0, 34.0, 35.0];
        let v4 = [41.0, 42.0, 43.0, 44.0, 45.0];
        ps.push_residual(&v1);
        ps.push_residual(&v2);
        ps.push_residual(&v3);
        ps.push_residual(&v4);
        let off = KPMAX - 1; // +139 base
                             // d(86..90) = v1
        for k in 0..IDIM {
            assert_eq!(ps.residual_buffer()[off + 86 + k], v1[k]);
        }
        // d(91..95) = v2
        for k in 0..IDIM {
            assert_eq!(ps.residual_buffer()[off + 91 + k], v2[k]);
        }
        // d(96..100) = v3
        for k in 0..IDIM {
            assert_eq!(ps.residual_buffer()[off + 96 + k], v3[k]);
        }
        // d(81..85) = v4
        for k in 0..IDIM {
            assert_eq!(ps.residual_buffer()[off + 81 + k], v4[k]);
        }
    }

    // ---------- Cold-start extract: defaults to KPMIN ----------------

    #[test]
    fn extract_on_all_zero_buffer_returns_kpmin() {
        // All correlations are zero → argmax is the first index of
        // the search range (KPMIN-derived) → p == KPMIN.
        let mut ps = PitchSearch::new();
        // Feed 3 zero vectors so the vec cursor is at extract-time
        // (the spec extracts at the third vector of each frame).
        ps.push_residual(&[0.0; IDIM]);
        ps.push_residual(&[0.0; IDIM]);
        ps.push_residual(&[0.0; IDIM]);
        let p = ps.extract();
        assert_eq!(p, KPMIN);
        // p̂ should now also be KPMIN.
        assert_eq!(ps.previous_pitch(), KPMIN);
    }

    #[test]
    fn extract_clamps_p_into_kpmin_kpmax() {
        let mut ps = PitchSearch::new();
        for _ in 0..3 {
            ps.push_residual(&[0.0; IDIM]);
        }
        let p = ps.extract();
        assert!((KPMIN..=KPMAX).contains(&p));
    }

    #[test]
    fn extract_sets_vec_cursor_to_three_for_fourth_vector_push() {
        // Per §4.7 prose the 4th vector of the just-extracted frame
        // is the NEXT residual the decoder pushes; it must land at
        // d(81..85) of the (post-slide) buffer. Cursor 3 in our
        // push_residual mapping writes to that slot, so extract
        // leaves the cursor at 3.
        let mut ps = PitchSearch::new();
        for _ in 0..3 {
            ps.push_residual(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        }
        let _ = ps.extract();
        assert_eq!(ps.vector_cursor(), 3);
    }

    // ---------- Real periodic signal → expected pitch ---------------

    /// Helper: fill both the residual buffer (period `period`) and
    /// the decimated buffer (period `period / 4`) with unit-impulse
    /// trains so that AFTER `extract` runs its step-1 left-shift of
    /// the decimated buffer (by `FRAME_DECIMATED = 5`), the impulses
    /// remain aligned at multiples of `dec_period` *as seen by the
    /// correlation search*.
    ///
    /// The decimated buffer needs explicit seeding because
    /// `extract` only refreshes its rightmost 5 entries from the
    /// freshly-decimated `d(81..100)`; the leftmost 55 entries hold
    /// previously-decimated samples that this helper synthesises.
    /// The +5 pre-shift offset on the decimated indices realigns
    /// pre-call slot `5 + j·dec_period` to post-shift slot
    /// `j·dec_period`, which is what the search reads.
    ///
    /// `period` must be a multiple of 4 ≤ KPMAX so the decimated
    /// period `period / 4` is an integer ≤ 35 (the coarse-search
    /// upper bound).
    fn fill_with_impulse_train(ps: &mut PitchSearch, period: usize) {
        assert!(period.is_multiple_of(4), "period must be a multiple of 4");
        for slot in 0..RESIDUAL_BUF_LEN {
            if slot.is_multiple_of(period) {
                ps.residual[slot] = 1.0;
            }
        }
        let dec_period = period / 4;
        for slot in 0..DECIMATED_BUF_LEN {
            // Offset by FRAME_DECIMATED so that AFTER the step-1
            // left-shift inside extract() the impulses re-land at
            // multiples of dec_period in the post-shift buffer.
            if slot >= FRAME_DECIMATED && (slot - FRAME_DECIMATED).is_multiple_of(dec_period) {
                ps.decimated[slot] = 1.0;
            }
        }
    }

    #[test]
    fn coarse_search_locks_to_period_on_clean_impulse_train() {
        // An impulse train with period P in the residual domain has
        // its autocorrelation maximised at lags that are multiples
        // of P. The coarse search runs over decimated lags 5..=35;
        // the full-resolution refinement then snaps to the true P.
        //
        // We pick P = 40 (mid-range, comfortably above KPMIN = 20 and
        // below KPMAX = 140). The decimated equivalent τ = 10 is in
        // [5, 35], and the refinement window 4·10−3..4·10+3 covers
        // exactly P = 40.
        let mut ps = PitchSearch::new();
        let period = 40usize;
        fill_with_impulse_train(&mut ps, period);
        // We have to skip the per-vector push (the test bypasses
        // push_residual and writes the buffer directly via the
        // module-private helper). The decimator pass inside extract
        // re-fills d̄(21..25) from d(81..100); for the impulse train
        // those last 20 samples are zero/one according to the
        // pattern, which is fine — the bulk of the signal lives at
        // d(-139..80).
        // Force vec_idx to 3 so extract treats this as a third-vector
        // call (any value works since we're bypassing push).
        ps.vec_idx = 3;
        let p = ps.extract();
        assert_eq!(p, period, "expected pitch {period}, got {p}");
    }

    #[test]
    fn coarse_search_locks_to_kpmin_period() {
        // Pitch period = KPMIN. The coarse search's lowest lag is 5;
        // the refinement window covers 4·5−3..4·5+3 = 17..23, which
        // includes 20. Resulting p should equal KPMIN = 20.
        let mut ps = PitchSearch::new();
        fill_with_impulse_train(&mut ps, KPMIN);
        ps.vec_idx = 3;
        let p = ps.extract();
        assert_eq!(p, KPMIN);
    }

    #[test]
    fn coarse_search_locks_to_kpmax_period() {
        // Pitch period = KPMAX = 140. The decimated equivalent is
        // τ = 35 (the upper coarse-search bound). Refinement window
        // 4·35−3..4·35+3 = 137..143 includes 140.
        let mut ps = PitchSearch::new();
        fill_with_impulse_train(&mut ps, KPMAX);
        ps.vec_idx = 3;
        let p = ps.extract();
        assert_eq!(p, KPMAX);
    }

    // ---------- Previous-pitch carry-over ----------------------------

    #[test]
    fn extract_stashes_p_into_p_prev() {
        let mut ps = PitchSearch::new();
        fill_with_impulse_train(&mut ps, 60);
        ps.vec_idx = 3;
        let p = ps.extract();
        // p_prev must now equal whatever extract returned.
        assert_eq!(ps.previous_pitch(), p);
        assert_eq!(ps.last_pitch(), p);
    }

    #[test]
    fn fundamental_vs_multiple_path_runs_without_panic() {
        // When p̂ is far away from p0 (i.e. > KPDELTA away), the
        // extractor evaluates the second-window correlation. Seed
        // p_prev to a value that's outside the coarse-search-derived
        // neighbourhood for a P = 40 impulse train, and confirm the
        // call completes and the result is in [KPMIN, KPMAX].
        let mut ps = PitchSearch::new();
        fill_with_impulse_train(&mut ps, 40);
        ps.p_prev = 100; // 100 is > 6 away from 40 → triggers the
                         // multiple-resolution path.
        ps.vec_idx = 3;
        let p = ps.extract();
        assert!((KPMIN..=KPMAX).contains(&p));
    }

    // ---------- Buffer slide after extract ---------------------------

    #[test]
    fn extract_slides_residual_buffer_left_by_nfrsz() {
        let mut ps = PitchSearch::new();
        // Mark d(81..100) with a known pattern (vec 4 of last frame
        // + 3 vectors of current frame). After extract, the last
        // NFRSZ=20 samples of the buffer should be zero (slid in),
        // and the slot that USED to hold d(81) should now hold what
        // was at d(101). Since d(101) doesn't exist in the buffer
        // (the buffer only covers d(−139..100)), the post-slide
        // sample at the d(81) slot will be zero.
        ps.push_residual(&[10.0, 20.0, 30.0, 40.0, 50.0]); // → d(86..90)
        ps.push_residual(&[60.0, 70.0, 80.0, 90.0, 100.0]); // → d(91..95)
        ps.push_residual(&[110.0, 120.0, 130.0, 140.0, 150.0]); // → d(96..100)
        let _ = ps.extract();
        // After the post-extract slide, the rightmost NFRSZ samples
        // are zeroed (no future data to fill them yet).
        for slot in (RESIDUAL_BUF_LEN - FRAME_RESIDUAL)..RESIDUAL_BUF_LEN {
            assert_eq!(ps.residual_buffer()[slot], 0.0, "slot {slot}");
        }
        // The pre-slide content of d(86..90) (slots 225..230) now
        // lives at slots 205..210 (shifted left by 20).
        let off = KPMAX - 1;
        for k in 0..IDIM {
            assert_eq!(
                ps.residual_buffer()[off + 86 - FRAME_RESIDUAL + k],
                [10.0, 20.0, 30.0, 40.0, 50.0][k]
            );
        }
    }

    // ---------- Lowpass + decimation behaviour -----------------------

    #[test]
    fn decimator_advances_decimated_buffer_by_five() {
        // Each extract shifts the decimated buffer left by 5 and
        // inserts 5 new samples on the right. Seed the buffer with
        // a tagged pattern and confirm the shift.
        let mut ps = PitchSearch::new();
        for slot in 0..DECIMATED_BUF_LEN {
            ps.decimated[slot] = slot as f64 + 1.0;
        }
        ps.vec_idx = 3;
        let _ = ps.extract();
        // Post-extract: slots [0..55] hold the pre-extract slots
        // [5..60] = values 6.0..60.0; slots [55..60] hold the 5 new
        // decimated outputs (whose value depends on the residual
        // input, which is zero here, so they are zero or driven by
        // residual-buffer churn from the slide).
        for slot in 0..(DECIMATED_BUF_LEN - FRAME_DECIMATED) {
            assert_eq!(
                ps.decimated_buffer()[slot],
                (slot + FRAME_DECIMATED) as f64 + 1.0
            );
        }
    }

    #[test]
    fn extract_produces_finite_pitch_on_random_inputs() {
        // Smoke: random-looking residuals shouldn't push the search
        // out of bounds or into a non-deterministic pitch.
        let mut ps = PitchSearch::new();
        for v in 0usize..16 {
            let mut sd = [0.0f64; IDIM];
            for k in 0..IDIM {
                sd[k] = ((v * 7 + k * 13) as f64 % 100.0) - 50.0;
            }
            ps.push_residual(&sd);
            // After every 3rd vector of a frame (spec ICOUNT = 3), run
            // an extract. With push_residual the per-frame cursor is
            // at 0 right after extract, so the 3rd vector of frame N
            // corresponds to push index 3·N + 2 → v % NUPDATE == 2.
            if v % NUPDATE == 2 {
                let p = ps.extract();
                assert!((KPMIN..=KPMAX).contains(&p));
            }
        }
    }

    // ---------- Lowpass-only smoke (filter state survives) ----------

    #[test]
    fn lowpass_filter_memory_advances_across_extracts() {
        // Two successive extracts on identical impulse-train data
        // should produce the same pitch decision. The lowpass
        // memory state advances but the periodicity of the input
        // keeps the correlation peak in place.
        let mut ps = PitchSearch::new();
        fill_with_impulse_train(&mut ps, 60);
        ps.vec_idx = 3;
        let p1 = ps.extract();
        // Re-fill the right half of the buffer (it was slid out)
        // with the same pattern so the second extract sees the same
        // periodic content.
        fill_with_impulse_train(&mut ps, 60);
        ps.vec_idx = 3;
        let p2 = ps.extract();
        assert_eq!(p1, p2, "successive same-input extracts must agree");
    }
}
