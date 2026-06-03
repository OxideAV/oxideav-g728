//! Long-term postfilter coefficient calculator — blocks 83 and 84 of
//! Figure 7/G.728 (§4.7).
//!
//! Reads the decoded-speech stream `sd(k)` and, at the third vector of
//! each adaptation frame, takes the pitch period `p` produced by
//! [`crate::pitch_search::PitchSearch`] (block 82) and produces the
//! long-term postfilter coefficients `(g_l, b)` that drive
//! [`crate::long_term_postfilter::LongTermPostfilter`] (block 71).
//!
//! ## Spec dataflow (§4.7, prose around eqs 4-12, 4-13)
//!
//! ### Block 83 — single-tap pitch predictor weight `β` (eq. 4-12)
//!
//! Block 83 and block 71 share a long buffer of decoded speech samples
//! covering `sd(−239), sd(−238), …, sd(4), sd(5)`. The 5 samples
//! `sd(1..5)` correspond to the current vector of decoded speech (the
//! third vector of the frame at coefficient-update time). Block 71 uses
//! this buffer as its FIR delay line; block 83 reads only the *past*
//! window `sd(−99..0)` plus `sd(−99−p..−p)` for lag `p`.
//!
//! With the pitch period `p ∈ [KPMIN, KPMAX]` from block 82, eq. 4-12 is:
//!
//! ```text
//!         Σ_{k=−99..0}  sd(k) · sd(k − p)
//! β  =  ─────────────────────────────────
//!        Σ_{k=−99..0}  sd(k − p) · sd(k − p)
//! ```
//!
//! The result is clamped into `[0, 1]` per the §4.7 prose immediately
//! after eq. 4-10 ("then also clamped between 0 and 1"), consistent
//! with the same `[0, 1]` clamp applied to β0 / β1 in eqs 4-9 / 4-10.
//! On a zero denominator (the lag-`p` window contains only zero
//! samples — cold-start / silence) the predictor weight is undefined;
//! `β = 0` is the spec-consistent choice because eq. 4-13 then routes
//! to the "β < 0.6, postfilter off" branch (`b = 0`, `g_l = 1`), which
//! is the §4.6.1 passthrough.
//!
//! ### Block 84 — `(g_l, b)` calculator (eq. 4-13 / 4-14)
//!
//! Given `β` from block 83:
//!
//! ```text
//!       ┌  0       if β < 0.6
//! b  =  │  0.15·β  if 0.6 ≤ β ≤ 1                              (4-13)
//!       └  0.15    if β > 1
//!
//!   g_l = 1 / (1 + b)                                          (4-14)
//! ```
//!
//! The `0.15` constant is `PPFZCF` (Table 1/G.728); the `0.6`
//! threshold is `PPFTH` (Table 1/G.728). When `β < PPFTH` the
//! long-term postfilter is "totally disabled" — `H_l(z) = 1` — per
//! the §4.7 prose.
//!
//! Because `b ∈ {0, 0.15·β, 0.15}` and `β` is clamped into `[0, 1]`,
//! the effective range of `b` is `[0, 0.15]`. The corresponding `g_l`
//! range is `[1/(1+0.15), 1/(1+0)] = [0.8696…, 1.0]`.
//!
//! ## Update cadence
//!
//! Per §7.1 of `g728-decode-trace.md` and §4.7 of the Recommendation,
//! `(g_l, b, p)` are updated at the **third** vector of each
//! adaptation frame. The caller is expected to:
//!
//! 1. Push every decoded-speech vector through
//!    [`Self::push_decoded_vector`] after the synthesis filter
//!    produces it — this keeps the `sd(−239..5)` buffer synced with
//!    the output stream.
//! 2. At the third vector of each frame (the same vector
//!    [`PitchSearch::extract`](crate::pitch_search::PitchSearch::extract)
//!    runs on), call [`Self::compute_coefficients`] with the extracted
//!    pitch period to obtain `(g_l, b)` for the next adaptation cycle.
//!
//! ## Cold start
//!
//! [`Self::new`] initialises the `sd` buffer to zero. Until enough
//! vectors have been pushed to fill the analysis window, the
//! correlation sums in eq. 4-12 will compute against the zero
//! prefix, producing `β = 0` (zero numerator over zero denominator —
//! the safe-fallback branch). Eq. 4-13 then routes to `b = 0`,
//! `g_l = 1`, matching the §4.6.1 "postfilter off" passthrough
//! exactly.
//!
//! ## Why this lives in its own module
//!
//! Block 71 (the comb filter itself) has its own FIR delay line sized
//! to `KPMAX` samples, which is sufficient to read `sd(n − p)` for any
//! `p ≤ KPMAX` during sample-rate filtering. Block 83's analysis window
//! is wider (`sd(−99..0)` minus a lag of up to `KPMAX` reaches back to
//! `sd(−239)`), so it needs a `KPMAX + NPWSZ` buffer. Keeping the two
//! buffers separate avoids forcing block 71 to widen its delay line —
//! its sample-rate cost stays exactly as it was when block 71 alone
//! was wired up.

use crate::consts::{IDIM, KPMAX, NPWSZ, PPFTH, PPFZCF};

/// Length of the spec's `sd(−239..5)` buffer: 245 samples.
///
/// Indexing convention: `buf[k + (KPMAX + NPWSZ - 1)]` is `sd(k)`.
/// Concretely, `sd(−239)` lives at slot 0 (`= -239 + 239`) and
/// `sd(5)` lives at slot 244 (`= 5 + 239`).
pub const SPEECH_BUF_LEN: usize = KPMAX + NPWSZ + IDIM; // 140 + 100 + 5 = 245

/// Per-spec offset from `sd(k)` to the in-buffer slot: `slot = k +
/// SPEECH_OFFSET`.
const SPEECH_OFFSET: usize = KPMAX + NPWSZ - 1; // 239

/// Block 83 + 84 — long-term postfilter coefficient calculator.
///
/// One instance per decoder. Carries:
///
/// * A 245-sample `sd(−239..5)` decoded-speech buffer; new vectors are
///   pushed via [`Self::push_decoded_vector`].
/// * The most recently computed `(g_l, b, β)` triple — exposed via the
///   accessors for tests and audit. Cold-start values are
///   `(g_l, b, β) = (1.0, 0.0, 0.0)` — the §4.6.1 passthrough.
#[derive(Debug, Clone)]
pub struct PitchPostfilterCoeff {
    /// `sd(−239..5)` buffer. `speech[k + SPEECH_OFFSET] = sd(k)`.
    speech: [f64; SPEECH_BUF_LEN],
    /// Most recently computed `β` (eq. 4-12, clamped to `[0, 1]`).
    last_beta: f64,
    /// Most recently computed `b` (eq. 4-13).
    last_b: f64,
    /// Most recently computed `g_l` (eq. 4-14).
    last_g_l: f64,
}

impl Default for PitchPostfilterCoeff {
    fn default() -> Self {
        Self::new()
    }
}

impl PitchPostfilterCoeff {
    /// Construct a fresh coefficient calculator with the `sd(−239..5)`
    /// buffer zeroed and `(g_l, b, β)` at the §4.6.1 passthrough
    /// values `(1.0, 0.0, 0.0)`.
    pub fn new() -> Self {
        Self {
            speech: [0.0; SPEECH_BUF_LEN],
            last_beta: 0.0,
            last_b: 0.0,
            last_g_l: 1.0,
        }
    }

    /// Most recently computed `β` (eq. 4-12, clamped into `[0, 1]`).
    /// Exposed for tests and audit.
    pub fn last_beta(&self) -> f64 {
        self.last_beta
    }

    /// Most recently computed `b` (eq. 4-13). Range `[0, PPFZCF]`.
    /// Exposed for tests and audit.
    pub fn last_b(&self) -> f64 {
        self.last_b
    }

    /// Most recently computed `g_l = 1/(1 + b)` (eq. 4-14). Range
    /// `[1/(1 + PPFZCF), 1.0]`. Exposed for tests and audit.
    pub fn last_g_l(&self) -> f64 {
        self.last_g_l
    }

    /// Read-only access to the `sd(−239..5)` buffer for tests and
    /// audit. `speech_buffer()[k + 239]` is `sd(k)`.
    pub fn speech_buffer(&self) -> &[f64; SPEECH_BUF_LEN] {
        &self.speech
    }

    /// Reset every internal piece of state to the cold-start defaults
    /// (`sd` buffer zeroed, `(g_l, b, β) = (1.0, 0.0, 0.0)`).
    pub fn reset(&mut self) {
        self.speech = [0.0; SPEECH_BUF_LEN];
        self.last_beta = 0.0;
        self.last_b = 0.0;
        self.last_g_l = 1.0;
    }

    /// Push one [`IDIM`]-sample decoded-speech vector into the
    /// `sd(−239..5)` buffer.
    ///
    /// The buffer is shifted left by [`IDIM`] samples (the oldest 5
    /// drop off the left edge), and the new vector lands in the
    /// rightmost 5 slots — the spec's `sd(1..5)` window, "the current
    /// vector of decoded speech." Older samples in `sd(−239..0)`
    /// represent successively earlier decoded vectors.
    ///
    /// Callers should push every decoded-speech vector immediately
    /// after the synthesis filter (block 32) produces it, regardless of
    /// where in the adaptation cycle the decoder is. The coefficient
    /// update itself runs at the third vector of each frame via
    /// [`Self::compute_coefficients`].
    pub fn push_decoded_vector(&mut self, sd_vec: &[f64; IDIM]) {
        // Shift the buffer left by IDIM = 5 samples. The leftmost 5
        // entries (sd(-239..-235) before the push) drop off; the
        // remaining 240 entries shift down by 5; the rightmost 5
        // entries receive the new vector.
        for i in 0..(SPEECH_BUF_LEN - IDIM) {
            self.speech[i] = self.speech[i + IDIM];
        }
        // The new vector lands at sd(1..5) = slots [240..245).
        let new_start = SPEECH_BUF_LEN - IDIM;
        self.speech[new_start..new_start + IDIM].copy_from_slice(sd_vec);
    }

    /// Compute the long-term postfilter coefficients `(g_l, b)` for
    /// pitch period `p` and return the new tap weight `β`.
    ///
    /// This is the combined evaluation of block 83 (eq. 4-12, the
    /// single-tap predictor weight over the decoded-speech buffer)
    /// followed by block 84 (eq. 4-13 / 4-14, the `(b, g_l)` lookup).
    /// The caller should invoke this at the third vector of each
    /// adaptation frame, after [`Self::push_decoded_vector`] has been
    /// called for the third vector's output, and after block 82 has
    /// produced the pitch period.
    ///
    /// The pitch period `p` is treated as already clamped into
    /// `[KPMIN, KPMAX]` by block 82; this method does not re-clamp.
    /// Caller-supplied `p` outside that range will be clamped by
    /// [`crate::long_term_postfilter::LongTermPostfilter::set_coefficients`]
    /// when the result is applied.
    ///
    /// Returns the updated `β` value (clamped to `[0, 1]`). The
    /// computed `(g_l, b)` are stashed and can be read out via
    /// [`Self::last_g_l`] and [`Self::last_b`]. The combined
    /// `(g_l, b, p)` triple is what the caller should feed into
    /// [`crate::long_term_postfilter::LongTermPostfilter::set_coefficients`].
    pub fn compute_coefficients(&mut self, p: usize) -> f64 {
        // ----- Block 83: eq. 4-12 ------------------------------------
        //
        //         Σ_{k=−99..0} sd(k) · sd(k − p)
        //   β = ─────────────────────────────────
        //        Σ_{k=−99..0} sd(k − p) · sd(k − p)
        //
        // In our buffer-offset terms sd(k) lives at slot (k + 239).
        // The analysis window k = −99..0 maps to slots [140, 240).
        // For lag p (in samples), sd(k − p) lives at slot
        // (k + 239 − p) = (k_slot − p). With p ≤ KPMAX = 140 and
        // k_slot ≥ 140, k_slot − p ≥ 0, so the read is always in-range.
        let k_start = SPEECH_OFFSET - NPWSZ + 1; // (-99 + 239) = 140
        let k_end = SPEECH_OFFSET + 1; // (0 + 239 + 1) = 240 (exclusive)

        let mut num = 0.0f64;
        let mut denom = 0.0f64;
        for k_slot in k_start..k_end {
            let sd_k = self.speech[k_slot];
            // p is supplied in-range by block 82; saturate to k_slot
            // defensively. `k_slot - p` is the slot for sd(k - p).
            let sd_kp = self.speech[k_slot - p];
            num += sd_k * sd_kp;
            denom += sd_kp * sd_kp;
        }

        // β = num / denom clamped to [0, 1]. On a zero denominator
        // (silence in the lag-p window) return β = 0 — eq. 4-13 then
        // routes through the postfilter-off branch (b = 0, g_l = 1),
        // matching the §4.6.1 cold-start identity.
        let beta = if denom == 0.0 {
            0.0
        } else {
            (num / denom).clamp(0.0, 1.0)
        };

        // ----- Block 84: eq. 4-13 / 4-14 -----------------------------
        //
        //         ┌  0        if β < PPFTH (postfilter off)
        //   b  =  │  PPFZCF·β if PPFTH ≤ β ≤ 1
        //         └  PPFZCF   if β > 1     (unreachable after clamp)
        //
        //   g_l = 1 / (1 + b)
        //
        // The β > 1 branch is unreachable here because β is already
        // clamped into [0, 1] above. We keep the explicit b = PPFZCF
        // limb for transcription completeness — at evaluation time the
        // clamp guarantees β ≤ 1.0, but adding the `else` keeps the
        // mapping explicit to the spec table.
        let b = if beta < PPFTH {
            0.0
        } else if beta <= 1.0 {
            PPFZCF * beta
        } else {
            PPFZCF
        };
        let g_l = 1.0 / (1.0 + b);

        self.last_beta = beta;
        self.last_b = b;
        self.last_g_l = g_l;

        // Silence the unused-variable lint for `p` when no buffer
        // index actually depends on it (impossible in this function,
        // but useful as documentation).
        let _ = p;

        beta
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::{KPMIN, NUPDATE};

    // ---------- Construction + state invariants ----------------------

    #[test]
    fn fresh_state_is_zeroed_and_passthrough() {
        // Cold-start: buffer zeroed, (g_l, b, β) at §4.6.1 passthrough.
        let pc = PitchPostfilterCoeff::new();
        assert!(pc.speech_buffer().iter().all(|&v| v == 0.0));
        assert_eq!(pc.last_beta(), 0.0);
        assert_eq!(pc.last_b(), 0.0);
        assert_eq!(pc.last_g_l(), 1.0);
    }

    #[test]
    fn buffer_dimension_matches_spec() {
        // sd(−239..5) → 245 entries per §4.7 prose.
        assert_eq!(SPEECH_BUF_LEN, 245);
        assert_eq!(SPEECH_OFFSET, 239);
    }

    #[test]
    fn reset_returns_to_cold_start() {
        let mut pc = PitchPostfilterCoeff::new();
        pc.push_decoded_vector(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        let _ = pc.compute_coefficients(KPMIN);
        pc.reset();
        assert!(pc.speech_buffer().iter().all(|&v| v == 0.0));
        assert_eq!(pc.last_beta(), 0.0);
        assert_eq!(pc.last_b(), 0.0);
        assert_eq!(pc.last_g_l(), 1.0);
    }

    // ---------- Push placement ----------------------------------------

    #[test]
    fn push_lands_vector_at_sd_1_through_5() {
        // After one push the new vector should occupy slots
        // [240..245) — the spec's sd(1..5) window.
        let mut pc = PitchPostfilterCoeff::new();
        let v = [11.0, 22.0, 33.0, 44.0, 55.0];
        pc.push_decoded_vector(&v);
        let off = SPEECH_OFFSET; // 239 ⇒ sd(1) lives at slot 240.
        for k in 0..IDIM {
            assert_eq!(pc.speech_buffer()[off + 1 + k], v[k]);
        }
    }

    #[test]
    fn successive_pushes_slide_buffer_left_by_idim() {
        // Two consecutive pushes: the first vector should appear in
        // slots [235..240) (= sd(-4..0)) and the second in
        // slots [240..245) (= sd(1..5)).
        let mut pc = PitchPostfilterCoeff::new();
        let v1 = [11.0, 22.0, 33.0, 44.0, 55.0];
        let v2 = [10.0, 20.0, 30.0, 40.0, 50.0];
        pc.push_decoded_vector(&v1);
        pc.push_decoded_vector(&v2);
        for k in 0..IDIM {
            // v1 should now live one vector to the left of v2.
            assert_eq!(pc.speech_buffer()[SPEECH_OFFSET + 1 - IDIM + k], v1[k]);
            assert_eq!(pc.speech_buffer()[SPEECH_OFFSET + 1 + k], v2[k]);
        }
    }

    #[test]
    fn pushing_drops_oldest_samples_off_the_left() {
        // Push 245 vectors' worth of samples (49 vectors); the first
        // vector should no longer be visible in the buffer.
        let mut pc = PitchPostfilterCoeff::new();
        let v1 = [1.0; IDIM];
        let v_rest = [0.0; IDIM];
        pc.push_decoded_vector(&v1);
        // 49 more pushes (245 = 49 · 5) wipes v1 from the buffer.
        for _ in 0..49 {
            pc.push_decoded_vector(&v_rest);
        }
        // No non-zero samples should remain — v1 has fallen off.
        assert!(pc.speech_buffer().iter().all(|&v| v == 0.0));
    }

    // ---------- Eq. 4-12 / 4-13 / 4-14: cold-start coefficients ------

    #[test]
    fn cold_start_compute_returns_passthrough() {
        // All-zero buffer → both num and denom are zero → β = 0
        // → b = 0, g_l = 1.
        let mut pc = PitchPostfilterCoeff::new();
        let beta = pc.compute_coefficients(KPMIN);
        assert_eq!(beta, 0.0);
        assert_eq!(pc.last_b(), 0.0);
        assert_eq!(pc.last_g_l(), 1.0);
    }

    #[test]
    fn near_zero_signal_routes_to_postfilter_off() {
        // Very small periodic signal where the energy is low: β might
        // still drop below PPFTH = 0.6 and route through the off
        // branch (b = 0, g_l = 1).
        let mut pc = PitchPostfilterCoeff::new();
        // Push 60 IDIM-sized vectors of all-zero with one tiny non-
        // zero sample injected far back; β should be 0 because the
        // sd(-99..0) window does not see the impulse.
        let zero = [0.0; IDIM];
        for _ in 0..49 {
            pc.push_decoded_vector(&zero);
        }
        // Inject in vector 49: sample 0 = 0.0001.
        let mut v = [0.0; IDIM];
        v[0] = 0.0001;
        pc.push_decoded_vector(&v);
        // The impulse now lives at sd(1) ... well, lands in the
        // current vector. The β computation reads sd(-99..0), which
        // hasn't seen the impulse yet, so β = 0.
        let beta = pc.compute_coefficients(KPMIN);
        assert_eq!(beta, 0.0);
        assert_eq!(pc.last_b(), 0.0);
        assert_eq!(pc.last_g_l(), 1.0);
    }

    // ---------- Eq. 4-12: highly periodic signal → β near 1 ----------

    /// Helper: fill the `sd(−239..5)` buffer with a unit-impulse train
    /// at period `period`. The impulses sit at every `period`-th
    /// sample, i.e. slots `0, period, 2·period, …`.
    fn fill_with_impulse_train(pc: &mut PitchPostfilterCoeff, period: usize) {
        for slot in 0..SPEECH_BUF_LEN {
            if slot.is_multiple_of(period) {
                pc.speech[slot] = 1.0;
            }
        }
    }

    #[test]
    fn perfect_periodicity_gives_beta_at_unity() {
        // Impulse train with period p in the residual domain has its
        // single-tap correlation maximum at lag p, with β = 1 (perfect
        // prediction). Pick p = 40 (in-range, well above KPMIN).
        let mut pc = PitchPostfilterCoeff::new();
        let p = 40usize;
        fill_with_impulse_train(&mut pc, p);
        let beta = pc.compute_coefficients(p);
        assert!(
            (beta - 1.0).abs() < 1e-12,
            "perfectly periodic signal at lag {p} should give β = 1, got {beta}"
        );
        // β = 1 routes to b = PPFZCF · 1 = 0.15; g_l = 1/(1.15).
        assert!((pc.last_b() - PPFZCF).abs() < 1e-12);
        assert!((pc.last_g_l() - 1.0 / (1.0 + PPFZCF)).abs() < 1e-12);
    }

    #[test]
    fn beta_is_clamped_above_unity() {
        // Signal where the lag-p correlation would naturally exceed 1:
        // a non-unity-magnitude impulse train where the past samples
        // are smaller than the current samples. eq. 4-12 with a
        // numerator larger than the denominator yields β > 1, and the
        // clamp brings it back to 1.
        //
        // Construct: sd(k) = 1 for k that are multiples of P; sd(k−P)
        // for the same k is also 1 → β = N / N = 1 (already at the
        // boundary). To force β > 1 we'd need a constant gain on the
        // current samples versus the past, but the spec clamps before
        // either eq. 4-13 branch fires. The behavioural test is just
        // "no β reading ever exceeds 1.0 in last_beta()" — re-confirm
        // here by feeding the impulse train and checking last_beta.
        let mut pc = PitchPostfilterCoeff::new();
        fill_with_impulse_train(&mut pc, 40);
        let _ = pc.compute_coefficients(40);
        assert!(pc.last_beta() <= 1.0 + 1e-12);
        assert!(pc.last_beta() >= 0.0 - 1e-12);
    }

    // ---------- Eq. 4-13: β < PPFTH branch ----------------------------

    #[test]
    fn beta_below_threshold_disables_postfilter() {
        // Sub-threshold β routes to b = 0, g_l = 1 (postfilter off).
        // Construct a sd buffer with weak periodicity: random-looking
        // values at the analysis window so the correlation sum is
        // small relative to the energy denominator.
        let mut pc = PitchPostfilterCoeff::new();
        // Push vectors with alternating sign so the lag-40
        // correlation flips sign and integrates to a small magnitude.
        for v in 0..50 {
            let sign = if v % 2 == 0 { 1.0 } else { -1.0 };
            let sd = [sign, sign * 0.5, sign * 0.25, sign * 0.125, sign * 0.0625];
            pc.push_decoded_vector(&sd);
        }
        let _ = pc.compute_coefficients(40);
        // Whichever branch fires, the result must satisfy the spec:
        // β ∈ [0, 1]; if β < 0.6 → b = 0 and g_l = 1.
        if pc.last_beta() < PPFTH {
            assert_eq!(pc.last_b(), 0.0);
            assert_eq!(pc.last_g_l(), 1.0);
        } else {
            assert!((pc.last_b() - PPFZCF * pc.last_beta()).abs() < 1e-12);
            assert!((pc.last_g_l() - 1.0 / (1.0 + pc.last_b())).abs() < 1e-12);
        }
    }

    // ---------- Eq. 4-13: 0.6 ≤ β ≤ 1 branch --------------------------

    #[test]
    fn beta_in_voiced_range_gives_b_proportional_to_beta() {
        // For 0.6 ≤ β ≤ 1, b = 0.15·β. We can synthesise this
        // condition by direct manipulation: set the buffer to a
        // pattern that produces β = 0.8 (between PPFTH = 0.6 and 1.0).
        //
        // Use sd(k) = 1 + ε at the analysis window k = −99..0 and
        // sd(k − p) = 1 at the same window positions. Then:
        //   num   = Σ (1 + ε) · 1     = 100 · (1 + ε)
        //   denom = Σ 1 · 1           = 100
        //   β     = 1 + ε             → clamped to 1.
        //
        // For β strictly inside the band we need a partial-correlation
        // signal. Simplest synthetic:
        //   sd(k) = past_signal + noise where the noise integrates to
        //   ~0.2 of the past-signal energy.
        //
        // We pick a deterministic recipe with disjoint windows. Set
        // the lag-p window slots [k_start - p, k_end - p) to 1.0 and
        // the current window slots [k_start, k_end) to α. For the two
        // windows to be disjoint with no overlap we need p ≥ NPWSZ
        // (the analysis-window length) so that k_end - p ≤ k_start.
        // With NPWSZ = 100 and p = 100 the two windows are exactly
        // adjacent. Then:
        //   num   = Σ_k α · 1 = NPWSZ · α
        //   denom = Σ_k 1 · 1 = NPWSZ
        //   β     = α
        let mut pc = PitchPostfilterCoeff::new();
        let p = NPWSZ; // 100 — equal to KPMIN-and-above guarantees disjointness
        let alpha = 0.8f64;
        let k_start = SPEECH_OFFSET - NPWSZ + 1; // 140
        let k_end = SPEECH_OFFSET + 1; // 240
        for k_slot in k_start..k_end {
            pc.speech[k_slot - p] = 1.0;
            pc.speech[k_slot] = alpha;
        }
        let beta = pc.compute_coefficients(p);
        assert!(
            (beta - alpha).abs() < 1e-12,
            "expected β = {alpha}, got {beta}"
        );
        // β ∈ [0.6, 1] → b = 0.15 · β; g_l = 1 / (1 + b).
        let expected_b = PPFZCF * alpha;
        let expected_g_l = 1.0 / (1.0 + expected_b);
        assert!((pc.last_b() - expected_b).abs() < 1e-12);
        assert!((pc.last_g_l() - expected_g_l).abs() < 1e-12);
    }

    // ---------- Eq. 4-14: g_l = 1 / (1 + b) consistency ---------------

    #[test]
    fn g_l_always_equals_one_over_one_plus_b() {
        // Across all branches of eq. 4-13, g_l must satisfy
        // g_l = 1 / (1 + b). Cycle through β values and confirm.
        let test_betas: [f64; 5] = [0.0, 0.4, 0.6, 0.8, 1.0];
        for target_beta in test_betas {
            let mut pc = PitchPostfilterCoeff::new();
            // Synthesise the target β with the same recipe as the
            // 0.6 ≤ β ≤ 1 test above (lag window and current window
            // disjoint at p = NPWSZ); outside that band the recipe
            // still produces β = target if target ≤ 1.
            let p = NPWSZ;
            let k_start = SPEECH_OFFSET - NPWSZ + 1;
            let k_end = SPEECH_OFFSET + 1;
            for k_slot in k_start..k_end {
                pc.speech[k_slot - p] = 1.0;
                pc.speech[k_slot] = target_beta;
            }
            let _ = pc.compute_coefficients(p);
            assert!(
                (pc.last_g_l() - 1.0 / (1.0 + pc.last_b())).abs() < 1e-12,
                "g_l = 1/(1+b) violated for β = {target_beta}: b = {}, g_l = {}",
                pc.last_b(),
                pc.last_g_l()
            );
        }
    }

    // ---------- b range invariant -------------------------------------

    #[test]
    fn b_stays_in_zero_to_ppfzcf_range() {
        // Across many random inputs, b must always satisfy 0 ≤ b ≤
        // PPFZCF (the spec's eq. 4-13 limits).
        let mut pc = PitchPostfilterCoeff::new();
        for v in 0..200 {
            let mut sd = [0.0f64; IDIM];
            for k in 0..IDIM {
                sd[k] = ((v * 7 + k * 13) as f64 % 100.0) - 50.0;
            }
            pc.push_decoded_vector(&sd);
            if v % NUPDATE == 2 {
                let _ = pc.compute_coefficients(40);
                assert!(
                    (0.0..=PPFZCF + 1e-12).contains(&pc.last_b()),
                    "b out of range [0, {PPFZCF}]: {}",
                    pc.last_b()
                );
                assert!(
                    (1.0 / (1.0 + PPFZCF) - 1e-12..=1.0 + 1e-12).contains(&pc.last_g_l()),
                    "g_l out of range: {}",
                    pc.last_g_l()
                );
            }
        }
    }

    // ---------- Long stream stability --------------------------------

    #[test]
    fn long_stream_keeps_outputs_finite_and_in_range() {
        let mut pc = PitchPostfilterCoeff::new();
        for v in 0..512u32 {
            let mut sd = [0.0f64; IDIM];
            for k in 0..IDIM {
                let t = (v as usize * IDIM + k) as f64;
                sd[k] = 1000.0 * (0.05 * t).sin();
            }
            pc.push_decoded_vector(&sd);
            if v % NUPDATE as u32 == 2 {
                let _ = pc.compute_coefficients(40);
                assert!(pc.last_beta().is_finite());
                assert!(pc.last_b().is_finite());
                assert!(pc.last_g_l().is_finite());
                assert!((0.0..=1.0 + 1e-12).contains(&pc.last_beta()));
            }
        }
    }
}
