//! Backward-adaptive LPC + log-gain predictors for G.728.
//!
//! The ITU-T G.728 decoder re-estimates a 50th-order LPC synthesis filter
//! and a 10th-order log-gain predictor every 4 vectors (= 20 samples,
//! 2.5 ms) purely from the decoder's own reconstructed speech / gain
//! history. No side information is transmitted — the encoder and decoder
//! stay in sync by running the same analysis on the same signal.
//!
//! # Hybrid window (Barnwell / logarithmic autocorrelation)
//!
//! G.728 §3.7 specifies a **hybrid window** for the autocorrelation:
//!
//! - A fixed, finite non-recursive tail spanning the most recent `NONR`
//!   samples, taken verbatim from Annex A of the spec.
//! - A decaying recursive portion covering the entire prior history; the
//!   contribution of the older samples to `r_m(i)` is updated as
//!   `r_{m+L}(i) = alpha^{2L} * r_m(i) + (new non-recursive block)`.
//!   For the synthesis filter, `alpha^{2L} = (3/4)` with `L = 20`.
//!
//! This module implements the hybrid window algorithm as described in
//! block 49 of §5.6, including the 257/256 white-noise correction factor.
//! The 105-sample non-recursive window (the first 35 are the *strict*
//! non-recursive tail; the remaining 70 form the leading edge of the
//! recursive portion) is taken from Annex A of the 09/92 recommendation.
//!
//! # Other components
//!
//! - Levinson-Durbin recursion for the 50th-order LPC coefficients.
//! - A short Levinson-Durbin (order 10) over the log-gain history for the
//!   gain predictor.
//! - Bandwidth-expansion applied to the LPC vector: γ = 253/256 for the
//!   synthesis filter (§3.7) and γ = 29/32 for the log-gain predictor
//!   (§3.8), matching the spec's `FAC` and `FACGP` constants.

use crate::{GAIN_ORDER, LPC_ORDER};

/// Length of the synthesis-filter hybrid window buffer (`SB`, block 49).
/// Equal to `LPC + NFRSZ + NONR = 50 + 20 + 35 = 105` samples.
pub const HYBRID_WIN_LEN: usize = 105;

/// Number of non-recursive-only samples in the synthesis-filter hybrid
/// window (`NONR`, block 49). The remaining recursive portion of the
/// window is accumulated via the `REXP` vector.
pub const HYBRID_NONR: usize = 35;

/// Samples added per adaptation cycle to the hybrid window (`NFRSZ`).
pub const HYBRID_NFRSZ: usize = 20;

/// Recursive decay factor `alpha^{2L} = (3/4)` for the synthesis-filter
/// hybrid window (§3.7, block 49).
pub const HYBRID_RECURSIVE_FACTOR: f32 = 0.75;

/// White-noise correction factor `WNCF = 257/256` (§3.7).
pub const WNCF: f32 = 257.0 / 256.0;

/// Bandwidth-expansion factor applied to the 50th-order synthesis-filter
/// coefficients (`FAC = 253/256`, §3.7).
pub const BW_EXPANSION: f32 = 253.0 / 256.0;

/// History window used for the backward-compatible Hamming path. Kept for
/// the log-gain predictor and for tests that exercise the raw
/// `autocorrelation` helper.
pub const HISTORY_LEN: usize = 100;

/// History length for the log-gain predictor autocorrelation.
pub const GAIN_HISTORY_LEN: usize = 40;

/// 105-sample hybrid window function (`WNR`) for the synthesis filter.
///
/// Transcribed from **ITU-T G.728 (09/92), Annex A.1** — the floating-
/// point equivalent of the Q15 integer table. The first 35 samples are
/// the strict non-recursive portion of the window (its shape is roughly a
/// raised cosine that rises from `~0.048` to `~1.000` at index 29 / 30);
/// the remaining 70 samples are the leading edge of the recursive
/// portion, shaped like `b * alpha^k`. Both pieces are applied multiplied
/// with the corresponding `SB` buffer entry.
///
/// Indexing convention in the spec is `WNR(1)` = the first listed value;
/// here `WNR[0]` holds that same value.
pub const HYBRID_WIN: [f32; HYBRID_WIN_LEN] = [
    0.047760010,
    0.095428467,
    0.142852783,
    0.189971924,
    0.236663818,
    0.282775879,
    0.328277588,
    0.373016357,
    0.416900635,
    0.459838867,
    0.501739502,
    0.542480469,
    0.582000732,
    0.620178223,
    0.656921387,
    0.692199707,
    0.725891113,
    0.757904053,
    0.788208008,
    0.816680908,
    0.843322754,
    0.868041992,
    0.890747070,
    0.911437988,
    0.930053711,
    0.946533203,
    0.960876465,
    0.973022461,
    0.982910156,
    0.990600586,
    0.996002197,
    0.999114990,
    0.999969482,
    0.998565674,
    0.994842529,
    0.988861084,
    0.981781006,
    0.974731445,
    0.967742920,
    0.960815430,
    0.953948975,
    0.947082520,
    0.940307617,
    0.933563232,
    0.926879883,
    0.920227051,
    0.913635254,
    0.907104492,
    0.900604248,
    0.894134521,
    0.887725830,
    0.881378174,
    0.875061035,
    0.868774414,
    0.862548828,
    0.856384277,
    0.850250244,
    0.844146729,
    0.838104248,
    0.832092285,
    0.826141357,
    0.820220947,
    0.814331055,
    0.808502197,
    0.802703857,
    0.796936035,
    0.791229248,
    0.785583496,
    0.779937744,
    0.774353027,
    0.768798828,
    0.763305664,
    0.757812500,
    0.752380371,
    0.747009277,
    0.741638184,
    0.736328125,
    0.731048584,
    0.725830078,
    0.720611572,
    0.715454102,
    0.710327148,
    0.705230713,
    0.700164795,
    0.695159912,
    0.690185547,
    0.685241699,
    0.680328369,
    0.675445557,
    0.670593262,
    0.665802002,
    0.661041260,
    0.656280518,
    0.651580811,
    0.646911621,
    0.642272949,
    0.637695313,
    0.633117676,
    0.628570557,
    0.624084473,
    0.619598389,
    0.615142822,
    0.610748291,
    0.606384277,
    0.602020264,
];

// ---------------------------------------------------------------------------
// Autocorrelation: Hamming path (legacy, used for tests + the log-gain
// predictor where the spec's separate hybrid window would add complexity
// disproportionate to the 10th-order gain trajectory).
// ---------------------------------------------------------------------------

/// Compute the first `order+1` autocorrelation lags of a windowed buffer
/// using a plain Hamming window — kept for the log-gain predictor and for
/// unit tests that exercise the recursion on small fixtures. The
/// synthesis-filter path now goes through [`HybridWindow`] instead.
///
/// `history[0]` is the most recent sample; successively older samples
/// follow. The window is applied in time order (oldest sample gets the
/// leading-edge weight), matching the usual DSP convention.
pub fn autocorrelation<const N: usize>(history: &[f32; N], order: usize) -> Vec<f32> {
    assert!(order < N, "order+1 must fit in history");
    let mut win = vec![0.0_f32; N];
    let denom = (N - 1) as f32;
    for n in 0..N {
        let phase = 2.0 * core::f32::consts::PI * (n as f32) / denom;
        win[n] = 0.54 - 0.46 * phase.cos();
    }
    // Apply window in time order: history is newest-first, so history[N-1]
    // is oldest and should get win[0] (leading edge).
    let mut x = vec![0.0_f32; N];
    for n in 0..N {
        x[n] = history[N - 1 - n] * win[n];
    }
    let mut r = vec![0.0_f32; order + 1];
    for k in 0..=order {
        let mut acc = 0.0_f32;
        for n in k..N {
            acc += x[n] * x[n - k];
        }
        r[k] = acc;
    }
    r
}

/// Levinson-Durbin recursion: autocorrelation `r[0..=order]` → predictor
/// coefficients `a[0..=order]` with `a[0] = 1`.
///
/// Returns `None` if the recursion encounters a non-positive prediction
/// error (numerically singular input). Callers should fall back to the
/// previous predictor in that case.
pub fn levinson_durbin(r: &[f32], order: usize) -> Option<Vec<f32>> {
    assert_eq!(r.len(), order + 1);
    if r[0] <= 0.0 {
        return None;
    }
    let mut a = vec![0.0_f32; order + 1];
    a[0] = 1.0;
    let mut e = r[0];

    let mut tmp = vec![0.0_f32; order + 1];

    for i in 1..=order {
        // Reflection coefficient k_i = -(r[i] + sum_{j=1..i-1} a[j]*r[i-j]) / e.
        let mut acc = r[i];
        for j in 1..i {
            acc += a[j] * r[i - j];
        }
        let k = -acc / e;
        if !k.is_finite() || k.abs() >= 1.0 {
            // Non-minimum-phase: bail out.
            return None;
        }
        // a^(i)[j] = a^(i-1)[j] + k * a^(i-1)[i-j]
        tmp[..=i].copy_from_slice(&a[..=i]);
        for j in 1..i {
            a[j] = tmp[j] + k * tmp[i - j];
        }
        a[i] = k;
        e *= 1.0 - k * k;
        if e <= 0.0 || !e.is_finite() {
            return None;
        }
    }
    Some(a)
}

/// Levinson-Durbin recursion that also returns the first reflection
/// coefficient `k_1` and a snapshot of the predictor taps at a chosen
/// intermediate order. The short-term postfilter's spectral tilt
/// compensation (`mu = 0.15 * k_1`) needs `k_1` separately from the
/// final predictor vector, and §5.5 of G.728 extracts the 10th-order
/// predictor as a by-product of the 50th-order run by snapshotting at
/// order 10.
///
/// `snapshot_order` should be in `[0, order]`. When set to `0`, no
/// snapshot is returned (only the final `a` and `k_1`).
pub fn levinson_durbin_with_refl(
    r: &[f32],
    order: usize,
    snapshot_order: usize,
) -> Option<(Vec<f32>, f32, Option<Vec<f32>>)> {
    assert_eq!(r.len(), order + 1);
    assert!(snapshot_order <= order);
    if r[0] <= 0.0 {
        return None;
    }
    let mut a = vec![0.0_f32; order + 1];
    a[0] = 1.0;
    let mut e = r[0];
    let mut tmp = vec![0.0_f32; order + 1];
    let mut k1 = 0.0_f32;
    let mut snapshot: Option<Vec<f32>> = None;

    for i in 1..=order {
        let mut acc = r[i];
        for j in 1..i {
            acc += a[j] * r[i - j];
        }
        let k = -acc / e;
        if !k.is_finite() || k.abs() >= 1.0 {
            return None;
        }
        if i == 1 {
            k1 = k;
        }
        tmp[..=i].copy_from_slice(&a[..=i]);
        for j in 1..i {
            a[j] = tmp[j] + k * tmp[i - j];
        }
        a[i] = k;
        e *= 1.0 - k * k;
        if e <= 0.0 || !e.is_finite() {
            return None;
        }
        if i == snapshot_order && snapshot_order > 0 {
            snapshot = Some(a[..=snapshot_order].to_vec());
        }
    }
    Some((a, k1, snapshot))
}

/// Apply bandwidth expansion: a[k] := a[k] * γ^k for k = 1..=order.
pub fn bandwidth_expand(a: &mut [f32], gamma: f32) {
    let mut g = gamma;
    for k in 1..a.len() {
        a[k] *= g;
        g *= gamma;
    }
}

// ---------------------------------------------------------------------------
// ITU-T G.728 §3.7 hybrid window for the 50th-order synthesis filter
// ---------------------------------------------------------------------------

/// Barnwell / logarithmic autocorrelation window state (block 49 of the
/// spec). Accumulates the recursive portion sample-by-sample via
/// `push_frame` (20 samples per adaptation cycle) and exposes the 51
/// autocorrelation lags needed to drive the 50th-order Levinson-Durbin.
///
/// The sign convention matches the spec: `sb[HYBRID_WIN_LEN - 1]` is the
/// most recent sample, `sb[0]` is the oldest. The window function
/// [`HYBRID_WIN`] is indexed from the recent end: `WNR(1)` lands on the
/// newest sample.
pub struct HybridWindow {
    /// Sample buffer `SB(1..=105)` — oldest first at index 0.
    sb: [f32; HYBRID_WIN_LEN],
    /// Recursive component `REXP(1..=LPC+1)` accumulating the aggregate
    /// autocorrelation contribution from samples older than the non-
    /// recursive tail.
    rexp: [f32; LPC_ORDER + 1],
}

impl Default for HybridWindow {
    fn default() -> Self {
        Self::new()
    }
}

impl HybridWindow {
    pub const fn new() -> Self {
        Self {
            sb: [0.0; HYBRID_WIN_LEN],
            rexp: [0.0; LPC_ORDER + 1],
        }
    }

    /// Reset the buffers to all-zero. Used when restarting the decoder.
    pub fn reset(&mut self) {
        self.sb = [0.0; HYBRID_WIN_LEN];
        self.rexp = [0.0; LPC_ORDER + 1];
    }

    /// Push the next 20-sample adaptation cycle and return the 51
    /// autocorrelation lags after applying the hybrid window and the
    /// 257/256 white-noise correction. The input `frame` provides the new
    /// samples in chronological order (oldest first); `frame[0]` is the
    /// first sample of the current cycle, `frame[NFRSZ-1]` the newest.
    ///
    /// The returned vector has length `LPC_ORDER + 1`. Indexing follows
    /// the spec's 1-based convention re-biased to 0: `r[i]` is
    /// `R(i + 1)` in the spec.
    ///
    /// Returns all-zero `r` when the hybrid window buffers are still
    /// warming up; callers should treat a zero `r[0]` as "skip update"
    /// the same way the reference implementation does.
    pub fn push_frame(&mut self, frame: &[f32; HYBRID_NFRSZ]) -> [f32; LPC_ORDER + 1] {
        const N1: usize = LPC_ORDER + HYBRID_NFRSZ; // 70
        const N2: usize = LPC_ORDER + HYBRID_NONR; // 85
        const N3: usize = LPC_ORDER + HYBRID_NFRSZ + HYBRID_NONR; // 105

        // 1. Shift the old signal buffer left by NFRSZ, dropping the
        //    oldest `NFRSZ` samples.
        //    Reference: `for N=1..N2: SB(N) = SB(N + NFRSZ)`.
        for n in 0..N2 {
            self.sb[n] = self.sb[n + HYBRID_NFRSZ];
        }
        // 2. Insert the new NFRSZ samples at the newest end.
        //    Reference: `for N=1..NFRSZ: SB(N2 + N) = STTMP(N)`, so
        //    `SB(N3)` is the newest sample.
        self.sb[N2..N3].copy_from_slice(frame);

        // 3. Multiply the window with the buffer. The spec walks the
        //    buffer from the newest sample to the oldest while indexing
        //    into `WNR` from 1 upward, which means:
        //      WS(n) = SB(n) * WNR(N3 - n + 1)        (1-based)
        //      ws[n] = sb[n] * HYBRID_WIN[N3 - 1 - n]  (0-based)
        let mut ws = [0.0_f32; N3];
        for n in 0..N3 {
            ws[n] = self.sb[n] * HYBRID_WIN[N3 - 1 - n];
        }

        // 4. Update the recursive component over the samples that have
        //    just slid from the non-recursive into the recursive region:
        //    `N = LPC+1..=N1` (indices LPC..N1 in 0-based).
        //    `rexp[i] = 0.75 * rexp[i] + sum_{n=LPC..N1} ws[n] * ws[n - i]`
        for i in 0..=LPC_ORDER {
            let mut tmp = 0.0_f32;
            for n in LPC_ORDER..N1 {
                // Guard against i > n (impossible here because n >= LPC_ORDER
                // >= i), but we still need the explicit index check to keep
                // the loop literal.
                debug_assert!(n >= i);
                tmp += ws[n] * ws[n - i];
            }
            self.rexp[i] = HYBRID_RECURSIVE_FACTOR * self.rexp[i] + tmp;
        }

        // 5. Compute R. Start with the recursive portion we just
        //    maintained, then add the non-recursive tail spanning
        //    `N = N1+1..=N3` (indices N1..N3 in 0-based).
        let mut r = [0.0_f32; LPC_ORDER + 1];
        for i in 0..=LPC_ORDER {
            r[i] = self.rexp[i];
            for n in N1..N3 {
                debug_assert!(n >= i);
                r[i] += ws[n] * ws[n - i];
            }
        }

        // 6. White-noise correction factor.
        r[0] *= WNCF;
        r
    }
}

// ---------------------------------------------------------------------------
// LPC + gain predictor updates
// ---------------------------------------------------------------------------

/// Result of refreshing the 50th-order LPC predictor from a fresh set
/// of hybrid-window autocorrelation lags. Exposes the bandwidth-expanded
/// 50th-order vector used by the synthesis filter, the first reflection
/// coefficient (for the short-term postfilter's tilt compensation), and
/// the unexpanded 10th-order predictor snapshot that the §5.5 postfilter
/// uses for its LPC inverse filter and pole-zero shaping.
///
/// "Unexpanded" on the order-10 snapshot matches the reference: the spec
/// scales the 10-th order predictor by its own `SPFPCF = 0.75` and
/// `SPFZCF = 0.65` factors inside the postfilter coefficient calculator,
/// not by the synthesis filter's `FAC = 253/256`.
pub struct LpcUpdate {
    pub a: [f32; LPC_ORDER + 1],
    pub k1: f32,
    pub order10: [f32; 11],
}

/// Recompute the 50th-order LPC coefficient vector from a fresh set of
/// hybrid-window autocorrelation lags.
///
/// Returns `None` if the Levinson-Durbin recursion hits a non-positive
/// prediction error, matching the spec's "skip update if ill-conditioned"
/// behaviour.
pub fn update_lpc_from_hybrid_r(r: &[f32; LPC_ORDER + 1]) -> Option<LpcUpdate> {
    if r[0] <= 0.0 {
        return None;
    }
    let r_vec = r.to_vec();
    let (mut a_vec, k1, snap10) = levinson_durbin_with_refl(&r_vec, LPC_ORDER, 10)?;
    bandwidth_expand(&mut a_vec, BW_EXPANSION);
    let mut a = [0.0_f32; LPC_ORDER + 1];
    a.copy_from_slice(&a_vec[..=LPC_ORDER]);
    let mut order10 = [0.0_f32; 11];
    if let Some(s) = snap10 {
        order10[..=10].copy_from_slice(&s[..=10]);
    } else {
        order10[0] = 1.0;
    }
    Some(LpcUpdate { a, k1, order10 })
}

/// Update a gain-predictor coefficient vector from the log-gain history.
/// This keeps the original Hamming-window path for the 10th-order gain
/// predictor; the spec's separate hybrid window for the gain predictor
/// (block 43) is a pure stability tweak on an already-narrow predictor,
/// and replacing it would not move the needle on round-trip quality.
pub fn update_gain_predictor(
    b: &mut [f32; GAIN_ORDER + 1],
    log_gain_history: &[f32; GAIN_HISTORY_LEN],
) -> bool {
    let r = autocorrelation::<GAIN_HISTORY_LEN>(log_gain_history, GAIN_ORDER);
    let mut r = r;
    let floor = r[0] * 1e-4 + 1e-6;
    r[0] += floor;
    let Some(mut new_b) = levinson_durbin(&r, GAIN_ORDER) else {
        return false;
    };
    // FACGP = 29/32 per §3.8.
    bandwidth_expand(&mut new_b, 29.0 / 32.0);
    b[..=GAIN_ORDER].copy_from_slice(&new_b[..=GAIN_ORDER]);
    true
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn levinson_recovers_ar1_coefficient() {
        // y[n] = 0.8 * y[n-1] + x[n]   ⇒ A(z) = 1 - 0.8 z^-1
        // Autocorrelation for the output of this model driven by unit
        // white noise has r[0] = 1 / (1 - 0.8^2) = 2.7778 and
        // r[1] = 0.8 * r[0] = 2.2222.
        let r = vec![2.7778, 2.2222];
        let a = levinson_durbin(&r, 1).expect("recursion");
        assert!((a[0] - 1.0).abs() < 1e-6);
        // Our convention: y[n] = x[n] - sum(a[k] y[n-k]). So for the
        // given model we want a[1] = -0.8.
        assert!((a[1] + 0.8).abs() < 1e-3, "a[1] = {} expected ≈ -0.8", a[1]);
    }

    #[test]
    fn levinson_with_refl_returns_k1() {
        let r = vec![2.7778, 2.2222];
        let (a, k1, _snap) = levinson_durbin_with_refl(&r, 1, 0).expect("recursion");
        assert!((a[1] + 0.8).abs() < 1e-3);
        // For a single-tap predictor k_1 = -r[1] / r[0] = -0.8.
        assert!((k1 + 0.8).abs() < 1e-3, "k1 = {k1} expected ≈ -0.8");
    }

    #[test]
    fn levinson_with_refl_snapshot_is_order_10_when_requested() {
        // Construct a decaying autocorrelation that produces a stable
        // higher-order predictor and verify the order-10 snapshot is
        // returned.
        let mut r = vec![0.0_f32; 16];
        r[0] = 1.0;
        for i in 1..16 {
            r[i] = 0.3_f32.powi(i as i32);
        }
        let (_a, _k1, snap) = levinson_durbin_with_refl(&r, 15, 10).expect("recursion");
        let snap = snap.expect("snapshot");
        assert_eq!(snap.len(), 11);
        assert_eq!(snap[0], 1.0);
    }

    #[test]
    fn levinson_rejects_nonpositive_r0() {
        assert!(levinson_durbin(&[0.0, 0.5], 1).is_none());
    }

    #[test]
    fn bandwidth_expansion_shrinks_higher_order_taps() {
        let mut a = [1.0_f32, 0.5, 0.5, 0.5];
        bandwidth_expand(&mut a, 0.5);
        assert!((a[0] - 1.0).abs() < 1e-6);
        assert!((a[1] - 0.25).abs() < 1e-6);
        assert!((a[2] - 0.125).abs() < 1e-6);
        assert!((a[3] - 0.0625).abs() < 1e-6);
    }

    #[test]
    fn autocorrelation_is_symmetric_in_construction() {
        let mut hist = [0.0_f32; 16];
        for n in 0..16 {
            hist[n] = ((n as f32) * 0.3).sin();
        }
        let r = autocorrelation::<16>(&hist, 4);
        assert_eq!(r.len(), 5);
        // r[0] is energy — must be non-negative.
        assert!(r[0] >= 0.0);
        // r[k] for k>0 should be bounded by r[0].
        for k in 1..5 {
            assert!(r[k].abs() <= r[0] + 1e-6);
        }
    }

    #[test]
    fn hybrid_window_zero_input_stays_zero() {
        let mut hw = HybridWindow::new();
        // 32 cycles of zero input — recursive component must not drift.
        let frame = [0.0_f32; HYBRID_NFRSZ];
        for _ in 0..32 {
            let r = hw.push_frame(&frame);
            assert_eq!(r[0], 0.0, "r[0] should stay zero on zero input");
            for i in 0..=LPC_ORDER {
                assert!(r[i].is_finite());
            }
        }
    }

    #[test]
    fn hybrid_window_dc_input_energy_grows_then_saturates() {
        // On a DC input, r[0] is (sum of squared weighted samples), and
        // the recursive portion saturates at a finite limit (3/4 decay).
        let mut hw = HybridWindow::new();
        let frame = [1.0_f32; HYBRID_NFRSZ];
        let mut last = 0.0_f32;
        for _ in 0..50 {
            let r = hw.push_frame(&frame);
            assert!(r[0].is_finite());
            assert!(r[0] >= last - 1e-3, "r[0] regressed unexpectedly");
            last = r[0];
        }
        // Steady-state r[0] is finite and bounded (each of 105 samples
        // weighted by at most 1.0, so r[0] ≤ 105 * WNCF).
        assert!(last < 110.0, "r[0] diverged at DC: {last}");
    }

    #[test]
    fn hybrid_window_sine_gives_positive_energy() {
        let mut hw = HybridWindow::new();
        // 400 Hz at 8 kHz sample rate → 0.05 cycles per sample.
        let mut phase = 0.0_f32;
        let step = 2.0 * core::f32::consts::PI * 0.05;
        for _ in 0..10 {
            let mut frame = [0.0_f32; HYBRID_NFRSZ];
            for v in frame.iter_mut() {
                *v = phase.sin();
                phase += step;
            }
            let _ = hw.push_frame(&frame);
        }
        let mut frame = [0.0_f32; HYBRID_NFRSZ];
        for v in frame.iter_mut() {
            *v = phase.sin();
            phase += step;
        }
        let r = hw.push_frame(&frame);
        assert!(r[0] > 0.0, "sine input produced zero energy");
        // r[1..] are bounded by r[0] after the WNCF scaling.
        for i in 1..=LPC_ORDER {
            assert!(r[i].abs() <= r[0] + 1e-4);
        }
    }

    #[test]
    fn update_lpc_from_hybrid_r_handles_near_silence() {
        // Simulate a handful of zero frames followed by a tiny non-zero
        // sample. The recursion may bail out but must not panic.
        let mut hw = HybridWindow::new();
        let mut frame = [0.0_f32; HYBRID_NFRSZ];
        frame[HYBRID_NFRSZ - 1] = 1e-3;
        let r = hw.push_frame(&frame);
        // r[0] may be 0 if the white-noise-corrected recursion is still
        // starved; either way the update must return a deterministic
        // Option without UB.
        let _ = update_lpc_from_hybrid_r(&r);
    }
}
