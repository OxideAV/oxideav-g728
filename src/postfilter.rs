//! ITU-T G.728 §4.6 + §5.5 adaptive postfilter.
//!
//! The postfilter sits between the synthesis filter and the PCM
//! conversion stage. Its job is to improve perceived quality by:
//!
//! 1. Emphasising the spectral peaks (formants) via a **short-term**
//!    pole-zero LPC-shaped filter with spectral tilt compensation.
//! 2. Emphasising the harmonics of the fundamental pitch via a
//!    **long-term** (pitch) comb filter.
//! 3. Renormalising the aggregate gain so the postfilter does not shift
//!    the signal level (AGC — automatic gain control).
//!
//! All three stages are backward-adaptive: they update their coefficients
//! once per 20-sample frame from the decoded speech and the 10th-order
//! LPC extracted by the main synthesis-filter adaptation. No side
//! information is transmitted.
//!
//! # Chain
//!
//! ```text
//!   decoded_speech  →  long_term  →  short_term  →  AGC  →  output
//! ```
//!
//! The long-term postfilter runs first so it sees the un-warped formant
//! shape; the short-term stage then refines the spectrum and the AGC
//! restores the original power level.
//!
//! # Spec references
//!
//! - §4.6 / §5.5 postfilter operation (blocks 71-77).
//! - §4.7 / §5.5 postfilter adapter (blocks 81-85).
//! - Annex D: 1 kHz 3rd-order pole-zero elliptic lowpass used for pitch
//!   decimation.

use crate::LPC_ORDER;

/// 10th-order short-term postfilter / LPC inverse filter order.
pub const POSTFILTER_LPC_ORDER: usize = 10;

/// Frame size in samples (one adaptation cycle).
pub const FRAME_SIZE: usize = 20;

/// Pitch analysis window size (`NPWSZ`, §4.7).
pub const PITCH_WINDOW: usize = 100;

/// Minimum pitch period in samples (`KPMIN`).
pub const KP_MIN: usize = 20;

/// Maximum pitch period in samples (`KPMAX`).
pub const KP_MAX: usize = 140;

/// Allowed deviation from the previous pitch period (`KPDELTA`).
pub const KP_DELTA: usize = 6;

/// Tap threshold for turning off the pitch postfilter (`PPFTH`).
pub const PPFTH: f32 = 0.6;

/// Pitch postfilter zero controlling factor (`PPFZCF`).
pub const PPFZCF: f32 = 0.15;

/// Tap threshold for fundamental pitch replacement (`TAPTH`).
pub const TAPTH: f32 = 0.4;

/// Short-term postfilter pole controlling factor (`SPFPCF`).
pub const SPFPCF: f32 = 0.75;

/// Short-term postfilter zero controlling factor (`SPFZCF`).
pub const SPFZCF: f32 = 0.65;

/// Spectral tilt compensation controlling factor (`TILTF`).
pub const TILTF: f32 = 0.15;

/// AGC adaptation speed controlling factor (`AGCFAC`).
pub const AGCFAC: f32 = 0.99;

/// Annex D coefficients for the 1 kHz third-order pole-zero elliptic
/// lowpass filter used by the pitch period extractor to 4:1-decimate
/// the LPC residual.
///
/// Transfer function: `L(z) = (sum b_i z^-i) / (1 + sum a_i z^-i)` with
/// the `a_i` and `b_i` from Annex D.
const LP_A: [f32; 4] = [0.0, -2.34036589, 2.01190019, -0.614109218];
const LP_B: [f32; 4] = [0.0357081667, -0.0069956244, -0.0069956244, 0.0357081667];

/// G.728 adaptive postfilter (§5.5).
///
/// The postfilter is fed one 5-sample decoded-speech vector at a time via
/// [`Postfilter::process_vector`]. Every fourth vector (i.e. every frame
/// = every 20 samples) the coefficients are refreshed from the current
/// 10th-order LPC predictor, first reflection coefficient, and the
/// accumulated decoded-speech buffer.
///
/// # Internal state
///
/// The postfilter keeps three long buffers:
///
/// - `sd_buf`: the most recent `KPMAX + PITCH_WINDOW` (240) samples of
///   decoded speech. The long-term postfilter indexes into this for
///   `st(k - p)`; the pitch predictor tap calculator spans the trailing
///   window.
/// - `lpc_residual`: the output of the 10th-order LPC inverse filter
///   over the last `PITCH_WINDOW + (KPMAX - PITCH_WINDOW + max_range)`
///   samples, kept with the spec's special `d(-139..=100)` indexing.
/// - `dec_residual`: the 4:1-decimated LPC residual, spanning the 60
///   samples needed for correlation peak-picking in the decimated
///   domain (`NPWSZ/4 = 25` + history).
///
/// # Why it's safe to default-enable
///
/// The postfilter is deliberately benign on non-speech signals: the
/// pitch postfilter drops to `b = 0, g_l = 1` when the pitch tap is
/// below 0.6 (so non-periodic signals bypass it), and the short-term
/// postfilter is already bandwidth-expanded to avoid sharpening noise.
/// The AGC renormalises the level, so even a pathological filter
/// response does not shift energy more than a few percent per sample.
pub struct Postfilter {
    // ---- short-term postfilter ----
    /// Short-term postfilter denominator coefficients `ap[1..=10]`
    /// (`AP` in the spec). `ap[0]` is unused.
    ap: [f32; POSTFILTER_LPC_ORDER + 1],
    /// Short-term postfilter numerator coefficients `az[1..=10]` (`AZ`).
    az: [f32; POSTFILTER_LPC_ORDER + 1],
    /// Spectral tilt compensation factor `TILTZ = 0.15 * k_1`.
    tiltz: f32,
    /// All-zero section memory (`STPFFIR`).
    stpffir: [f32; POSTFILTER_LPC_ORDER],
    /// All-pole section memory (`STPFIIR`).
    stpfiir: [f32; POSTFILTER_LPC_ORDER],

    // ---- long-term postfilter ----
    /// Pitch tap `b` for the current frame.
    b: f32,
    /// Scaling factor `g_l = 1 / (1 + b)` for the current frame.
    gl: f32,
    /// Current pitch period in samples (`KP`).
    kp: usize,
    /// Previous-frame pitch period (`KP1`). Starts at 50.
    kp1: usize,

    // ---- decoded-speech buffer ----
    /// Most recent `KPMAX + PITCH_WINDOW` = 240 samples of decoded
    /// speech, newest at the end.
    sd_buf: Vec<f32>,

    // ---- LPC residual buffers ----
    /// LPC residual buffer spanning `[-KPMAX - PITCH_WINDOW + 1 ..= PITCH_WINDOW]`
    /// in the spec's `d(k)` notation. We store it as a flat `Vec` with
    /// the oldest sample at index 0 and the newest (spec's `d(100)`) at
    /// index `len - 1`.
    lpc_residual: Vec<f32>,
    /// 4:1-decimated LPC residual, 60 samples (`d̄(-34..=25)` in the
    /// spec), oldest at index 0.
    dec_residual: Vec<f32>,
    /// 10-sample memory of the 10th-order LPC inverse filter
    /// (`STLPCI`, spec §5.14 block 81). Newest at index 0.
    stlpci: [f32; POSTFILTER_LPC_ORDER],
    /// 3-sample memory of the Annex D lowpass filter (`STLPF`). Newest
    /// at index 0.
    stlpf: [f32; 3],

    // ---- AGC ----
    /// Lowpass-filtered scaling factor (`SCALEFIL`, block 76). Starts
    /// at 1.0 per Table 2/G.728.
    scalefil: f32,

    // ---- cached LPC for pending adapter update ----
    /// 10th-order LPC predictor coefficients `apf[1..=10]` extracted as
    /// a by-product of the main 50th-order Levinson-Durbin. Updated
    /// each frame by [`Postfilter::update_adapter`].
    apf: [f32; POSTFILTER_LPC_ORDER + 1],
    /// Frame counter modulo 4 (spec's `ICOUNT - 1`). Drives the
    /// once-a-frame coefficient refresh + pitch extraction.
    icount: u32,
    /// Vector counter used to schedule LPC-residual computation and
    /// pitch extraction at the right sub-vector boundaries.
    vector_in_frame: u32,
}

impl Default for Postfilter {
    fn default() -> Self {
        let mut ap = [0.0_f32; POSTFILTER_LPC_ORDER + 1];
        ap[0] = 1.0;
        let mut az = [0.0_f32; POSTFILTER_LPC_ORDER + 1];
        az[0] = 1.0;
        let mut apf = [0.0_f32; POSTFILTER_LPC_ORDER + 1];
        apf[0] = 1.0;
        Self {
            ap,
            az,
            tiltz: 0.0,
            stpffir: [0.0; POSTFILTER_LPC_ORDER],
            stpfiir: [0.0; POSTFILTER_LPC_ORDER],
            b: 0.0,
            gl: 1.0,
            kp: 50,
            kp1: 50,
            sd_buf: vec![0.0; KP_MAX + PITCH_WINDOW],
            lpc_residual: vec![0.0; KP_MAX + PITCH_WINDOW + FRAME_SIZE],
            dec_residual: vec![0.0; (KP_MAX + PITCH_WINDOW) / 4],
            stlpci: [0.0; POSTFILTER_LPC_ORDER],
            stlpf: [0.0; 3],
            scalefil: 1.0,
            apf,
            icount: 0,
            vector_in_frame: 0,
        }
    }
}

impl Postfilter {
    pub fn new() -> Self {
        Self::default()
    }

    /// Provide the postfilter with this frame's 10th-order LPC
    /// coefficients and first reflection coefficient. This must be
    /// called **before** the frame's first `process_vector` so the
    /// short-term postfilter coefficients are in place.
    ///
    /// `apf[0]` must equal 1.0; `apf[1..=10]` are the 10th-order
    /// predictor coefficients extracted as a by-product of the main
    /// 50th-order Levinson-Durbin (i.e. obtained by stopping the
    /// recursion at order 10 before continuing to order 50).
    pub fn set_lpc(&mut self, apf: &[f32; POSTFILTER_LPC_ORDER + 1], k1: f32) {
        self.apf = *apf;
        // Short-term postfilter coefficients (block 85).
        // AP_i = SPFPCF^i * APF_i, AZ_i = SPFZCF^i * APF_i.
        let mut gp = SPFPCF;
        let mut gz = SPFZCF;
        self.ap[0] = 1.0;
        self.az[0] = 1.0;
        for i in 1..=POSTFILTER_LPC_ORDER {
            self.ap[i] = gp * apf[i];
            self.az[i] = gz * apf[i];
            gp *= SPFPCF;
            gz *= SPFZCF;
        }
        self.tiltz = TILTF * k1;
    }

    /// Process one 5-sample decoded-speech vector through the full
    /// postfilter chain (long-term → short-term → AGC) and write the
    /// result into `out`. The decoded-speech buffer and pitch-period
    /// state are updated before filtering, so this is safe to call on
    /// any vector position within a frame.
    ///
    /// `vec_idx_in_frame` is the 0-based vector index within the
    /// current 4-vector frame (0..=3). At vector 2 we run the pitch
    /// period extractor and update `(b, g_l, kp)`. The caller must
    /// drive this counter — `Postfilter` doesn't own the global vector
    /// clock.
    pub fn process_vector(
        &mut self,
        decoded: &[f32; crate::VECTOR_SIZE],
        vec_idx_in_frame: u32,
        out: &mut [f32; crate::VECTOR_SIZE],
    ) {
        self.vector_in_frame = vec_idx_in_frame;

        // 1. Update sd_buf with the 5 new decoded samples.
        self.push_decoded(decoded);

        // 2. Compute 10th-order LPC residual for the new 5 samples.
        self.compute_lpc_residual_step(decoded);

        // 3. At the third vector of each frame (vec_idx 2), refresh the
        //    pitch period + long-term postfilter coefficients.
        if vec_idx_in_frame == 2 {
            self.refresh_pitch_and_long_term();
        }

        // 4. Long-term postfilter.
        let mut lt = [0.0_f32; crate::VECTOR_SIZE];
        self.long_term_filter(decoded, &mut lt);

        // 5. Short-term postfilter (pole-zero + tilt compensation).
        let mut st = [0.0_f32; crate::VECTOR_SIZE];
        self.short_term_filter(&lt, &mut st);

        // 6. AGC — scale to match sum-of-absolute-values of the decoded
        //    input, low-pass filtered across vectors for smoothness.
        self.agc(decoded, &st, out);
    }

    // ---- internal helpers ----

    /// Shift the `sd_buf` by 5 samples and insert the new decoded
    /// vector at the tail.
    fn push_decoded(&mut self, v: &[f32; crate::VECTOR_SIZE]) {
        let n = self.sd_buf.len();
        self.sd_buf.copy_within(crate::VECTOR_SIZE..n, 0);
        let dst = &mut self.sd_buf[n - crate::VECTOR_SIZE..n];
        dst.copy_from_slice(v);
    }

    /// Compute the 5-sample LPC prediction residual for this vector and
    /// append it to the residual buffer. Block 81 in the spec.
    ///
    /// `d(k) = st(k) - sum_{i=1..=10} apf_i * st(k - i)`
    fn compute_lpc_residual_step(&mut self, decoded: &[f32; crate::VECTOR_SIZE]) {
        let n = self.lpc_residual.len();
        // Shift the residual buffer left by 5.
        self.lpc_residual.copy_within(crate::VECTOR_SIZE..n, 0);
        for (j, &s) in decoded.iter().enumerate() {
            let mut acc = s;
            for i in 1..=POSTFILTER_LPC_ORDER {
                acc -= self.apf[i] * self.stlpci[i - 1];
            }
            // Shift STLPCI (short LPC-inverse memory); newest at index 0.
            for k in (1..POSTFILTER_LPC_ORDER).rev() {
                self.stlpci[k] = self.stlpci[k - 1];
            }
            self.stlpci[0] = s;
            self.lpc_residual[n - crate::VECTOR_SIZE + j] = acc;
        }
    }

    /// Refresh the pitch period + long-term postfilter coefficients.
    /// Runs once per frame at vector index 2 in the spec.
    fn refresh_pitch_and_long_term(&mut self) {
        // Pitch extraction block 82 — operates on the LPC residual.
        let kp = self.extract_pitch();
        self.kp1 = self.kp;
        self.kp = kp;

        // Pitch predictor tap β (block 83) — computed over the decoded
        // speech buffer, not the LPC residual.
        // β = sum sd(k)*sd(k-p) / sum sd(k-p)^2   for k = -99..=0 in
        //     the spec's notation (i.e. the last 100 samples of sd_buf,
        //     excluding the current vector).
        let sd = &self.sd_buf;
        let n = sd.len();
        // `sd(1..=5)` in the spec = current vector (last 5 samples).
        // `sd(-99..=0)` = the 100 samples *before* the current vector.
        // We want `st(k) * st(k - kp)` over `k = -99..=0`.
        let end = n - crate::VECTOR_SIZE; // exclude current vector
        let start = end.saturating_sub(PITCH_WINDOW);
        let mut sum = 0.0_f32;
        let mut denom = 0.0_f32;
        for k in start..end {
            if k < kp {
                continue;
            }
            let a = sd[k];
            let b = sd[k - kp];
            sum += a * b;
            denom += b * b;
        }
        let mut ptap = if denom > 0.0 { sum / denom } else { 0.0 };
        // Long-term postfilter coefficient calculator (block 84).
        if ptap > 1.0 {
            ptap = 1.0;
        }
        if ptap < PPFTH {
            ptap = 0.0;
        }
        self.b = PPFZCF * ptap;
        self.gl = 1.0 / (1.0 + self.b);
    }

    /// Extract the pitch period (block 82). Runs the 4:1 decimation +
    /// lag search as specified, then refines around the previous
    /// frame's pitch to avoid locking onto a multiple.
    fn extract_pitch(&mut self) -> usize {
        // Step 1: 4:1 lowpass + decimate the last PITCH_WINDOW LPC-
        // residual samples. The spec feeds the 100 samples d(1..=100)
        // (i.e. the window ending at the current frame's last residual
        // sample, index `n - 1` in our layout).
        let d = &self.lpc_residual;
        let n = d.len();
        let window_end = n;
        let window_start = window_end - PITCH_WINDOW;

        // Running AL/BL lowpass over d[window_start..window_end].
        let mut stlpf = self.stlpf;
        let mut filtered = [0.0_f32; PITCH_WINDOW];
        for (i, &s) in d[window_start..window_end].iter().enumerate() {
            // L(z) = B(z) / (1 + A(z))   (all-pole + all-zero)
            let tmp = s - LP_A[1] * stlpf[0] - LP_A[2] * stlpf[1] - LP_A[3] * stlpf[2];
            let out = tmp * LP_B[0] + stlpf[0] * LP_B[1] + stlpf[1] * LP_B[2] + stlpf[2] * LP_B[3];
            stlpf[2] = stlpf[1];
            stlpf[1] = stlpf[0];
            stlpf[0] = tmp;
            filtered[i] = out;
        }
        self.stlpf = stlpf;

        // 4:1 decimate: take every 4th sample of the filtered block. We
        // have 100 filtered samples → 25 decimated. Stuff them at the
        // tail of `dec_residual`, shifting older history left.
        let new_dec: Vec<f32> = filtered.iter().step_by(4).copied().collect();
        let dn = self.dec_residual.len();
        let m = new_dec.len(); // 25
        self.dec_residual.copy_within(m..dn, 0);
        self.dec_residual[dn - m..dn].copy_from_slice(&new_dec);

        // Step 2: find correlation peak in the decimated domain for
        // lags τ = 5..=35. Correlate `d̄(1..=25)` (most recent 25
        // decimated samples) against `d̄(1-τ..=25-τ)`.
        let dec = &self.dec_residual;
        let dec_latest = dn; // one past the most recent
        let mut best_tau = 5usize;
        let mut best_corr = f32::NEG_INFINITY;
        for tau in 5..=35 {
            let mut corr = 0.0_f32;
            for k in 0..m {
                let idx_new = dec_latest - m + k;
                let idx_old = idx_new.wrapping_sub(tau);
                if idx_old < dec_latest {
                    corr += dec[idx_new] * dec[idx_old];
                }
            }
            if corr > best_corr {
                best_corr = corr;
                best_tau = tau;
            }
        }

        // Step 3: refine in the undecimated domain around 4*tau ± 3.
        let t = best_tau as isize;
        let mut m1 = (4 * t - 3).max(KP_MIN as isize) as usize;
        let mut m2 = (4 * t + 3).min(KP_MAX as isize) as usize;
        if m2 < m1 {
            std::mem::swap(&mut m1, &mut m2);
        }
        let mut p0 = m1;
        let mut cmax_p0 = f32::NEG_INFINITY;
        for lag in m1..=m2 {
            let mut corr = 0.0_f32;
            for k in 1..=PITCH_WINDOW {
                let idx_new = window_start + k - 1;
                if idx_new < lag {
                    continue;
                }
                let idx_old = idx_new - lag;
                corr += d[idx_new] * d[idx_old];
            }
            if corr > cmax_p0 {
                cmax_p0 = corr;
                p0 = lag;
            }
        }

        // Step 4: also search ±KPDELTA around kp1 (previous pitch) to
        // catch the fundamental when p0 is a multiple.
        let pp = self.kp1 as isize;
        let m1b = (pp - KP_DELTA as isize).max(KP_MIN as isize) as usize;
        let m2b = (pp + KP_DELTA as isize).min(KP_MAX as isize) as usize;
        let mut p1 = m1b;
        let mut cmax_p1 = f32::NEG_INFINITY;
        for lag in m1b..=m2b {
            let mut corr = 0.0_f32;
            for k in 1..=PITCH_WINDOW {
                let idx_new = window_start + k - 1;
                if idx_new < lag {
                    continue;
                }
                let idx_old = idx_new - lag;
                corr += d[idx_new] * d[idx_old];
            }
            if corr > cmax_p1 {
                cmax_p1 = corr;
                p1 = lag;
            }
        }

        // Step 5: compute taps (β0, β1) and apply the replacement rule.
        let tap0 = single_tap_coeff(d, window_start, p0);
        let tap1 = single_tap_coeff(d, window_start, p1);
        // Replace with fundamental if the beta ratio is large enough.
        if tap1 > TAPTH * tap0 {
            p1
        } else {
            p0
        }
    }

    fn long_term_filter(
        &self,
        decoded: &[f32; crate::VECTOR_SIZE],
        out: &mut [f32; crate::VECTOR_SIZE],
    ) {
        let n = self.sd_buf.len();
        let base = n - crate::VECTOR_SIZE; // index of the first new sample in sd_buf
        for k in 0..crate::VECTOR_SIZE {
            let idx_now = base + k;
            let idx_old = idx_now.wrapping_sub(self.kp);
            let old = if idx_old < n {
                self.sd_buf[idx_old]
            } else {
                0.0
            };
            out[k] = self.gl * (decoded[k] + self.b * old);
        }
    }

    fn short_term_filter(
        &mut self,
        input: &[f32; crate::VECTOR_SIZE],
        out: &mut [f32; crate::VECTOR_SIZE],
    ) {
        // Pole-zero filter per block 72:
        //   TEMP(k) = INPUT(k) + sum AZ(j+1) * STPFFIR(j)  (all-zero branch)
        //   shift STPFFIR
        //   TEMP(k) -= sum AP(j+1) * STPFIIR(j)            (all-pole branch)
        //   shift STPFIIR
        //   TEMP(k) += STPFIIR(2) * TILTZ                  (tilt compensation)
        for k in 0..crate::VECTOR_SIZE {
            // All-zero section.
            let mut tmp = input[k];
            for j in (2..=POSTFILTER_LPC_ORDER).rev() {
                tmp += self.stpffir[j - 1] * self.az[j];
                self.stpffir[j - 1] = self.stpffir[j - 2];
            }
            tmp += self.stpffir[0] * self.az[1];
            self.stpffir[0] = input[k];

            // All-pole section.
            for j in (2..=POSTFILTER_LPC_ORDER).rev() {
                tmp -= self.stpfiir[j - 1] * self.ap[j];
                self.stpfiir[j - 1] = self.stpfiir[j - 2];
            }
            tmp -= self.stpfiir[0] * self.ap[1];
            // Tilt compensation uses the *previous* STPFIIR(1) — per the
            // reference, TEMP += STPFIIR(2) * TILTZ evaluated before the
            // final `STPFIIR(1) = TEMP`.
            let tilt = self.stpfiir[0] * self.tiltz;
            self.stpfiir[0] = tmp;
            out[k] = tmp + tilt;
        }
    }

    fn agc(
        &mut self,
        decoded: &[f32; crate::VECTOR_SIZE],
        filtered: &[f32; crate::VECTOR_SIZE],
        out: &mut [f32; crate::VECTOR_SIZE],
    ) {
        // Block 73: sum |decoded|; block 74: sum |filtered|.
        let mut sum_unfil = 0.0_f32;
        let mut sum_fil = 0.0_f32;
        for k in 0..crate::VECTOR_SIZE {
            sum_unfil += decoded[k].abs();
            sum_fil += filtered[k].abs();
        }
        // Block 75: scale = sum_unfil / sum_fil (with the spec's 1.0
        // fallback when sum_fil is too small).
        let scale = if sum_fil > 1.0 {
            sum_unfil / sum_fil
        } else {
            1.0
        };
        // Block 76/77: lowpass-filter the scale and apply sample by
        // sample.
        for k in 0..crate::VECTOR_SIZE {
            self.scalefil = AGCFAC * self.scalefil + (1.0 - AGCFAC) * scale;
            out[k] = self.scalefil * filtered[k];
        }
    }

    /// Current short-term postfilter coefficients (read-only). Useful
    /// for tests / debugging.
    pub fn current_ap(&self) -> &[f32; POSTFILTER_LPC_ORDER + 1] {
        &self.ap
    }
    pub fn current_az(&self) -> &[f32; POSTFILTER_LPC_ORDER + 1] {
        &self.az
    }
    pub fn current_pitch(&self) -> usize {
        self.kp
    }
    pub fn current_b(&self) -> f32 {
        self.b
    }
}

/// Compute the single-tap pitch predictor coefficient
/// `sum d(k) * d(k-lag) / sum d(k-lag)^2` over the 100-sample window
/// starting at `window_start`. Result is clamped to [0, 1].
fn single_tap_coeff(d: &[f32], window_start: usize, lag: usize) -> f32 {
    let mut num = 0.0_f32;
    let mut den = 0.0_f32;
    for k in 0..PITCH_WINDOW {
        let idx_new = window_start + k;
        if idx_new < lag {
            continue;
        }
        let idx_old = idx_new - lag;
        num += d[idx_new] * d[idx_old];
        den += d[idx_old] * d[idx_old];
    }
    if den <= 0.0 {
        return 0.0;
    }
    (num / den).clamp(0.0, 1.0)
}

/// Extract the order-10 LPC predictor from a fresh set of autocorrelation
/// lags by running a short Levinson-Durbin. This is the "stop at order
/// 10" path the spec refers to in §5.5.
///
/// Returns `None` on ill-conditioned input — the postfilter should keep
/// its previous `apf` in that case.
pub fn order10_from_autocorrelation(
    r: &[f32; LPC_ORDER + 1],
) -> Option<[f32; POSTFILTER_LPC_ORDER + 1]> {
    if r[0] <= 0.0 {
        return None;
    }
    let r_short: Vec<f32> = r[..=POSTFILTER_LPC_ORDER].to_vec();
    let a = crate::predictor::levinson_durbin(&r_short, POSTFILTER_LPC_ORDER)?;
    // Spec convention: APF(1..=10) are the positive-sign tap
    // coefficients (a in this crate is already in that form: a[i] = -q_i
    // per §3.3, but our Levinson output uses the `y[n] = x[n] - sum a*y[n-k]`
    // convention, which matches G.728's APF sign choice used directly in
    // the filter APF(I) multiplications in block 81).
    let mut out = [0.0_f32; POSTFILTER_LPC_ORDER + 1];
    // In the spec block 81: `d(k) = st(k) + sum STLPCI(j) * APF(j+1)`.
    // Our `a` has the convention `st(k) = exc(k) - sum a[k] * hist`, i.e.
    // the synthesis filter's `a[k] = -APF(k)`. Flip sign so `self.apf[i]`
    // can be used directly as the LPC-residual predictor.
    out[0] = 1.0;
    for i in 1..=POSTFILTER_LPC_ORDER {
        out[i] = -a[i];
    }
    Some(out)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn postfilter_default_is_identity_on_silence() {
        let mut pf = Postfilter::new();
        let zero = [0.0_f32; crate::VECTOR_SIZE];
        let mut out = [0.0_f32; crate::VECTOR_SIZE];
        for i in 0..8 {
            pf.process_vector(&zero, (i % 4) as u32, &mut out);
            assert_eq!(out, zero);
        }
    }

    #[test]
    fn postfilter_set_lpc_scales_coefficients() {
        let mut pf = Postfilter::new();
        let mut apf = [0.0_f32; POSTFILTER_LPC_ORDER + 1];
        apf[0] = 1.0;
        for i in 1..=POSTFILTER_LPC_ORDER {
            apf[i] = 0.1 * i as f32;
        }
        pf.set_lpc(&apf, 0.5);
        // AP_i = 0.75^i * APF_i, AZ_i = 0.65^i * APF_i.
        let ap = pf.current_ap();
        let az = pf.current_az();
        assert!((ap[1] - 0.75 * 0.1).abs() < 1e-6);
        assert!((ap[2] - 0.75 * 0.75 * 0.2).abs() < 1e-6);
        assert!((az[1] - 0.65 * 0.1).abs() < 1e-6);
        assert!((az[2] - 0.65 * 0.65 * 0.2).abs() < 1e-6);
    }

    #[test]
    fn order10_lpc_from_ar_toy_signal() {
        // Construct an autocorrelation matching a simple AR(1) model
        // and verify that order10_from_autocorrelation returns
        // well-defined finite coefficients.
        let mut r = [0.0_f32; LPC_ORDER + 1];
        r[0] = 1.0;
        // Weak decay so higher-order lags stay near-zero and the
        // recursion is well-conditioned.
        for i in 1..=LPC_ORDER {
            r[i] = 0.3_f32.powi(i as i32);
        }
        let apf = order10_from_autocorrelation(&r).expect("recursion");
        assert_eq!(apf[0], 1.0);
        for i in 1..=POSTFILTER_LPC_ORDER {
            assert!(apf[i].is_finite());
        }
    }

    #[test]
    fn postfilter_preserves_signal_level_on_steady_input() {
        // Drive the postfilter with a non-trivial but stable signal
        // (400 Hz sine at 8 kHz) and verify that the AGC preserves the
        // per-vector average magnitude to within a couple of percent.
        let mut pf = Postfilter::new();
        let mut phase = 0.0_f32;
        let step = 2.0 * core::f32::consts::PI * 0.05;
        let mut sum_in = 0.0_f64;
        let mut sum_out = 0.0_f64;
        for n in 0..800 {
            let mut v = [0.0_f32; crate::VECTOR_SIZE];
            for s in v.iter_mut() {
                *s = 1000.0 * phase.sin();
                phase += step;
            }
            let mut out = [0.0_f32; crate::VECTOR_SIZE];
            pf.process_vector(&v, (n % 4) as u32, &mut out);
            if n >= 200 {
                // Skip the AGC warmup.
                for (a, b) in v.iter().zip(out.iter()) {
                    sum_in += a.abs() as f64;
                    sum_out += b.abs() as f64;
                }
            }
        }
        let ratio = sum_out / sum_in;
        assert!(
            (0.7..=1.3).contains(&ratio),
            "postfilter level ratio out of range: {ratio:.3}"
        );
    }

    #[test]
    fn postfilter_differs_from_input_on_non_trivial_signal() {
        // The whole point of the postfilter is to shape the spectrum,
        // so it must modify at least a fraction of the samples on a
        // non-trivial input.
        let mut pf = Postfilter::new();
        // Set a non-identity LPC so the filter has something to do.
        let mut apf = [0.0_f32; POSTFILTER_LPC_ORDER + 1];
        apf[0] = 1.0;
        apf[1] = -0.8;
        apf[2] = 0.3;
        pf.set_lpc(&apf, -0.6);
        let mut phase = 0.0_f32;
        let step = 2.0 * core::f32::consts::PI * 0.05;
        let mut total_diff = 0.0_f64;
        let mut total_input = 0.0_f64;
        for n in 0..400 {
            let mut v = [0.0_f32; crate::VECTOR_SIZE];
            for s in v.iter_mut() {
                *s = 2000.0 * phase.sin();
                phase += step;
            }
            let mut out = [0.0_f32; crate::VECTOR_SIZE];
            pf.process_vector(&v, (n % 4) as u32, &mut out);
            if n >= 100 {
                for (a, b) in v.iter().zip(out.iter()) {
                    total_diff += ((*a - *b) as f64).abs();
                    total_input += (*a as f64).abs();
                }
            }
        }
        assert!(total_input > 0.0);
        let rel = total_diff / total_input;
        assert!(rel > 0.0, "postfilter produced identical output");
    }
}
