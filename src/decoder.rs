//! G.728 LD-CELP decoder core state + first-cut synthesis loop.
//!
//! This module implements a first-cut G.728 decoder:
//!
//! 1. Bit-unpack: 10-bit codebook indices from the packed stream.
//! 2. Excitation: 128-entry shape codebook × 8-entry gain codebook
//!    (sign comes from the extra bit in the index).
//! 3. Backward-adaptive LPC (50th order): autocorrelation +
//!    Levinson-Durbin over the recent synthesis history, refreshed every
//!    4 vectors (2.5 ms).
//! 4. Synthesis: all-pole IIR filter fed by the per-vector excitation.
//! 5. Backward-adaptive log-gain prediction (10th order): tracks the
//!    excitation energy trajectory in the log domain.
//!
//! The autocorrelation windowing follows the spec's §3.7 hybrid
//! (Barnwell / logarithmic) window: a fixed 35-sample non-recursive
//! tail plus a decaying recursive portion with `alpha^{2L} = 3/4`, using
//! the 105-sample window table from Annex A. Bandwidth expansion uses
//! the spec's `FAC = 253/256` factor.
//!
//! The adaptive long-term (pitch) + short-term postfilters from §5.5
//! are wired through [`crate::postfilter`] and are enabled by default.
//! Callers that want the raw synthesis output can disable the postfilter
//! via [`G728Decoder::set_postfilter_enabled`].

use oxideav_codec::Decoder;
use oxideav_core::{
    AudioFrame, CodecId, CodecParameters, Error, Frame, Packet, Result, SampleFormat, TimeBase,
};

use crate::bitreader::{BitReader, UnpackedIndex};
use crate::postfilter::Postfilter;
use crate::predictor::{
    update_gain_predictor, update_lpc_from_hybrid_r, HybridWindow, GAIN_HISTORY_LEN, HYBRID_NFRSZ,
};
use crate::tables::{GAIN_CB, SHAPE_CB};
use crate::{CODEC_ID_STR, GAIN_ORDER, INDEX_BITS, LPC_ORDER, SAMPLE_RATE, VECTOR_SIZE};

/// Number of vectors between backward-adaptation refreshes (G.728 §3.7:
/// LPC re-estimation every 4 vectors = 20 samples = 2.5 ms).
pub const VECTORS_PER_BLOCK: u32 = 4;

// ---------------------------------------------------------------------------
// Core decoder state
// ---------------------------------------------------------------------------

/// Backward-adaptive 50th-order LPC synthesis filter state.
///
/// `a[0] ≡ 1.0` and `a[1..=LPC_ORDER]` are the AR predictor taps. The
/// synthesis filter realises `1 / A(z)` as:
///
/// ```text
///   y[n] = x[n] - sum_{k=1..=50} a[k] * y[n-k]
/// ```
///
/// The tap vector is refreshed every `VECTORS_PER_BLOCK` vectors by
/// pushing the new adaptation-cycle samples into [`HybridWindow`], then
/// running Levinson-Durbin on the resulting autocorrelation lags.
pub struct LpcPredictor {
    /// Current LPC synthesis coefficients `a[1..=50]` (a[0] ≡ 1.0).
    pub a: [f32; LPC_ORDER + 1],
    /// Delay line of past synthesised samples (most recent at index 0).
    pub history: [f32; LPC_ORDER],
    /// G.728 §3.7 hybrid-window state. Samples are pushed into it in
    /// 20-sample adaptation cycles (`HYBRID_NFRSZ`); each push returns
    /// the 51 autocorrelation lags used to refresh `a`.
    pub hybrid: HybridWindow,
    /// Pending 20-sample cycle being accumulated; flushed into `hybrid`
    /// every `VECTORS_PER_BLOCK` vectors.
    pub frame_buf: [f32; HYBRID_NFRSZ],
    /// Number of samples already written into `frame_buf`.
    pub frame_fill: usize,
    /// Number of vectors processed since the last coefficient update.
    pub vectors_since_update: u32,
    /// First reflection coefficient `k_1` from the latest successful
    /// Levinson-Durbin update. Consumed by the short-term postfilter's
    /// spectral-tilt compensation (§4.6, `mu = 0.15 * k_1`).
    pub k1: f32,
    /// 10th-order predictor snapshot taken during the 50th-order
    /// recursion (§5.5 "by-product of the 50th-order Levinson-Durbin").
    /// Consumed by the postfilter's LPC inverse filter and pole-zero
    /// coefficients. `apf[0] ≡ 1.0`.
    pub apf: [f32; 11],
    /// Whether the latest refresh produced usable coefficients. Lets
    /// downstream consumers (postfilter) skip updates on the first
    /// few frames, where the recursion may be starved.
    pub last_update_ok: bool,
}

impl Default for LpcPredictor {
    fn default() -> Self {
        let mut a = [0.0_f32; LPC_ORDER + 1];
        a[0] = 1.0;
        let mut apf = [0.0_f32; 11];
        apf[0] = 1.0;
        Self {
            a,
            history: [0.0; LPC_ORDER],
            hybrid: HybridWindow::new(),
            frame_buf: [0.0; HYBRID_NFRSZ],
            frame_fill: 0,
            vectors_since_update: 0,
            k1: 0.0,
            apf,
            last_update_ok: false,
        }
    }
}

impl LpcPredictor {
    pub fn new() -> Self {
        Self::default()
    }

    /// Apply the all-pole synthesis filter to one 5-sample excitation
    /// vector, producing 5 reconstructed speech samples. Both
    /// `history` (short delay line) and the hybrid-window cycle buffer
    /// are advanced by the output.
    pub fn synthesise(&mut self, excitation: &[f32; VECTOR_SIZE], out: &mut [f32; VECTOR_SIZE]) {
        for n in 0..VECTOR_SIZE {
            let mut acc = excitation[n];
            for k in 1..=LPC_ORDER {
                acc -= self.a[k] * self.history[k - 1];
            }
            // Hard clip to prevent runaway if the filter is briefly
            // ill-conditioned (can only happen if update_lpc_from_hybrid_r
            // produces a marginally-stable filter and is hit with large
            // excitation).
            let y = acc.clamp(-1.0e4, 1.0e4);
            out[n] = y;
            // Shift short history: newest sample at index 0.
            for k in (1..LPC_ORDER).rev() {
                self.history[k] = self.history[k - 1];
            }
            self.history[0] = y;
            // Accumulate the 20-sample hybrid-window cycle in
            // chronological order so `push_frame` sees the oldest
            // sample at index 0.
            if self.frame_fill < HYBRID_NFRSZ {
                self.frame_buf[self.frame_fill] = y;
                self.frame_fill += 1;
            }
        }
        self.vectors_since_update = self.vectors_since_update.wrapping_add(1);
    }

    /// Feed the pending 20-sample adaptation-cycle buffer through the
    /// hybrid window and re-estimate `a` from the resulting
    /// autocorrelation lags. Returns `true` if the update succeeded; the
    /// filter is left unchanged on failure (matching the spec's
    /// "ill-conditioned → keep previous coefficients" behaviour).
    pub fn refresh_coefficients(&mut self) -> bool {
        // `refresh_coefficients` is called only at cycle boundaries, so
        // `frame_buf` should always hold exactly NFRSZ samples. If for
        // some reason fewer are present (e.g. tests that call
        // `refresh` manually), zero-pad the trailing samples rather
        // than short-cycling the hybrid window.
        for i in self.frame_fill..HYBRID_NFRSZ {
            self.frame_buf[i] = 0.0;
        }
        let r = self.hybrid.push_frame(&self.frame_buf);
        self.frame_fill = 0;
        self.vectors_since_update = 0;
        if let Some(update) = update_lpc_from_hybrid_r(&r) {
            self.a = update.a;
            self.k1 = update.k1;
            self.apf = update.order10;
            self.last_update_ok = true;
            true
        } else {
            self.last_update_ok = false;
            false
        }
    }
}

/// Backward-adaptive 10th-order log-gain predictor (§3.9 of G.728).
///
/// Predicts the log-domain excitation gain from the 10 most recent
/// log-gains. Like the LPC predictor it is updated every
/// `VECTORS_PER_BLOCK` vectors, but from the gain trajectory rather
/// than the synthesis signal.
pub struct GainPredictor {
    /// Prediction coefficients `b[1..=GAIN_ORDER]` (b[0] ≡ 1.0).
    pub b: [f32; GAIN_ORDER + 1],
    /// Short delay line for prediction (newest at index 0).
    pub history: [f32; GAIN_ORDER],
    /// Longer history window for the Levinson-Durbin update.
    pub analysis_history: [f32; GAIN_HISTORY_LEN],
    /// Most recently predicted log gain (linear dB-ish units).
    pub last_log_gain: f32,
    /// Vectors since the last coefficient update.
    pub vectors_since_update: u32,
}

impl Default for GainPredictor {
    fn default() -> Self {
        let mut b = [0.0_f32; GAIN_ORDER + 1];
        b[0] = 1.0;
        Self {
            b,
            history: [0.0; GAIN_ORDER],
            analysis_history: [0.0; GAIN_HISTORY_LEN],
            last_log_gain: 0.0,
            vectors_since_update: 0,
        }
    }
}

impl GainPredictor {
    pub fn new() -> Self {
        Self::default()
    }

    /// Produce the next predicted log-gain (log base e). The actual
    /// excitation gain is recovered as `exp(log_gain)`.
    ///
    /// The spec's predictor is an AR model of the *mean-removed* log
    /// gain — that way a slowly-varying signal (near-DC log-gain) is
    /// tracked by the mean term and the AR coefficients only need to
    /// model deviations. We mirror that here: subtract the running
    /// mean from the history, apply `b` to the residual, add the mean
    /// back. Without this split the `b` vector from Levinson-Durbin
    /// on a near-constant history degenerates (reflection coefficient
    /// hits ±1 and the recursion bails out, leaving `b[1..]` at zero),
    /// which would pin the predicted log gain to zero forever.
    pub fn predict(&mut self) -> f32 {
        let mean = {
            let mut s = 0.0_f32;
            for k in 0..GAIN_ORDER {
                s += self.history[k];
            }
            s / GAIN_ORDER as f32
        };
        let mut acc = 0.0_f32;
        for k in 1..=GAIN_ORDER {
            acc -= self.b[k] * (self.history[k - 1] - mean);
        }
        let predicted = mean + acc;
        // Clamp against wild predictions (e.g., ±20 ≈ 9 dec e-folds).
        self.last_log_gain = predicted.clamp(-6.0, 6.0);
        self.last_log_gain
    }

    /// Slide histories, inserting the latest observed log-gain.
    pub fn push(&mut self, log_gain: f32) {
        let g = log_gain.clamp(-6.0, 6.0);
        for k in (1..GAIN_ORDER).rev() {
            self.history[k] = self.history[k - 1];
        }
        self.history[0] = g;
        for k in (1..GAIN_HISTORY_LEN).rev() {
            self.analysis_history[k] = self.analysis_history[k - 1];
        }
        self.analysis_history[0] = g;
        self.vectors_since_update = self.vectors_since_update.wrapping_add(1);
    }

    /// Re-estimate `b` from `analysis_history`.
    pub fn refresh_coefficients(&mut self) -> bool {
        let ok = update_gain_predictor(&mut self.b, &self.analysis_history);
        self.vectors_since_update = 0;
        ok
    }
}

/// Aggregate decoder state — LPC + gain predictor + vector counters.
pub struct G728State {
    pub lpc: LpcPredictor,
    pub gain: GainPredictor,
    /// Counter of decoded vectors since the start of the stream.
    pub vector_count: u64,
}

impl Default for G728State {
    fn default() -> Self {
        Self {
            lpc: LpcPredictor::new(),
            gain: GainPredictor::new(),
            vector_count: 0,
        }
    }
}

impl G728State {
    pub fn new() -> Self {
        Self::default()
    }

    /// Decode one raw 10-bit index into an excitation vector.
    ///
    /// This looks up the shape codebook row indexed by the top 7 bits,
    /// applies the sign bit, and scales by the 3-bit gain codebook entry
    /// combined with the current backward-predicted log-gain. The
    /// returned vector is *pre*-synthesis; feed it through
    /// `LpcPredictor::synthesise` to obtain PCM.
    pub fn excitation_from_index(&self, raw: u16) -> [f32; VECTOR_SIZE] {
        let idx = UnpackedIndex::from_raw(raw);
        let shape = &SHAPE_CB[idx.shape_index as usize];
        let mag = GAIN_CB[idx.gain_mag as usize];
        let sign: f32 = if idx.sign != 0 { -1.0 } else { 1.0 };
        // Backward-adaptive gain: exp(last_log_gain) scales the codebook
        // magnitude. last_log_gain is clamped to ±6 nats ≈ ±52 dB.
        let adaptive = self.gain.last_log_gain.exp();
        let scale = sign * mag * adaptive;
        [
            shape[0] * scale,
            shape[1] * scale,
            shape[2] * scale,
            shape[3] * scale,
            shape[4] * scale,
        ]
    }

    /// Decode a single 10-bit index into 5 PCM f32 samples and advance
    /// all state (LPC history, gain predictor, adaptation counters).
    ///
    /// Per §3.9 of the spec, the predicted log gain is refreshed **once
    /// per vector** (using the 10th-order gain predictor's current `b`
    /// coefficients over the running log-gain history) *before* the
    /// excitation is scaled. The coefficients themselves are re-estimated
    /// every 4 vectors via Levinson-Durbin. This split — fast per-vector
    /// prediction + slow coefficient adaptation — is what lets the log-
    /// gain track the signal envelope in real time without transmitting
    /// any side information.
    pub fn decode_vector(&mut self, raw: u16, out: &mut [f32; VECTOR_SIZE]) {
        // Per-vector log-gain prediction. This updates `last_log_gain`,
        // which `excitation_from_index` reads for the adaptive scale.
        self.gain.predict();

        let excitation = self.excitation_from_index(raw);

        // Synthesise: run excitation through 1 / A(z).
        self.lpc.synthesise(&excitation, out);

        // Update log-gain history from the magnitude of this excitation.
        // We use the log of the excitation RMS so the predictor tracks
        // the envelope in the log domain.
        let mut ss = 0.0_f32;
        for n in 0..VECTOR_SIZE {
            ss += excitation[n] * excitation[n];
        }
        let rms = (ss / VECTOR_SIZE as f32).sqrt();
        // ln(max(rms, eps))
        let log_g = rms.max(1.0e-6).ln();
        self.gain.push(log_g);

        self.vector_count = self.vector_count.wrapping_add(1);

        // Backward-adaptive coefficient refresh every 4 vectors.
        if self.lpc.vectors_since_update >= VECTORS_PER_BLOCK {
            self.lpc.refresh_coefficients();
        }
        if self.gain.vectors_since_update >= VECTORS_PER_BLOCK {
            self.gain.refresh_coefficients();
        }
    }
}

// ---------------------------------------------------------------------------
// Decoder trait wiring
// ---------------------------------------------------------------------------

pub fn make_decoder(params: &CodecParameters) -> Result<Box<dyn Decoder>> {
    make_decoder_with_options(params, true)
}

/// Build a G.728 decoder with explicit control over the §5.5 postfilter.
/// Set `postfilter_enabled` to `false` for raw synthesis output (useful
/// for bit-exact testing and for non-speech signals where the pitch
/// postfilter would hurt SNR more than it helps perception).
pub fn make_decoder_with_options(
    params: &CodecParameters,
    postfilter_enabled: bool,
) -> Result<Box<dyn Decoder>> {
    let sample_rate = params.sample_rate.unwrap_or(SAMPLE_RATE);
    if sample_rate != SAMPLE_RATE {
        return Err(Error::unsupported(format!(
            "G.728 decoder: only 8000 Hz is supported (got {sample_rate})"
        )));
    }
    let channels = params.channels.unwrap_or(1);
    if channels != 1 {
        return Err(Error::unsupported(format!(
            "G.728 decoder: only mono is supported (got {channels} channels)"
        )));
    }
    if params.codec_id.as_str() != CODEC_ID_STR {
        return Err(Error::unsupported(format!(
            "G.728 decoder: unexpected codec id {:?}",
            params.codec_id
        )));
    }
    let mut dec = G728Decoder::new();
    dec.set_postfilter_enabled(postfilter_enabled);
    Ok(Box::new(dec))
}

struct G728Decoder {
    codec_id: CodecId,
    state: G728State,
    /// Optional §5.5 adaptive postfilter — enabled by default. Disable
    /// via [`G728Decoder::set_postfilter_enabled`] when an application
    /// wants the raw synthesis output (e.g. non-speech test tones where
    /// perceptual shaping gets in the way).
    postfilter: Postfilter,
    postfilter_enabled: bool,
    /// Vector index within the current 4-vector frame (0..=3). Used to
    /// drive the once-a-frame pitch / long-term postfilter refresh.
    vec_in_frame: u32,
    pending: Option<Packet>,
    eof: bool,
    time_base: TimeBase,
}

impl G728Decoder {
    fn new() -> Self {
        Self {
            codec_id: CodecId::new(CODEC_ID_STR),
            state: G728State::new(),
            postfilter: Postfilter::new(),
            postfilter_enabled: true,
            vec_in_frame: 0,
            pending: None,
            eof: false,
            time_base: TimeBase::new(1, SAMPLE_RATE as i64),
        }
    }

    /// Enable or disable the §5.5 adaptive postfilter. Defaults to
    /// enabled. Disabling yields the raw synthesis output (without
    /// pitch emphasis, formant sharpening, or AGC).
    pub fn set_postfilter_enabled(&mut self, enabled: bool) {
        self.postfilter_enabled = enabled;
    }

    /// Decode a packet's worth of 10-bit indices into an f32 PCM buffer.
    fn decode_packet(&mut self, data: &[u8]) -> Result<Vec<f32>> {
        // Each index is 10 bits = 1.25 bytes. The packet length must be
        // enough to hold at least one index. We accept any byte length
        // that corresponds to a whole number of 5-sample vectors (i.e.
        // the bit count is a multiple of 10), or 10-bit-rounded-up-to-
        // bytes framing.
        let total_bits = (data.len() as u64) * 8;
        let vectors = total_bits / (INDEX_BITS as u64);
        if vectors == 0 {
            return Err(Error::invalid(format!(
                "G.728: packet too short ({} bytes; need at least 2 bytes for 1 index)",
                data.len()
            )));
        }
        let mut br = BitReader::new(data);
        let mut pcm = Vec::with_capacity((vectors as usize) * VECTOR_SIZE);
        let mut vec_out = [0.0_f32; VECTOR_SIZE];
        let mut post_out = [0.0_f32; VECTOR_SIZE];
        for _ in 0..vectors {
            let raw = br.read_index10()?;
            // Track whether the main LPC predictor *just* refreshed, so
            // we can hand its by-product 10th-order coefficients to the
            // postfilter. Check before decode_vector (which may trigger
            // the refresh).
            let will_refresh = self.state.lpc.vectors_since_update + 1 >= VECTORS_PER_BLOCK;
            self.state.decode_vector(raw, &mut vec_out);
            if will_refresh && self.state.lpc.last_update_ok {
                self.postfilter
                    .set_lpc(&self.state.lpc.apf, self.state.lpc.k1);
            }
            if self.postfilter_enabled {
                self.postfilter
                    .process_vector(&vec_out, self.vec_in_frame, &mut post_out);
                pcm.extend_from_slice(&post_out);
            } else {
                pcm.extend_from_slice(&vec_out);
            }
            self.vec_in_frame = (self.vec_in_frame + 1) % VECTORS_PER_BLOCK;
        }
        Ok(pcm)
    }
}

impl Decoder for G728Decoder {
    fn codec_id(&self) -> &CodecId {
        &self.codec_id
    }

    fn send_packet(&mut self, packet: &Packet) -> Result<()> {
        if self.pending.is_some() {
            return Err(Error::other(
                "G.728 decoder: receive_frame must be called before sending another packet",
            ));
        }
        self.pending = Some(packet.clone());
        Ok(())
    }

    fn receive_frame(&mut self) -> Result<Frame> {
        let Some(pkt) = self.pending.take() else {
            return if self.eof {
                Err(Error::Eof)
            } else {
                Err(Error::NeedMore)
            };
        };
        let samples = self.decode_packet(&pkt.data)?;
        // Convert f32 -> S16 LE.
        let mut bytes = Vec::with_capacity(samples.len() * 2);
        for &s in &samples {
            let v = s.round().clamp(-32768.0, 32767.0) as i16;
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        Ok(Frame::Audio(AudioFrame {
            format: SampleFormat::S16,
            channels: 1,
            sample_rate: SAMPLE_RATE,
            samples: samples.len() as u32,
            pts: pkt.pts,
            time_base: self.time_base,
            data: vec![bytes],
        }))
    }

    fn flush(&mut self) -> Result<()> {
        self.eof = true;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lpc_synthesise_zero_excitation_stays_silent() {
        let mut lpc = LpcPredictor::new();
        let exc = [0.0_f32; VECTOR_SIZE];
        let mut out = [0.0_f32; VECTOR_SIZE];
        lpc.synthesise(&exc, &mut out);
        assert!(out.iter().all(|&s| s == 0.0));
    }

    #[test]
    fn lpc_with_unit_a0_passes_impulse_through() {
        // With a[1..] == 0 the filter is just y[n] = x[n], so an impulse
        // excitation comes out unchanged.
        let mut lpc = LpcPredictor::new();
        let exc = [1.0, 0.0, 0.0, 0.0, 0.0];
        let mut out = [0.0; 5];
        lpc.synthesise(&exc, &mut out);
        assert_eq!(out, exc);
        // History is newest-first; after emitting [1, 0, 0, 0, 0] the impulse
        // is 4 slots deep.
        assert_eq!(lpc.history[0], 0.0);
        assert_eq!(lpc.history[4], 1.0);
    }

    #[test]
    fn gain_predictor_defaults_to_zero() {
        let mut gp = GainPredictor::new();
        assert_eq!(gp.predict(), 0.0);
    }

    #[test]
    fn excitation_from_index_is_nonzero_for_itu_tables() {
        // Any non-zero gain index should produce output.
        // Layout: shape(7) | sign(1) | mag(2). The 2-bit mag field selects
        // GAIN_CB[0..=3] directly; the sign bit flips polarity.
        let st = G728State::new();
        let idx: u16 = (5 << 3) | 2; // shape=5, sign=+, mag=2 → GAIN_CB[2]
        let v = st.excitation_from_index(idx);
        // Gain predictor starts at 0 ⇒ exp(0) = 1; GAIN_CB[2] ≠ 0.
        assert!(v.iter().any(|&s| s.abs() > 0.0));
    }

    #[test]
    fn state_defaults_are_neutral() {
        let st = G728State::new();
        assert_eq!(st.vector_count, 0);
        assert_eq!(st.lpc.a[0], 1.0);
        assert_eq!(st.gain.b[0], 1.0);
    }

    #[test]
    fn make_decoder_returns_working_decoder() {
        let mut params = CodecParameters::audio(CodecId::new(CODEC_ID_STR));
        params.sample_rate = Some(SAMPLE_RATE);
        params.channels = Some(1);
        assert!(make_decoder(&params).is_ok());
    }

    #[test]
    fn make_decoder_rejects_wrong_sample_rate() {
        let mut params = CodecParameters::audio(CodecId::new(CODEC_ID_STR));
        params.sample_rate = Some(16_000);
        assert!(make_decoder(&params).is_err());
    }

    #[test]
    fn decode_vector_produces_bounded_output() {
        let mut st = G728State::new();
        let mut out = [0.0_f32; VECTOR_SIZE];
        // Drive a fixed non-trivial index through 64 vectors.
        let raw: u16 = 0b00_1101_0010; // shape=26, sign=+, mag=2
        for _ in 0..64 {
            st.decode_vector(raw, &mut out);
            for &s in &out {
                assert!(s.is_finite(), "synthesis went non-finite: {s}");
                // LpcPredictor clamps to ±1e4 internally; the sample must
                // sit within that range plus a small i16-scale margin.
                assert!(s.abs() <= 1.0e4 + 1.0, "synthesis exploded: {s}");
            }
        }
    }

    #[test]
    fn zero_excitation_keeps_filter_stable() {
        // Feeding all-zero indices for a while must not cause growth.
        let mut st = G728State::new();
        let mut out = [0.0_f32; VECTOR_SIZE];
        // Shape row 0 is not the zero vector (ITU Annex B has 668,-2950,
        // ... / 2048), so the lowest `GAIN_CB` level still produces
        // non-zero excitation. The assertion is that the filter stays
        // bounded as the gain predictor ratchets the excitation scale
        // down over many vectors.
        let raw: u16 = 0; // shape=0, sign=0, mag=0
        let mut max_abs = 0.0_f32;
        for _ in 0..200 {
            st.decode_vector(raw, &mut out);
            for &s in &out {
                assert!(s.is_finite());
                if s.abs() > max_abs {
                    max_abs = s.abs();
                }
            }
        }
        // Output should stay bounded — much less than i16 full-scale.
        assert!(
            max_abs < 1.0e4,
            "constant-excitation output grew to {max_abs}"
        );
    }
}
