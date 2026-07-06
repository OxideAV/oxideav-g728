//! `oxideav-core` registry wiring — the packet/frame face of the codec.
//!
//! Registers one `g728` [`CodecInfo`] carrying both factories:
//!
//! * **decoder** — §5.11 serial byte stream in (10-bit channel indices
//!   packed MSB-first; the natural framing unit is one adaptation cycle
//!   = 4 indices = 5 bytes), 16-bit mono 8 kHz PCM out. Runs the
//!   **Annex G fixed-point** decoder ([`DecoderFixed`]) with the
//!   adaptive postfilter — the chain that reproduces the official ITU-T
//!   conformance vectors bit-exactly — and renders the `Q2` postfiltered
//!   speech to 16-bit linear (`×2`, the same conversion the conformance
//!   `outb4g` reference uses).
//! * **encoder** — 16-bit mono 8 kHz PCM in (`S = sample >> 1` §G.3.1
//!   input conversion), §5.11 byte stream out via the bit-exact
//!   [`EncoderFixed`]. Input is buffered to whole adaptation cycles
//!   (20 samples → 5 bytes) so every emitted packet is byte-aligned;
//!   `flush` zero-pads a trailing partial cycle.
//!
//! The claimed container tag is the registered WAVE format id `0x0041`
//! (ITU-T G.728 CELP).
//!
//! The direct (registry-free) factories remain [`crate::make_decoder`] /
//! [`crate::make_encoder`] per the workspace dual-API convention.

use std::collections::VecDeque;

use oxideav_core::{
    AudioFrame, CodecCapabilities, CodecId, CodecInfo, CodecParameters, CodecTag, Error, Frame,
    Packet, Result, RuntimeContext, SampleFormat, TimeBase,
};

use crate::annex_g_coder::{DecoderFixed, EncoderFixed};
use crate::consts::IDIM;

/// G.728 fixed line rate: 8 kHz sampling.
const SAMPLE_RATE: u32 = 8000;
/// One adaptation cycle: 4 vectors × 5 samples.
const CYCLE_SAMPLES: usize = 20;
/// One adaptation cycle on the wire: 4 × 10 bits = 5 bytes.
const CYCLE_BYTES: usize = 5;

/// Registered WAVE format id for ITU-T G.728 CELP.
const WAVE_FORMAT_G728_CELP: u16 = 0x0041;

/// Register the `g728` codec (decoder + encoder) into the runtime
/// context's codec registry.
pub fn register(ctx: &mut RuntimeContext) {
    let mut caps = CodecCapabilities::audio("g728_sw");
    caps.lossy = true;
    ctx.codecs.register(
        CodecInfo::new(CodecId::new("g728"))
            .capabilities(caps)
            .decoder(new_decoder)
            .encoder(new_encoder)
            .tag(CodecTag::wave_format(WAVE_FORMAT_G728_CELP)),
    );
}

/// Validate the audio parameters shared by both factories: G.728 is
/// mono 8 kHz, 16-bit linear PCM on the uncompressed side.
fn check_params(params: &CodecParameters) -> Result<()> {
    if let Some(rate) = params.sample_rate {
        if rate != SAMPLE_RATE {
            return Err(Error::unsupported(format!(
                "g728: sample rate {rate} unsupported (G.728 is {SAMPLE_RATE} Hz only)"
            )));
        }
    }
    if let Some(ch) = params.channels {
        if ch != 1 {
            return Err(Error::unsupported(format!(
                "g728: {ch} channels unsupported (G.728 is mono)"
            )));
        }
    }
    if let Some(fmt) = params.sample_format {
        if fmt != SampleFormat::S16 {
            return Err(Error::unsupported(format!(
                "g728: sample format {fmt:?} unsupported (S16 only)"
            )));
        }
    }
    Ok(())
}

// ───────────────────────── decoder ─────────────────────────

struct G728RegistryDecoder {
    id: CodecId,
    inner: DecoderFixed,
    /// Carry-over compressed bytes (< [`CYCLE_BYTES`]) between packets.
    pending: Vec<u8>,
    /// Decoded frames awaiting `receive_frame`.
    out: VecDeque<Frame>,
    /// End-of-stream flag set by `flush`.
    eof: bool,
    /// Running sample counter for synthesized PTS (1/8000 time base).
    next_pts: i64,
}

fn new_decoder(params: &CodecParameters) -> Result<Box<dyn oxideav_core::Decoder>> {
    check_params(params)?;
    Ok(Box::new(G728RegistryDecoder {
        id: CodecId::new("g728"),
        inner: DecoderFixed::new(),
        pending: Vec::new(),
        out: VecDeque::new(),
        eof: false,
        next_pts: 0,
    }))
}

impl G728RegistryDecoder {
    /// Decode every whole adaptation cycle in `self.pending` into one
    /// [`AudioFrame`] (S16 interleaved mono).
    fn drain_cycles(&mut self, pts: Option<i64>) {
        let cycles = self.pending.len() / CYCLE_BYTES;
        if cycles == 0 {
            return;
        }
        let consumed = cycles * CYCLE_BYTES;
        let bytes: Vec<u8> = self.pending.drain(..consumed).collect();
        let vectors = self
            .inner
            .decode_bytes(&bytes)
            .expect("whole cycles are whole 10-bit indices");
        let mut data = Vec::with_capacity(vectors.len() * IDIM * 2);
        for v in &vectors {
            for &q2 in v.iter() {
                // Q2 → 16-bit linear: ×2, saturated (the conformance
                // reference conversion).
                let s = (i64::from(q2) * 2).clamp(i64::from(i16::MIN), i64::from(i16::MAX)) as i16;
                data.extend_from_slice(&s.to_le_bytes());
            }
        }
        let samples = (vectors.len() * IDIM) as u32;
        let pts = pts.or(Some(self.next_pts));
        self.next_pts += i64::from(samples);
        self.out.push_back(Frame::Audio(AudioFrame {
            samples,
            pts,
            data: vec![data],
        }));
    }
}

impl oxideav_core::Decoder for G728RegistryDecoder {
    fn codec_id(&self) -> &CodecId {
        &self.id
    }

    fn send_packet(&mut self, packet: &Packet) -> Result<()> {
        if self.eof {
            return Err(Error::invalid("g728: send_packet after flush"));
        }
        self.pending.extend_from_slice(&packet.data);
        self.drain_cycles(packet.pts);
        Ok(())
    }

    fn receive_frame(&mut self) -> Result<Frame> {
        match self.out.pop_front() {
            Some(f) => Ok(f),
            None if self.eof => Err(Error::Eof),
            None => Err(Error::NeedMore),
        }
    }

    fn flush(&mut self) -> Result<()> {
        // A conforming §5.11 stream is 5-byte-framed; trailing bytes
        // cannot hold a whole adaptation cycle and are dropped.
        self.pending.clear();
        self.eof = true;
        Ok(())
    }

    fn reset(&mut self) -> Result<()> {
        self.inner = DecoderFixed::new();
        self.pending.clear();
        self.out.clear();
        self.eof = false;
        self.next_pts = 0;
        Ok(())
    }
}

// ───────────────────────── encoder ─────────────────────────

struct G728RegistryEncoder {
    id: CodecId,
    params: CodecParameters,
    inner: EncoderFixed,
    /// Carry-over PCM samples (< [`CYCLE_SAMPLES`]) between frames.
    pending: Vec<i16>,
    out: VecDeque<Packet>,
    /// Running compressed-sample counter for packet PTS (1/8000 base).
    next_pts: i64,
    flushed: bool,
}

fn new_encoder(params: &CodecParameters) -> Result<Box<dyn oxideav_core::Encoder>> {
    check_params(params)?;
    let mut out_params = CodecParameters::audio(CodecId::new("g728"));
    out_params.sample_rate = Some(SAMPLE_RATE);
    out_params.channels = Some(1);
    out_params.sample_format = Some(SampleFormat::S16);
    Ok(Box::new(G728RegistryEncoder {
        id: CodecId::new("g728"),
        params: out_params,
        inner: EncoderFixed::new(),
        pending: Vec::new(),
        out: VecDeque::new(),
        next_pts: 0,
        flushed: false,
    }))
}

impl G728RegistryEncoder {
    /// Encode every whole adaptation cycle in `self.pending` into one
    /// packet (5 bytes per cycle).
    fn drain_cycles(&mut self) {
        let cycles = self.pending.len() / CYCLE_SAMPLES;
        if cycles == 0 {
            return;
        }
        let consumed = cycles * CYCLE_SAMPLES;
        let vectors: Vec<[i16; IDIM]> = self.pending[..consumed]
            .chunks_exact(IDIM)
            // §G.3.1: 16-bit linear → Q2 is a right shift of 1 bit.
            .map(|c| std::array::from_fn(|k| c[k] >> 1))
            .collect();
        self.pending.drain(..consumed);
        let bytes = self.inner.encode_buffer(&vectors);
        let samples = (cycles * CYCLE_SAMPLES) as i64;
        let packet = Packet::new(0, TimeBase::new(1, i64::from(SAMPLE_RATE)), bytes)
            .with_pts(self.next_pts)
            .with_duration(samples)
            .with_keyframe(self.next_pts == 0);
        self.next_pts += samples;
        self.out.push_back(packet);
    }
}

impl oxideav_core::Encoder for G728RegistryEncoder {
    fn codec_id(&self) -> &CodecId {
        &self.id
    }

    fn output_params(&self) -> &CodecParameters {
        &self.params
    }

    fn send_frame(&mut self, frame: &Frame) -> Result<()> {
        if self.flushed {
            return Err(Error::invalid("g728: send_frame after flush"));
        }
        let audio = match frame {
            Frame::Audio(a) => a,
            _ => return Err(Error::invalid("g728: encoder accepts audio frames only")),
        };
        let plane = audio
            .data
            .first()
            .ok_or_else(|| Error::invalid("g728: audio frame has no data plane"))?;
        if plane.len() % 2 != 0 {
            return Err(Error::invalid("g728: S16 plane has odd byte length"));
        }
        self.pending.extend(
            plane
                .chunks_exact(2)
                .map(|b| i16::from_le_bytes([b[0], b[1]])),
        );
        self.drain_cycles();
        Ok(())
    }

    fn receive_packet(&mut self) -> Result<Packet> {
        match self.out.pop_front() {
            Some(p) => Ok(p),
            None if self.flushed => Err(Error::Eof),
            None => Err(Error::NeedMore),
        }
    }

    fn flush(&mut self) -> Result<()> {
        if !self.pending.is_empty() {
            // Zero-pad the trailing partial cycle so the last samples
            // are not lost (the decoder will render the pad as
            // near-silence).
            self.pending.resize(
                self.pending.len().div_ceil(CYCLE_SAMPLES) * CYCLE_SAMPLES,
                0,
            );
            self.drain_cycles();
        }
        self.flushed = true;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn audio_params() -> CodecParameters {
        let mut p = CodecParameters::audio(CodecId::new("g728"));
        p.sample_rate = Some(SAMPLE_RATE);
        p.channels = Some(1);
        p.sample_format = Some(SampleFormat::S16);
        p
    }

    fn s16_frame(samples: &[i16]) -> Frame {
        let mut data = Vec::with_capacity(samples.len() * 2);
        for s in samples {
            data.extend_from_slice(&s.to_le_bytes());
        }
        Frame::Audio(AudioFrame {
            samples: samples.len() as u32,
            pts: None,
            data: vec![data],
        })
    }

    fn tone(n: usize) -> Vec<i16> {
        (0..n)
            .map(|i| ((i as f64 * 0.11).sin() * 11000.0 + (i as f64 * 0.041).sin() * 4000.0) as i16)
            .collect()
    }

    #[test]
    fn register_exposes_decoder_and_encoder() {
        let mut ctx = RuntimeContext::new();
        register(&mut ctx);
        let id = CodecId::new("g728");
        assert!(ctx.codecs.has_decoder(&id));
        assert!(ctx.codecs.has_encoder(&id));
        // Tag lookup: WAVE format 0x0041 resolves to g728.
        let tag = CodecTag::wave_format(WAVE_FORMAT_G728_CELP);
        let probe = oxideav_core::ProbeContext::new(&tag);
        assert_eq!(
            ctx.codecs.resolve_tag_ref(&probe).map(|c| c.as_str()),
            Some("g728")
        );
    }

    #[test]
    fn factories_reject_wrong_parameters() {
        let mut p = audio_params();
        p.sample_rate = Some(16000);
        assert!(new_decoder(&p).is_err());
        assert!(new_encoder(&p).is_err());
        let mut p = audio_params();
        p.channels = Some(2);
        assert!(new_decoder(&p).is_err());
        let mut p = audio_params();
        p.sample_format = Some(SampleFormat::U8);
        assert!(new_encoder(&p).is_err());
    }

    #[test]
    fn encode_decode_roundtrip_through_registry() {
        // 16 kbit/s LD-CELP: 200 samples → 10 packets-worth of cycles →
        // decoded speech tracks the input closely after convergence.
        let mut ctx = RuntimeContext::new();
        register(&mut ctx);
        let params = audio_params();
        let mut enc = ctx.codecs.first_encoder(&params).expect("encoder");
        let mut dec = ctx.codecs.first_decoder(&params).expect("decoder");

        let pcm = tone(400);
        enc.send_frame(&s16_frame(&pcm)).expect("send_frame");
        enc.flush().expect("flush");

        let mut decoded = Vec::new();
        loop {
            let packet = match enc.receive_packet() {
                Ok(p) => p,
                Err(Error::Eof) => break,
                Err(e) => panic!("receive_packet: {e}"),
            };
            // 16 kbit/s: 5 bytes per 20 samples.
            assert_eq!(packet.data.len() % CYCLE_BYTES, 0);
            dec.send_packet(&packet).expect("send_packet");
            loop {
                match dec.receive_frame() {
                    Ok(Frame::Audio(a)) => {
                        assert_eq!(a.data.len(), 1);
                        decoded.extend(
                            a.data[0]
                                .chunks_exact(2)
                                .map(|b| i16::from_le_bytes([b[0], b[1]])),
                        );
                    }
                    Ok(_) => panic!("non-audio frame"),
                    Err(Error::NeedMore) => break,
                    Err(e) => panic!("receive_frame: {e}"),
                }
            }
        }
        assert_eq!(decoded.len(), pcm.len());

        // Correlation over the converged tail.
        let (mut num, mut den_a, mut den_b) = (0.0f64, 0.0f64, 0.0f64);
        for i in 100..pcm.len() {
            let a = f64::from(pcm[i]);
            let b = f64::from(decoded[i]);
            num += a * b;
            den_a += a * a;
            den_b += b * b;
        }
        let corr = num / (den_a.sqrt() * den_b.sqrt());
        assert!(corr > 0.9, "decoded speech should track input ({corr})");
    }

    #[test]
    fn decoder_buffers_partial_cycles_across_packets() {
        let mut ctx = RuntimeContext::new();
        register(&mut ctx);
        let params = audio_params();
        let mut enc = ctx.codecs.first_encoder(&params).expect("encoder");
        let mut dec = ctx.codecs.first_decoder(&params).expect("decoder");

        enc.send_frame(&s16_frame(&tone(40))).expect("send");
        enc.flush().expect("flush");
        let packet = enc.receive_packet().expect("packet");
        assert_eq!(packet.data.len(), 10, "40 samples = 2 cycles = 10 bytes");

        // Split mid-cycle: 3 bytes now, 7 bytes later.
        let head = Packet::new(0, packet.time_base, packet.data[..3].to_vec());
        let tail = Packet::new(0, packet.time_base, packet.data[3..].to_vec());
        dec.send_packet(&head).expect("head");
        assert!(matches!(dec.receive_frame(), Err(Error::NeedMore)));
        dec.send_packet(&tail).expect("tail");
        let mut total = 0u32;
        while let Ok(Frame::Audio(a)) = dec.receive_frame() {
            total += a.samples;
        }
        assert_eq!(total, 40);
    }

    #[test]
    fn flush_pads_partial_cycle_and_signals_eof() {
        let mut ctx = RuntimeContext::new();
        register(&mut ctx);
        let params = audio_params();
        let mut enc = ctx.codecs.first_encoder(&params).expect("encoder");
        // 27 samples: one whole cycle + 7 residual → flush pads to 40.
        enc.send_frame(&s16_frame(&tone(27))).expect("send");
        enc.flush().expect("flush");
        let mut bytes = 0usize;
        loop {
            match enc.receive_packet() {
                Ok(p) => bytes += p.data.len(),
                Err(Error::Eof) => break,
                Err(e) => panic!("{e}"),
            }
        }
        assert_eq!(bytes, 10, "27 samples pad to 2 cycles = 10 bytes");
        // Decoder-side flush → Eof.
        let mut dec = ctx.codecs.first_decoder(&params).expect("decoder");
        dec.flush().expect("flush");
        assert!(matches!(dec.receive_frame(), Err(Error::Eof)));
    }

    #[test]
    fn decoder_reset_restores_cold_state() {
        let mut ctx = RuntimeContext::new();
        register(&mut ctx);
        let params = audio_params();
        let mut enc = ctx.codecs.first_encoder(&params).expect("encoder");
        enc.send_frame(&s16_frame(&tone(40))).expect("send");
        enc.flush().expect("flush");
        let packet = enc.receive_packet().expect("packet");

        let mut dec1 = ctx.codecs.first_decoder(&params).expect("decoder");
        dec1.send_packet(&packet).expect("send");
        let f1 = dec1.receive_frame().expect("frame");
        // Reset, decode the same packet again: identical output (cold
        // backward-adapter state both times).
        dec1.reset().expect("reset");
        dec1.send_packet(&packet).expect("send again");
        let f2 = dec1.receive_frame().expect("frame again");
        match (f1, f2) {
            (Frame::Audio(a), Frame::Audio(b)) => {
                assert_eq!(a.data, b.data, "reset must restore cold-start state")
            }
            _ => panic!("audio frames expected"),
        }
    }
}
