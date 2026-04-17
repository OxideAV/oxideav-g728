//! Round-trip SNR characterisation with the ITU-T G.728 Annex B tables.
//!
//! Encodes a 400 Hz sine and a 3-tone synthetic voiced-speech mix, then
//! decodes and measures the per-sample reconstruction error in dB after
//! the first 100 ms of adaptation transient. With the Annex B `CODEBK` /
//! `GQ` tables shipped in this crate the sine lands around 44 dB SNR and
//! the voiced mix around 30 dB; swapping the placeholder codebooks that
//! preceded Annex B left these at ~38 dB and ~25 dB respectively.
//!
//! The numeric floors are intentionally loose — they protect against
//! regressions in the analysis-by-synthesis loop, not against deviations
//! from the ITU reference decoder (this crate still ships a Hamming
//! autocorrelation window rather than the spec's recursive Barnwell
//! window, and has no postfilter).

use oxideav_codec::{Decoder, Encoder};
use oxideav_core::{
    AudioFrame, CodecId, CodecParameters, Error, Frame, Packet, SampleFormat, TimeBase,
};
use oxideav_g728::encoder::{PACKET_BYTES, PACKET_SAMPLES};
use oxideav_g728::{CODEC_ID_STR, SAMPLE_RATE};

const PACKETS: usize = 400;
const TOTAL_SAMPLES: usize = PACKETS * PACKET_SAMPLES;

fn make_params() -> CodecParameters {
    let mut p = CodecParameters::audio(CodecId::new(CODEC_ID_STR));
    p.sample_rate = Some(SAMPLE_RATE);
    p.channels = Some(1);
    p.sample_format = Some(SampleFormat::S16);
    p
}

fn build_sine(n: usize) -> Vec<i16> {
    let sr = SAMPLE_RATE as f32;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f32 / sr;
        let s = 4000.0 * (2.0 * core::f32::consts::PI * 400.0 * t).sin();
        out.push(s.round().clamp(-32768.0, 32767.0) as i16);
    }
    out
}

fn pack_audio_frame(samples: &[i16]) -> AudioFrame {
    let mut bytes = Vec::with_capacity(samples.len() * 2);
    for &s in samples {
        bytes.extend_from_slice(&s.to_le_bytes());
    }
    AudioFrame {
        format: SampleFormat::S16,
        channels: 1,
        sample_rate: SAMPLE_RATE,
        samples: samples.len() as u32,
        pts: None,
        time_base: TimeBase::new(1, SAMPLE_RATE as i64),
        data: vec![bytes],
    }
}

fn encode_all(enc: &mut Box<dyn Encoder>, input: &[i16]) -> Vec<Packet> {
    let mut packets = Vec::new();
    for block in input.chunks(PACKET_SAMPLES) {
        let af = pack_audio_frame(block);
        enc.send_frame(&Frame::Audio(af)).unwrap();
        loop {
            match enc.receive_packet() {
                Ok(p) => packets.push(p),
                Err(Error::NeedMore) => break,
                Err(e) => panic!("recv: {e}"),
            }
        }
    }
    enc.flush().unwrap();
    loop {
        match enc.receive_packet() {
            Ok(p) => packets.push(p),
            Err(Error::NeedMore) | Err(Error::Eof) => break,
            Err(e) => panic!("recv post-flush: {e}"),
        }
    }
    packets
}

fn decode_all(dec: &mut Box<dyn Decoder>, packets: &[Packet]) -> Vec<i16> {
    let mut pcm = Vec::new();
    for p in packets {
        dec.send_packet(p).unwrap();
        loop {
            match dec.receive_frame() {
                Ok(Frame::Audio(af)) => {
                    for chunk in af.data[0].chunks_exact(2) {
                        pcm.push(i16::from_le_bytes([chunk[0], chunk[1]]));
                    }
                }
                Ok(_) => {}
                Err(Error::NeedMore) | Err(Error::Eof) => break,
                Err(e) => panic!("rf: {e}"),
            }
        }
    }
    let _ = dec.flush();
    pcm
}

fn snr_db(input: &[i16], output: &[i16]) -> f64 {
    let n = input.len().min(output.len());
    // Skip leading transient (100 ms) so predictors can settle.
    let skip = (SAMPLE_RATE as usize) / 10;
    if n <= skip {
        return f64::NAN;
    }
    let mut sig_e = 0.0f64;
    let mut err_e = 0.0f64;
    for i in skip..n {
        let s = input[i] as f64;
        let o = output[i] as f64;
        sig_e += s * s;
        let d = s - o;
        err_e += d * d;
    }
    if err_e == 0.0 {
        f64::INFINITY
    } else {
        10.0 * (sig_e / err_e).log10()
    }
}

/// Mix three speech-band tones (male pitch + vowel formants) to stand in
/// for narrowband voiced speech without needing an external WAV file.
fn build_voiced(n: usize) -> Vec<i16> {
    let sr = SAMPLE_RATE as f32;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f32 / sr;
        let s = 2000.0 * (2.0 * core::f32::consts::PI * 150.0 * t).sin()
            + 1500.0 * (2.0 * core::f32::consts::PI * 700.0 * t).sin()
            + 800.0 * (2.0 * core::f32::consts::PI * 1600.0 * t).sin();
        out.push(s.round().clamp(-32768.0, 32767.0) as i16);
    }
    out
}

#[test]
fn print_roundtrip_snr() {
    let params = make_params();

    // 400 Hz pure tone.
    {
        let input = build_sine(TOTAL_SAMPLES);
        let mut enc = oxideav_g728::encoder::make_encoder(&params).unwrap();
        let packets = encode_all(&mut enc, &input);
        assert_eq!(packets.len(), PACKETS);
        for p in &packets {
            assert_eq!(p.data.len(), PACKET_BYTES);
        }
        let mut dec = oxideav_g728::decoder::make_decoder(&params).unwrap();
        let out = decode_all(&mut dec, &packets);
        let snr = snr_db(&input, &out);
        println!("roundtrip SNR (400 Hz sine): {snr:.2} dB");
        assert!(snr.is_finite(), "SNR is not finite: {snr}");
        assert!(
            snr >= 40.0,
            "sine SNR regressed below 40 dB floor: {snr:.2}"
        );
    }

    // Synthetic voiced-speech-like mix.
    {
        let input = build_voiced(TOTAL_SAMPLES);
        let mut enc = oxideav_g728::encoder::make_encoder(&params).unwrap();
        let packets = encode_all(&mut enc, &input);
        let mut dec = oxideav_g728::decoder::make_decoder(&params).unwrap();
        let out = decode_all(&mut dec, &packets);
        let snr = snr_db(&input, &out);
        println!("roundtrip SNR (voiced-speech mix): {snr:.2} dB");
        assert!(snr.is_finite(), "voiced SNR is not finite: {snr}");
        assert!(
            snr >= 27.0,
            "voiced SNR regressed below 27 dB floor: {snr:.2}"
        );
    }
}
