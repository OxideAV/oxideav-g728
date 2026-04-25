//! Integration tests for G.728 frame-erasure concealment (Annex A.3 /
//! §5.8).
//!
//! A `Packet` whose `flags.corrupt` bit is set is treated as a lost
//! frame: the decoder doesn't try to unpack its payload and instead
//! emits an extrapolation based on the last good excitation, attenuated
//! across successive erasures. These tests exercise that path via the
//! standard `Decoder` trait so no private APIs are touched.

use oxideav_core::{
    AudioFrame, CodecId, CodecParameters, Error, Frame, Packet, SampleFormat, TimeBase,
};
use oxideav_core::{Decoder, Encoder};
use oxideav_g728::encoder::{PACKET_BYTES, PACKET_SAMPLES};
use oxideav_g728::{CODEC_ID_STR, SAMPLE_RATE, VECTOR_SIZE};

const PACKETS: usize = 80; // 200 ms of audio
const TOTAL_SAMPLES: usize = PACKETS * PACKET_SAMPLES;

fn make_params() -> CodecParameters {
    let mut p = CodecParameters::audio(CodecId::new(CODEC_ID_STR));
    p.sample_rate = Some(SAMPLE_RATE);
    p.channels = Some(1);
    p.sample_format = Some(SampleFormat::S16);
    p
}

fn make_encoder() -> Box<dyn Encoder> {
    oxideav_g728::encoder::make_encoder(&make_params()).expect("encoder ctor")
}

fn make_decoder() -> Box<dyn Decoder> {
    oxideav_g728::decoder::make_decoder(&make_params()).expect("decoder ctor")
}

fn build_sine(n: usize) -> Vec<i16> {
    let sr = SAMPLE_RATE as f32;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f32 / sr;
        let s = 2000.0 * (2.0 * core::f32::consts::PI * 400.0 * t).sin();
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

#[test]
fn single_corrupt_packet_emits_plausible_audio() {
    // Encode a sine, then replay with one packet flagged corrupt and
    // confirm the decoder emits bounded audio rather than an error.
    let input = build_sine(TOTAL_SAMPLES);
    let mut enc = make_encoder();
    let packets = encode_all(&mut enc, &input);
    assert!(packets.len() >= 30);

    let mut dec = make_decoder();

    // Warm up on 20 clean packets.
    for p in packets.iter().take(20) {
        dec.send_packet(p).unwrap();
        let _ = dec.receive_frame().unwrap();
    }

    // Erase one packet.
    let mut bad = packets[20].clone();
    bad.flags.corrupt = true;
    dec.send_packet(&bad).unwrap();
    let Frame::Audio(af) = dec.receive_frame().unwrap() else {
        panic!("expected audio frame");
    };
    assert_eq!(af.samples as usize, PACKET_SAMPLES);
    // First-erasure attenuation is 0.9, so output should still carry
    // most of the pre-erasure energy.
    let mut peak: i32 = 0;
    for chunk in af.data[0].chunks_exact(2) {
        let s = i16::from_le_bytes([chunk[0], chunk[1]]) as i32;
        peak = peak.max(s.abs());
    }
    assert!(peak > 50, "concealed output too quiet: peak = {peak}");
    assert!(peak < 20_000, "concealed output too loud: peak = {peak}");

    // Resume clean decoding without error.
    for p in packets.iter().skip(21).take(5) {
        dec.send_packet(p).unwrap();
        let _ = dec.receive_frame().unwrap();
    }
}

#[test]
fn long_erasure_fades_to_silence() {
    // Many erased packets in a row should ramp the output down.
    let input = build_sine(TOTAL_SAMPLES);
    let mut enc = make_encoder();
    let packets = encode_all(&mut enc, &input);

    let mut dec = make_decoder();
    for p in packets.iter().take(20) {
        dec.send_packet(p).unwrap();
        let _ = dec.receive_frame().unwrap();
    }

    // Build 20 corrupt packets and feed them back-to-back. Measure the
    // peak on each erasure; by the end the output must be very quiet.
    let template = packets[20].clone();
    let mut peaks: Vec<i32> = Vec::with_capacity(20);
    for _ in 0..20 {
        let mut bad = template.clone();
        bad.flags.corrupt = true;
        dec.send_packet(&bad).unwrap();
        let Frame::Audio(af) = dec.receive_frame().unwrap() else {
            panic!("expected audio");
        };
        let mut peak: i32 = 0;
        for chunk in af.data[0].chunks_exact(2) {
            let s = i16::from_le_bytes([chunk[0], chunk[1]]) as i32;
            peak = peak.max(s.abs());
        }
        peaks.push(peak);
    }
    // Final erasure: must be near-silent.
    assert!(
        *peaks.last().unwrap() < 500,
        "tail of long erasure not silent: {:?}",
        peaks
    );
    // Overall trend: the 5th erasure should already be substantially
    // quieter than the 1st (both absolute and per the 0.9/0.7/0.5/0.3/0
    // attenuation ladder).
    assert!(
        peaks[4] * 2 < peaks[0] + 100,
        "fifth-erasure peak not decaying: {:?}",
        peaks
    );
}

#[test]
fn corrupt_packet_without_duration_falls_back_to_4_vectors() {
    // A packet flagged corrupt with no duration hint should still
    // emit PCM: the decoder uses the standard 4-vector packet size.
    let mut dec = make_decoder();
    let mut pkt = Packet::new(0, TimeBase::new(1, SAMPLE_RATE as i64), Vec::new());
    pkt.flags.corrupt = true;
    dec.send_packet(&pkt).unwrap();
    let Frame::Audio(af) = dec.receive_frame().unwrap() else {
        panic!("expected audio");
    };
    assert_eq!(af.samples as usize, 4 * VECTOR_SIZE);
}

#[test]
fn corrupt_packet_with_custom_duration_sizes_output() {
    // 3 vectors = 15 samples.
    let mut dec = make_decoder();
    let mut pkt = Packet::new(0, TimeBase::new(1, SAMPLE_RATE as i64), Vec::new());
    pkt.flags.corrupt = true;
    pkt.duration = Some(15);
    dec.send_packet(&pkt).unwrap();
    let Frame::Audio(af) = dec.receive_frame().unwrap() else {
        panic!("expected audio");
    };
    assert_eq!(af.samples as usize, 15);
}

#[test]
fn concealment_recovers_cleanly_after_erasure_burst() {
    // After a 3-packet burst loss, the decoder should resume on the
    // next clean packet and not throw or produce silence forever.
    let input = build_sine(TOTAL_SAMPLES);
    let mut enc = make_encoder();
    let packets = encode_all(&mut enc, &input);

    let mut dec = make_decoder();
    for p in packets.iter().take(20) {
        dec.send_packet(p).unwrap();
        let _ = dec.receive_frame().unwrap();
    }

    // 3-packet burst erasure.
    for p in packets.iter().skip(20).take(3) {
        let mut bad = p.clone();
        bad.flags.corrupt = true;
        dec.send_packet(&bad).unwrap();
        let _ = dec.receive_frame().unwrap();
    }

    // Feed 10 clean packets and expect nontrivial output by the end.
    let mut post_loss_peak: i32 = 0;
    for p in packets.iter().skip(23).take(10) {
        dec.send_packet(p).unwrap();
        let Frame::Audio(af) = dec.receive_frame().unwrap() else {
            panic!("expected audio");
        };
        for chunk in af.data[0].chunks_exact(2) {
            let s = i16::from_le_bytes([chunk[0], chunk[1]]) as i32;
            post_loss_peak = post_loss_peak.max(s.abs());
        }
    }
    assert!(
        post_loss_peak > 200,
        "decoder stayed near-silent after burst loss recovery: peak = {post_loss_peak}"
    );
}

#[test]
fn first_packet_erasure_is_silent_not_a_crash() {
    // If the very first packet is marked corrupt, the decoder has no
    // history to extrapolate from. It must still not panic and must
    // emit PCM of the right length.
    let mut dec = make_decoder();
    let mut pkt = Packet::new(
        0,
        TimeBase::new(1, SAMPLE_RATE as i64),
        vec![0u8; PACKET_BYTES],
    );
    pkt.flags.corrupt = true;
    pkt.duration = Some(PACKET_SAMPLES as i64);
    dec.send_packet(&pkt).unwrap();
    let Frame::Audio(af) = dec.receive_frame().unwrap() else {
        panic!("expected audio");
    };
    assert_eq!(af.samples as usize, PACKET_SAMPLES);
    // Cold-start concealment: excitation history is zero, so output
    // must be silent.
    let peak: i32 = af.data[0]
        .chunks_exact(2)
        .map(|c| (i16::from_le_bytes([c[0], c[1]]) as i32).abs())
        .max()
        .unwrap_or(0);
    assert!(
        peak < 100,
        "cold-start erasure emitted audible content: {peak}"
    );
}

#[test]
fn clean_packet_clears_erasure_run() {
    // Two erasures, one clean, then two more erasures. The attenuation
    // on the 4th and 5th erased packets should be as if they were the
    // first and second after the clean packet — i.e. 0.9 then 0.7 —
    // not 0.3 then 0.0.
    let input = build_sine(TOTAL_SAMPLES);
    let mut enc = make_encoder();
    let packets = encode_all(&mut enc, &input);

    let mut dec = make_decoder();
    for p in packets.iter().take(20) {
        dec.send_packet(p).unwrap();
        let _ = dec.receive_frame().unwrap();
    }

    // Two corrupt packets.
    for p in packets.iter().skip(20).take(2) {
        let mut bad = p.clone();
        bad.flags.corrupt = true;
        dec.send_packet(&bad).unwrap();
        let _ = dec.receive_frame().unwrap();
    }

    // One clean packet to reset the run counter.
    dec.send_packet(&packets[22]).unwrap();
    let _ = dec.receive_frame().unwrap();

    // Next erasure should be treated as a fresh one — high peak.
    let mut bad = packets[23].clone();
    bad.flags.corrupt = true;
    dec.send_packet(&bad).unwrap();
    let Frame::Audio(af) = dec.receive_frame().unwrap() else {
        panic!("expected audio");
    };
    let peak: i32 = af.data[0]
        .chunks_exact(2)
        .map(|c| (i16::from_le_bytes([c[0], c[1]]) as i32).abs())
        .max()
        .unwrap_or(0);
    // Attenuation 0.9 over a fresh excitation should land well above
    // the silence floor produced by a fully-faded 5-packet burst.
    assert!(
        peak > 100,
        "erasure run didn't reset after clean packet: peak = {peak}"
    );
}
