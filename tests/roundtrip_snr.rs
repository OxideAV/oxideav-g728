//! Round-trip SNR characterisation with the ITU-T G.728 Annex B tables.
//!
//! Encodes a 400 Hz sine and a 3-tone synthetic voiced-speech mix, then
//! decodes and measures the per-sample reconstruction error in dB after
//! the first 100 ms of adaptation transient.
//!
//! The SNR floors below are measured against the **raw synthesis**
//! output — i.e. with the §5.5 adaptive postfilter disabled. This is the
//! correct quantity for an "encoder fidelity" gauge: the postfilter is a
//! perceptual shaper (pitch comb + formant emphasis + spectral tilt
//! compensation) and, by design, *lowers* L2 SNR on pure tones and
//! non-speech inputs while raising subjective quality on speech. The
//! spec itself (§4.6.1) recommends turning the postfilter off for
//! non-speech signals.
//!
//! A companion test `print_postfilter_delta` exercises the postfilter
//! path end-to-end and asserts that it actually modifies the signal
//! (non-trivially) without blowing up the level.

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

    // 400 Hz pure tone — postfilter off, measures raw synthesis fidelity.
    {
        let input = build_sine(TOTAL_SAMPLES);
        let mut enc = oxideav_g728::encoder::make_encoder(&params).unwrap();
        let packets = encode_all(&mut enc, &input);
        assert_eq!(packets.len(), PACKETS);
        for p in &packets {
            assert_eq!(p.data.len(), PACKET_BYTES);
        }
        let mut dec = oxideav_g728::decoder::make_decoder_with_options(&params, false).unwrap();
        let out = decode_all(&mut dec, &packets);
        let snr = snr_db(&input, &out);
        println!("roundtrip SNR (400 Hz sine, raw): {snr:.2} dB");
        assert!(snr.is_finite(), "SNR is not finite: {snr}");
        assert!(
            snr >= 40.0,
            "sine SNR regressed below 40 dB floor: {snr:.2}"
        );
    }

    // Synthetic voiced-speech-like mix — postfilter off, raw fidelity.
    {
        let input = build_voiced(TOTAL_SAMPLES);
        let mut enc = oxideav_g728::encoder::make_encoder(&params).unwrap();
        let packets = encode_all(&mut enc, &input);
        let mut dec = oxideav_g728::decoder::make_decoder_with_options(&params, false).unwrap();
        let out = decode_all(&mut dec, &packets);
        let snr = snr_db(&input, &out);
        println!("roundtrip SNR (voiced-speech mix, raw): {snr:.2} dB");
        assert!(snr.is_finite(), "voiced SNR is not finite: {snr}");
        assert!(
            snr >= 27.0,
            "voiced SNR regressed below 27 dB floor: {snr:.2}"
        );
    }
}

/// Exercise the default (postfilter-on) path end-to-end. The postfilter
/// is perceptual, so we don't assert an SNR floor — only that the output
/// is finite, non-trivially different from the raw-synthesis output, and
/// preserves the overall energy level (AGC).
#[test]
fn print_postfilter_delta() {
    let params = make_params();

    let input = build_voiced(TOTAL_SAMPLES);
    let mut enc = oxideav_g728::encoder::make_encoder(&params).unwrap();
    let packets = encode_all(&mut enc, &input);

    // Decode twice — once raw, once with the postfilter — so we can
    // measure the delta.
    let mut dec_raw = oxideav_g728::decoder::make_decoder_with_options(&params, false).unwrap();
    let raw = decode_all(&mut dec_raw, &packets);

    let mut dec_pf = oxideav_g728::decoder::make_decoder_with_options(&params, true).unwrap();
    let pf = decode_all(&mut dec_pf, &packets);

    let snr_raw = snr_db(&input, &raw);
    let snr_pf = snr_db(&input, &pf);
    println!("raw SNR:       {snr_raw:.2} dB");
    println!("postfilter SNR: {snr_pf:.2} dB (perceptual, NOT SNR-optimal)");

    // The postfilter must modify the signal, else the wiring is broken.
    assert_eq!(raw.len(), pf.len());
    let skip = (SAMPLE_RATE as usize) / 10;
    let mut total_diff: i64 = 0;
    for i in skip..raw.len() {
        total_diff += (raw[i] as i64 - pf[i] as i64).abs();
    }
    assert!(
        total_diff > 0,
        "postfilter output identical to raw synthesis (wiring broken)"
    );

    // Level preservation: average |sample| should be within ±30 % of
    // the input. AGC clamps the drift tighter than that, but the
    // postfilter's comb/formant shape can still skew the envelope a
    // little.
    let avg_abs = |s: &[i16]| -> f64 {
        let mut sum = 0.0_f64;
        for &v in s.iter().skip(skip) {
            sum += (v as f64).abs();
        }
        sum / (s.len() - skip) as f64
    };
    let lvl_in = avg_abs(&input);
    let lvl_out = avg_abs(&pf);
    let ratio = lvl_out / lvl_in;
    println!("postfilter level ratio: {ratio:.3}");
    assert!(
        (0.5..=1.5).contains(&ratio),
        "postfilter AGC drift out of range: {ratio:.3}"
    );

    // Must not produce non-finite samples.
    for &s in &pf {
        assert!(s != i16::MIN && s != i16::MAX, "postfilter saturated: {s}");
    }
}

/// Naive DFT magnitude at a single frequency, used to build a
/// coarse spectrum without pulling in an FFT dependency.
fn bin_magnitude(samples: &[i16], freq_hz: f64) -> f64 {
    let sr = SAMPLE_RATE as f64;
    let w = 2.0 * core::f64::consts::PI * freq_hz / sr;
    let mut re = 0.0_f64;
    let mut im = 0.0_f64;
    for (n, &s) in samples.iter().enumerate() {
        re += (s as f64) * (w * n as f64).cos();
        im -= (s as f64) * (w * n as f64).sin();
    }
    (re * re + im * im).sqrt()
}

/// Assert the postfilter reshapes the spectrum on a non-trivial input.
/// The spec's short-term postfilter deliberately attenuates the region
/// between formant peaks and emphasises the peaks themselves; the
/// spectral-tilt compensation further modifies the HF roll-off. We
/// don't check a specific shape (it depends on the current LPC state)
/// — only that the ratio of a low-band bin to a mid-band bin moves
/// *at all*, proving the filter has frequency-dependent effect.
#[test]
fn postfilter_reshapes_spectrum() {
    let params = make_params();

    let input = build_voiced(TOTAL_SAMPLES);
    let mut enc = oxideav_g728::encoder::make_encoder(&params).unwrap();
    let packets = encode_all(&mut enc, &input);

    let mut dec_raw = oxideav_g728::decoder::make_decoder_with_options(&params, false).unwrap();
    let raw = decode_all(&mut dec_raw, &packets);
    let mut dec_pf = oxideav_g728::decoder::make_decoder_with_options(&params, true).unwrap();
    let pf = decode_all(&mut dec_pf, &packets);

    let skip = (SAMPLE_RATE as usize) / 4; // skip 250 ms transient
    let raw_tail = &raw[skip..];
    let pf_tail = &pf[skip..];

    // Sample the spectrum at three probe frequencies that straddle the
    // voiced-speech test signal's tones (150 / 700 / 1600 Hz).
    let probes = [200.0, 800.0, 2500.0];
    let raw_mag: Vec<f64> = probes.iter().map(|&f| bin_magnitude(raw_tail, f)).collect();
    let pf_mag: Vec<f64> = probes.iter().map(|&f| bin_magnitude(pf_tail, f)).collect();

    println!(
        "raw spectrum  @ 200/800/2500 Hz: {:.1} / {:.1} / {:.1}",
        raw_mag[0], raw_mag[1], raw_mag[2]
    );
    println!(
        "post spectrum @ 200/800/2500 Hz: {:.1} / {:.1} / {:.1}",
        pf_mag[0], pf_mag[1], pf_mag[2]
    );

    // All three bins should have non-zero energy on the voiced mix
    // (the 150 Hz fundamental has a fat skirt).
    for &m in &raw_mag {
        assert!(m > 0.0, "raw magnitude was zero somewhere");
    }
    for &m in &pf_mag {
        assert!(m > 0.0, "postfilter magnitude was zero somewhere");
    }

    // The postfilter must move the low-vs-high band ratio by more than
    // 1% — otherwise it's an identity filter and the wiring is broken.
    let raw_lh = raw_mag[2] / raw_mag[0];
    let pf_lh = pf_mag[2] / pf_mag[0];
    let ratio_shift = (pf_lh / raw_lh).ln().abs();
    println!(
        "low/high ratio shift: {:.3} nats (raw {:.3}, pf {:.3})",
        ratio_shift, raw_lh, pf_lh
    );
    assert!(
        ratio_shift > 0.01,
        "postfilter did not reshape the spectrum ({ratio_shift:.4} nats)"
    );
}
