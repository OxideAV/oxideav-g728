# oxideav-g728

Pure-Rust **ITU-T G.728 Low-Delay CELP** (LD-CELP, 16 kbit/s) codec —
decoder + encoder, exhaustive 128 x 8 analysis-by-synthesis shape/gain
search. G.728 specifies a single 16 kbit/s operating point; this crate
covers both directions at that point. Zero C dependencies.

Part of the [oxideav](https://github.com/OxideAV/oxideav-workspace)
framework but usable standalone.

## Installation

```toml
[dependencies]
oxideav-core = "0.0"
oxideav-codec = "0.0"
oxideav-g728 = "0.0"
```

## Quick use

G.728 is narrowband (8 kHz) mono, fixed at 16 kbit/s. Each codeword
covers 5 samples (0.625 ms) and is 10 bits wide. The encoder groups four
consecutive codewords into a 40-bit / 5-byte packet spanning 20 samples
(2.5 ms) — that is the smallest packet unit the encoder emits, but the
decoder accepts any byte buffer that holds one or more whole 10-bit
codewords (trailing partial bits are dropped).

```rust
use oxideav_codec::CodecRegistry;
use oxideav_core::{CodecId, CodecParameters, Frame, Packet, SampleFormat, TimeBase};

let mut reg = CodecRegistry::new();
oxideav_g728::register(&mut reg);

let mut params = CodecParameters::audio(CodecId::new("g728"));
params.sample_rate = Some(8_000);
params.channels = Some(1);
params.sample_format = Some(SampleFormat::S16);

let mut dec = reg.make_decoder(&params)?;
dec.send_packet(&Packet::new(0, TimeBase::new(1, 8_000), g728_bytes))?;
let Frame::Audio(a) = dec.receive_frame()? else { unreachable!() };
// a.data[0] is S16LE PCM at 8 kHz mono; 5 samples per 10-bit input index.
# Ok::<(), oxideav_core::Error>(())
```

### Encoder

```rust
let mut enc = reg.make_encoder(&params)?;
enc.send_frame(&Frame::Audio(s16_mono_8khz_frame))?;
let pkt = enc.receive_packet()?;
// Every encoder packet is exactly 5 bytes = 4 x 10-bit indices = 2.5 ms.
```

The encoder only accepts 8 kHz mono S16 input; any other sample rate,
channel count, or sample format is rejected at construction. This is a
fixed-rate codec, not a configurable one — there are no bitrate or
channel knobs to turn.

### Codec IDs

- `"g728"`; single capability entry `"g728_sw"`.

## Status and caveats

What is implemented:

- Backward-adaptive 50th-order LPC synthesis predictor (windowed
  autocorrelation + Levinson-Durbin with bandwidth expansion), refreshed
  every four vectors (2.5 ms).
- Backward-adaptive 10th-order log-gain predictor.
- 10-bit MSB-first bit reader / packer: 7-bit shape index, 1-bit sign,
  2-bit gain magnitude.
- **ITU-T G.728 Annex B codebooks**: the 128 x 5 `CODEBK` shape table
  (Q11 integer values divided by 2048) and the 4-entry `GQ` gain
  magnitude table (`33/64`, `231/256`, `1617/1024`, `11319/4096`),
  transcribed verbatim from the 09/92 recommendation.
- Exhaustive 128 x 8 analysis-by-synthesis search in the encoder, sharing
  its LPC / gain machinery verbatim with the decoder so a bitstream
  produced here round-trips cleanly through the in-tree decoder.
- Encode -> decode round-trip coverage (sine carries non-trivial energy,
  silence decays after the initial transient, PTS rises monotonically).
- Round-trip SNR sanity floor: 400 Hz sine reconstructs at ~44 dB SNR,
  3-tone voiced-speech mix at ~30 dB, both measured after a 100 ms
  adaptation transient.

What is deliberately not shipped yet:

- The spec's recursive Barnwell / logarithmic autocorrelation window;
  this crate uses a fixed 100-sample Hamming window instead.
- The adaptive long-term (pitch) and short-term postfilters from the
  2012 edition, ITU-T G.728 section 5.5.

Practical consequence: the excitation codebooks are spec-accurate, but
the surrounding LPC / log-gain adaptation still uses a simplified window
and no postfilter, so the decoder is not yet bit-compatible with ITU
reference streams. It is deterministic, stable on long zero-excitation
runs, and suitable for pipeline testing and round-trip transcodes within
oxideav.

## License

MIT — see [LICENSE](LICENSE).
