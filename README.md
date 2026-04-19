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
oxideav-core = "0.1"
oxideav-codec = "0.1"
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

- Backward-adaptive 50th-order LPC synthesis predictor using the spec's
  **§3.7 hybrid (Barnwell / logarithmic) autocorrelation window** — a
  35-sample non-recursive tail plus a decaying recursive portion with
  `alpha^{2L} = 3/4`, using the 105-sample window table from Annex A.1
  and `FAC = 253/256` bandwidth expansion. Refreshed every four vectors
  (2.5 ms) with Levinson-Durbin + the spec's 257/256 white-noise
  correction.
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
- **§5.5 adaptive postfilter**: long-term (pitch comb) + short-term
  (pole-zero formant emphasis with spectral-tilt compensation) + AGC
  level renormalisation. The 10th-order postfilter LPC is extracted as
  a by-product of the main 50th-order Levinson-Durbin. Enabled by
  default; pass `make_decoder_with_options(&params, false)` for raw
  synthesis output (recommended by the spec §4.6.1 for non-speech
  signals like modems).
- Encode -> decode round-trip coverage (sine carries non-trivial energy,
  silence decays after the initial transient, PTS rises monotonically).
- Round-trip SNR (raw-synthesis floor, postfilter off): 400 Hz sine at
  ~44 dB, 3-tone voiced-speech mix at ~34 dB, both measured after a
  100 ms adaptation transient. The postfilter is a perceptual shaper:
  it lowers L2 SNR on pure tones (~19 dB on the voiced mix) while
  raising subjective quality on speech.

The remaining deviation from strict ITU bit-compat is the log-gain
predictor's window — the spec's separate hybrid window for the
10th-order gain trajectory (block 43) is still driven by a Hamming-
windowed accumulator in this crate. That is a perceptually invisible
stability tweak that does not affect round-trip quality in practice,
but it does mean bitstream-level comparisons against the ITU reference
decoder will drift on long speech runs.

Pipeline use is otherwise clean: deterministic, stable on long zero-
excitation runs, and suitable for round-trip transcodes within
oxideav.

## License

MIT — see [LICENSE](LICENSE).
