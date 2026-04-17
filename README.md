# oxideav-g728

Pure-Rust **ITU-T G.728 Low-Delay CELP** (LD-CELP, 16 kbit/s) codec —
decoder + encoder, exhaustive 128×8 analysis-by-synthesis shape/gain
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

G.728 is 8 kHz mono · 5 bytes per 2.5 ms vector (5 × 10-bit codebook
indices packed) · 16 kbit/s · low-delay (< 2 ms algorithmic latency).

```rust
use oxideav_codec::CodecRegistry;
use oxideav_core::{CodecId, CodecParameters, Frame, Packet, TimeBase};

let mut reg = CodecRegistry::new();
oxideav_g728::register(&mut reg);

let mut params = CodecParameters::audio(CodecId::new("g728"));
params.sample_rate = Some(8_000);
params.channels = Some(1);

let mut dec = reg.make_decoder(&params)?;
// A typical packet bundles multiple 5-byte vectors.
dec.send_packet(&Packet::new(0, TimeBase::new(1, 8_000), g728_bytes))?;
let Frame::Audio(a) = dec.receive_frame()? else { unreachable!() };
// `a.data[0]` is S16 PCM at 8 kHz mono (5 samples per 5-byte input
// vector = 1.0 expansion in samples-per-bit terms).
# Ok::<(), oxideav_core::Error>(())
```

Encoder is symmetric via `reg.make_encoder(&params)`.

### Codec IDs

- `"g728"`

## License

MIT — see [LICENSE](LICENSE).
