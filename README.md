# oxideav-g728

Pure-Rust **ITU-T G.728 Low-Delay CELP** (LD-CELP, 16 kbit/s) codec —
decoder + encoder, exhaustive 128×8 analysis-by-synthesis shape/gain
search.

G.728 specifies a single 16 kbit/s operating point; this crate covers
both directions at that point. Zero C dependencies, no FFI, no `*-sys`
crates.

Originally part of the [oxideav](https://github.com/KarpelesLab/oxideav)
framework; extracted to its own crate for independent publication.

## Usage

```toml
[dependencies]
oxideav-g728 = "0.0.3"
```

Plugs into [`oxideav-codec`](https://crates.io/crates/oxideav-codec):

```rust
let mut reg = oxideav_codec::CodecRegistry::new();
oxideav_g728::register(&mut reg);
```

Decoder id: `"g728"` — input packets are 10-bit indices packed into
bytes (five 10-bit indices per 20 ms frame = 5 bytes at 16 kbit/s);
output is 16-bit PCM at 8 kHz mono.

## License

MIT — see [LICENSE](LICENSE).
