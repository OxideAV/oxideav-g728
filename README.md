# oxideav-g728

Pure-Rust ITU-T G.728 Low-Delay CELP (LD-CELP, 16 kbit/s) speech codec.
Zero C dependencies, no FFI, no `*-sys` crates.

## Status

Clean-room rebuild from the published ITU-T G.728 (1992-09)
Recommendation prose, built one spec-cited unit at a time. The codec
is implemented as a **standalone `Encoder` / `Decoder` API**; it does
not yet register a trait-surface codec in the `oxideav-core` runtime
(the registry hook is a no-op pending A-law / µ-law PCM I/O wiring).

Implemented end-to-end:

- **Autonomous decoder** (`Decoder::decode_vector`) — runs the full
  Figure 3/G.728 per-vector chain (excitation lookup, gain scaling,
  synthesis filter, the two backward adapters) off the raw 10-bit
  channel-index stream with all backward adaptation handled internally.
- **§4.7 postfilter chain** (`Decoder::decode_vector_postfiltered`) —
  the complete long-term (pitch) + short-term (spectral) + AGC
  pipeline: LPC inverse filter (block 81) → pitch extractor (82) →
  coefficient calculator (83/84) → long-term comb (71) → short-term
  pole-zero + tilt (72) → AGC (73–77).
- **Complete encoder loop** (`Encoder::encode_vector`) — the full
  analysis-by-synthesis Figure 2/G.728 dataflow: σ(n) prediction,
  perceptual weighting filter (blocks 4 / 36 / 37 / 38), zero-input
  response + VQ target (9 / 10 / 11), §3.9 codebook search (12–18),
  excitation + gain scaling (19 / 21), and §3.10 memory update. The
  encoder runs in **bit-exact lockstep** with the decoder (proven by a
  200-vector integration test asserting the decoder reproduces the
  encoder's quantized speech `sq(n)` bit for bit).
- **ICOUNT-faithful update stagger** — the backward synthesis filter
  and weighting coefficients swap exactly at the third vector of each
  adaptation cycle per the spec's ICOUNT = 3 rule.
- **§3.11 in-band signalling** — `encode_vector_with_sync_bit` /
  `extract_sync_bit` rob the leftmost codeword bit of a chosen vector
  for a synchronization / side-channel bit via a half-codebook search.
- **§5.11 serial byte-stream framing** — `bitstream::pack_indices` /
  `unpack_indices` (MSB-first 10-bit-per-vector packing) and the
  whole-stream wrappers `Encoder::encode_buffer` /
  `Decoder::decode_bytes`. The natural framing unit is one adaptation
  cycle = 4 indices = 40 bits = 5 bytes.
- **Annex I §I.4.1 frame-erasure excitation extrapolation** —
  `FrameErasureExcitation` (blocks 31SF / 31FE / 31E) maintains the
  `ETPAST()` gain-scaled-excitation history and synthesises one `ET()`
  vector per erased vector: at onset the last good `PTAP` (vs lowered
  `VTH = PPFTH/1.4`) picks **voiced** (periodic decaying repeat of the
  last `KP` samples, `FESCALE` stepping `0.8 → 0` over 50 ms) or
  **unvoiced** (random `5..=140`-sample slide-back, magnitude-matched to
  the pre-erasure `AVMAG`, attenuated over 60 ms). The "random"
  slide-back is a pluggable `SlideSource` (default dependency-free LCG)
  since the spec licenses any aperiodic `5..=140` sequence during
  erasures. Floating-point only.
- **Annex I §I.4.2 frame-erasure LPC softening** — `soften_predictor`
  + `FrameErasureLpc` step bookkeeping (progressive `(0.97)^{k·i}`
  bandwidth expansion of the last good predictor during erased frames).
- **Annex I §I.4.5 gain-growth limit** — `limit_log_gain` (Block 47AF,
  §I.5.6) clamps the backward-adapted log-gain to at most `FEGAINMAX =
  +2 dB`/vector above the last log-gain for the first few good frames
  after an erasure, plus `GainGrowthLimiter` driving the §I.5.1
  `AFTERFE` / `FECOUNT` / `OGAINDB` bookkeeping (clamp lasts the erasure
  length, saturated at `AFTERFEMAX = 16` cycles = 40 ms).

### Not yet implemented

- Registry-side `decoder` / `encoder` factory wiring (waits on
  A-law / µ-law PCM I/O per §5.3 / §3.1, handled by `oxideav-g711`).
- The remaining Annex I concealment mechanisms (§I.4.3 continued
  backward adaptation, §I.4.4 floating post-filter) and the decoder
  erasure-flag drive path that wires §I.4.1/§I.4.2/§I.4.5 together.
- Annex G fixed-point variant — the §G.1.2 numerical representations
  and the §G.1.3 arithmetic primitives (`shr_sat` / `shl_sat`, `rnd`,
  `findnls` / `vscale`, `ScalarFloat`, `divide`) are landed in the
  `annex_g_arith` module as the bit-exact foundation, and the §G.2.1
  reformulated backward vector gain adapter — the two log-gain dB tables
  (`gain_log_db` = `20·log10|g_i|`, `shape_log_db` = `10·log10·P[y_j]`,
  Figure G.1 blocks 93/94) plus the adder-96/limiter-97
  `offset_removed_log_gain` (eq. G-14) — is landed in the `annex_g_gain`
  module. The §G.2.2 variable-precision Levinson-Durbin changes and the
  §G.3 per-block fixed-point pseudo-code that build on this foundation
  remain deferred behind the floating-point build.

## Usage

```rust
use oxideav_g728::{make_decoder, make_encoder};

let mut enc = make_encoder();
let mut dec = make_decoder();

// One IDIM = 5-sample vector per call.
let vector = [0.0f64; 5];
let ichan = enc.encode_vector(&vector).unwrap();   // 10-bit channel index
let pcm = dec.decode_vector(ichan);                // 5 reconstructed samples
```

Or drive the whole serial stream at once:

```rust
use oxideav_g728::make_encoder;

let mut enc = make_encoder();
let vectors = vec![[0.0f64; 5]; 8];          // 8 vectors = 2 adaptation cycles
let bytes = enc.encode_buffer(&vectors).unwrap();  // §5.11 byte stream
```

## Clean-room provenance

Every numeric value lives in `src/tables.rs`, transcribed directly
from the integer columns of Annexes A–D of the staged ITU-T G.728
(1992-09) PDF under `docs/audio/g728/`. Derived quantities (the
`Y_ENERGY` shape-codevector energy table, the `WeightingFilterCoeff`
bandwidth-broadened coefficients) are computed at compile time from
the transcribed columns per their spec equations, each with a per-test
cross-check. Every control-flow line in the per-block modules carries
a comment pointing at its §3 / §4 / §5 pseudo-code, and the per-test
cross-checks act as in-repo audit anchors against transcription typos.

## License

MIT — see [LICENSE](LICENSE).
