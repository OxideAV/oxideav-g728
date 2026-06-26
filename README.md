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
- **Annex G §G.3 fixed-point coder** (`annex_g_codebook`) — the
  bit-exact fixed-point analysis-by-synthesis codebook search of
  blocks 11 – 21, transcribed from the §G.3.4 – §G.3.10 pseudo-code one
  sub-clause at a time: block 11 VQ target (`Q2`), block 12 cascade
  impulse response `H` (`Q13`, combined-accumulator `>> 14`), block 13
  time-reversed convolution `PN` (`NLSPN = 7`), block 14/15
  filtered-shape energies `Y2` (`NLSY2 = 5`), block 16 target
  normalization via `DIVIDE` + `VSCALE`, block 17/18 `|correlation|`
  gain decision + double-precision distortion search + end-of-loop
  sign re-derivation, and block 19/21 gain-scaled excitation `ET` with
  the `NNGQ`-normalized `GQ·GAIN` product. The `encode_vector` driver
  chains all blocks for one vector. The §G.5/§G.3.9 Q-format gain
  tables (`GQ_Q13`, `GB_Q13`, `G2_Q12`, `GSQ_Q11`, `NNGQ`) are
  const-derived and cross-checked against the float tables; the search
  is proven to track the floating-point `codebook_search` shape
  decisions over a 256-target sweep on the identity cascade.
- **Annex G §G.3.11 fixed-point decoder synthesis filter**
  (`annex_g_decoder`) — block 32, the 50th-order LPC synthesis filter
  `1/A(z)` reconstructing the quantized-speech vector `ST` from the
  gain-scaled excitation `ET`. `SynthesisFilterFixed` holds `STATELPC`
  in segmented block-floating-point form (one shared NLS per
  `IDIM`-sample segment, `NLSSTATE[1..=10]` init 16 + scratch min slot)
  and runs the four §G.3.11 phases in lockstep — the zero-input
  "ringing" response with per-segment memory shift, the
  `VSCALE`-to-15-bit re-normalization, the zero-state response of `ET`
  with the `AA0 << 3` overflow probe + `LABEL1` retry, and the
  three-case `LABEL2` alignment that sums ZIR + ZSR, clips to `±4095` at
  segment scale, `VSCALE`s to 14 bits, and reverses the top `IDIM` taps
  into `ST`. Cross-checked against an in-test transcription of the
  §G.3.11 floating-point pseudo-code.
- **Annex G §G.3.20–§G.3.23 fixed-point adaptive postfilter**
  (`annex_g_postfilter`) — the decoder's bit-exact postfilter chain
  (blocks 71 – 77). `PostfilterFixed::filter_vector` runs the long-term
  pitch postfilter (`GL` Q14 / `GLB` Q16 / lag `KP` over the Q2/Q0 `SST`
  buffer) cascaded with the short-term all-zero (`AZ` Q14), all-pole
  (`AP` Q14) and spectral-tilt (`TILTZ` Q14) postfilter (blocks 71/72,
  combined per §G.3.20 with Q2 `STPFFIR`/`STPFIIR` memory), the
  sum-of-absolute-value calculators (73/74), the `SCALE = SUMUNFIL /
  SUMFIL` ratio via the §G.1.3.4 `DIVIDE` (75), and the first-order
  low-pass + output gain scaling (76/77, `SCALEFIL` Q14 init 16384,
  `AGCFAC = 16220` Q14 / `AGCFAC1 = 20972` Q21) yielding the
  postfiltered vector `SPF` (Q2). Cross-checked against an in-test
  transcription of the §G.3.20 floating-point short-term postfilter.
- **Annex G §G.3.17–§G.3.19 fixed-point backward synthesis-filter adapter**
  (`annex_g_synth_adapter`) — the decoder's bit-exact backward LPC
  adaptation chain (blocks 49 / HWMCORE / 51), the fixed-point analogue of
  the floating-point `synthesis_adapter`. `SynthAdapterFixed::adapt` runs
  block 49 (the segmented-block-floating-point hybrid window on past
  quantized speech `SB`, with the `NRSH` per-segment alignment to the
  common `NLSTMP` and the Q15-window-multiply `−1` compensation), HWMCORE
  (the recursive 3/4-decay + non-recursive autocorrelation accumulation
  across the three `NLSRREC`-vs-`NLSAA0` alignment cases, the `>> 8`
  white-noise correction, the `VSCALE`/`MLS = 30` normalization, and the
  32-bit `R(LPC+1) = 0` `ILLCOND` interoperability test), block 50 (the
  §G.2.2 `levinson_durbin_fixed`), and block 51 (`FACV` bandwidth
  expansion to the live Q14 predictor `A`, with the `NLSATMP`-driven
  `<< 3/2/1 → Q30` shift and the keep-old-`A`-on-overflow/`ILLCOND`
  path). The emitted Q14 `A` is exactly the coefficient set the §G.3.11
  block-32 fixed-point synthesis filter consumes — closing the fixed-point
  decoder's backward-adaptation loop. Cross-checked against the §G.3.17
  floating-point hybrid window (`HybridWindowState`) on identical
  requantized speech. The block-51 expansion is factored into a shared
  `bandwidth_expand_q14` core that also backs the §G.3.15 block-45
  log-gain bandwidth expansion (`log_gain_bandwidth_expand`, `FACGPV`,
  `LPCLG = 10`, `ILLCONDG`-gated) — the first fixed-point piece of the
  backward *gain* adapter.
- **Annex G §G.3.17/§G.3.18 shared fixed-point hybrid-window core**
  (`annex_g_hybrid`) — the §G.3.17 windowing step + §G.3.18 HWMCORE,
  factored out of the LPC-hard-coded synthesis adapter into a
  dimension-agnostic `HybridWindowFixed` / `HybridWindowFixedState` (the
  fixed-point analogue of the floating-point `HybridWindowState`). It
  carries the permanent SBFL signal-history `SB` / `NLSSB` and the BFL
  recursive-autocorrelation `REXP` / `NLSREXP`, and runs the
  segment-shift / window-multiply / three-NLS-case
  recursive+non-recursive autocorrelation / white-noise-correction chain
  for any `(order, l, n, window)` — emitting `RTMP(1..=order+1)` plus the
  §G.2.2 `ILLCOND` verdict. `SynthAdapterFixed` (block 49) and the new
  `WeightAdapterFixed` (block 36) both drive it.
- **Annex G §G.3 fixed-point perceptual-weighting-filter backward adapter**
  (`annex_g_weight_adapter`) — the encoder's `W(z) = Q̃(z/γ₁) / Q̃(z/γ₂)`
  backward adaptation in bit-exact fixed point, the analogue of the
  floating-point `weighting_filter_adapter` (blocks 36/37) +
  `weighting_filter_coeff` (block 38). `WeightAdapterFixed::adapt` chains
  block 36 (the shared hybrid-window core at `LPCW = 10` / `NONRW = 30` /
  `WNRW`), block 37 (the §G.2.2 `levinson_durbin_fixed` at order `LPCW`),
  and block 38 (two `bandwidth_expand_q14` expansions: the numerator with
  `WZCFV_Q14`, the denominator with `WPCFV_Q14`, eq. 3-4b/3-4c), emitting
  the Q14 `WeightCoeffFixed` pair. Keeps the previous cycle's pair on
  ill-conditioning / Q14 overflow; `disabled` collapses to the §3.4.1
  non-speech `W(z) = 1` all-pass. Cross-checked against the
  floating-point chain on identical requantized speech.

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
  module. The §G.2.2 variable-precision Levinson-Durbin recursion is
  landed in the `annex_g_levinson` module — the bit-exact reformulation
  of the three recursion blocks (37 / 44 / 50): the `SIMPDIV` 16-bit
  restoring division, the `Q15→Q14→Q13` `NRS` overflow-rescale strategy
  for the predictor coefficients, the saved-high-word `ALPHATMP`, the
  `NLSATMP = 15 − NRS` precision signal to the bandwidth-expansion
  module, the `ILLCOND` / `ILLCONDP` ill-conditioning flags, and the
  decoder `MINC0 = 10` restart — proven against the floating-point
  `levinson` reference and with a bit-exact resume-vs-one-shot test.
  The §G.3 fixed-point **coder** (analysis-by-synthesis codebook
  search), the §G.3.11 decoder **synthesis filter** (block 32), and the
  §G.3.20–§G.3.23 adaptive **postfilter** (blocks 71–77) are now landed
  in the `annex_g_codebook` / `annex_g_decoder` / `annex_g_postfilter`
  modules (see the Implemented section). The synthesis-filter backward
  adaptation chain (block 49 hybrid window, HWMCORE, block 51 bandwidth
  expansion) is now landed in the `annex_g_synth_adapter` module (see the
  Implemented section), closing the fixed-point decoder's backward
  adaptation loop with the already-landed §G.3.11 block-32 synthesis
  filter. Block 45 (the §G.3.15 log-gain bandwidth expansion) is also
  landed via the shared `bandwidth_expand_q14` core. The block 49 / block
  36 hybrid window + HWMCORE is now factored into the shared
  `annex_g_hybrid` core, and the *perceptual-weighting-filter* backward
  adapter (blocks 36 / 37 / 38) is landed in `annex_g_weight_adapter`
  (see the Implemented section). The remaining §G.3 per-block modules —
  the *log-gain* hybrid window + Levinson wiring (block 43, which reuses
  the same hybrid core at the `LPCLG` dimensions; note its 4-scalar
  log-gain input differs from the segmented-speech input of blocks 36/49),
  the §G.3.16 block-46 log-gain linear prediction, and the postfilter
  *coefficient* calculators (block 81 LPC inverse, 82 pitch extraction,
  83 pitch tap, 84 long-term coeff, 85 short-term coeff) — stay deferred
  behind the floating-point build. The full fixed-point log-gain adapter
  (blocks 43–48) additionally requires the §G.3.12–§G.3.16 log-gain
  Q-format scaling, which the staged §G.3 trace flags as a documented gap
  (per-module Q scaling / rounding / saturation not reproduced).

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
