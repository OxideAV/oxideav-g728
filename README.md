# oxideav-g728

[![CI](https://github.com/OxideAV/oxideav-g728/actions/workflows/ci.yml/badge.svg)](https://github.com/OxideAV/oxideav-g728/actions/workflows/ci.yml) [![crates.io](https://img.shields.io/crates/v/oxideav-g728.svg)](https://crates.io/crates/oxideav-g728) [![docs.rs](https://docs.rs/oxideav-g728/badge.svg)](https://docs.rs/oxideav-g728) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Pure-Rust ITU-T G.728 Low-Delay CELP (LD-CELP, 16 kbit/s) speech codec.
Zero C dependencies, no FFI, no `*-sys` crates.

## Status

Clean-room rebuild from the published ITU-T G.728 (1992-09)
Recommendation prose + Annex G (1994-11, fixed point) + Annex I
(1999-05, PLC), built one spec-cited unit at a time.

**Conformance — the Annex G fixed-point coder is BIT-EXACT on the
official ITU-T G.728 test vectors** (the Appendix I sequences,
input → reference-output black-box pairs, gated by
`tests/conformance.rs` whenever the private docs tree is present):

| Leg | Vectors | Result |
|-----|---------|--------|
| fixed encoder → `incw1g…incw6g` | 98 048 indices | **bit-exact** |
| fixed decoder (postfilter off) → `outa1g…outa6g` | 491 520 samples | **bit-exact** |
| fixed decoder (postfilter on) → `outb4g` | 51 200 samples | **bit-exact** |
| float decoder (postfilter off) → `outa1…outa6` | 491 520 samples | **bit-exact** |
| float encoder → `incw1…incw4` | 14 336 indices | **bit-exact** |
| float encoder → `incw5` / `incw6` | 84 736 indices | near-tie prefix pinned (float references are reproducible only up to the reference codec's own floating-point tie behaviour) |
| float decoder (postfilter on) → `outb4` | 51 200 samples | 85 % exact, max ±3 LSB (same float precision bound) |

The codec is exposed both as **direct `Encoder` / `Decoder` /
`EncoderFixed` / `DecoderFixed` APIs** and as a registered
`oxideav-core` codec: `register(ctx)` wires `g728` decoder + encoder
factories (§5.11 byte-stream packets ↔ S16 mono 8 kHz frames through
the conformance-bit-exact fixed-point chains; WAVE tag `0x0041`).

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
- **Annex G §G.3 fixed-point short-term postfilter coefficient calculator**
  (`annex_g_pf_coeff`) — the §4.6 block-85 calculator (eq. 4-3/4-4/4-5) in
  bit-exact fixed point, producing the Q14 `ShortTermCoeff` (`AZ` / `AP` /
  `TILTZ`) the §G.3.20–§G.3.23 postfilter consumes from the order-10
  synthesis-filter LPC by-product. `AZ` (`ã_i·0.65^i`, `SPFZCFV_Q14`) and
  `AP` (`ã_i·0.75^i`, `SPFPCFV_Q14`) are two `bandwidth_expand_q14`
  expansions of the same predictor (the postfilter analogue of weighting
  block 38); `TILTZ = 0.15·k1` is the Q14 tilt term (`TILTF_Q14 = 2458`,
  `k1` the §G.2.2 Q15 `RC1`). Keeps the previous coefficients on Q14
  overflow / out-of-range `NLSATMP`. Cross-checked against the
  floating-point `set_from_synthesis_byproduct`.
- **Annex G fixed-point backward vector gain adapter** (`annex_g_gain_adapter`)
  — the complete Figure G.1 chain: the §G.3.14 block-43 flat-Q9 log-gain
  hybrid window (`LogGainWindowFixed`, driving the shared HWMCORE), the
  §G.3.16 per-vector blocks — 46 (`log_gain_predict`), 98
  (`limit_log_gain_q9`, −32…+28 dB), 99 + 48 (`inverse_log_gain`: `GOFF`
  add + the Q15 Horner Taylor `2^x` antilog → scalar-floating `(GAIN,
  NLSGAIN)`), 93/94/96/97 (`gstate1_update` via the §G.5 Q11 dB tables
  `GCBLG_Q11` / `SHAPELG_Q11`) — and the `GainAdapterFixed` driver with
  the §G.7 commit timing. σ tracks the floating-point adapter to
  ≈ 0.008 dB steady-state on an encoder-style index stream.
- **Annex G §G.3.24–§G.3.27 fixed-point long-term-postfilter adaptation**
  (`annex_g_pitch`) — blocks 81 (BFL `ST` → Q2 `SST` + Q13 LPC inverse
  filter into the Q1 residual `D`), 82 (1 kHz low-pass + 4:1 decimation,
  two-stage correlation peak-picking, cross-multiplied `ITAPTH`
  sub-multiple decision), 83 (`PTAP` Q14 via `DIVIDE`) and 84 (`GL` Q14 /
  `GLB` Q16), plus the §G.3.28 `LABEL` `apf_to_q13` conversion.
- **Annex G §G.3.1–§G.3.3 fixed-point encoder filter loop**
  (`annex_g_encoder`) — block 4 (Q2 perceptual weighting filter),
  blockzir (segmented-BFL synthesis ZIR + Q2 weighting ZIR) and the
  blocks-9/10 memory update (the `LABEL1` `ET >> 1` overflow-retry
  probe, the ±4095 segment-scale clip, `VSCALE` re-normalizations and
  the reversed `ST` / `NLSST` emission), cross-checked against in-test
  float transcriptions and a direct-form cascade.
- **Annex G §G.7 fixed-point encoder + decoder main programs**
  (`annex_g_coder`) — `EncoderFixed` / `DecoderFixed` chain every §G.3
  block in the §G.7 execution order (block-51/38 commits at
  `ICOUNT = 3`, block 45 at 2, blocks 49/50 at 4, blocks 36/37 at 2,
  blocks 43/44 at 1; the decoder interrupts block 50 at order 10 for
  the `APF`/`RC1` postfilter by-product and resumes). The decoder
  reproduces the encoder's quantized speech `ST` **bit-exactly**
  (mantissas + `NLSST`) from the channel indices alone over a
  200-vector lockstep test; the postfiltered output tracks the input at
  correlation > 0.9 within 3 dB.

### Errata

The Annex-G fixed-point pseudo-code carries several print typos that a
literal transcription would inherit; this crate follows the
conformance-verified corrections catalogued in
`docs/audio/g728/g728-errata.md`. Each correction is anchored at its
source line in `src/` with an `E1…E6` / `N1` / `N2` reference:

| Id | §        | Correction |
|----|----------|------------|
| E1 | G.3.18 `HWMCORE` case 1 | `RREC` recursion is a `−(RREC<<NLSATT)+(RREC<<16)` decay, not the printed `+` growth |
| E2 | G.3.12 block 36 | hybrid-window decay is `1/2` (`NLSATTW=15`), not the sibling blocks' `3/4` |
| E3 | G.3.24 block 81 | `IP` reset target is `NPWSZ − NFRSZ`, not `NFRSZ` |
| E5 | G.3.20 block 72 | short-term IIR store sources `AA1` (all-pole result), not `AA0` |
| E6 | G.3.2 `Blockzir` block 9 | middle MAC loop accumulates (`AA0 = AA0 − …`), not resets (`AA0 = 0 − …`) |
| N1 | G.3.5 block 12 | impulse response consumes `A(2..5)`, never `A(6)` |
| N2 | G.3.28 block 85 | `TILTF = 4915` in Q15, not `2458` in Q14 |

(E4 — the §G.3.25 refined-pitch first-correlation-factor ambiguity — is
resolved here in the equivalent `D`-autocorrelation form `Σ D(K)·D(K−J)`,
bit-exact against the vectors.)

### Not yet implemented

- The remaining Annex I concealment mechanisms (§I.4.3 continued
  backward adaptation, §I.4.4 floating post-filter) and the §I.5.1
  decoder erasure-flag drive path (`VEC_LOOP`) that wires the existing
  §I.4.1 excitation extrapolation, §I.4.2 LPC softening and §I.4.5
  gain-growth limiter modules together. **Docs gap:** the Annex I PLC
  path cannot be conformance-adjudicated from the in-repo corpus — the
  staged `conformance/mask1` supplies a packet-loss mask but no paired
  concealed-PCM reference output, so a PLC drive path can only be
  behaviourally validated, not proven bit-exact, until the ITU Annex I
  concealed-PCM reference is staged (see the outstanding-gap note in
  `docs/audio/g728/g728-errata.md`).

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
