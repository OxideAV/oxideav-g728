# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

- **Annex I §I.5.1 `VEC_LOOP` frame-erasure driver (float decoder).**
  `Decoder::conceal_vector()` / `Decoder::conceal_vector_postfiltered()`
  conceal one lost vector each through the full Annex I loop: block
  31SF flags/FESCALE at erasure onset and every 10 ms, block 31FE
  voiced/unvoiced excitation extrapolation, block 31E `ETPAST()`
  update (now also maintained on every good vector, as §I.5.4
  requires), §I.4.2 LPC softening via block 51FE (`0.97^i` over the
  last good `ATMP`, compounding on the live `A` every further 10 ms),
  blocks 97FE + 43FE gain-adapter vital operations, block 49FE +
  order-10-only block 50 at cycle end (the postfilter by-product keeps
  "floating" per §I.4.4), the first-good-`ICOUNT = 3` block-51 resync
  on the current mixed `ATMP`, and the §I.4.5 / 47AF gain-growth clamp
  (`FEGAINMAX = 2 dB`/vector for `min(FECOUNT, AFTERFEMAX)` good
  cycles). Erasure granularity is one adaptation cycle (`FESIZE = 1`
  semantics; conceal 4 calls per lost 2.5 ms cycle). New audit
  accessors `Decoder::fe_counters()` /
  `Decoder::frame_erasure_excitation()`; new
  `SynthesisAdapter::{adapt_erased, raw_atmp, illcond,
  expand_current}`, `GainAdapter::{advance_erased,
  predict_next_limited}` and
  `HybridWindowState::{advance_erased, run_erased}` building blocks.
  The never-erased path is bit-identical (all six ITU conformance legs
  still pass bit-exactly).

  **Conformance caveat:** the concealed output itself is
  *non-adjudicable in-repo* — the ITU corpus stages `mask1` but no
  paired concealed-PCM reference, upstream builds no Annex-I decoder
  variant, and the Annex I PDF ships no I/O vectors (docs
  `g728-errata.md`, issue #232). `tests/annex_i.rs` pins the
  structural properties instead (prefix bit-exactness up to the first
  erasure of a `mask1`-driven decode of the official `cw4` bitstream,
  FESCALE silencing schedule, clamp arm/release/saturation, recovery,
  hostile-mask robustness).

- **Anchored regression tests for errata E3–E6 (docs answer, r408).**
  The four Annex-G fixed-point pseudo-code errata are now each pinned by
  a standalone unit test that runs everywhere (no conformance corpus
  needed): the block-81 `IP` rewind trajectory (E3), the block-82
  refined undecimated peak-pick argmax identity (E4), the block-71/72
  `STPFIIR(1)` all-pole-source recursion (E5), and the Blockzir block-9
  middle-loop accumulation against an independent flat-scale ZIR
  transcription (E6). E4 additionally gains its source-line anchor
  comment (the misprinted factor is the *array name* — the printed
  `DEC(K − J)` cannot even be indexed at `K − J`; the paired
  floating-point listing prints the corrected `D(K)·D(K − J)` form
  literally), plus a `PitchAdapterFixed::ip()` audit accessor.

### Fixed

- **Block-31FE slide-back was off by one vector (§I.5.3)** — the
  frame-erasure excitation extrapolator read its 5-sample segment
  relative to the *previous* vector's start instead of the current
  vector position, so a slide-back of `d` effectively reached `d + IDIM`
  samples into the past: voiced concealment repeated the excitation with
  period `KP + 5` instead of `KP`. The index base is now the spec's
  `ETPAST(I − FEDELAY)` mapping (spec slot 0 = the newest stored
  sample), and a new regression test drives the full
  extrapolate-then-push loop to pin the generated stream's period at
  exactly `KP`.

- **Block-83 PTAP shift-overflow panic guard (§G.3.27)** — the
  fixed-point pitch-tap calculator's final `PTAP >> (NLSPTAP − 14)` used
  a raw `i16` shift amount that reaches the type width and panics
  ("shift with overflow") on a degenerate or corrupted codeword stream.
  `NLSPTAP` stays in `[14, 29]` on the conformance path, so the in-range
  arithmetic is unchanged (all six ITU legs remain bit-exact); only the
  out-of-range shift amounts are guarded (`≥16 → 0`; `≤−2 → 16384`, the
  Q14 unit ceiling — both unreachable on valid input). Surfaced by the
  new robustness suite.

- **Block 12 impulse-response off-by-one (§G.3.5) — fixed encoder now
  BIT-EXACT against the official ITU-T conformance vectors (r392)** —
  `annex_g_codebook::impulse_response` read the caller's 0-based
  coefficient arrays (`a[0] = A(1)` unity head) with the spec's 1-based
  index, so the cascade filter used `A(3..6)`/`AWZ(3..6)`/`AWP(3..6)`
  in place of `A(2..5)`… — wrong `H`, hence wrong `Y2`/`PN` and ~85 %
  wrong codebook decisions once the backward adapters went live. With
  the fix the full fixed-point encoder free-runs **byte-exactly** onto
  `incw1g…incw6g` (98 048 consecutive channel indices, all six official
  test cases).
- **`TILTF` is 4915 in Q15, not 2458 in Q14 (§G.3.28) — fixed decoder
  postfilter now BIT-EXACT (r392)** — block 85's spectral-tilt term is
  `TILTZ = RND(TILTF·RC1)` with the §G.3.28 inline constant
  `TILTF = 4915` (Q15, the *truncation* of 0.15·2¹⁵ = 4915.2); the
  previous Q14 `2458` rendering (≡ 4916 Q15, with a different rounding
  point) left ~19.5 % of postfiltered samples ±1 Q2 LSB off. With the
  fix the postfiltered decode of `cw4` matches `outb4g` byte-exactly
  (51 200 samples). Exported constant renamed `TILTF_Q14` →
  `TILTF_Q15`.

- **Float pitch extractor locked onto pitch doubles — §4.7 block-82
  fundamental-vs-multiple gate is one-sided (r392)** — the block-82
  pseudo-code only runs the second (fundamental-candidate) search when
  the raw correlation winner `p0` *exceeds* `p̂ + KPDELTA` ("If KP <
  M2 + 1, go to LABEL"); the float `PitchSearch` used a symmetric
  neighbourhood test that also fired when `p0 ≪ p̂` and — since
  β1 ≈ β0 on periodic signals — replaced a correct fundamental with a
  stale-neighbourhood multiple. On the official `cw4` decode the float
  chain tracked exactly 2× the Annex-G reference pitch (until 2·KP left
  `[KPMIN, KPMAX]`), pushing the postfiltered output up to 4 041 16-bit
  LSBs off `outb4`. With the one-sided gate the float postfiltered
  decode agrees with `outb4` on 43 664 / 51 200 samples with max
  |Δ| ≤ 3 — the reference codec's own floating-point precision bound.

### Documentation

- **Anchored the conformance-verified Annex-G errata at their source
  lines.** The staged `docs/audio/g728/g728-errata.md` formalizes the
  fixed-point pseudo-code print typos this crate already corrects; each
  correction now carries an `E1`/`E2`/`E3`/`E5`/`E6`/`N1`/`N2` reference
  at its `src/` site, a stale `TILTF_Q14 = 2458` module doc-comment was
  corrected to `TILTF_Q15 = 4915` (Q15), and the README's stale "Not yet
  implemented" section was replaced with an errata table and the
  Annex-I PLC-reference docs-gap note.

### Added

- **`tests/robustness.rs` — self-consistency hardening gates.** No ITU
  reference material; asserts invariants the private vectors never
  exercise (the full 1024-codeword alphabet + long adversarial LCG
  streams): the encoder never emits an out-of-10-bit-range index;
  neither decoder path panics, emits `NaN`/`Inf` or lets the
  backward-adapted state run away; fixed + float encode→decode
  round-trips stay bounded.

- **`oxideav-core` registry wiring (r392)** — `register(ctx)` now
  registers a real `g728` `CodecInfo` (decoder + encoder factories,
  WAVE-format tag `0x0041`): the decoder consumes the §5.11 serial byte
  stream (5 bytes per adaptation cycle, buffered across packet
  boundaries) through the conformance-bit-exact Annex-G fixed-point
  chain with the adaptive postfilter and emits S16 mono 8 kHz frames;
  the encoder consumes S16 frames (whole-cycle buffering, zero-padded
  final cycle on flush) through the bit-exact `EncoderFixed`. `reset`
  restores cold-start state for seeks. The direct `make_decoder` /
  `make_encoder` factories remain per the dual-API convention.
- **Official ITU-T G.728 conformance gate (r392)** — new
  `tests/conformance.rs` drives the complete official test-vector
  corpus (`docs/audio/g728/conformance/`, the ITU-T G.728 Appendix I
  sequences; auto-skips when the private docs tree is absent):
  - **Annex G fixed encoder**: `in1…in6` → `incw1g…incw6g` bit-exact
    (98 048 indices);
  - **Annex G fixed decoder, postfilter off**: `cw1…cw6` →
    `outa1g…outa6g` bit-exact (491 520 samples, block-floating `ST`
    rendered at `st·2^(3−NLSST)` with §G.1.3.1 rounding);
  - **Annex G fixed decoder, postfilter on**: `cw4` → `outb4g`
    bit-exact (51 200 samples, `SPF` Q2 doubled);
  - **float decoder, postfilter off**: `cw1…cw6` → `outa1…outa6`
    bit-exact (491 520 samples, `round(x·8)` output conversion);
  - **float encoder**: `in1…in4` → `incw1…incw4` bit-exact (14 336
    indices); case 6 pinned at ≤ 1 near-tie flip (index 29 of 256) and
    case 5 pinned to track ≥ 200 decisions before its near-tie
    free-run drift — the float references are only reproducible up to
    the reference codec's own floating-point tie behaviour.
- **Conformance-adjudication probes (r392)** —
  `EncoderFixed::encode_vector_resynced` (the search runs normally
  while the state update follows a caller-forced channel index,
  isolating per-vector decisions from divergence cascades when driving
  the reference bitstream through the encoder) and read-only accessors
  `EncoderFixed::{weight_coeff, predictor}`,
  `DecoderFixed::postfilter`.

- **Annex G full fixed-point encoder + decoder main programs (§G.7)
  (r389)** — new `annex_g_coder` module lands `EncoderFixed` /
  `DecoderFixed`, chaining every §G.3 fixed-point block in the §G.7
  main-program execution order and closing the fixed-point coder end to
  end. The decoder reproduces the encoder's quantized speech `ST`
  **bit-exactly** (mantissas + `NLSST`) from the 10-bit channel-index
  stream alone (200-vector lockstep test), the postfiltered output
  tracks the input at correlation > 0.9 / ±3 dB, and the decoder's
  block 50 uses the §G.2.2 `MINC0 = 10` order-10 interrupt to capture
  the postfilter by-product (`APF` / `RC1` / `NLSAPF`, with the §G.3.28
  `LABEL` Q13 conversion) before resuming to order 50.
- **Annex G fixed-point backward vector gain adapter (Figure G.1,
  blocks 43–48 + 91–99) (r389)** — new `annex_g_gain_adapter` module:
  the §G.3.14 block-43 flat-Q9 log-gain hybrid window
  (`LogGainWindowFixed`), the §G.3.16 per-vector blocks — 46
  (`log_gain_predict`), 98 (`limit_log_gain_q9`), 99+48
  (`inverse_log_gain`, the 10·2⁻⁶ + 20649·2⁻²¹ split of
  0.05·log₂10 and the Q15 Horner Taylor 2^x) and 93/94/96/97
  (`gstate1_update` via the §G.5 Q11 dB tables `GCBLG_Q11` /
  `SHAPELG_Q11`, Tables G.3/G.4, each entry cross-checked against the
  float dB derivations) — and the `GainAdapterFixed` driver with the
  §G.7 commit timing (43/44 at `ICOUNT = 1`, 45 at `ICOUNT = 2`).
  Fixed-vs-float σ trajectories agree to ≈ 0.008 dB steady-state on an
  encoder-style index stream.
- **Annex G §G.3.24–§G.3.27 fixed-point long-term-postfilter adaptation
  (blocks 81–84) (r389)** — new `annex_g_pitch` module
  (`PitchAdapterFixed`): block 81 (BFL `ST` → Q2 `SST` + the Q13 LPC
  inverse filter into the Q1 residual `D`), block 82 (Q19/Q13 1 kHz
  low-pass + 4:1 decimation, decimated + refined correlation
  peak-picking, and the cross-multiplied `CMAX·SUM >
  ((CORMAX·TMP) >> 16)·ITAPTH` sub-multiple decision), block 83
  (`PTAP` Q14 via `FINDNLS`/`RND`/`DIVIDE`), block 84 (`GL` Q14 /
  `GLB` Q16 via `DIVIDE`), plus `apf_to_q13` (§G.3.28 `LABEL`).
- **Annex G §G.3.1–§G.3.3 fixed-point encoder filter loop (r389)** —
  new `annex_g_encoder` module (`EncoderFiltersFixed`): block 4 (the
  Q2 perceptual weighting filter), blockzir (the segmented-BFL
  synthesis-filter zero-input response + the Q2 weighting-filter ZIR)
  and the blocks-9/10 memory update with the `LABEL1` `ET >> 1`
  overflow-retry probe, the ±4095-at-segment-scale clip and the
  reversed `ST`/`NLSST` emission.

### Fixed

- **Block 36's recursive decay is 1/2, not 3/4 (r389)** — the base
  Recommendation's block-36 pseudo-code spells `REXPW(I) =
  (1/2)·REXPW(I) + TMP` and §G.3.12 passes `NLSATTW = 15` into
  HWMCORE; both the float hybrid window (which hardcoded 0.75 for all
  blocks) and `WeightAdapterFixed` (which reused `NLSATT50 = 14`) were
  decaying the weighting window too slowly. `HybridWindow` gains a
  `decay` field, the fixed core a `nlsatt` parameter. The fixed
  block-36 front end was also rewritten to the §G.3.12 flat 15-bit Q2
  `SBW` form (with the `NLSREXPW ≤ 41` clamp); `WeightAdapterFixed::
  adapt` now takes the cycle's 20 Q2 `STMP` samples directly.
- **§G.3.18 HWMCORE case-1 recursive-loop sign typo (r389)** — the
  published pseudo-code's `AA1 = AA1 + (RREC(I+1) << 16)` line (a 5/4
  growth) contradicts its own "Scale RREC by 3/4" margin comment, the
  R(1) line above it and all four case-2/case-3 instances; the
  transcription now follows the five consistent `−AA1 + (RREC << 16)`
  instances. Block 49's tests rarely reached case 1 with the old
  `NLSREXP = 0` init; block 43 hits it every cycle.
- **Table G.2 initial values (r389)** — `NLSSB` starts at 16 (was 0)
  and `NLSREXP`/`NLSREXPW`/`NLSREXPLG` at 31 (was 0), per the table's
  explicit "initial value" entries and the §G.3.12/§G.3.14 prose.
- **Float gain adapter GP seed (r389)** — `GainAdapter::new` seeded
  `GP = (1, 0, …)`, predicting δ̂ = 0 dB (σ ≈ 39.81) on the first
  cycle; Table 2/G.728 lists `1, −1, 0, 0, …`, which predicts
  δ̂(n) = δ(n−1) and yields the correct unity first σ. The two tests
  that had codified the wrong first-vector σ were rewritten.

### Changed

- `LevinsonStatus` gains an `rc1` field (the Q15 first reflection
  coefficient, Table G.2's `RC1`), captured at the first-order step for
  block 85's `k1` input.

- **Annex G §G.3 fixed-point short-term postfilter coefficient calculator
  (block 85) (r374)** — new `annex_g_pf_coeff` module lands the fixed-point
  reformulation of the §4.6 adaptive postfilter's short-term coefficient
  calculator (eq. 4-3 / 4-4 / 4-5), the analogue of the floating-point
  `ShortTermPostfilter::set_from_synthesis_byproduct`.
  `short_term_coeff_fixed(atil, k1, nlsatmp)` derives the Q14
  `ShortTermCoeff` the §G.3.20–§G.3.23 postfilter (`PostfilterFixed`)
  consumes from the order-10 synthesis-filter LPC by-product: **AZ** (the
  all-zero numerator, `b̄_i = ã_i·SPFZCF^i`, `SPFZCF = 0.65`) and **AP**
  (the all-pole denominator, `ā_i = ã_i·SPFPCF^i`, `SPFPCF = 0.75`) are the
  two §G.3.15/§G.3.19 `bandwidth_expand_q14` expansions of the same `ATMP`
  predictor with the staged Q14 broadening tables `SPFZCFV_Q14` /
  `SPFPCFV_Q14` (postfilter analogue of weighting-filter block 38), and
  **TILTZ** (the spectral-tilt term `µ = TILTF·k1`, `TILTF = 0.15`) is the
  one-multiply Q14 round of `TILTF_Q14 · k1` over `2¹⁵` (new
  `TILTF_Q14 = 2458` const). Returns `None` (keep the previous cycle's
  postfilter coefficients) on a Q14 overflow / out-of-range `NLSATMP`. No
  new dB-domain Q-format is introduced — AZ/AP reuse the already-specified
  Table-G.2 Q14 coefficient format and `k1` is the §G.2.2 recursion's Q15
  `RC1`. Public surface: `short_term_coeff_fixed`, `TILTF_Q14`. 6 module
  tests: `TILTF_Q14` constant, unity heads + orders, the `SPFZCF < SPFPCF`
  zero-vs-pole radius ordering, the tilt-term half-LSB match vs float, an
  AZ/AP/TILTZ cross-check against the floating-point block-85 calculator on
  an identical order-10 predictor, and the bad-`NLSATMP` decline. Crate
  tests: 414 → 420.
- **Annex G §G.3 fixed-point perceptual-weighting-filter backward adapter
  (blocks 36/37/38) (r374)** — new `annex_g_weight_adapter` module lands the
  fixed-point reformulation of the encoder's perceptual weighting filter
  `W(z) = Q̃(z/γ₁) / Q̃(z/γ₂)` backward adaptation, the fixed-point analogue
  of the floating-point `weighting_filter_adapter` (blocks 36/37) +
  `weighting_filter_coeff` (block 38). `WeightAdapterFixed::adapt` runs one
  cycle: **block 36** (§G.3.17 windowing structure) drives the shared
  §G.3.17/§G.3.18 hybrid-window core (`annex_g_hybrid`) with the `LPCW = 10`
  / `NONRW = 30` / `WNRW` dimensions on past input speech, producing the 11
  autocorrelation taps `R(1..=LPCW+1)`; **block 37** is the §G.2.2
  `levinson_durbin_fixed` at order `LPCW`, yielding the order-10 predictor
  `QTMP` in its `Q13/Q14/Q15` block-floating format; **block 38** (eq. 3-4b
  / 3-4c) bandwidth-expands `QTMP` twice via the shared
  `bandwidth_expand_q14` core — the numerator `Q̃(z/γ₁)` with `γ₁ = WZCF =
  0.9` (`WZCFV_Q14`) and the denominator `Q̃(z/γ₂)` with `γ₂ = WPCF = 0.6`
  (`WPCFV_Q14`) — emitting the Q14 `WeightCoeffFixed` pair. On
  ill-conditioning or a Q14 overflow the previous cycle's pair is kept (the
  §3.4 "keep the previous coefficients" rule); the §3.4.1 non-speech `W(z) =
  1` disabling collapses both vectors to the all-pass `(16384, 0, …, 0)` via
  `WeightCoeffFixed::disabled` / `WeightAdapterFixed::disabled_coeffs`. The
  emitted Q14 numerator/denominator pair is exactly the coefficient set a
  fixed-point encoder weighting filter (blocks 4/10) would consume. Public
  surface: `WeightAdapterFixed`, `WeightCoeffFixed`, `N5W`. 8 module tests:
  disabled cold-start, default-vs-new, spec dimensions (`N3 = 60`,
  `N5W = 4`), zero-input ⇒ ILLCOND ⇒ disabled pair preserved, block-38
  unity heads + the `γ₂ < γ₁` pole-vs-zero radius ordering, a broadband
  end-to-end drive that commits a real `W(z)`, the disabled-helper
  cross-check, and a normalized-predictor shape cross-check against the
  floating-point `WeightingFilterAdapter` + `WeightingFilterCoeff` on
  identical requantized speech. Crate tests: 406 → 414.
- **Annex G §G.3.17/§G.3.18 shared fixed-point hybrid-window core (r374)** —
  new `annex_g_hybrid` module extracts the §G.3.17 windowing step (block 49)
  + §G.3.18 HWMCORE from the LPC-hard-coded `SynthAdapterFixed` into a
  parameterized, dimension-agnostic core `HybridWindowFixed` /
  `HybridWindowFixedState` — the fixed-point analogue of the floating-point
  `HybridWindowState` that already factors the three G.728 hybrid windows
  (blocks 36/43/49) into one parameter object. The core carries the
  permanent SBFL signal-history buffer `SB` with its per-segment `NLSSB` and
  the BFL recursive-autocorrelation `REXP` / `NLSREXP`, and runs the
  segment-shift / window-multiply / three-NLS-case recursive+non-recursive
  autocorrelation / white-noise-correction chain for any `(order, l, n,
  window)`, emitting `RTMP(1..=order+1)` plus the §G.2.2 `ILLCOND` verdict;
  the `3/4` recursive decay (`NLSATT50 = 14`) is shared across blocks.
  `SynthAdapterFixed` now delegates blocks 49 + HWMCORE to the shared core
  (via a `BflSegment` conversion of its `StSegment` input), keeping only its
  block-50 Levinson + block-51 expansion wiring; its public `HwmcoreOut` /
  `NLSATT50` are re-exported from the new module unchanged. The
  autocorrelation-vs-float-window cross-check moves into the core and now
  also covers the `LPCW` (block 36) window. New public surface:
  `HybridWindowFixed`, `HybridWindowFixedState`, `BflSegment`, plus the
  relocated `HwmcoreOut` / `NLSATT50`. 11 module tests. Crate tests:
  400 → 406.
- **Annex G §G.3.17–§G.3.19 fixed-point backward synthesis-filter adapter
  (r367)** — new `annex_g_synth_adapter` module lands the fixed-point
  reformulation of the decoder's backward synthesis-filter adaptation
  chain (blocks 49 / HWMCORE / 51), the fixed-point analogue of the
  floating-point `synthesis_adapter`. `SynthAdapterFixed::adapt` runs one
  adaptation cycle: **block 49** (§G.3.17) applies the Q15 hybrid window
  `WNR` to the segmented-block-floating-point (SBFL) signal buffer `SB` of
  past quantized speech — shifting in the four newest `ST` segments with
  their per-segment `NLSSB`, finding the common minimum shift `NLSTMP`,
  aligning each segment via the `NRSH = NLSSB(J) − NLSTMP − 1` per-segment
  shift (with the `NRSH = −1 ⇒ P << 1` Q15-multiplication compensation),
  and rounding the products into the BFL windowed signal `WS`. **HWMCORE**
  (§G.3.18) accumulates the recursive autocorrelation component (3/4 decay
  realised as `−RREC << NLSATT50(=14) + RREC << 16`) and the non-recursive
  component across the three `NLSRREC`-vs-`NLSAA0` alignment cases, applies
  the `>> 8` white-noise correction (`257/256`) to `R(1)`, `VSCALE`-
  normalizes against `MLS = 30`, and emits `RTMP(1..=LPC+1)` (BFL, 16-bit)
  plus the `ILLCOND` verdict (the full 32-bit `R(LPC+1)` accumulator being
  zero — the §G.2.2 interoperability fix that moves the zero-test upstream
  of the 16-bit-rounded Levinson input). **Block 50** is the already-
  landed §G.2.2 `levinson_durbin_fixed`. **Block 51** (§G.3.19) bandwidth-
  expands `ATMP` (Q13/Q14/Q15 per `NLSATMP`) by `FACV` into the live Q14
  predictor `A`, shifting to `Q30` (`<< 3 / 2 / 1` for `NLSATMP =
  13 / 14 / 15`) before `RND`, and keeps the previous cycle's `A` on
  `ILLCOND` or a Q14 overflow (the §G.3.19 `LABEL` "do not update" path).
  The produced Q14 `A` array is exactly the coefficient set the §G.3.11
  block-32 fixed-point synthesis filter (`SynthesisFilterFixed`) consumes,
  closing the fixed-point decoder's backward-adaptation loop. New
  `consts`: `N1` / `N2` / `N3` / `N4` / `N5` / `N6` (the §G.3.17 buffer
  dimensions 70 / 85 / 105 / 21 / 4 / 17), `NLSATT50 = 14`, `NLSSB_INIT`.
  Public surface: `SynthAdapterFixed`, `StSegment` (one 14-bit BFL `ST`
  segment + its `NLS`), `HwmcoreOut`. 11 module tests: fresh all-pass
  predictor, the §G.3.17 buffer dimensions, the SBFL shift-in of new
  segments at the `SB` tail, zero-input ⇒ `ILLCOND` ⇒ all-pass preserved,
  a broadband-input end-to-end drive that completes the §G.2.2 Levinson
  and commits a non-trivial Q14 predictor once the `N3 = 105`-sample
  hybrid-window history fills, block-51 `FACV·ATMP` scaling + `ILLCOND`
  decline, and an autocorrelation-shape cross-check against the §G.3.17
  *floating-point* hybrid window (`hybrid_window::HybridWindowState`) on
  the same requantized speech (normalized `R(k)/R(1)` tracked within the
  BFL-arithmetic tolerance). Floating-point build still drives the live
  codec; this is the §G.3 fixed-point variant. **Docs note:** the
  §G.3.18 Case-2/3 recursive-tap lines read `AA1 = RREC << NLSATT`, then
  `AA1 = −AA1 + RREC << 16` (negated), matching the `I = 0` line; only
  Case-1's recursive *taps* use the un-negated `AA1 = AA1 + RREC << 16`
  (the `I = 0` Case-1 line stays negated). Total crate tests: 385 → 396.
- **Annex G §G.3.15 fixed-point log-gain bandwidth expansion (block 45)
  + shared expansion core (r367)** — the §G.3.15 and §G.3.19 bandwidth-
  expansion modules (blocks 45 and 51) have identical fixed-point
  structure, so the block-51 body is factored into a shared
  `bandwidth_expand_q14(coeff, fac, nlsatmp)` function: scale each tap
  `2..=order+1` by the Q14 broadening table, shift the `Q27/Q28/Q29`
  product up to `Q30` per the Levinson `NLSATMP` (`<< 3 / 2 / 1` for
  `Q13/Q14/Q15`), `RND` the high word to `Q14`, and return `None` (keep
  the previous coefficients) on a `Q14` overflow or out-of-range
  `NLSATMP`. `log_gain_bandwidth_expand(gptmp, nlsgptmp, illcondg)` is the
  block-45 wrapper — the `LPCLG = 10`-order log-gain predictor expansion
  with the `FACGPV` table, gated by the log-gain ill-conditioning flag
  `ILLCONDG` — the gain-adapter analogue of block 51 (committed at
  `ICOUNT = 2` vs block 51's `ICOUNT = 3`). 4 new tests: the shared helper
  reproduces block 51 tap-for-tap, declines on bad `NLSATMP`, the
  `FACGPV·GPTMP` block-45 scaling, and the `ILLCONDG` decline. Total crate
  tests: 396 → 400.
- **Annex G §G.3.20–§G.3.23 fixed-point adaptive postfilter (r358)** —
  new `annex_g_postfilter` module lands the decoder's bit-exact
  fixed-point postfilter chain (blocks 71 – 77), the final §4.7/G.728
  stage applied to the reconstructed speech. `PostfilterFixed::
  filter_vector` runs, per `IDIM`-sample vector: blocks 71/72 (combined
  per §G.3.20) — the long-term pitch postfilter (`GL` Q14, `GLB = GL·B`
  Q16, lag `KP`, reading the Q2/Q0 decoded-speech buffer `SST`) cascaded
  with the short-term all-zero FIR (`AZ` Q14), all-pole IIR (`AP` Q14),
  and spectral-tilt (`TILTZ` Q14) postfilter (`STPFFIR`/`STPFIIR` memory
  Q2); blocks 73/74 — the sums of absolute values of decoded and
  filtered speech; block 75 (`scale_factor`) — their ratio `SCALE` in
  scalar-float form via the §G.1.3.4 `DIVIDE` (with the `SUMFIL ≤ 4 ⇒
  SCALE = 1` guard); blocks 76/77 — the first-order low-pass of `SCALE`
  (`SCALEFIL` Q14 init 16384, `AGCFAC = 16220` Q14 / `AGCFAC1 = 20972`
  Q21) and the output gain scaling to the postfiltered vector `SPF`
  (Q2). The §G.4 Table G.1/G.728 fixed-point constants are transcribed
  verbatim. An in-test transcription of the §G.3.20 *floating-point*
  short-term postfilter is the cross-check oracle (tracked within the
  block's stated precision across a multi-vector sweep). **Docs note:**
  the §G.3.20 fixed-point all-pole memory store renders literally as
  `AA0 = AA0 >> 14; STPFIIR(1) = AA0`, but `AA0` there holds
  `longterm << 2`, not the all-pole output; the §G.3.20 floating-point
  pseudo-code on the same page (`STPFIIR(1) = TEMP(K)`, the all-pole
  output before tilt) and the stated Q-formats force the source to be
  the Q16 accumulator `AA1` — implemented as `STPFIIR(1) = AA1 >> 14`.
- **Annex G §G.3.11 fixed-point decoder synthesis filter (r358)** — new
  `annex_g_decoder` module lands block 32, the 50th-order LPC decoder
  synthesis filter `1/A(z)`, as the first fixed-point **decoder** block
  (the prior Annex G work covered the §G.2 gain adapter / Levinson and
  the §G.3.4–§G.3.10 coder search). `SynthesisFilterFixed::synthesize`
  drives one gain-scaled excitation vector `ET` (15-bit block-floating,
  `NLSET`) through the filter to reconstruct the quantized-speech vector
  `ST` (14-bit block-floating, `NLSST`). The `STATELPC` 50-tap memory is
  held in segmented block-floating-point form with one shared left-shift
  per `IDIM`-sample segment (`NLSSTATE[1..=10]`, init 16, plus the
  scratch running-minimum slot 11). The four §G.3.11 phases are
  transcribed in lockstep: (A) the zero-input "ringing" response with the
  per-segment memory shift and `NLSSTATE`-aligned accumulation, (B) the
  `VSCALE`-to-15-bit re-normalization that refreshes `NLSSTATE`, (C) the
  zero-state response of `ET` with the `AA0 << 3` overflow probe and the
  "halve `ET`, drop `NLSET`, restart `LABEL1`" retry, and (D) the
  three-case `LABEL2` alignment that sums ZIR + ZSR, clips the synthesis
  memory to `±4095` at segment scale, `VSCALE`s to 14 bits, and reverses
  the top `IDIM` taps into `ST`. Built on the §G.1.3 `vscale` /
  `shl_sat` / `shr_sat` primitives in `annex_g_arith`; an in-test
  transcription of the §G.3.11 *floating-point* pseudo-code is the
  cross-check oracle (identity-filter pass-through and a first-order
  all-pole filter tracked within the annex's stated precision).
- **Annex G §G.2.2 fixed-point Levinson-Durbin recursion (r337)** — new
  `annex_g_levinson` module lands the §G.2.2 bit-exact reformulation of
  the three Levinson-Durbin modules (block 50 synthesis filter, block 37
  weighting filter, block 44 log-gain predictor — identical code, renamed
  variables per the §G.2.2 translation table). Built on the
  `annex_g_arith` `RND(.)` primitive. `simpdiv()` is the §G.2.2 `SIMPDIV`
  subroutine — the "different and simpler division algorithm" the
  recursion uses for its two divisions (the first-order reflection
  coefficient and each `RC = −SUM/ALPHATMP`): a plain 16-iteration
  restoring long division of two 16-bit integers, output in the low 17
  bits of an accumulator (distinct from the §G.1.3.4 `DIVIDE` used
  elsewhere). `levinson_durbin_fixed()` realises the full fresh-run
  pseudo-code: the autocorrelation `RTMP` arrives `Q15` block-floating-
  point (largest magnitude in `[0.5, 1)`); `ATMP` starts `Q15` and
  down-shifts to `Q14`/`Q13` on overflow (the `NRS` strategy — on an
  `AA0`/`AA1` overflow, `NRS += 1`, halve every `ATMP(2..=MINC+1)`, and
  recompute the overflowed term, restarting the half-step with both terms
  re-derived when the second product overflows); `ALPHATMP` is carried as
  the rounded 16-bit accumulator high word; the output `NLSATMP = 15 −
  NRS` signals the `Q` format to the downstream bandwidth-expansion
  module. The first floating-point `RTMP(LPC+1) = 0` test is taken as the
  upstream `ILLCOND` input (block 49 tests the full 32-bit accumulator,
  per the §G.2.2 interoperability fix), and failures set `ILLCOND` /
  `ILLCONDP` (postfilter coeffs also invalid when the stop order ≤ 10)
  plus the `NLSATMP < 13` post-check. `levinson_durbin_fixed_resume()`
  realises the decoder restart at `MINC0 = 10` (the recursion is
  interrupted after the order-10 postfilter coefficients, then resumed
  carrying `ATMP` / `NRS` / `ALPHATMP`). 10 module tests: SIMPDIV vs
  scaled integer division, both pre-recursion failure paths, AR(1)/AR(2)
  coefficient agreement with the floating-point `levinson` reference,
  `NLSATMP ∈ {13,14,15}`, unity-tap reconstruction, a 50th-order
  well-conditioned completion, and a decoder-resume run that reproduces
  the one-shot 1→50 run's orders 11..50 bit for bit. Floating-point build
  still drives the live codec; the §G.5/§G.6 Q-format packing stays
  deferred.

- **Annex G §G.2.1 reformulated backward vector gain adapter (r332)** —
  new `annex_g_gain` module lands the §G.2.1 mathematically-equivalent
  log-gain measurement that replaces the floating-point adapter's
  per-vector logarithm (blocks 39 / 40) with two precomputed dB
  table-lookups (Figure G.1 blocks 93 / 94). `gain_log_db()` is the
  block-93 gain-codebook log-gain table `20·log10|g_i|` (8 entries,
  sign-symmetric over the four distinct `|GQ[i]|` magnitudes);
  `shape_log_db()` is the block-94 shape-codebook log-gain table
  `10·log10·P[y_j]` (128 entries) where `P[y_j] = E_j / IDIM =
  Y_ENERGY[j]·DIMINV`; both are derived at runtime from the already-
  transcribed `GQ` / `Y_ENERGY` tables (the dB conversion needs
  `f64::log10`, so unlike `Y_ENERGY` they cannot be `const`).
  `offset_removed_log_gain()` realises adder 96 + limiter 97:
  `δ(n−1) = δ̂(n−1) + 20·log10|g_i| + 10·log10·P[y_j]` (eq. G-14,
  corrected), floored to `−32` dB (eq. G-9 / eq. G-7's `P[e(n)] ≥ 1`
  clip) — the value mathematically equal to the adder-42 output of
  Figure 6/G.728. The printed eq. G-14 shows the shape term with
  coefficient 1, but eq. G-8 / G-13 and the prose both carry the factor
  10; the module follows the factor-10 prose, and a per-test proves the
  factor-10 form reproduces the direct `10·log10·P[e(n)]−32` while the
  factor-1 form diverges. 12 module tests including the full
  reformulation-vs-direct-log equivalence sweep across σ / gain index /
  shape index. Floating-point only; the §G.5 / §G.6 Q-format packing of
  these dB tables stays deferred behind the floating-point build.
- **Annex G fixed-point arithmetic foundation §G.1.2 / §G.1.3 (r326)** —
  new `annex_g_arith` module lands the bottom layer the Annex G
  (1994-11) bit-exact fixed-point variant is built on. Numerical
  representations of §G.1.2 (single/double-precision fixed point, scalar
  floating point as a normalized `[16384, 32767]` mantissa + `NLS`,
  block floating point) plus the §G.1.3 arithmetic primitives: `shr_sat`
  / `shl_sat` (the §G.1.3.1 sign-extending right shift with the
  documented `3>>1=1` / `−3>>1=−2` magnitude anomaly, and the
  over-wide-shift-collapses-to-sign-fill rule), `rnd` (the §G.1.3.1
  `RND(.)` round-to-nearest-with-saturation on the 32-bit accumulator,
  matching the worked 1.5→2 / −1.5→−1 examples), `findnls` / `vscale`
  (the §G.1.3 `FINDNLS` / `VSCALE` block-floating-point normalization,
  with the three sign/magnitude cases and the zero-vector `NLS = MLS+1`
  rule, `MLS = 14` single / `30` double), `ScalarFloat::from_i16`
  (`FLOAT(.)`), and `divide` (the §G.1.3.4 `DIVIDE` 16-bit-precision
  scalar-float long division used inside Durbin's recursion, with the
  15-iteration restoring loop + 16th-bit rounding). 18 module tests
  including float-reference cross-checks for `DIVIDE` and the spec's
  worked `RND` / shift examples. Arithmetic-only foundation; the Annex G
  per-block coder is still deferred.
- **Annex I §I.4.1 frame-erasure excitation extrapolation (r322)** — new
  `excitation_extrapolation` module realises blocks 31SF / 31FE / 31E
  (§I.5.2–§I.5.4 floating-point pseudo-code). `FrameErasureExcitation`
  maintains the `ETPAST()` gain-scaled-excitation history (`KPMAX + IDIM`
  samples) and extrapolates one `ET()` vector per erased vector. Block
  31SF (`on_erased_cycle`) latches the voiced/unvoiced decision at
  erasure onset from the last good `PTAP` against the lowered threshold
  `VTH = PPFTH/1.4`, saving `FEDELAY = KP` (voiced) or the magnitude
  reference `AVMAG = (1/8)·Σ|last 40 ETPAST|` (unvoiced), and refreshes
  the scale `FESCALE` from the `VOICEDFEGAIN`/`UNVOICEDFEGAIN` schedules
  per 10 ms of erasure (forced to 0 past 50/60 ms). Block 31FE
  (`extrapolate`) produces the vector: voiced = `FESCALE · ETPAST(i −
  FEDELAY)` (periodic decaying repeat), unvoiced = a random 5-sample
  slide-back magnitude-matched via `RATIO = FESCALE / MAG`. Block 31E
  (`push`) slides `ETPAST()` by `IDIM` and appends, run every frame so
  the history is ready at erasure onset. The unvoiced "random" slide-back
  is a `SlideSource` trait (spec licenses any aperiodic `5..=140`
  sequence); a dependency-free `LcgSlideSource` is the default, tests
  drive a fixed sequence. Floating-point only — fixed-point Q2 `ETPAST`
  stays with the deferred Annex G build. 13 module tests + a `VTH`
  cross-check.
- **Annex I §I.4.5 gain-growth limitation after frame erasure (r316)** —
  new `gain_growth_limiter` module realises Block 47AF (the §I.5.6
  floating-point replacement for the normal Block 47 log-gain limiter)
  plus the §I.5.1 main-loop `AFTERFE` / `FECOUNT` / `OGAINDB`
  bookkeeping. `limit_log_gain(gain, ogaindb, afterfe) -> f64` is the
  pure transform: the normal `[0, 60]` dB range clamp followed by the
  Annex I addition `if AFTERFE > 0 and GAIN - OGAINDB > FEGAINMAX, GAIN =
  OGAINDB + FEGAINMAX` (caps backward-adapted log-gain growth to
  `FEGAINMAX = +2 dB`/vector for the first few good frames after an
  erasure, preventing the post-erasure gain "pop"). `GainGrowthLimiter`
  drives the cycle bookkeeping: a good cycle decrements an active
  `AFTERFE`; the first good cycle after an erasure loads `AFTERFE` with
  the erasure length `FECOUNT` saturated at `AFTERFEMAX = 16` cycles
  (= 40 ms); an erased cycle increments `FECOUNT`; `OGAINDB` starts at
  the `-32` dB log-gain floor. New `consts` `FEGAINMAX`, `AFTERFEMAX`,
  `OGAINDB_INIT`. Twelve per-tests pin the range-clamp-before-growth
  ordering, the `+2 dB`/vector cap, decrease-never-clamped, the
  `AFTERFE = FECOUNT` load + saturation + decrement cadence, `OGAINDB`
  threading through `limit`, and the Table I.1 literals (incl. the
  `16 × 2.5 ms = 40 ms` cross-check).
- **Annex I §I.4.2 frame-erasure LPC filter softening (r312)** — new
  `frame_erasure_lpc` module realises the first Annex I (frame/packet
  loss concealment) mechanism: during an erased frame the decoder reuses
  the last good 50th-order LPC predictor but "softens" it by an extra
  bandwidth expansion `a′_i = (0.97)^{k·i}·a_i` (factor `FEFAC = 0.97`
  rather than the normal `FAC = 253/256`), progressively re-softened by
  another factor of `0.97` every additional 10 ms of erasure (step `k`
  incrementing every `FE_LPC_CYCLES_PER_STEP = 4` adaptation cycles).
  `soften_predictor(&[f64; 51], step) -> [f64; 51]` is the pure block-51-
  style transform (leaving the implicit `a_0 = 1` tap untouched);
  `FrameErasureLpc` tracks the §I.4.2 step `k` across runs of erased /
  good adaptation cycles (step 1 on the first erased cycle, +1 per 10 ms,
  reset to 0 on the next good frame). New `consts` `FEFAC` and
  `FE_LPC_CYCLES_PER_STEP`. Nine per-tests pin the closed-form
  `(0.97)^{k·i}` scaling, the `a_0` invariant, the 4-cycle step cadence
  (`1,1,1,1,2,2,2,2,3,…`), good-cycle reset, and the `FEFAC < FAC`
  guard.
- **§5.11 serial byte-stream framing (r303)** — new `bitstream` module
  serialises the 16 kbit/s channel per the block-18 transmit epilogue
  of the §5.11 pseudo-code: each 10-bit `ICHAN` is packed
  most-significant-bit-first (`b9` first, … `b0` last), indices
  concatenated with no padding, packed into bytes MSB-first.
  `pack_indices(&[u16]) -> Vec<u8>` and `unpack_indices(&[u8]) ->
  Result<Vec<u16>>` are exact inverses; `unpack_indices` returns
  `Error::InvalidInputLength` when the buffer's bit length is not a
  whole multiple of `CHANNEL_INDEX_BITS = 10`. Whole-stream wrappers
  `Encoder::encode_buffer(&[[f64; 5]]) -> Result<Vec<u8>>` and
  `Decoder::decode_bytes(&[u8]) -> Result<Vec<f64>>` bridge the
  per-vector encode / decode path to the serial wire format. Eleven
  per-tests pin the MSB-first layout, byte-boundary straddles, masking,
  rejection of non-multiple-of-10 lengths, and an end-to-end
  `encode_buffer` → `decode_bytes` drive.

### Changed

- `Error::InvalidInputLength` now documents its spec-faithful §5.11
  meaning (bit length not a multiple of 10) instead of the earlier
  placeholder 16-bit-little-endian-word layout.

### Added (earlier)

- **§3.11 synchronization / in-band-signalling bit robbing (r297)** —
  the encoder can now rob the leftmost codeword bit of a chosen vector
  to carry a synchronization or in-band-signalling bit, the in-stream
  side-channel of §3.11. `CodebookSearch::search_with_sync_bit(target,
  gain, bit)` restricts the §3.9 shape scan to one half of the codebook
  — shape indices `0..=63` for a `0`, `64..=127` for a `1` — so the
  emitted index's most-significant shape bit (which, because the seven
  shape bits precede the three sign-and-gain bits, is the leftmost bit
  of the whole codeword) equals the desired bit. `Encoder::
  encode_vector_with_sync_bit(input, sync_bit)` runs the full
  per-vector analysis-by-synthesis loop with that half-codebook search;
  the new crate-root `extract_sync_bit(ichan)` helper recovers the bit
  on the decoder side (bit 9 of the channel index). Because the robbed
  vector runs the identical backward-adaptation dataflow (only the
  shape range narrows) and the half-codebook codevector is a perfectly
  valid 10-bit index, the encoder and decoder stay in bit-exact
  lockstep on a robbed stream — both ends must agree out of band on the
  every-`N`-th robbing schedule (§3.11 recommends `N` a multiple of 4,
  e.g. `N = 16` ≈ 100 bit/s, robbed from the last vector of a cycle).
  Four new per-tests: the half-codebook search confines the shape index
  to the requested half + sets the top codeword bit; it still finds the
  half-restricted brute-force eq. 3-16 minimum; the sync bit round-trips
  through encode + `extract_sync_bit` on an `N = 16` schedule; and the
  encoder ↔ decoder `sq(n)` lockstep is preserved bit for bit across a
  robbed stream.
- **ICOUNT = 3 synthesis-filter update stagger — spec-faithful on
  both encode and decode paths (Table E-1/G.728 + block 51, r287)** —
  `Decoder::decode_vector` and `Encoder::encode_vector` now carry a
  1-based `sf_icount` cycle counter (`1 → 2 → 3 → 4 → 1` per vector)
  and a `pending` / `active` synthesis-predictor pair. The backward
  synthesis-filter adapter (blocks 49/50/51) still runs at the cycle
  boundary on the previous cycle's quantised speech, but the
  bandwidth-expanded predictor is stashed and only swapped into the
  live filter when `ICOUNT` reaches 3 — the block-51 `Wait until
  ICOUNT = 3, then A = ATMP` instruction, and the "first use at
  decoding/encoding vector 3" row of Table E-1/G.728. Previously both
  ends applied the freshly-adapted predictor on the very next vector
  (vector 1), a two-vector phase error relative to the spec; that
  error is now removed, so the crate interoperates with a third-party
  G.728 endpoint that follows the normative cadence. The encoder
  additionally defers its block-38 weighting-filter coefficient commit
  and the block-12/14/15 impulse-response + filtered-shape-energy
  refresh to the same `ICOUNT = 3` boundary (the spec runs blocks
  12/14/15 at ICOUNT = 3 once `A`, `AWZ`, `AWP` are ready). The
  backward vector gain adapter already updated its log-gain predictor
  at `ICOUNT = 2` internally, matching Table E-1's "first use at
  vector 2", so no change was needed there. New accessors
  `Decoder::active_predictor` / `Decoder::sf_icount` and
  `Encoder::active_predictor` / `Encoder::sf_icount` expose the
  staggered state; five new tests pin the cadence (ICOUNT walk,
  swap-only-at-ICOUNT-3 on both paths, and encoder↔decoder active-
  predictor lockstep), while the pre-existing 200-vector bit-exact
  `sq(n)` lockstep test confirms the two ends still reconstruct
  identical quantised speech.
- **Codebook search + memory update — encoder loop CLOSED (§3.9
  blocks 12–18 + §3.10 / §5.13, r276)** — new `codebook_search`
  module transcribes the §5.11 pseudo-code: block 12 (five-sample
  impulse response of the cascade `F(z)·W(z)` from the zero-memory
  `{1,0,0,0,0}` excitation), blocks 14 + 15 (per-cycle convolution
  of all 128 shape codevectors with `h` and the filtered energy
  table `E_j = ‖H·y_j‖²`, eq. 3-20; Table 2 initials `H = 1,0,0,0,0`
  and `Y2 =` energy of `y_j`), block 16 (target normalisation
  `x̂ = x/σ(n)`), block 13 (time-reversed convolution
  `p(n) = Hᵀ·x̂(n)`, eq. 3-19), and blocks 17 + 18 (the §3.9.2
  division-free gain decision tree over the `GB` cell boundaries
  plus the eq. 3-23 distortion `D̂ = −b_i·P_j + c_i·E_j` minimum and
  the `ICHAN = (IS−1)·NG + (IG−1)` concatenation).
  `ZeroInputResponse::update_memory` adds the §5.13 filter-memory
  update: zero-state responses of the chosen `e(n)` added onto the
  post-ring memory, `STATELPC` clipped at the ±4095 MAX/MIN
  envelope, `ZIRWFIR` re-anchored to the top `STATELPC` taps, and
  `sq(n)` read out of the top five taps (block 22 omitted per
  §3.10). `Encoder::encode_vector` now runs the full Figure 2
  per-vector loop — block 20 → 4 → 9/10/11 → 12–18 → 19/21 → §5.13
  → end-of-cycle adapter bookkeeping (same cadence as
  `Decoder::decode_vector`) — and emits real 10-bit channel
  indices.
- **Public surface additions (r276)**: `CodebookSearch`,
  `SearchResult`, the `codebook_search` module,
  `ZeroInputResponse::update_memory`, `Encoder::codebook_search`,
  `Encoder::quantized_speech`.
- 21 new tests (r276), including: a brute-force cross-check of the
  §3.9.2 decision tree against the raw eq. 3-16 MSE over all 1024
  (gain, shape) pairs; exact-codevector recovery across both gain
  halves; block-16 σ(n) invariance; the COR = 0 tie route of the
  §5.11 pseudo-code; bit-for-bit equality of ring + `update_memory`
  against the decoder's block-32 `Synthesizer` (output and all 50
  memory taps, 32 vectors); a 200-vector encoder ↔ decoder
  lockstep test (decoder output == encoder `sq(n)` bit for bit);
  and a coding-error-energy ≪ signal-energy tracking property
  after adapter convergence.

- **Perceptual weighting filter adapter — upstream blocks 36 + 37
  (§3.3, r258)** — new `weighting_filter_adapter` module wires the
  hybrid window on input speech (block 36, Annex A.3 `wnrw`, 60-tap
  window with dimensions `LPCW + NFRSZ + NONRW = 10 + 20 + 30`)
  into the order-10 Levinson-Durbin recursion (block 37) and
  produces the predictor `q_i = a_i^{(10)}` (eq. 3-2f) that block 38
  (landed r248) consumes. The module reuses the existing
  `HybridWindow` / `HybridWindowState` machinery (already shared
  with blocks 43 and 49) and the existing `levinson_durbin`
  routine; on Levinson failure the cached predictor is kept and the
  error is propagated, mirroring `SynthesisAdapter::adapt`'s policy
  for block 33 of §3.7. The adapter is owned by the `Encoder` and
  exposed via three new methods: `adapt_weighting_filter` runs one
  cycle of `NFRSZ = 20` input-speech samples through blocks 36 +
  37, `weighting_filter_adapter` borrows the cached predictor for
  inspection, and `commit_weighting_filter_coefficients` performs
  the §3.3 "third vector of the cycle" coefficient swap by pushing
  the predictor through block 38 into the live block-4 filter
  (`set_weighting_filter_coeff_from_lpc`). Per §3.4 spec rule the
  per-sample memory of block 4 is preserved across the swap.
- **Public surface additions**: `WeightingFilterAdapter`, the
  `weighting_filter_adapter` module, `Encoder::adapt_weighting_filter`,
  `Encoder::commit_weighting_filter_coefficients`,
  `Encoder::weighting_filter_adapter`.
- 15 new tests (`weighting_filter_adapter` module: fresh adapter
  predictor is all-pass `q_i = (1, 0, ..., 0)`, `Default` matches
  `new`, zero input drives Levinson to `ZeroSignal` /
  `TrailingZero` and keeps the all-pass predictor, speech-like
  input produces `q_0 = 1` and at least one non-zero tap, accessor
  exposes `LPCW + 1` taps, predictor round-trips through
  `WeightingFilterCoeff::from_lpc` with both leading unity taps
  preserved, predictor drives `Encoder::set_weighting_filter_coeff_from_lpc`
  to a non-disabled state, `q_0 = 1` invariant survives mixed
  success/failure cycles, hybrid-window dimensions match
  `LPCW + NFRSZ + NONRW = 60`; `encoder` module: fresh encoder's
  upstream adapter is all-pass, `adapt_weighting_filter` surfaces
  Levinson errors verbatim on zero input, `adapt_weighting_filter`
  on speech-like input succeeds at least once over 6 cycles,
  `commit_weighting_filter_coefficients` propagates the predictor
  through block 38 and into block 4, block 4's per-sample memory
  survives the commit-driven swap, adapter accessor returns the
  same state as a standalone `WeightingFilterAdapter::new`).

- **Perceptual weighting filter — block 4 application path
  (§3.4, r249)** — new `weighting_filter` module transcribes the
  spec's one-paragraph §3.4: the current input speech vector `s(n)`
  is passed through `W(z) = Q̃(z/γ₁) / Q̃(z/γ₂)` (eq. 3-4a) to produce
  the weighted speech vector `v(n)`. Realised as a direct-form I
  pole-zero filter `v(n) = s(n) + Σ_{i=1..10} qγ₁_i · s(n-i)
  − Σ_{i=1..10} qγ₂_i · v(n-i)`, derived by cross-multiplying
  `V(z) · Q̃(z/γ₂) = S(z) · Q̃(z/γ₁)`. The two implicit `q_0·γ_k^0 = 1`
  leading taps of `Q̃(z)` show up as the standalone `+ s(n)` /
  `v(n)` terms in the difference equation; only the broadened taps
  `qγ_k_i` for `i = 1..10` interact with the delay lines. Per §3.4
  spec note ("filter memory should not be reset to zero at any time
  except during initialization"), the two 10-tap memories
  (`s(n-1..n-10)`, `v(n-1..n-10)`) are zeroed only by
  `PerceptualWeightingFilter::new`; coefficient swaps, the §3.4.1
  disable switch, and per-frame freeze-and-swap updates all leave
  the memories untouched. Block 10 (the §3.5 zero-input response of
  `F(z)·W(z)`) is intentionally NOT wired up — it requires the
  §3.10 pre-/post-save memory dance with the synthesis filter and
  is queued for a later round.
- **`Encoder::apply_weighting_filter` + `Encoder::weighting_filter`
  accessor.** The encoder now carries a live block-4 filter
  initialised at construction to the §3.4 / §3.4.1 all-pass
  `W(z) = 1` state. `set_weighting_filter_coeff_from_lpc` and
  `disable_weighting_filter` now also mirror the new coefficient set
  into the live filter — the per-sample memory is preserved across
  the swap. `apply_weighting_filter` consumes one
  `FRAME_LEN`-sample input vector and emits the corresponding
  weighted vector `v(n)`.
- **Public surface additions**: `PerceptualWeightingFilter`, the
  `weighting_filter` module, `Encoder::apply_weighting_filter`, and
  `Encoder::weighting_filter`.
- 14 new tests (`weighting_filter` module: fresh filter is the
  all-pass `W(z) = 1`, fresh memory is zeroed per §3.4
  initialisation, `Default` equals `new`, memory advances one
  sample per emitted output and carries across vector boundaries,
  coefficient swap leaves memory untouched, first-sample-after-init
  always equals `s[0]` regardless of `W(z)` shape, δ-input matches
  the hand-traced first five outputs of the difference equation
  term-by-term, `γ₁ = γ₂` collapse drives the filter to the identity
  on 64 vectors of sinusoidal input, 1024-vector finite-output
  stability; `encoder` module: fresh encoder filter is the
  pass-through, setter mirrors block-38 output into the live block-4
  taps, memory survives the coefficient swap, the §3.4.1 disable
  switch zeroes the coefficients but leaves the memory). Total
  crate tests: 193 → 207.

- **Perceptual weighting filter coefficient calculator (block 38,
  Figure 4a/G.728, §3.3 / §3.4, r248)** — new
  `weighting_filter_coeff` module transcribes block 38 of the encoder's
  perceptual weighting filter adapter. §3.3 defines the filter as
  `W(z) = Q̃(z/γ₁) / Q̃(z/γ₂)` (eq. 3-4a) where
  `Q̃(z) = 1 + Σ q_i · z⁻ⁱ` (eq. 3-3a) is the 10th-order prediction-
  error filter built from the order-10 LPC predictor `q_i` output of
  block 37, and `γ₁ = WZCF = 0.9`, `γ₂ = WPCF = 0.6` are the
  Table 1/G.728 zero-control / pole-control bandwidth factors. Block
  38's transform is the mechanical substitution `z ← z/γ_k` of
  eq. 3-4b / 3-4c, producing the bandwidth-broadened sequences
  `q_i · γ₁ⁱ` (numerator) and `q_i · γ₂ⁱ` (denominator) for
  `i = 1..=LPCW`. Both sequences carry the spec's 1-based layout
  (`q_gamma_k[0] = 1.0` implicit leading tap; `q_gamma_k[1..=LPCW]`
  the broadened taps), matching the `SynthesisAdapter::coefficients`
  convention already in use elsewhere in the crate.
- **`WeightingFilterCoeff::disabled` — §3.4.1 non-speech-mode
  constructor** — produces the `W(z) = 1` all-pass state for modem /
  other non-speech signals (the spec's "set γ₁ = γ₂ = 0" shortcut,
  which drops every non-unity tap of both sequences). Also used as
  the start-up value before the first adaptation cycle has produced
  `q_i`.
- **`Encoder::weighting_filter_coeff` accessor +
  `Encoder::set_weighting_filter_coeff_from_lpc` +
  `Encoder::disable_weighting_filter`** — block 38 driven through
  the encoder's typed surface. The fresh encoder starts at the
  §3.4.1 `W(z) = 1` state; downstream encoder rounds wire block 37's
  Levinson-Durbin output into the calculator via the typed setter.
  The disable hook flips the §3.4.1 switch without re-running the
  transform.
- **Public surface additions**: `WeightingFilterCoeff`, the
  `weighting_filter_coeff` module, and the three new `Encoder`
  methods above.
- 10 new tests (`weighting_filter_coeff` module: disabled state is
  `W(z) = 1` shape, `Default` equals `disabled`, unity-q input
  reproduces the `γ_k^i` geometric progression to within `1e-14`,
  zero-γ inputs disable the filter, hand-traced `q_i = (−0.5)^i`
  matches the eq. 3-4b / 3-4c termwise expected values to within
  `1e-14`, `γ₂ < γ₁ < 1` geometric-consequence invariant at the
  longest tap, non-unity `q[0]` correctly panics; `encoder` module:
  fresh encoder matches `WeightingFilterCoeff::disabled`,
  `set_weighting_filter_coeff_from_lpc` is the spec-direct
  dataflow, `disable_weighting_filter` returns to the all-pass).
  Total crate tests: 183 → 193.

- **Typed encoder scaffold + `E_j` shape-energy table (§3.9, r235)** —
  new `encoder` module exposes a stable [`Encoder`] type and a
  [`make_encoder`] factory mirroring the decoder's dual-API
  convention. The encoder carries the two backward adapters
  (`SynthesisAdapter` and `GainAdapter`, §3.7 / §3.8) plus the
  previous gain-scaled excitation `ET(1..IDIM)` slot, all initialised
  to Table 2/G.728 values. §4.4 / §4.5 of the Recommendation require
  the two backward adapters to be bit-for-bit identical at both ends
  of the channel, so the encoder reuses the exact types the decoder
  already drives.
- **`tables::Y_ENERGY` — precomputed shape codevector energy
  `E_j = Σ_k y_j(k)²`** — new `[f64; NCWD]` constant derived at
  compile time from the already-transcribed Annex B `Y_Q11` shape
  codebook via `Σ_k Y_Q11[j][k]² / 2²²` (the `Y_Q11 / 2¹¹` Q-shift
  squared). §3.9 of the Recommendation rearranges the
  analysis-by-synthesis distortion (eqs. 3-14..3-23) into
  `D_{i,j} = b_i · <x̃, ỹ_j> + c_i · E_j` (with `b_i = 2·g_i`,
  `c_i = g_i²`), and `E_j` is a constant of the codebook that the
  search loop will look up once per shape scan. Exposed via the new
  `Encoder::shape_energy()` accessor and the runtime `y_energy_f64()`
  view function (same convention as `facv_f64` / `facgpv_f64`).
- **`Encoder::encode_vector` typed signature** returning
  `Result<u16, Error>` — `Error::NotImplemented` in round 235 because
  the §3.9 analysis-by-synthesis search (blocks 1..28 + 67..70) is
  not yet implemented. Signature is final for the future pipeline:
  one `FRAME_LEN`-sample input vector in, one 10-bit channel index
  out.
- **Public surface additions**: `Encoder`, `make_encoder`,
  `CHANNEL_INDEX_BITS = 10` (§2.1), `tables::Y_ENERGY`,
  `tables::y_energy_f64`.
- 13 new tests (`tables` module: `Y_ENERGY` dimension matches NCWD,
  every entry equals the direct `y_f64()` dot product to machine
  precision, every entry is finite + strictly positive, row 0 hits
  the hand-computed value 20 443 149 / 2²², runtime accessor mirrors
  the const; `encoder` module: fresh state matches Table 2 `ET = 0`,
  `Default` == `new`, shape-energy accessor returns the right
  table, `encode_vector` returns `NotImplemented` in r235,
  `CHANNEL_INDEX_BITS == 10`, `make_encoder` produces fresh state,
  the synthesis adapter exposes `A(1) = 1` at fresh construction,
  the gain adapter starts at `ICOUNT = 1`). Total crate tests: 170 →
  183.

- **Long-term postfilter coefficient calculator (blocks 83 + 84,
  Figure 7/G.728, §4.7, r229)** — new `PitchPostfilterCoeff` type
  transcribes the §4.7 single-tap pitch-predictor weight and the
  long-term-postfilter coefficient calculator. 245-sample
  `sd(−239..5)` decoded-speech buffer (kept separate from block 71's
  `KPMAX`-sample comb-filter delay line so block 71's sample-rate
  cost is unchanged). At the third vector of each frame, given the
  pitch period `p` from block 82, eq. 4-12 evaluates
  `β = Σ_{k=−99..0} sd(k)·sd(k−p) / Σ_{k=−99..0} sd(k−p)²` clamped
  into `[0, 1]`. Eq. 4-13 then maps `β` to `b`:
  `0` if `β < PPFTH = 0.6` (postfilter off — `H_l(z) = 1`);
  `PPFZCF · β = 0.15 · β` if `0.6 ≤ β ≤ 1`;
  `PPFZCF = 0.15` if `β > 1` (unreachable after the `[0, 1]` clamp).
  Eq. 4-14 derives `g_l = 1 / (1 + b)`. The `(g_l, b)` range is
  `[1/(1+PPFZCF), 1] × [0, PPFZCF]`. New public API:
  `PitchPostfilterCoeff::{new, last_beta, last_b, last_g_l,
  speech_buffer, push_decoded_vector, compute_coefficients, reset}`;
  public constant `SPEECH_BUF_LEN = 245`.
- **`Decoder::decode_vector_postfiltered` drives blocks 83 + 84.**
  Every decoded-speech vector is pushed into the coefficient
  calculator's `sd(−239..5)` buffer; at the third vector of every
  adaptation cycle (post-increment `icount == 3`), the pitch
  extractor's output `p` is fed to
  `PitchPostfilterCoeff::compute_coefficients`, and the resulting
  `(g_l, b, p)` triple is applied to
  `LongTermPostfilter::set_coefficients`. From the second adaptation
  cycle onward the long-term comb filter operates at the
  spec-prescribed coefficients (rather than the previous §4.6.1
  passthrough fallback).
- **`Decoder::pitch_pf_coeff()` accessor** for tests and audit.
- 22 new tests (`pitch_postfilter_coeff` module: cold-start zeroed
  buffer + `(g_l, b, β) = (1, 0, 0)` invariants, spec 245-sample
  buffer dimension, reset round-trip, push lands the new vector at
  `sd(1..5)`, successive pushes slide buffer left by `IDIM`, pushing
  enough vectors drops the oldest off the left edge, cold-start
  compute returns passthrough, near-zero signal routes to
  postfilter-off, perfectly periodic impulse-train signal gives
  `β = 1` at the analysis lag with `b = PPFZCF` and
  `g_l = 1/(1+PPFZCF)`, β clamped above unity, sub-threshold β
  disables the postfilter, voiced-range β gives `b = PPFZCF·β`,
  `g_l = 1/(1+b)` consistency across all eq. 4-13 branches,
  `b ∈ [0, PPFZCF]` and `g_l ∈ [1/(1+PPFZCF), 1]` range invariant on
  random inputs, finite-output stability over a 512-vector sinusoidal
  drive; `Decoder`-level: accessor exposes cold-start state, buffer
  advances with decoded speech, coefficient calculator updates at
  the third vector of each frame, propagation of `(g_l, b)` from the
  coefficient calculator into the long-term comb at every extract,
  block-83/84 wiring does not perturb the register-only
  `decode_vector` path, cold-start `sf == sd` identity over the
  first frame is preserved, long-term coefficients satisfy the
  spec's `(g_l, b, p)` invariants after many frames, eq. 4-13 / 4-14
  branch invariant after a long decode). Total crate tests: 148 →
  170.

- **Pitch period extractor (block 82, Figure 7/G.728, §4.7, r223)** —
  new `PitchSearch` type transcribes the §4.7 pitch-period dataflow.
  240-sample LPC-residual buffer covering `d(−139..100)` with the
  §4.7-prescribed `d(81..85) / d(86..90) / d(91..95) / d(96..100)`
  vector slot ordering (4th vector of the previous frame stored at
  `d(81..85)`, then current frame's three vectors at the following
  blocks). 60-sample decimated `d̄(−34..25)` buffer fed by the Annex
  D third-order elliptic 1 kHz lowpass + 4:1 decimator (using `AL` /
  `BL` coefficients already in `tables.rs`). Coarse correlation
  search `ρ(i) = Σ d̄(n)·d̄(n−i)` for `n ∈ [1, 25]` over decimated
  lags `i ∈ [5, 35]` (eq. 4-7) selects τ; full-resolution refinement
  `C(i) = Σ d(k)·d(k−i)` for `k ∈ [1, 100]` over `[4τ−3, 4τ+3]`
  clamped into `[KPMIN, KPMAX]` (eq. 4-8) selects `p0`. When `p̂` is
  outside `[p0 − KPDELTA, p0 + KPDELTA]` a second `C(i)` search
  over `[p̂−KPDELTA, p̂+KPDELTA]` selects `p1`; the single-tap
  predictor weights `β0 = C(p0)/Σ d(k−p0)²` (eq. 4-9) and
  `β1 = C(p1)/Σ d(k−p1)²` (eq. 4-10) clamped into `[0, 1]` drive
  the eq. 4-11 selection `p = p0` if `β1 ≤ TAPTH·β0` else `p1`
  (`TAPTH = 0.4`). The chosen `p ∈ [KPMIN, KPMAX] = [20, 140]` is
  stashed as `p̂` for the next frame. New public API:
  `PitchSearch::{new, previous_pitch, last_pitch, vector_cursor,
  residual_buffer, decimated_buffer, push_residual, extract,
  reset}`; public constants `RESIDUAL_BUF_LEN = 240`,
  `DECIMATED_BUF_LEN = 60`, `FRAME_RESIDUAL = 20`,
  `FRAME_DECIMATED = 5`.
- **`Decoder::decode_vector_postfiltered` runs block 82 each call.**
  After the inverse filter (block 81) produces the residual vector,
  it is pushed into `PitchSearch`'s buffer; at the third vector of
  every adaptation cycle (spec ICOUNT = 3 / our post-increment
  `icount == 3`) the extractor runs the full lowpass + decimate +
  correlation + refinement + resolution pipeline. The extracted
  pitch is not yet driving the long-term postfilter — that stays at
  the §4.6.1 passthrough until blocks 83 (single-tap `β` over the
  decoded-speech buffer, eq. 4-12) and 84 (the `(g_l, b)`
  calculator of eq. 4-13 / 4-14) land in a later round. Cold-start
  `sf = sd` invariant is preserved bit-for-bit (regression test
  `pitch_search_wiring_does_not_perturb_decode_vector_path`).
- **`Decoder::pitch_search()` accessor** for tests and audit, and
  as the read surface for blocks 83 / 84 when they land.
- 23 new tests (`pitch_search` module: cold-start zeroed-buffer +
  `p̂ = KPMIN` invariants, spec buffer dimensions, reset round-trip,
  vector-cursor walks `0..3` then wraps, push-residual lands the
  four per-frame vectors at the §4.7 prose's `d(81..85) / d(86..90) /
  d(91..95) / d(96..100)` slots, extract on all-zero buffer returns
  KPMIN with the `[KPMIN, KPMAX]` clamp surviving, extract sets
  the cursor to 3 for the upcoming 4th-vector push, coarse search
  locks to the spec-correct pitch at periods KPMIN / KPMAX / 40 on
  a unit-impulse train, `p_prev` carry-over after extract,
  fundamental-vs-multiple path runs without panic when `p̂` falls
  outside the neighbourhood, post-extract residual buffer slide by
  `NFRSZ = 20`, decimator advances by `FRAME_DECIMATED = 5` each
  extract, finite-pitch output on random inputs, lowpass-filter
  memory survives across extracts; `Decoder`-level: accessor
  exposes cold-start state, vector cursor advances with each
  postfiltered call, extract runs at the third vector of each
  frame, output is finite over 64 vectors, wiring preserves cold-
  start `sf == sd`, wiring leaves long-term postfilter at
  `(g_l, b, p) = (1, 0, KPMIN)`). Total crate tests: 125 → 148.

- **Pitch-extractor LPC inverse filter (block 81, Figure 7/G.728, §4.7,
  r220)** — new `PitchInverseFilter` type transcribes the spec's
  10th-order LPC inverse filter `Ã(z) = 1 − Σ ã_i · z^{-i}` (eq. 4-6)
  applied to the decoded-speech stream `sd(k)`. Sign convention
  `ã_i = -a_i` is applied inside `set_from_synthesis_byproduct`; the
  10-tap FIR memory carries across vector boundaries. The filter is
  driven by the same order-10 by-product already published by
  `SynthesisAdapter::order10_predictor()` and refreshes at the same
  first-vector-of-cycle boundary as the short-term postfilter. New
  public API: `PitchInverseFilter::{new, set_from_synthesis_byproduct,
  coefficients, memory, reset, filter_vector}`.
- **`Decoder::decode_vector_postfiltered` runs block 81 each call.**
  After the synthesis filter produces `sd(n)`, the inverse filter
  advances on it to keep the residual state synced with the decoded-
  speech stream. The residual `d(k)` is computed but not yet consumed
  — it is the input that block 82 (pitch search) will read once that
  module lands. The long-term (block 71) comb filter therefore stays
  at the §4.6.1 passthrough; the cold-start `sf = sd` invariant is
  preserved bit-for-bit (regression test
  `pitch_inverse_filter_wiring_preserves_cold_start_passthrough_equality`).
- **`Decoder::pitch_inverse_filter()` accessor** for tests and audit
  and as the read surface for the upcoming block-82 search.
- 12 new tests (`pitch_inverse_filter` module: fresh-filter pass-
  through, cold-start zero-coefficient invariants, all-pass predictor
  drives identity across many vectors, sign-flip `ã_i = −a_i` math,
  eq. 4-6 impulse response matches `−ã_i`, AR(1) roundtrip recovers
  unit excitation exactly, 10-tap memory carries across vector
  boundaries, reset returns to identity, sinusoidal-drive finite-
  output stability; `Decoder`-level: accessor exposes cold-start
  state, memory advances per postfiltered call, coefficients refresh
  each cycle, wiring preserves cold-start `sf = sd`, wiring leaves
  long-term postfilter at `(g_l, b, p) = (1, 0, KPMIN)`). Total
  crate tests: 111 → 125.

- **Long-term (pitch) postfilter comb-filter machinery (block 71,
  Figure 7/G.728, r213)** — new `LongTermPostfilter` type transcribes
  the spec's `H_l(z) = g_l · (1 + b · z^{-p})` cascade (eq. 4-1) as
  a `KPMAX`-sample circular FIR delay line. Pitch period `p` is
  clamped into `[KPMIN, KPMAX] = [20, 140]` per Table 1/G.728.
  Cold-start `(g_l, b, p) = (1, 0, KPMIN)` realises the §4.6.1
  "postfilter off" identity for unvoiced / weakly-voiced frames
  (decode-trace §7.1, equations 4-13/4-14). New public API:
  `LongTermPostfilter::{new, g_l, b, p, set_coefficients, reset,
  filter_vector}`.
- **`Decoder::decode_vector_postfiltered` chain expanded to block 71.**
  The long-term comb filter now slots between the synthesis output
  `sd(n)` and the short-term postfilter (block 72) per Figure 7/G.728.
  Per the same figure, the AGC numerator (block 73) reads the raw
  decoded `sd` directly — only the `sf` denominator (block 74) sees
  the long-term + short-term cascade. While block 71 stays in the
  §4.6.1 passthrough (`g_l = 1, b = 0`) this is observationally
  equivalent to the r207 chain; the cold-start `sf = sd` invariant
  is preserved bit-for-bit (regression test
  `long_term_passthrough_preserves_short_term_postfilter_behaviour`).
- **`Decoder::long_term_postfilter()` accessor** for tests and audit.
- 15 new tests (`long_term_postfilter` module: fresh-filter identity,
  passthrough cold-start coefficients, 64-vector identity hold,
  pitch clamp into `[KPMIN, KPMAX]`, gain+tap storage round-trip,
  eq. 4-1 comb verification on a unit impulse at `p = KPMIN`, pure-
  gain branch when `b = 0`, `b = 0` independence from `p`, echo
  exactly at lag `KPMAX`, FIR no-recursion check (exactly two non-
  zero samples), reset restores passthrough, sinusoidal-drive
  finite-output stability, circular delay-line wrap-around past
  `KPMAX`; `Decoder`-level: passthrough state still reported after
  16 vectors, cold-start `sf = sd` invariant preserved). Total
  crate tests: 96 → 111.

- **Short-term (spectral) postfilter (block 72, Figure 7/G.728, r207)**
  — new `ShortTermPostfilter` type transcribes the spec's
  `H_s(z) = (1 − Σ b̄_i z⁻ⁱ) / (1 − Σ ā_i z⁻ⁱ) · (1 + µ·z⁻¹)` cascade
  (eq. 4-2..4-5). Coefficient bandwidth expansion: `b̄_i = ã_i·SPFZCF^i`
  (`SPFZCF = 0.65`), `ā_i = ã_i·SPFPCF^i` (`SPFPCF = 0.75`),
  `µ = TILTF · k1` (`TILTF = 0.15`). The spec's `ã_i = -a_i` sign
  flip is applied inside `set_from_synthesis_byproduct`; the
  pole-zero + tilt filter cascade carries per-tap memory across
  vector boundaries.
- **Order-10 by-product extraction from the synthesis adapter** —
  `SynthesisAdapter::adapt` runs a second Levinson on the first 11
  entries of RTMP to produce `ã_1..ã_10`, and computes
  `k1 = -R(2)/R(1)` directly. New accessors:
  `SynthesisAdapter::order10_predictor()` and
  `SynthesisAdapter::k1()`. On order-10 Levinson failure (R(11) = 0
  case) the previous cache is retained, matching the spec's "skip
  block 51" semantics at the postfilter scale. New public constant
  `synthesis_adapter::PF_LPC_ORDER = 10`.
- **`Decoder::decode_vector_postfiltered(ichan)` is now live.** Runs
  the full block 29 → 33 chain, refreshes the short-term postfilter
  coefficients at the first vector of each adaptation cycle (spec
  ICOUNT = 1), filters `sd → sf` through `H_s(z)`, then runs the
  AGC against the real `sd, sf` pair. The §4.6.1 "postfilter off"
  path is the cold-start fallback (all coefficients zero → identity)
  until the first cycle commits a non-trivial Levinson result.
- **`Decoder::short_term_postfilter()` accessor** for tests and audit.
- 15 new tests (`synthesis_adapter`: order-10 by-product all-pass
  cold start, non-trivial-input order-10 with A(1) = 1 invariant,
  `|k1| < 1` reflection-coefficient bound; `short_term_postfilter`:
  fresh-filter identity, zero-coefficient cold-start invariants,
  bandwidth-expansion `ã_i · λ^i` match for both pole and zero
  arrays, `µ = TILTF · k1` invariant across signs, non-zero
  coefficients change output, memory state advances per sample,
  finite-output stability on synthetic speech, reset restores
  identity, DC steady-state stays finite, zero-coefficient
  passthrough for every input; `Decoder`-level: cold-start cycle
  matches `decode_vector` exactly, post-cycle output diverges from
  raw, SCALEFIL pinned at 1.0 during cold-start cycle, coefficients
  update each cycle, finite-output stability over 64 vectors). Total
  crate tests: 81 → 96.

### Fixed

- **`Decoder::decode_vector` dropped the gain-codebook level (r276)**
  — the block-29 lookup scaled the shape codevector by the backward-
  adapted σ(n) only, omitting the `GQ(IG)` factor that §5.14
  prescribes (`YN(K) = GQ(IG)·Y(NN + K)`); all eight gain levels in
  the same shape row decoded identically. The lookup now routes
  through `ExcitationVector::from_channel_index` (which implements
  the spec's block 29 exactly) and a new test pins the `GQ` ratio
  and the sign-mirrored negative half on cold-start vectors.

### Changed

- `Decoder::decode_vector_postfiltered` semantics: previously a strict
  §4.6.1 passthrough (sf = sd, GAINSF = 1, SCALEFIL pinned at 1.0).
  Now: the cold-start cycle (first NUPDATE vectors before the
  synthesis adapter commits) still matches `decode_vector` bit-for-bit,
  but once a non-trivial Levinson result lands the short-term
  postfilter coefficients refresh and `sf` diverges from `sd`. The
  long-term (block 71) stage is still on §4.6.1.

- **Output gain control / postfilter AGC tail (blocks 73 / 74 / 75 /
  76 / 77, Figure 7/G.728, r201)** — new `Agc` type transcribes the
  per-vector Σ|sd| / Σ|sf| ratio (blocks 73 / 74 / 75), the first-
  order lowpass `H(z) = 0.01 / (1 − 0.99·z⁻¹)` with `AGCFAC = 0.99`
  (block 76), and the per-sample multiplication into the post-
  filtered output (block 77). Includes the §4.6.1 "postfilter off"
  passthrough path that pins `SCALEFIL` at `1.0` when `sf = sd`.
- **`Decoder::decode_vector_postfiltered(ichan)`** — new entry that
  runs the full block 29 → 33 chain followed by the AGC pass. Until
  blocks 71 / 72 land it follows the §4.6.1 passthrough rule and is
  provably equal to `decode_vector` (regression test
  `decode_vector_postfiltered_matches_decode_vector_in_pf_off_mode`).
- **`Decoder::agc()` accessor** for tests and audit.
- 13 new tests (`agc` module: cold-start invariants, lowpass DC =
  unity convergence, AGCFAC geometric decay, sf = sd passthrough,
  zero-sf safe fallback, finite-output stability, reset, IIR
  trajectory cross-check, drift-free steady-state; `Decoder`-level:
  AGC pass equals raw decode in §4.6.1 mode, SCALEFIL stays at 1.0,
  output finite). Total crate tests: 68 → 81.

- **Backward synthesis-filter adapter (block 33)** — new
  `SynthesisAdapter` type wires the spec's three §5.6 sub-blocks
  (block 49 hybrid window → block 50 Levinson-Durbin → block 51
  bandwidth expansion via `FACV`) so the decoder can derive its
  50th-order LPC predictor autonomously from quantised speech, one
  NFRSZ-sample adaptation cycle at a time.
- **Backward vector gain adapter (block 30)** — new `GainAdapter`
  type wires the per-vector chain of §5.7 blocks
  67/39/40/42/46/47/48 plus the once-per-cycle (ICOUNT=2)
  43/44/45 hybrid-window → Levinson → `FACGPV` expansion path.
  Produces one σ(n) per call; preserves the previous-cycle
  log-gain predictor on Levinson failure.
- **`hybrid_window` module** — single transcription of the
  block-36 / block-43 / block-49 pseudocode shape, parameterised
  by order / per-call input length / non-recursive tail length /
  window data. Both adapters share the same code path; their state
  (`SB` / `SBLG` buffers + `REXP` / `REXPLG` recursive
  autocorrelation) lives in per-adapter `HybridWindowState`
  instances.
- **`Decoder::decode_vector(ichan)`** — autonomous decode path
  that walks block 29 → 30 → 31 → 32 → 33 per Figure 3/G.728 on
  every call, with the gain adapter producing σ(n) and the
  synthesis adapter producing A() internally. The earlier
  `Decoder::decode_index` / `set_synthesis_predictor` hooks are
  preserved for the register-only path.
- New tests (r195): 12 added (hybrid-window dimensions / WNCF / zero
  input, synthesis-adapter A(1)=1 invariant / zero-input
  ill-conditioning / nonzero-input convergence / decoder
  integration, gain-adapter Table 2 initial state / ICOUNT cycle /
  σ(n) limiter / first-vector 10^(GOFF/20), full-chain
  decode_vector finiteness / state propagation / cycle update).
  Total crate tests: 56 → 68.

## [0.0.7](https://github.com/OxideAV/oxideav-g728/releases/tag/v0.0.7) - 2026-05-29

### Other

- r189 decoder front end + Annex A/B/C/D tables
- orphan rebuild: reset to scaffold under clean-room policy (round 171)

### Added

- **Annex A.1 / A.2 / A.3 hybrid windows** — synthesis-filter (105
  samples), log-gain predictor (34 samples) and perceptual-weighting
  filter (60 samples) windows, transcribed from the Q15 integer
  columns of Recommendation G.728 (1992-09) Annex A. Float views are
  derived by the prose-stated `value / 2¹⁵` rule; per-table peak-value
  tests cross-check the integer against the float column the prose
  prints alongside it.
- **Annex B excitation codebook** — the 128 × 5-sample Q11 shape
  codebook (`Y_Q11`) and the 8-level Q13 gain codebook (`GQ` / `GB`)
  with precomputed `G2` and `GSQ` arrays from equations 3-21 and 3-22.
  Row 0, row 127 and the sign-symmetry of the gain codebook back half
  are pinned by spot-tests.
- **Annex C bandwidth-broadening tables** — `FACV` (51 entries),
  `FACGPV` / `WPCFV` / `WZCFV` / `SPFPCFV` / `SPFZCFV` (11 entries
  each), all Q14. Per-table tests cross-check against the prose-stated
  geometric progression `λⁱ` / `γⁱ` to one Q14 quantum.
- **Annex D 1 kHz lowpass filter coefficients** — the 4-numerator +
  3-denominator floating-point taps used by the pitch-extraction
  prefilter.
- **Table 1/G.728 parameters** — every named codec constant (IDIM,
  LPC, LPCLG, LPCW, NCWD, NG, NFRSZ, NPWSZ, NUPDATE, GOFF, FAC, FACGP,
  WNCF, perceptual / pitch / postfilter controls) lives in
  `src/consts.rs`.
- **Levinson-Durbin recursion** — line-for-line transcription of the
  block 50 pseudocode, parameterised by order (LPC = 50 / LPCLG = 10 /
  LPCW = 10). The three ill-conditioning exits surface as a typed
  `LevinsonError`. An analytical AR(1) round trip and a 2-pole AR(2)
  round trip pin the recursion's sign and indexing conventions.
- **Decoder front end** — `Decoder::decode_index(ichan)` splits a
  10-bit channel index into the 7-bit shape and 3-bit gain components
  per §3.11, looks up the gain-scaled codevector (block 29), and runs
  the 50th-order all-pole synthesis filter (block 32) including the
  spec's saturation-on-write to filter memory (§5.13).
- **Direct API** — `make_decoder()` factory entry mirrors the
  workspace's dual-API convention.

### Changed

- Crate state moves from orphan-rebuild scaffold to clean-room rebuild
  in progress; `Error::NotImplemented` is still the response from the
  encoder-side entry points (encoder will land in a later round), but
  the decoder side is now functional behind `Decoder` /
  `make_decoder`.
