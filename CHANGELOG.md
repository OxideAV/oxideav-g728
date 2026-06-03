# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

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
- New tests: 12 added (hybrid-window dimensions / WNCF / zero
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
