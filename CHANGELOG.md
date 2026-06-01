# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

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
