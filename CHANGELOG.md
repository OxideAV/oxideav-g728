# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

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
