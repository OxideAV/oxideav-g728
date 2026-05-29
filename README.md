# oxideav-g728

Pure-Rust ITU-T G.728 Low-Delay CELP (LD-CELP, 16 kbit/s) speech codec.

## Status: clean-room rebuild — decoder front end landed (round 189)

The crate was reset to a register-only scaffold under the workspace
clean-room policy (round 171, master `14e3bad`): the previous
implementation had numeric tables extracted from an external reference
C distribution, which the policy forbids regardless of the source's
licence. The crate is being rebuilt one spec-cited unit at a time from
the published ITU-T G.728 (1992-09) Recommendation prose alone.

Round 189 lands:

- **Table 1/G.728 parameters** — the complete set of LD-CELP coder
  constants (IDIM, LPC, LPCLG, LPCW, NCWD, NG, NFRSZ, NPWSZ, NUPDATE,
  GOFF, FAC, FACGP, WNCF, perceptual-weighting / postfilter / pitch
  controls), with sanity tests cross-checking the values that the
  prose lists twice.
- **Annex A.1 / A.2 / A.3 hybrid windows** for the synthesis filter
  (105 samples), log-gain predictor (34 samples) and perceptual
  weighting filter (60 samples). Stored as their normative Q15 i16
  integers; the prose-stated `value / 2¹⁵` rule produces the matching
  float view. Per-table peak-value tests cross-check the integer
  column against the float column the Annex prints alongside it.
- **Annex B excitation codebook** — the full 128-codevector × 5-sample
  shape codebook (Q11) and the 8-level gain codebook with the prose's
  closed-form `GQ(1) = 33/64`, `GQ(i) = (7/4)·GQ(i-1)` recurrence. The
  derived `G2 = 2·GQ` and `GSQ = GQ²` arrays from equations 3-21 /
  3-22 are precomputed at compile time.
- **Annex C bandwidth-broadening vectors** — `FACV` (51 entries, Q14),
  `FACGPV` / `WPCFV` / `WZCFV` / `SPFPCFV` / `SPFZCFV` (11 entries
  each, Q14), with per-table geometric-progression tests against the
  prose-stated `λⁱ` / `γⁱ` formulas.
- **Annex D 1 kHz lowpass coefficients** — the four numerator + three
  denominator taps of the pitch-extraction prefilter.
- **Levinson-Durbin recursion** — line-for-line transcription of
  block 50's pseudocode, dispatchable by order (LPC = 50, LPCLG = 10,
  LPCW = 10), with the three ill-conditioning exits surfaced as a
  typed error. AR(1) / AR(2) round-trip tests pin the sign and index
  conventions; a 50th-order smoke test exercises the spec's
  block-50 / block-37 / block-44 hot path.
- **Decoder front end** — `Decoder::decode_index(ichan)` splits a
  10-bit channel index into the 7-bit shape index and 3-bit gain
  index per §3.11, looks up the codebook entry (block 29), and runs
  the 50th-order all-pole synthesis filter (block 32) with the
  saturation prescribed by §5.13. `Synthesizer::set_predictor` is the
  hook for the future backward-adapter wiring (block 33).

## What is NOT yet wired up

- The backward adapters — block 30 (vector gain) and block 33 (LPC
  predictor). The synthesiser exposes a stable predictor-injection
  hook so the next round can plug the hybrid-window → Levinson-Durbin
  → bandwidth-expansion chain in without churning the API.
- The adaptive postfilter (blocks 71–77). The postfilter is a
  perceptual quality stage and does not affect bit-exact agreement
  against the standard's reference decoded vectors.
- A-law / µ-law PCM I/O — `oxideav-g711` handles this per §5.3 / §3.1.
- The encoder side, Annex G fixed-point variant, and Annex I frame-
  loss concealment all come later.

## Clean-room provenance

Every numeric value lives in `src/tables.rs` and is transcribed
directly from the integer columns of Annexes A, B, C, D in the
ITU-T G.728 1992-09 PDF that lives under `docs/audio/g728/`. No
external implementation source has been opened or consulted during
this rebuild — the per-test cross-checks (peak Q15 values, AR(1)
predictor round-trip, `λⁱ` geometric progression, codebook row
spot-checks) act as in-repo audit anchors against transcription
typos.

## License

MIT — see [LICENSE](LICENSE).
