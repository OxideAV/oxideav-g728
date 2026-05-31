# oxideav-g728

Pure-Rust ITU-T G.728 Low-Delay CELP (LD-CELP, 16 kbit/s) speech codec.

## Status: clean-room rebuild — autonomous decoder pipeline (round 195)

The crate was reset to a register-only scaffold under the workspace
clean-room policy (round 171, master `14e3bad`): the previous
implementation had numeric tables extracted from an external reference
C distribution, which the policy forbids regardless of the source's
licence. The crate is being rebuilt one spec-cited unit at a time from
the published ITU-T G.728 (1992-09) Recommendation prose alone.

Round 195 lands the two **backward adapters**, which together let the
decoder run autonomously off the raw 10-bit channel-index stream
without a caller-supplied predictor:

- **`SynthesisAdapter` — block 33** wires the spec's three sub-blocks
  end to end:
  - block 49 (hybrid window on the quantised-speech buffer, producing
    the 51-tap autocorrelation `RTMP(1..LPC+1)`),
  - block 50 (the Levinson-Durbin recursion landed in r189), and
  - block 51 (bandwidth expansion via the Annex C `FACV` Q14 vector).
  One call per NFRSZ-sample adaptation cycle returns the cycle's
  bandwidth-expanded predictor in the spec's `A(1..LPC+1)` layout.
  Levinson-Durbin failures (zero-signal, trailing-zero, or
  ill-conditioned autocorrelation) keep the previous cycle's
  predictor per the spec's "skip block 51" rule, and surface a typed
  `LevinsonError` for trace.
- **`GainAdapter` — block 30** wires the per-vector dataflow of
  blocks 67 / 39 / 40 (1-vector delay + RMS + 10·log10 with
  ETRMS<1 clip), 42 (`GSTATE(1) = ETRMS − GOFF`), 43 / 44 / 45 (the
  ICOUNT=2 log-gain hybrid window → Levinson → `FACGPV` expansion),
  46 (the log-gain linear predictor), 47 (limiter clamp `[0, 60]` dB)
  and 48 (`σ(n) = 10^(GAIN/20)`). One call produces one σ(n) for the
  vector that is about to be decoded.
- **`hybrid_window` module** — the spec's block-36 / block-43 /
  block-49 pseudocode shares an identical shape (only the LPC order,
  per-call input length, non-recursive tail length, and window data
  differ). The transcription lives in one place and is dispatched by
  a parameter object; each adapter brings its own `HybridWindowState`
  carrying the `SB`/`SBLG` signal buffer and the `REXP`/`REXPLG`
  recursive autocorrelation memory.
- **`Decoder::decode_vector(ichan)`** — autonomous decode path that
  walks block 29 → 30 → 31 → 32 → 33 per Figure 3/G.728 on every
  call. Caller passes one 10-bit channel index; the decoder returns
  one `FRAME_LEN`-sample PCM vector with all backward adaptation
  handled internally. The earlier `Decoder::decode_index` /
  `set_synthesis_predictor` hooks are preserved for the register-
  only path.

Round 189 (preserved) provides the Annex A.1/A.2/A.3 hybrid windows,
the full Annex B excitation codebook (128 × 5 shape + 8-level gain),
the Annex C bandwidth-broadening vectors (Q14), the Annex D 1 kHz
lowpass coefficients, Table 1/G.728 parameter constants, the
order-parameterised Levinson-Durbin recursion, and the
`Decoder::decode_index` caller-driven entry point.

## What is NOT yet wired up

- The adaptive postfilter (blocks 71–77). Long-term + short-term +
  AGC stages are perceptual; they do not affect bit-exact agreement
  against the standard's reference decoded vectors. The Annex C
  postfilter bandwidth vectors and the Annex D 1 kHz prefilter
  coefficients are already in `tables/` for that round.
- A-law / µ-law PCM I/O — `oxideav-g711` handles this per §5.3 / §3.1.
- The encoder side (blocks 1..28 + 67..70 + the codebook search of
  §3.9 / blocks 12..18). The decoder-side synthesis-filter adapter
  (block 23 in the encoder = block 33 in the decoder) and the
  decoder-side gain adapter (block 20 in the encoder = block 30 in
  the decoder) are now shareable via the `SynthesisAdapter` /
  `GainAdapter` types.
- Annex G fixed-point variant and Annex I frame-loss concealment
  remain deferred behind the floating-point decoder.

## Clean-room provenance

Every numeric value lives in `src/tables.rs` and is transcribed
directly from the integer columns of Annexes A, B, C, D in the
ITU-T G.728 1992-09 PDF that lives under `docs/audio/g728/`. Every
control-flow line in the new `hybrid_window`, `synthesis_adapter`
and `gain_adapter` modules carries a comment pointing at the
spec's §5.6 / §5.7 pseudocode for blocks 36/43/49 (hybrid window),
50 (Levinson — already in `levinson.rs`), 51 / 45 (bandwidth
expansion), and 67/39/40/42/46/47/48 (the per-vector gain chain).
No external implementation source has been opened or consulted
during this rebuild — the per-test cross-checks (peak Q15 values,
AR(1) predictor round-trip, `λⁱ` geometric progression, codebook
row spot-checks, ICOUNT cycle ordering, limiter clamp bounds,
first-vector σ(n) = 10^(GOFF/20)) act as in-repo audit anchors
against transcription typos.

## License

MIT — see [LICENSE](LICENSE).
