# oxideav-g728

Pure-Rust ITU-T G.728 Low-Delay CELP (LD-CELP, 16 kbit/s) speech codec.

## Status: clean-room rebuild — autonomous decoder + long-term + short-term postfilter + AGC + block-81 pitch inverse filter (round 220)

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

Round 220 lands the **first stage of the §4.7 pitch-extraction
pipeline** — block 81, the 10th-order LPC inverse filter:

- **`PitchInverseFilter` — block 81.** New module transcribes the
  10th-order LPC inverse filter `Ã(z) = 1 − Σ ã_i · z^{-i}` (eq. 4-6)
  applied to the decoded-speech stream `sd(k)`. The output residual
  `d(k) = sd(k) − Σ ã_i · sd(k-i)` is the input that block 82's pitch
  search will consume. Spec sign convention `ã_i = -a_i` is applied
  inside `set_from_synthesis_byproduct`; the 10-tap FIR memory carries
  across vector boundaries.
- **Coefficient source.** The inverse filter shares the order-10
  by-product surface that already feeds the short-term postfilter
  (block 72): per `SynthesisAdapter::order10_predictor()`, refreshed
  at the first vector of each adaptation cycle (§7.1 — same boundary
  as §7.2).
- **`Decoder::decode_vector_postfiltered` advances block 81 each call.**
  After `sd(n)` is produced by the synthesis filter, the inverse
  filter runs on it to keep the residual state synced with the
  decoded-speech stream. The residual is computed but not yet
  consumed downstream — block 82 (1 kHz lowpass + 4:1 decimate +
  lag 5..35 coarse search + refinement) is the next chunk. The
  long-term (block 71) postfilter therefore stays at the §4.6.1
  passthrough; the cold-start `sf = sd` invariant is preserved
  bit-for-bit.
- **`Decoder::pitch_inverse_filter()` accessor** for tests and audit.

Round 213 (preserved) adds the **long-term (pitch) postfilter**
comb-filter machinery (block 71 of Figure 7/G.728):

- **`LongTermPostfilter` — block 71.** New module transcribes the
  spec's `H_l(z) = g_l · (1 + b · z^{-p})` cascade (eq. 4-1) as a
  `KPMAX`-sample FIR delay line, with `(g_l, b, p)` set externally
  via `LongTermPostfilter::set_coefficients`. Pitch period `p` is
  clamped into `[KPMIN, KPMAX] = [20, 140]` per Table 1/G.728. The
  filter is purely all-zero (no recursion); per-sample memory is
  carried across vector boundaries via the circular delay line.
- **§4.6.1 passthrough cold start.** Cold-start coefficients
  `(g_l, b, p) = (1, 0, KPMIN)` make the comb filter the identity —
  exactly the spec's "postfilter off" rule for unvoiced / weakly-
  voiced frames (decode-trace §7.1, equations 4-13/4-14). Until the
  §4.7 block-81..84 pitch-extraction / coefficient pipeline lands,
  the [`Decoder`] keeps the filter at this state.
- **`Decoder::decode_vector_postfiltered` chain expanded.** Block 71
  now slots in between the synthesis-filter output `sd(n)` and the
  short-term postfilter (block 72), per Figure 7/G.728. Per the
  same figure, the AGC numerator (block 73) reads the *raw* decoded
  `sd` directly — the long-term and short-term filters only feed
  into block 74's `sf` denominator. While block 71 is in passthrough
  this is observationally equivalent to the r207 chain, and the
  cold-start `sf = sd` invariant is preserved bit-for-bit.
- **`Decoder::long_term_postfilter()` accessor** for tests and audit.

Round 207 adds the **short-term (spectral) postfilter** (block 72 of
Figure 7/G.728):

- **`ShortTermPostfilter` — block 72.** New module transcribes the
  spec's `H_s(z) = (1 − Σ b̄_i z⁻ⁱ) / (1 − Σ ā_i z⁻ⁱ) · (1 + µ·z⁻¹)`
  cascade (eq. 4-2..4-5). The bandwidth-expanded coefficients are
  derived from the synthesis-filter Levinson recursion's **order-10
  by-products** `ã_1..ã_10` (eq. 4-3..4-4 with `SPFZCF = 0.65`,
  `SPFPCF = 0.75`) and the **first reflection coefficient `k1`** for
  the tilt-compensation factor `µ = TILTF · k1 = 0.15·k1` (eq. 4-5).
- **Order-10 by-product extraction** — `SynthesisAdapter::adapt`
  now also runs an order-10 Levinson on the same RTMP autocorrelation
  used for the order-50 synthesis predictor and caches `ã_1..ã_10`
  + `k1` for the postfilter to consume. Exposed via
  `SynthesisAdapter::order10_predictor()` and `::k1()`.
- **`Decoder::decode_vector_postfiltered(ichan)` is now live.** It
  runs the full decode chain, refreshes the short-term postfilter
  coefficients at the first vector of each adaptation cycle, filters
  `sd → sf` through `H_s(z)`, then runs the AGC against the real
  `sd, sf` pair. The §4.6.1 "postfilter off" path is still the
  cold-start fallback (all coefficients zero → identity; AGC stays at
  SCALEFIL = 1) until the first cycle commits a non-trivial Levinson
  result.
- The long-term (block 71) pitch postfilter still follows §4.6.1
  (`b = 0`, `g_l = 1`) because the block-81..85 pitch-extraction
  front end is not yet implemented. When it lands, the long-term
  filter slots in between `decode_vector` and `short_term_pf.filter_vector`
  with no changes to the AGC stage.

Round 201 (preserved) adds the **AGC tail of the postfilter** (blocks
73 / 74 / 75 / 76 / 77 of Figure 7/G.728):

- **`Agc` — blocks 73–77.** New standalone type transcribes the
  per-vector Σ|sd| / Σ|sf| ratio (blocks 73 / 74 / 75), the
  first-order lowpass `H(z) = 0.01 / (1 − 0.99·z⁻¹)` with the spec's
  `AGCFAC = 0.99` (block 76), and the per-sample multiplication into
  the post-filtered output (block 77). The lowpass has unity DC gain
  by construction, so a constant input ratio settles to the same
  ratio at steady state — the property the AGC relies on.
- **`Decoder::decode_vector_postfiltered(ichan)`** — new entry point
  that runs the full block 29 → 33 chain and then the block 73–77
  AGC pass. Until the long-term (block 71) and short-term (block 72)
  postfilter coefficient calculators land, this entry follows the
  §4.6.1 "postfilter disabled, raw synthesis output" path and feeds
  `sf = sd` to the AGC; the AGC's IIR pins SCALEFIL at `1.0` and the
  pass is provably a no-op (the regression test
  `decode_vector_postfiltered_matches_decode_vector_in_pf_off_mode`
  asserts exact equality with `decode_vector` across 32 vectors).
- The Decoder field layout (`agc: Agc` always present) and the
  `Agc::apply(&sd, &sf)` signature are already correct for blocks
  71 / 72 to slot in: a future round drops the §4.6.1 passthrough,
  feeds the real `sf` stream into `Agc::apply`, and the AGC stage
  needs no further changes.

Round 189 (preserved) provides the Annex A.1/A.2/A.3 hybrid windows,
the full Annex B excitation codebook (128 × 5 shape + 8-level gain),
the Annex C bandwidth-broadening vectors (Q14), the Annex D 1 kHz
lowpass coefficients, Table 1/G.728 parameter constants, the
order-parameterised Levinson-Durbin recursion, and the
`Decoder::decode_index` caller-driven entry point.

## What is NOT yet wired up

- §4.7 blocks 82 / 83 / 84 of the pitch-extraction pipeline driving
  block 71's `(g_l, b, p)`. Block 81 (10th-order LPC inverse filter
  producing the residual `d(k)`) is wired up (r220). Remaining are
  block 82 (Annex D 1 kHz lowpass + 4:1 decimate, coarse correlation
  search over lags 5..35, full-resolution refinement 4τ−3..4τ+3,
  fundamental-vs-multiple resolution against the previous frame's
  pitch), 83 (single-tap predictor weight β) and 84
  (`b = PPFZCF·β` / `g_l = 1/(1+b)` clamped per eq. 4-13/4-14). The
  comb-filter machinery itself (block 71) is wired up; until the
  block-82..84 sub-pipeline drives non-trivial coefficients the comb
  follows the §4.6.1 "long-term postfilter off" rule — `b = 0`,
  `g_l = 1`, the filter is the identity. The short-term postfilter
  (block 72) and AGC tail (blocks 73–77) are wired up.
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
control-flow line in the `hybrid_window`, `synthesis_adapter`,
`gain_adapter`, `long_term_postfilter`, `short_term_postfilter`,
`pitch_inverse_filter` and `agc` modules carries a comment pointing
at the spec's §5.6 / §5.7 / §4.6 / §4.7 pseudocode for blocks
36/43/49 (hybrid window), 50 (Levinson — already in `levinson.rs`),
51 / 45 (bandwidth expansion), 67/39/40/42/46/47/48 (the per-vector
gain chain), 71 (long-term comb filter `g_l · (1 + b·z^{-p})` from
eq. 4-1, KPMIN/KPMAX clamp), 72 (short-term postfilter pole-zero +
tilt cascade with `SPFZCF / SPFPCF / TILTF` and the `ã_i = -a_i`
sign-convention flip spelled out in `short_term_postfilter.rs`),
73/74/75/76/77 (the AGC tail with the lowpass form spelled out in
§4.6 immediately after eq. 4-5), and 81 (10th-order LPC inverse
filter `Ã(z) = 1 − Σ ã_i · z^{-i}` from eq. 4-6, same `ã_i = -a_i`
sign flip spelled out in `pitch_inverse_filter.rs`). No external
implementation has been opened or consulted during this rebuild —
the per-test cross-checks (peak Q15 values, AR(1) predictor round-
trip, `λⁱ` geometric progression, codebook row spot-checks, ICOUNT
cycle ordering, limiter clamp bounds, first-vector σ(n) = 10^(GOFF/20),
AGC unity-DC convergence + AGCFAC geometric decay, postfilter
identity at cold start + non-identity divergence after the first
adaptation cycle + bandwidth-expansion `ã_i · λ^i` invariants +
tilt `µ = TILTF · k1` invariant + block-81 impulse response matching
`−ã_i` + AR(1) synthesis/inverse round-trip to zero excitation
residual) act as in-repo audit anchors against transcription typos.

## License

MIT — see [LICENSE](LICENSE).
