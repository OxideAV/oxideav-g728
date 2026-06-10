# oxideav-g728

Pure-Rust ITU-T G.728 Low-Delay CELP (LD-CELP, 16 kbit/s) speech codec.

## Status: clean-room rebuild — autonomous decoder + full §4.7 pitch postfilter (blocks 81/82/83/84) + long-term comb + short-term postfilter + AGC + typed encoder scaffold + §3.9 `E_j` shape-energy table + §3.3 full perceptual-weighting adapter (blocks 36 + 37 + 38) + §3.4 weighting filter applied to input speech (block 4) + §3.5/§3.6 zero-input response + VQ target (blocks 9 + 10 + 11) (round 267)

The crate was reset to a register-only scaffold under the workspace
clean-room policy (round 171, master `14e3bad`): the previous
implementation had numeric tables extracted from an external reference
C distribution, which the policy forbids regardless of the source's
licence. The crate is being rebuilt one spec-cited unit at a time from
the published ITU-T G.728 (1992-09) Recommendation prose alone.

Round 267 lands the **zero-input response + VQ target vector** —
blocks 9, 10 and 11 of Figure 2/G.728 (§3.5 / §3.6; pseudo-code
§5.9 / §5.10) — the analysis-by-synthesis **search target** `x(n)`
that the §3.9 codebook search of equations 3-14..3-23 minimises
against:

- **`ZeroInputResponse` — blocks 9 + 10 + 11.** New
  `zero_input_response` module transcribes the §5.9 "ring for five
  samples" recurrence: with switch 5 open (zero input at node 7),
  the synthesis filter (block 9, `STATELPC(1..LPC)` memory) and the
  perceptual weighting filter (block 10, `ZIRWFIR(1..LPCW)` /
  `ZIRWIIR(1..LPCW)` memories) are let "ring" for one `IDIM = 5`-
  sample vector, producing the zero-input response `r(n)` from their
  carried-over memory of past gain-scaled excitation. Block 11 then
  subtracts: `x(n) = v(n) − r(n)` (§5.10, `TARGET(K) = SW(K) −
  ZIR(K)`). Block 9 reuses the exact synthesis-filter memory-shift
  idiom the decoder's `Synthesizer` already runs (the §5.9 block-9
  recurrence is the same zero-input recurrence the decoder computes
  on a zero excitation); a per-test cross-checks the two paths to
  machine precision. The §3.5 initialisation state (all three memory
  arrays zero) makes the first vector's `r(n) = 0` and `x(n) =
  v(n)` exactly, matching the spec note "except for the vector right
  after initialization, the memory of the filters 9 and 10 is in
  general non-zero".
- **`Encoder::compute_vq_target(v)` + `Encoder::zero_input_response_unit()`
  accessor.** The encoder now owns the ring unit. `compute_vq_target`
  drives blocks 9 + 10 + 11 with the encoder's own order-50 synthesis
  predictor (block 23, shared with the decoder per §4.5) and the live
  block-38 weighting coefficients (`AWZ = q_gamma1`, `AWP =
  q_gamma2`), consuming the weighted speech vector `v(n)` from
  `apply_weighting_filter` (block 4) and emitting the target `x(n)`.
- **The §3.10 memory-update phase is deliberately NOT in this round.**
  Saving the post-ring memory, zeroing it, refiltering the chosen
  excitation `e(n)` to add the zero-state response on top, and reading
  `sq(n)` back out of the top five `STATELPC` taps all depend on the
  §3.9 codebook-search output `e(n)` that is not yet wired up. This
  module is the half of §3.5 that runs **before** the search and only
  reads the current filter memory; the complementary memory-update
  half runs **after** the search and lands in a later round.

Round 258 lands the **upstream half of the perceptual weighting
filter adapter** (blocks 36 + 37 of §3.3) — the missing link
between the input-speech buffer and the block-38 coefficient
calculator already in place:

- **`WeightingFilterAdapter` — blocks 36 + 37.** New
  `weighting_filter_adapter` module wires the §3.3 hybrid window
  on past input speech (Annex A.3 `wnrw`, 60-tap window with
  dimensions `LPCW + NFRSZ + NONRW = 10 + 20 + 30`) into the
  order-10 Levinson-Durbin recursion. The result is the predictor
  `q_i = a_i^{(10)}` (eq. 3-2f) in canonical Levinson layout
  `[1.0, a_1, …, a_10]`. The module reuses the shared
  `HybridWindow` / `HybridWindowState` machinery already powering
  blocks 43 and 49 and the existing `levinson_durbin` routine. On
  Levinson failure the cached predictor is kept and the error is
  propagated — same "keep previous cycle" policy
  `SynthesisAdapter::adapt` honours for block 33.
- **`Encoder::adapt_weighting_filter(speech)` +
  `Encoder::commit_weighting_filter_coefficients()` +
  `Encoder::weighting_filter_adapter()` accessor.** The encoder
  now owns the upstream adapter. `adapt_weighting_filter` runs one
  cycle of `NFRSZ = 20` input-speech samples through blocks 36 +
  37; `commit_weighting_filter_coefficients` performs the §3.3
  "third vector of the cycle" coefficient swap by pushing the
  predictor through block 38 into the live block-4 filter. Per
  §3.4 spec rule the per-sample memory of block 4 is preserved
  across the swap. The cycle-timing gate (commit only at
  `ICOUNT = 3`) is still the caller's responsibility — the same
  way the synthesis-filter adapter's commit timing is the
  caller's responsibility for block 33.
- **The encoder's analysis-by-synthesis search loop is still NOT
  wired up.** Blocks 1..28 + 67..70 of §3.9 (the per-vector VQ
  cost expression, impulse-response convolution, gain-quantiser
  decision tree) remain to be lifted off the precomputed `E_j` /
  `G2` / `GSQ` / `GB` table set already in `tables.rs`.

Round 249 (preserved) lands the **block-4 application path** of
the perceptual weighting filter (§3.4) — given the current input
speech vector `s(n)`, push it through `W(z) = Q̃(z/γ₁)/Q̃(z/γ₂)`
(eq. 3-4a) to produce the weighted speech vector `v(n)`:

- **`PerceptualWeightingFilter` — block 4.** New `weighting_filter`
  module transcribes the spec's one-paragraph §3.4 application. The
  filter is realised as a direct-form I pole-zero stage,
  `v(n) = s(n) + Σ_{i=1..10} qγ₁_i · s(n-i) − Σ_{i=1..10} qγ₂_i · v(n-i)`,
  cross-derived from `V(z) · Q̃(z/γ₂) = S(z) · Q̃(z/γ₁)`. The two
  implicit `q_0·γ_k^0 = 1` leading taps of `Q̃(z)` (eq. 3-3a) appear
  explicitly as the standalone `+ s(n)` and `v(n)` terms — only the
  broadened taps `qγ_k_i` interact with the delay lines. Per §3.4
  the per-sample memory `(s(n-1..n-10), v(n-1..n-10))` is initialised
  to zero at construction and never reset thereafter; coefficient
  swaps, the §3.4.1 non-speech disable switch, and per-frame
  freeze-and-swap updates all leave the memory untouched.
- **`Encoder::apply_weighting_filter(s)` +
  `Encoder::weighting_filter()` accessor.** The encoder now carries a
  live block-4 filter initialised to the §3.4 / §3.4.1 all-pass
  state. `set_weighting_filter_coeff_from_lpc` and
  `disable_weighting_filter` propagate the new coefficient set into
  the live filter (preserving its per-sample memory across the
  swap), and `apply_weighting_filter` consumes one
  `FRAME_LEN`-sample input vector to emit `v(n)`.
- **Block 10 is still NOT wired up.** Its §3.5 zero-input-response
  job requires the cascade `F(z)·W(z)` with the synthesis filter's
  pre-/post-save memory dance from §3.10 and is queued for a later
  round.

Round 248 (preserved) lands the **perceptual weighting filter
coefficient calculator** — block 38 of the encoder's perceptual
weighting filter adapter (Figure 4a/G.728), the third and last
sub-block of §3.3 (after hybrid window 36 + Levinson 37 transcribed
earlier in the shared `hybrid_window` / `levinson` modules):

- **`WeightingFilterCoeff` — block 38.** New `weighting_filter_coeff`
  module transcribes the mechanical substitution `z ← z/γ_k` of
  eq. 3-4b / 3-4c. Given the order-10 LPC predictor `q_i` from block
  37 (eq. 3-3a with `q_0 = 1`), block 38 produces the bandwidth-
  broadened numerator `qγ₁_i = q_i · γ₁ⁱ` and denominator
  `qγ₂_i = q_i · γ₂ⁱ` for `i = 1..=LPCW`. The `(γ₁, γ₂) = (WZCF,
  WPCF) = (0.9, 0.6)` paragraph values of Table 1/G.728 are used by
  default; the same routine plumbed through `from_lpc_with_gammas`
  realises the §3.4.1 "non-speech, `γ₁ = γ₂ = 0`" disabled path
  without a branch in the production transform. Both output sequences
  carry the spec's 1-based layout (`q_gamma_k[0] = 1.0` implicit
  leading tap), matching the `SynthesisAdapter::coefficients`
  convention already used elsewhere in the crate so future
  apply-the-weighting-filter rounds (blocks 4 / 10) can reuse the
  same "skip element 0, dot product over `1..=LPCW`" idiom.
- **`Encoder::weighting_filter_coeff()` accessor +
  `Encoder::set_weighting_filter_coeff_from_lpc()` setter +
  `Encoder::disable_weighting_filter()` switch.** The encoder now
  carries a `WeightingFilterCoeff` field initialised to the
  §3.4 / §3.4.1 all-pass `W(z) = 1`. The setter is the spec-direct
  block-38 entry point (the future block-37 wiring round feeds its
  Levinson output straight through); the disable switch flips the
  §3.4.1 non-speech mode without re-running the transform.
- **Block 4 / block 10 — applying the filter to input speech /
  zero-input response — is NOT wired up.** Those carry their own
  per-vector delay-line state plus the §3.4 special initialisation
  rule ("filter memory should not be reset to zero at any time
  except during initialisation") and are queued for a later round.

Round 235 lands the **typed encoder scaffold** plus the **§3.9
precomputed `E_j` shape-codevector energy table** that the
analysis-by-synthesis search of equation 3-23 will consume:

- **`Encoder` — typed front-end scaffold.** New `encoder` module exposes
  a stable [`Encoder`] type and a [`make_encoder`] factory mirroring the
  decoder's dual-API convention (`make_decoder`). The encoder carries
  the two backward adapters (synthesis-filter adapter §3.7 = block 23,
  log-gain adapter §3.8 = block 20) plus the `ET(1..IDIM)` 1-vector
  delay slot, all initialised to Table 2/G.728 values. §4.4 / §4.5 of
  the Recommendation require the backward adapters to be bit-for-bit
  identical at both ends of the channel, so the encoder reuses the
  exact `SynthesisAdapter` / `GainAdapter` types the decoder already
  drives.
- **`Y_ENERGY` — `E_j = Σ_k y_j(k)²`.** New `[f64; NCWD]` constant in
  `tables.rs`, derived at compile time from the already-transcribed
  Annex B `Y_Q11` shape codebook. §3.9 of the Recommendation rearranges
  the analysis-by-synthesis distortion (eqs. 3-14..3-23) into
  `D_{i,j} = b_i · <x̃, ỹ_j> + c_i · E_j` (with `b_i = 2·g_i`,
  `c_i = g_i²`), where `E_j` is a constant of the codebook and
  therefore precomputable once at table-load time. Exposed via
  `Encoder::shape_energy()` for the future analysis-by-synthesis
  search; the per-test in `tables.rs` cross-checks every entry against
  a direct dot-product computation in `y_f64()` space to machine
  precision (and against a hand-computed value for row 0).
- **`Encoder::encode_vector` returns `Error::NotImplemented`.** The
  signature is final (one input vector → one 10-bit channel index)
  but the §3.9 search loop (perceptual-weighting filter + impulse
  response + zero-state response correlation + shape-codebook scan +
  gain-quantiser decision tree) is intentionally absent in round 235.
  Future rounds wire blocks 1..28 + 67..70 against the already-typed
  adapter surfaces this round establishes.

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

Round 229 lands the **third and fourth stages of the §4.7
pitch-extraction pipeline** — blocks 83 and 84:

- **`PitchPostfilterCoeff` — blocks 83 + 84.** New module transcribes
  the §4.7 single-tap pitch-predictor weight and the long-term
  postfilter coefficient calculator. Block 83 maintains a 245-sample
  `sd(−239..5)` decoded-speech buffer (separate from block 71's
  `KPMAX`-sample comb-filter delay line so block 71's sample-rate
  cost is unchanged) and at the third vector of each frame evaluates
  the eq. 4-12 correlation `β = Σ_{k=−99..0} sd(k)·sd(k−p) /
  Σ_{k=−99..0} sd(k−p)²` clamped into `[0, 1]`, where `p` is the
  pitch period from block 82. Block 84 maps `β` to `(b, g_l)` via
  the eq. 4-13 / 4-14 piecewise table:
  - `β < PPFTH = 0.6` → `b = 0, g_l = 1` (postfilter off, `H_l(z) = 1`).
  - `0.6 ≤ β ≤ 1` → `b = PPFZCF · β = 0.15·β, g_l = 1/(1+b)`.
  - `β > 1` → `b = PPFZCF = 0.15, g_l = 1/(1+0.15)` (clamped form
    handled inside; the `[0, 1]` clamp on β makes this branch
    unreachable from compute_coefficients in practice but it is kept
    explicit to mirror the spec table).
- **`Decoder::decode_vector_postfiltered` drives blocks 83 + 84 each
  third-vector extract.** Every decoded-speech vector is pushed into
  the coefficient calculator's `sd(−239..5)` buffer; at the third
  vector of each frame the pitch period from `PitchSearch::extract`
  is fed to `PitchPostfilterCoeff::compute_coefficients`, and the
  resulting `(g_l, b, p)` triple is applied to
  `LongTermPostfilter::set_coefficients`. From that point on the comb
  filter operates at the spec-prescribed coefficients for the
  upcoming adaptation cycle.
- **`Decoder::pitch_pf_coeff()` accessor** for tests and audit.

With blocks 81 + 82 + 83 + 84 + 71 in place the §4.7 long-term
postfilter pipeline is end-to-end complete: synthesis (block 32) →
LPC inverse filter (block 81) → pitch extractor (block 82) →
coefficient calculator (blocks 83 / 84) → long-term comb (block 71)
→ short-term postfilter (block 72) → AGC (blocks 73–77).

Round 223 (preserved) lands the **second stage of the §4.7
pitch-extraction pipeline** — block 82, the pitch period extractor:

- **`PitchSearch` — block 82.** New module transcribes the spec's
  §4.7 pitch-period dataflow: a 240-sample LPC-residual buffer
  covering `d(−139..100)`, a 60-sample decimated-residual buffer
  covering `d̄(−34..25)`, the third-order Annex D 1 kHz elliptic
  lowpass + 4:1 decimator (`AL` / `BL` coefficients from
  `tables.rs`), the coarse correlation search `ρ(i) = Σ d̄(n)·d̄(n−i)`
  over decimated lags `5..=35` (eq. 4-7), the full-resolution
  refinement `C(i) = Σ d(k)·d(k−i)` over `4τ−3..=4τ+3` (eq. 4-8)
  clamped into `[KPMIN, KPMAX]`, and the fundamental-vs-multiple
  resolution: when `p̂` is outside `[p0 − KPDELTA, p0 + KPDELTA]`,
  a second `C(i)` search over `p̂−KPDELTA..=p̂+KPDELTA` produces
  `p1`, then the single-tap predictor weights `β0 = C(p0)/Σ d(k−p0)²`
  and `β1 = C(p1)/Σ d(k−p1)²` clamped into `[0, 1]` (eqs. 4-9,
  4-10) select `p = p0` if `β1 ≤ TAPTH·β0` else `p1` (eq. 4-11).
  The output `p ∈ [KPMIN, KPMAX] = [20, 140]` is stashed as `p̂`
  for the next frame.
- **Spec slot ordering for `d(81..85)` / `d(86..90)` / `d(91..95)` /
  `d(96..100)`.** The §4.7 prose explicitly assigns the 4th vector
  of the previous frame to `d(81..85)`, current frame's 1st/2nd/3rd
  vectors to `d(86..90)` / `d(91..95)` / `d(96..100)`. `PitchSearch`
  realises this via a per-frame cursor that's reset to 3 after each
  extract (so the next push — the 4th vector — lands at `d(81..85)`
  of the post-slide buffer) and then wraps to 0 at the start of the
  next frame.
- **`Decoder::decode_vector_postfiltered` advances block 82 each call.**
  After the inverse filter (block 81) produces the residual vector,
  it is pushed into `PitchSearch`'s residual buffer; at the third
  vector of every adaptation cycle (spec ICOUNT = 3) the extractor
  runs the full lowpass + decimate + correlation + refinement +
  resolution pipeline and stashes `p` as `p̂`. The output is not
  yet consumed downstream — the long-term postfilter (block 71)
  still follows the §4.6.1 passthrough until blocks 83 (single-tap
  `β` over the decoded-speech buffer, eq. 4-12) and 84 (the
  `(g_l, b)` calculator of eq. 4-13 / 4-14) land in a later round.
- **`Decoder::pitch_search()` accessor** for tests and audit.

Round 220 (preserved) lands the **first stage of the §4.7
pitch-extraction pipeline** — block 81, the 10th-order LPC inverse
filter:

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

- A-law / µ-law PCM I/O — `oxideav-g711` handles this per §5.3 / §3.1.
- The encoder pipeline (blocks 1..28 + 67..70 + the codebook search
  of §3.9 / blocks 12..18). Round 235 lands the typed encoder front
  end (`Encoder` / `make_encoder` / `Encoder::encode_vector`) so the
  symbol is available now; round 248 adds block 38 (the perceptual-
  weighting coefficient calculator) as a typed transform consumable
  via `Encoder::set_weighting_filter_coeff_from_lpc`; round 249 wires
  block 4 (the §3.4 application of the filter to input speech) as
  `Encoder::apply_weighting_filter`; round 258 wires blocks 36 + 37
  (the §3.3 hybrid window + Levinson on the input-speech buffer) as
  `Encoder::adapt_weighting_filter` /
  `commit_weighting_filter_coefficients`; round 267 wires blocks 9 +
  10 + 11 (the §3.5 / §3.6 zero-input response + VQ target) as
  `Encoder::compute_vq_target`. The remaining work — the §3.10
  memory-update phase (adding the zero-state response of the chosen
  excitation `e(n)` back onto the saved ring memory and reading
  `sq(n)` out of the top `STATELPC` taps) and the full §3.9 codebook
  search loop (impulse-response convolution, time-reversed
  correlation, shape-codebook scan, gain-quantiser decision tree) —
  still surface `Error::NotImplemented` from `encode_vector`. The
  search and the memory update both depend on the chosen excitation
  `e(n)` that §3.9 selects, so they land together in a later round.
  The shared backward adapters (block 23 in the
  encoder = block 33 in the decoder; block 20 in the encoder = block
  30 in the decoder) are reused unchanged via the existing
  `SynthesisAdapter` / `GainAdapter` types per §4.4 / §4.5.
- Annex G fixed-point variant and Annex I frame-loss concealment
  remain deferred behind the floating-point decoder.

## Clean-room provenance

Every numeric value lives in `src/tables.rs` and is transcribed
directly from the integer columns of Annexes A, B, C, D in the
ITU-T G.728 1992-09 PDF that lives under `docs/audio/g728/`. Round
235's `Y_ENERGY` table is a **derived** quantity from the Annex B
shape codebook — not a separately printed column — computed at
compile time as `Σ_k (Y_Q11[j][k] / 2¹¹)²` per §3.9 equation 3-23,
with a per-test cross-check against the same `y_f64()` dot product
the rest of the crate already consumes. Round 248's
`WeightingFilterCoeff` is similarly a derived quantity: the spec-
paragraph γ₁ = 0.9, γ₂ = 0.6 of Table 1/G.728 power-multiplied
against the caller-supplied order-10 `q_i` per eq. 3-4b / 3-4c, with
per-test cross-checks against (i) the unity-q geometric progression
of `γ_k^i`, (ii) a hand-traced `q_i = (−1/2)^i` term-by-term, and
(iii) the `γ₂ < γ₁ < 1` ordering of the Table 1 paragraph values.
Every control-flow line in the `hybrid_window`, `synthesis_adapter`,
`gain_adapter`, `long_term_postfilter`, `short_term_postfilter`,
`pitch_inverse_filter`, `pitch_search`, `pitch_postfilter_coeff`,
`agc`, `encoder`, `weighting_filter_coeff`, `weighting_filter` and
`zero_input_response`
modules carries a
comment pointing at the spec's §3.3 / §3.4 / §3.5 / §3.6 / §3.7 /
§3.8 / §3.9 / §4.4 / §4.5 / §5.6 / §5.7 / §5.9 / §5.10 / §4.6 / §4.7
pseudocode
for blocks 36/43/49 (hybrid window), 50 (Levinson — already in
`levinson.rs`), 51 / 45 (bandwidth expansion), 67/39/40/42/46/47/48
(the per-vector gain chain), 71 (long-term comb filter
`g_l · (1 + b·z^{-p})` from eq. 4-1, KPMIN/KPMAX clamp), 72
(short-term postfilter pole-zero + tilt cascade with
`SPFZCF / SPFPCF / TILTF` and the `ã_i = -a_i` sign-convention flip
spelled out in `short_term_postfilter.rs`), 73/74/75/76/77 (the AGC
tail with the lowpass form spelled out in §4.6 immediately after
eq. 4-5), 81 (10th-order LPC inverse filter
`Ã(z) = 1 − Σ ã_i · z^{-i}` from eq. 4-6, same `ã_i = -a_i` sign
flip spelled out in `pitch_inverse_filter.rs`), and 82 (the §4.7
pitch period extractor: 240-sample `d(−139..100)` residual buffer
with the §4.7-prescribed `d(81..85) / d(86..90) / d(91..95) /
d(96..100)` vector slot ordering, 60-sample `d̄(−34..25)` decimated
buffer fed by the Annex D third-order elliptic lowpass + 4:1
decimator, coarse-search `ρ(i)` over decimated lags `5..=35` from
eq. 4-7, full-resolution refinement `C(i)` over `4τ−3..=4τ+3`
clamped into `[KPMIN, KPMAX]` from eq. 4-8, and the
fundamental-vs-multiple `β0`/`β1` resolution of eqs 4-9..4-11 with
`TAPTH = 0.4`), and 83 + 84 (the §4.7 long-term postfilter
coefficient calculator: 245-sample `sd(−239..5)` decoded-speech
buffer, single-tap predictor weight `β = Σ_{k=−99..0} sd(k)·sd(k−p)
/ Σ_{k=−99..0} sd(k−p)²` from eq. 4-12 clamped into `[0, 1]`, and
the piecewise `(b, g_l)` map of eqs. 4-13 / 4-14 with `PPFTH = 0.6`
and `PPFZCF = 0.15` driving `LongTermPostfilter::set_coefficients`
at the third vector of each frame). No external implementation has been opened or
consulted during this rebuild — the per-test cross-checks (peak Q15
values, AR(1) predictor round-trip, `λⁱ` geometric progression,
codebook row spot-checks, ICOUNT cycle ordering, limiter clamp
bounds, first-vector σ(n) = 10^(GOFF/20), AGC unity-DC convergence
+ AGCFAC geometric decay, postfilter identity at cold start +
non-identity divergence after the first adaptation cycle +
bandwidth-expansion `ã_i · λ^i` invariants + tilt `µ = TILTF · k1`
invariant + block-81 impulse response matching `−ã_i` + AR(1)
synthesis/inverse round-trip to zero excitation residual + block-82
pitch lock at periods `KPMIN`, `KPMAX` and `40` on a unit-impulse
train + spec-slot vector placement at `d(81..85) / d(86..90) /
d(91..95) / d(96..100)` + post-extract residual buffer slide by
`NFRSZ = 20` + cold-start `p̂ = KPMIN` invariant + block-83/84
unity-β invariant on a perfectly periodic signal at the analysis
lag + spec-table `b = PPFZCF·β`, `g_l = 1/(1+b)` mapping over the
voiced range + sub-threshold β routing to `(b, g_l) = (0, 1)`
postfilter-off + buffer slide by `IDIM` per push + cold-start
identity `sf = sd` over the first frame + propagation of
`(g_l, b)` from `PitchPostfilterCoeff::last_*` into
`LongTermPostfilter` at every third-vector extract) act as in-repo
audit anchors against transcription typos.

## License

MIT — see [LICENSE](LICENSE).
