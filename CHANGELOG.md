# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

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
