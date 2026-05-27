# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.0.8](https://github.com/OxideAV/oxideav-g728/compare/v0.0.7...v0.0.8) - 2026-05-27

### Other

- port ITU-T §3.10 block-43 hybrid window for the log-gain predictor

### Other

- port ITU-T §3.10 block-43 hybrid window for the 10th-order log-gain
  predictor: 34-sample `SBLG` buffer (`LPCLG=10`, `NUPDATE=4`, `NONRLG=20`),
  `3/4` recursive decay, `257/256` white-noise correction, 34-sample
  `WNRLG` table transcribed from the staged `loggain-hybrid-window.csv`.
  Wires `GainHybridWindow` into `GainPredictor::push` (flushes a cycle
  every 4 samples) and `update_gain_predictor_from_hybrid_r` into the
  refresh path with the spec's "skip blocks 44+45 when `R(LPCLG+1)=0`"
  guard. The legacy Hamming-window path stays as the cold-start fallback.

## [0.0.7](https://github.com/OxideAV/oxideav-g728/compare/v0.0.6...v0.0.7) - 2026-05-06

### Other

- drop dead `linkme` dep
- auto-register via oxideav_core::register! macro (linkme distributed slice)
- unify entry point on register(&mut RuntimeContext) ([#502](https://github.com/OxideAV/oxideav-g728/pull/502))

## [0.0.6](https://github.com/OxideAV/oxideav-g728/compare/v0.0.5...v0.0.6) - 2026-05-03

### Other

- cargo fmt: pending rustfmt cleanup
- replace never-match regex with semver_check = false
- migrate to centralized OxideAV/.github reusable workflows
- adopt slim VideoFrame/AudioFrame shape
- pin release-plz to patch-only bumps

## [0.0.5](https://github.com/OxideAV/oxideav-g728/compare/v0.0.4...v0.0.5) - 2026-04-25

### Other

- drop oxideav-codec/oxideav-container shims, import from oxideav-core
- frame-erasure concealment (Annex A.3 / §5.8)

## [0.0.4](https://github.com/OxideAV/oxideav-g728/compare/v0.0.3...v0.0.4) - 2026-04-19

### Other

- update README + lib.rs docs to reflect §3.7 window + §5.5 postfilter
- add spectral-shape test for the postfilter
- port §5.5 adaptive postfilter (long-term + short-term + AGC)
- port ITU-T §3.7 hybrid (Barnwell) autocorrelation window
- drop Cargo.lock — this crate is a library
- bump oxideav-core / oxideav-codec dep examples to "0.1"
- bump to oxideav-core 0.1.1 + codec 0.1.1
- migrate register() to CodecInfo builder
- bump oxideav-core + oxideav-codec deps to "0.1"
- ship ITU-T G.728 Annex B CODEBK + GQ values
- correct packet layout and document status + limits
- add 'Quick use' example for standalone decode/encode
- loosen oxideav-* pins to '0.0' (accept any 0.0.x)
