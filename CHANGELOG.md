# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.0.6](https://github.com/OxideAV/oxideav-g728/compare/v0.0.5...v0.0.6) - 2026-04-26

### Other

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
