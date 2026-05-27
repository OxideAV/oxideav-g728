# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Changed

- **Reset to orphan-rebuild scaffold (2026-05-27).** The prior
  implementation was retired under the workspace clean-room policy: its
  numeric tables had a provenance the policy does not permit — they
  were extracted from an external reference C distribution rather than
  transcribed from the Recommendation prose. All public APIs now return
  `Error::NotImplemented` pending a clean-room rebuild against a staged
  ITU-T G.728 Recommendation prose plus a future behavioural trace
  document.
