# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.0.3](https://github.com/OxideAV/oxideav-g728/releases/tag/v0.0.3) - 2026-04-17

### Other

- add GitHub Actions (tests + release-plz auto-publish)
- make crate standalone (pin deps to crates.io, update README + LICENSE)
- add publish metadata (readme/homepage/keywords/categories)
- address workspace-wide lints to unblock CI
- cargo fmt across the workspace
- add LD-CELP encoder (exhaustive 128×8 analysis-by-synthesis)
- first-cut LD-CELP decoder from scaffold — machinery real, tables placeholder
- harvest 5 more agent scaffolds (g7231, g728, g729, gsm, speex)
