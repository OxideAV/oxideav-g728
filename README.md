# oxideav-g728

Pure-Rust ITU-T G.728 Low-Delay CELP (LD-CELP, 16 kbit/s) speech codec.

## Status: orphan-rebuild scaffold (reset 2026-05-27)

The previous implementation was retired under the OxideAV clean-room
policy. Its numeric tables had a provenance the policy does not permit:
they were extracted from an external reference C distribution rather
than transcribed from the Recommendation prose. The policy forbids
consulting any external implementation's source for any reason —
including a normative reference shipped alongside a spec — so the
provenance could not be defended.

The crate will be re-built from scratch against the staged ITU-T G.728
Recommendation prose plus a future clean-room behavioural trace document
in a subsequent round, once that trace is staged under
`docs/audio/g728/`. Until then every public API returns
`Error::NotImplemented`.

## License

MIT — see [LICENSE](LICENSE).
