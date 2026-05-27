//! # oxideav-g728
//!
//! **Status:** orphan-rebuild scaffold (reset 2026-05-27).
//!
//! The prior implementation's numeric tables had a provenance that the
//! workspace clean-room policy does not permit: they were extracted from
//! an external reference C distribution rather than transcribed from the
//! Recommendation prose. The policy forbids consulting any external
//! implementation's source for any reason (including a normative
//! reference shipped alongside a spec), so the provenance could not be
//! defended.
//!
//! The crate will be re-implemented from scratch against the staged
//! ITU-T G.728 Recommendation prose + a future clean-room behavioural
//! trace document in a subsequent round, once that trace is staged under
//! `docs/audio/g728/`. Until then every public API returns
//! [`Error::NotImplemented`].

#![warn(missing_debug_implementations)]

use oxideav_core::RuntimeContext;

/// Crate-local error type. Until the clean-room rebuild lands every
/// public API path returns [`Error::NotImplemented`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// The crate has been reset to a scaffold pending clean-room
    /// rebuild; no decoder or encoder functionality is wired up yet.
    NotImplemented,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "oxideav-g728: orphan-rebuild scaffold — no codec wired up"
        )
    }
}

impl std::error::Error for Error {}

/// No-op codec registration — the orphan-rebuild scaffold registers
/// nothing into the runtime context.
pub fn register(_ctx: &mut RuntimeContext) {}

oxideav_core::register!("g728", register);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_displays() {
        assert!(Error::NotImplemented.to_string().contains("scaffold"));
    }

    #[test]
    fn register_is_no_op() {
        let mut ctx = RuntimeContext::default();
        register(&mut ctx);
    }
}
