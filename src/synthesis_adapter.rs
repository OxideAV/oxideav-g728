//! Backward synthesis-filter adapter — block 33 of Figure 5/G.728.
//!
//! Chains the three sub-blocks the spec lists in §5.6:
//!
//! * **Block 49** — hybrid window on the quantised-speech buffer
//!   `STTMP(1..NFRSZ)`, producing autocorrelation taps `RTMP(1..LPC+1)`.
//! * **Block 50** — Levinson-Durbin recursion on `RTMP` producing the
//!   bandwidth-unexpanded predictor `ATMP(1..LPC+1)`.
//! * **Block 51** — bandwidth expansion: `ATMP(I) = FACV(I) · ATMP(I)`;
//!   the resulting `ATMP` is copied into the live predictor `A` at
//!   the third vector of the adaptation cycle (`ICOUNT = 3`).
//!
//! ## Cycle timing (§3.7, §5.6 block 50 prose, block 51 prose)
//!
//! G.728's adaptation cycle is `NUPDATE = 4` vectors (`NFRSZ = 20`
//! samples). The block-50 prose spells out the staggered update:
//!
//! > Although the autocorrelation `RTMP` array is available at the
//! > first vector of each adaptation cycle, the actual updates of
//! > synthesis filter coefficients will not take place until the
//! > third vector. This intentional delay of updates allows the
//! > real-time hardware to spread the computation of this module
//! > over the first three vectors of each adaptation cycle.
//!
//! In a one-shot software implementation we don't need to spread
//! across vectors — but the *timing semantics* still matter, because
//! the adapter must observe the full previous cycle's quantised
//! speech before updating, not the current cycle's. We therefore
//! accept one buffered call per cycle (the cycle-boundary delivery
//! of `NFRSZ` samples) and produce one expanded predictor per cycle.

use crate::consts::{LPC, NFRSZ, NONR};
use crate::hybrid_window::{HybridWindow, HybridWindowState};
use crate::levinson::{levinson_durbin, LevinsonError};
use crate::tables::{facv_f64, wnr_f64, FACV_Q14, Q14};

/// Static (table-derived) parameters for the synthesis-filter
/// hybrid window. Constructed once per adapter; passed by reference to
/// the per-cycle [`SynthesisAdapter::adapt`] call.
#[derive(Debug)]
struct SynthCfg {
    wnr: [f64; 105],
    facv: [f64; LPC + 1],
}

impl SynthCfg {
    fn new() -> Self {
        Self {
            wnr: wnr_f64(),
            facv: facv_f64(),
        }
    }
}

/// Backward synthesis-filter adapter (block 33) state.
///
/// Carries the hybrid-window state (`SB` and `REXP` per Table 2) and
/// the most recently produced bandwidth-expanded predictor. The
/// predictor is in spec `A` layout: `coefficients()[0] = 1.0`,
/// `coefficients()[i] = -aᵢ` for `1 ≤ i ≤ LPC`.
#[derive(Debug)]
pub struct SynthesisAdapter {
    cfg: SynthCfg,
    hw_state: HybridWindowState,
    /// Most recently produced predictor (in spec `A` layout: leading
    /// `1.0`, then `-aᵢ` × LPC). Defaults to the all-pass predictor
    /// at construction time so the adapter is safe to read from before
    /// any speech has been pushed.
    last_predictor: [f64; LPC + 1],
}

impl Default for SynthesisAdapter {
    fn default() -> Self {
        Self::new()
    }
}

impl SynthesisAdapter {
    /// Construct a fresh adapter with all internal buffers initialised
    /// to the Table 2/G.728 values (`SB`, `REXP` zeroed; predictor
    /// `A(1) = 1`, `A(i ≥ 2) = 0`).
    pub fn new() -> Self {
        let cfg = SynthCfg::new();
        let hw = Self::window_descriptor(&cfg.wnr);
        let hw_state = HybridWindowState::new(&hw);

        let mut last_predictor = [0.0f64; LPC + 1];
        last_predictor[0] = 1.0;

        Self {
            cfg,
            hw_state,
            last_predictor,
        }
    }

    fn window_descriptor(wnr: &[f64]) -> HybridWindow<'_> {
        HybridWindow {
            m: LPC,
            l: NFRSZ,
            n: NONR,
            window: wnr,
        }
    }

    /// Most recently produced bandwidth-expanded predictor, in spec
    /// `A` layout. Borrowed read-only so the synthesiser can call
    /// `set_predictor` with an owned copy.
    pub fn coefficients(&self) -> &[f64; LPC + 1] {
        &self.last_predictor
    }

    /// Process one adaptation cycle of quantised-speech samples
    /// (`NFRSZ = 20` samples) through blocks 49, 50 and 51, and update
    /// the cached predictor.
    ///
    /// Returns the up-to-date predictor on success. On Levinson-Durbin
    /// failure (any of the three ill-conditioning exits) the cached
    /// predictor from the previous cycle is left in place per the
    /// spec's "keep using the previous adaptation cycle" rule, and the
    /// error is propagated so the caller can log/trace it.
    pub fn adapt(&mut self, sttmp: &[f64]) -> Result<&[f64; LPC + 1], LevinsonError> {
        assert_eq!(sttmp.len(), NFRSZ, "adapter input must be NFRSZ samples");

        let hw = Self::window_descriptor(&self.cfg.wnr);

        // ----- Block 49: hybrid window → RTMP -------------------------
        let mut rtmp = vec![0.0f64; LPC + 1];
        self.hw_state.run(&hw, sttmp, &mut rtmp);

        // ----- Block 50: Levinson-Durbin → ATMP -----------------------
        // The recursion writes ATMP into a fresh buffer; on failure
        // (any of the three ill-conditioning exits) levinson_durbin
        // populates the buffer with the canonical all-pass predictor,
        // but per the spec we must KEEP the previous adaptation
        // cycle's predictor on failure, not adopt the all-pass reset.
        // So we run into a scratch buffer and only commit on Ok.
        let mut atmp = vec![0.0f64; LPC + 1];
        levinson_durbin(&rtmp, &mut atmp, LPC)?;

        // ----- Block 51: bandwidth expansion --------------------------
        // | For I = 2,3,...,LPC+1: ATMP(I) = FACV(I) * ATMP(I)
        //
        // Note FACV(1) = 16384 (Q14) = 1.0; the spec writes the loop
        // starting at I = 2 because the leading tap is implicitly 1.0.
        // We multiply every tap; the I=1 multiplication is a no-op
        // and folds in cleanly.
        for i in 0..=LPC {
            atmp[i] *= self.cfg.facv[i];
        }

        // Commit ATMP into the cached predictor (this is the spec's
        // `Wait until ICOUNT = 3, then A(I) = ATMP(I)` step — in our
        // batched cycle model the wait collapses to "now").
        self.last_predictor[..=LPC].copy_from_slice(&atmp[..=LPC]);

        Ok(&self.last_predictor)
    }
}

/// Compile-time guard: FACV(1) in Q14 must be exactly 16384
/// (i.e. 1.0), so the implicit `A(1) = 1` invariant survives the
/// bandwidth-expansion multiplication.
const _: () = {
    let scale = 1i32 << Q14;
    assert!(FACV_Q14[0] as i32 == scale);
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::IDIM;
    use crate::tables::FACV_Q14;

    #[test]
    fn fresh_adapter_has_allpass_predictor() {
        // Table 2 initial: A(1) = 1, A(i ≥ 2) = 0 (the synthesis
        // filter is the identity before any adaptation has happened).
        let adapter = SynthesisAdapter::new();
        let a = adapter.coefficients();
        assert_eq!(a[0], 1.0);
        assert!(a[1..].iter().all(|&v| v == 0.0));
    }

    #[test]
    fn facv_first_tap_is_unity_q14() {
        // Block 51 multiplies every predictor tap by FACV(I); for the
        // I=1 tap (which must remain 1.0 to preserve the A(1)=1
        // invariant) FACV(1) must be exactly 16384 in Q14.
        assert_eq!(FACV_Q14[0], 16_384);
    }

    #[test]
    fn zero_input_leaves_predictor_as_allpass_after_levinson_failure() {
        // Driving zero speech vectors through the adapter produces a
        // zero autocorrelation. Block 50's R(1) ≤ 0 guard fires, the
        // spec says "keep using the previous cycle's predictor", and
        // the previous predictor is the all-pass identity. Confirm
        // both branches: the error is surfaced AND the cached
        // predictor stays at all-pass.
        let mut adapter = SynthesisAdapter::new();
        let input = [0.0f64; NFRSZ];
        let result = adapter.adapt(&input);
        assert!(
            matches!(
                result,
                Err(LevinsonError::ZeroSignal | LevinsonError::TrailingZero)
            ),
            "expected Levinson to flag zero signal, got {:?}",
            result
        );
        let a = adapter.coefficients();
        assert_eq!(a[0], 1.0);
        assert!(a[1..].iter().all(|&v| v == 0.0));
    }

    #[test]
    fn nonzero_input_produces_a1_eq_unity_and_nonzero_taps() {
        // A non-trivial speech signal should drive the autocorrelation
        // away from zero and let Levinson produce a non-trivial
        // predictor. The A(1) = 1 invariant must be preserved by
        // block 51 (FACV(1) = 1 in Q14, hand-verified above).
        let mut adapter = SynthesisAdapter::new();
        // Use a frequency-rich signal (mix of two sines + small noise)
        // so Levinson has enough excitation to land away from the
        // ill-conditioning boundary.
        let mut input = [0.0f64; NFRSZ];
        for k in 0..NFRSZ {
            let t = k as f64;
            input[k] = 100.0 * (0.4 * t).sin() + 50.0 * (0.9 * t).cos() + (t * 13.0).fract();
        }
        // Run enough cycles to populate the recursive part of the
        // hybrid window. A single cycle on a fresh adapter might still
        // produce a degenerate Levinson on some signals; this test
        // exercises the typical, steady-state behaviour.
        let mut last = None;
        for _ in 0..6 {
            if let Ok(a) = adapter.adapt(&input) {
                last = Some(*a);
            }
        }
        let a = last.expect("at least one adapt() should succeed");
        assert_eq!(a[0], 1.0, "A(1) = 1 invariant must hold");
        assert!(
            a[1..].iter().any(|&v| v != 0.0),
            "expected at least one non-zero predictor tap after adaptation"
        );
    }

    #[test]
    fn adapter_drives_decoder_without_panicking() {
        // Smoke test: wire the adapter into the existing decoder hook
        // and confirm the predictor swap path is exercised cleanly.
        // We round-trip three adaptation cycles' worth of vectors.
        use crate::decoder::Synthesizer;

        let mut adapter = SynthesisAdapter::new();
        let mut synth = Synthesizer::new();

        for cycle in 0..3 {
            // Run NUPDATE=4 vectors per cycle, IDIM samples each.
            let mut sttmp = [0.0f64; NFRSZ];
            for k in 0..NFRSZ {
                let t = (cycle * NFRSZ + k) as f64;
                sttmp[k] = 200.0 * (0.3 * t).sin();
            }
            let _ = adapter.adapt(&sttmp); // may be Err on the first cycle
            synth.set_predictor(*adapter.coefficients());

            for v in 0..(NFRSZ / IDIM) {
                let mut exc = [0.0f64; IDIM];
                for k in 0..IDIM {
                    exc[k] = 10.0 * ((cycle * NFRSZ + v * IDIM + k) as f64).sin();
                }
                let out = synth.filter_vector(&crate::decoder::ExcitationVector(exc));
                for &v in &out {
                    assert!(v.is_finite(), "decoder produced NaN/Inf with adapter");
                }
            }
        }
    }
}
