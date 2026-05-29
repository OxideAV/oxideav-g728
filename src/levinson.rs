//! Levinson-Durbin recursion (G.728 block 37 / block 44 / block 50).
//!
//! Transcribed line-for-line from the Fortran-like pseudocode in §5.6
//! of the Recommendation (block 50). Block 37 (perceptual weighting
//! adapter) and block 44 (log-gain predictor) use the same recursion
//! with different orders and input variables — implemented here once
//! and dispatched by passing the correct order in `lpc_order`.
//!
//! ## Indexing convention
//!
//! The spec writes the autocorrelation array `R(1..LPC+1)` and the
//! predictor array `A(1..LPC+1)` with a leading "1-based, A(1) is
//! always 1" convention (Table 1/G.728 closing paragraph). The Rust
//! transcription uses ordinary 0-based slices: `r[0]` is the spec's
//! `R(1)`, `r[k]` is the spec's `R(k+1)`, and the output `a[0]` is
//! the spec's `A(1) = 1`. The transcription is mechanical; see the
//! per-line comments tying each Rust statement to its source line in
//! block 50.
//!
//! ## Ill-conditioning
//!
//! The pseudocode aborts on three conditions:
//!
//! 1. `R(LPC+1) = 0` — last autocorrelation tap is exactly zero.
//! 2. `R(1) ≤ 0` — zero-signal (or negative-energy ill-conditioning).
//! 3. `ALPHA ≤ 0` during any iteration — the residual energy is non-
//!    positive, i.e. the autocorrelation matrix is not positive
//!    definite (most often a fixed-point underflow or a numerically
//!    constant input).
//!
//! When any of these fire, the spec says "skip block 51, do not
//! update the synthesis filter coefficients" — i.e. the caller keeps
//! the previous adaptation cycle's predictor. In this transcription
//! we surface that as a [`LevinsonError`] return; the caller maps it
//! to the same "keep previous predictor" behaviour.

/// Outcome of a Levinson-Durbin call.
///
/// The variants correspond directly to the three "go to LABEL" exits
/// in block 50. Callers treat all three as "keep the predictor from
/// the previous adaptation cycle" per the spec.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LevinsonError {
    /// `R(LPC + 1) = 0` — "skip if zero" branch (block 50 first test).
    TrailingZero,
    /// `R(1) ≤ 0` — "skip if zero signal" branch (block 50 second test).
    ZeroSignal,
    /// `ALPHA ≤ 0` during recursion — "abort if ill-conditioned"
    /// (block 50, after the residual-energy update).
    IllConditioned,
}

/// Run the Levinson-Durbin recursion for the given autocorrelation
/// sequence and LPC order.
///
/// * `r` — autocorrelation coefficients in spec 1-based convention
///   shifted to 0-based: `r[k]` carries `R(k+1)`. Must have length
///   at least `lpc_order + 1`.
/// * `lpc_order` — predictor order (block 37 uses LPCW = 10, block 44
///   uses LPCLG = 10, block 50 uses LPC = 50).
///
/// On success the predictor coefficients are written into `a` with
/// `a[0] = 1.0` (the spec's `A(1) = 1`) and the trailing
/// `a[1..=lpc_order]` filled in. The buffer must be at least
/// `lpc_order + 1` long; entries past `lpc_order` are left untouched.
///
/// On failure the predictor buffer is filled with the canonical
/// "all-pass" predictor `[1.0, 0.0, 0.0, ...]` so a careless caller
/// that ignores the error still produces a well-defined synthesis
/// filter.
pub fn levinson_durbin(r: &[f64], a: &mut [f64], lpc_order: usize) -> Result<(), LevinsonError> {
    assert!(
        r.len() > lpc_order,
        "autocorrelation buffer too short: need {} entries, got {}",
        lpc_order + 1,
        r.len()
    );
    assert!(
        a.len() > lpc_order,
        "predictor buffer too short: need {} entries, got {}",
        lpc_order + 1,
        a.len()
    );

    // Default-fill so that a propagated error still yields a well-
    // defined predictor (matches the "keep using the previous
    // adaptation cycle" intent of the spec).
    a[0] = 1.0;
    for slot in &mut a[1..=lpc_order] {
        *slot = 0.0;
    }

    // | If R(LPC + 1) = 0, go to LABEL              | Skip if zero
    if r[lpc_order] == 0.0 {
        return Err(LevinsonError::TrailingZero);
    }
    // | If R(1) ≤ 0, go to LABEL                    | Skip if zero signal
    if r[0] <= 0.0 {
        return Err(LevinsonError::ZeroSignal);
    }

    // We need a temporary reflection-coefficient buffer; the spec
    // uses an RC array of length LPC. Allocating on the stack via a
    // fixed maximum (LPC = 50) would force an unrelated dependency,
    // so we use a small heap vector — Levinson runs at most once per
    // 20-sample adaptation cycle, the allocation cost is negligible.
    let mut rc = vec![0.0f64; lpc_order + 1];

    // | RCTMP(1) = -RTMP(2) / RTMP(1)               |
    // | ATMP(1) = 1                                 | (already set above)
    // | ATMP(2) = RCTMP(1)                          | First-order predictor
    // | ALPHATMP = RTMP(1) + RTMP(2) * RCTMP(1)     |
    // | If ALPHATMP ≤ 0, go to LABEL                | Abort if ill-conditioned
    rc[0] = -r[1] / r[0];
    a[1] = rc[0];
    let mut alpha = r[0] + r[1] * rc[0];
    if alpha <= 0.0 {
        return Err(LevinsonError::IllConditioned);
    }

    // For MINC = 2,3,4,...,LPC, do the following:
    for minc in 2..=lpc_order {
        // | SUM = 0                                     |
        // | For IP = 1,2,3,...,MINC, do the next two lines
        // |     N1 = MINC - IP + 2                      |
        // |     SUM = SUM + RTMP(N1) * ATMP(IP)         |
        let mut sum = 0.0;
        for ip in 1..=minc {
            let n1 = minc - ip + 2; // spec 1-based
                                    // r is stored 0-based, so spec R(N1) is r[N1 - 1].
                                    // a is stored 0-based, so spec A(IP) is a[IP - 1].
            sum += r[n1 - 1] * a[ip - 1];
        }

        // | RCTMP(MINC) = -SUM / ALPHATMP               | Reflection coeff
        // | MH = MINC / 2 + 1                           |
        rc[minc - 1] = -sum / alpha;
        let mh = minc / 2 + 1;

        // | For IP = 2,3,4,...,MH, do the next four lines
        // |     IB = MINC - IP + 2                      |
        // |     AT = ATMP(IP) + RCTMP(MINC) * ATMP(IB)  |
        // |     ATMP(IB) = ATMP(IB) + RCTMP(MINC) * ATMP(IP)
        // |     ATMP(IP) = AT                           | Update predictor
        for ip in 2..=mh {
            let ib = minc - ip + 2;
            // Avoid aliasing when ip == ib (the midpoint of an odd-order
            // step) — split the read first, then assign both slots.
            let a_ip = a[ip - 1];
            let a_ib = a[ib - 1];
            let at = a_ip + rc[minc - 1] * a_ib;
            a[ib - 1] = a_ib + rc[minc - 1] * a_ip;
            a[ip - 1] = at;
        }

        // | ATMP(MINC + 1) = RCTMP(MINC)                |
        // | ALPHATMP = ALPHATMP + RCTMP(MINC) * SUM     | Residual energy
        // | If ALPHATMP ≤ 0, go to LABEL                | Abort if ill-conditioned
        a[minc] = rc[minc - 1];
        alpha += rc[minc - 1] * sum;
        if alpha <= 0.0 {
            return Err(LevinsonError::IllConditioned);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_tap_is_unity_after_run() {
        // Block 50 invariant: A(1) = 1 throughout.
        let r = [4.0, 1.5, 0.5, 0.1, 0.01];
        let mut a = [0.0; 5];
        levinson_durbin(&r, &mut a, 4).unwrap();
        assert_eq!(a[0], 1.0);
    }

    #[test]
    fn trailing_zero_returns_error_and_resets_predictor() {
        // R(LPC + 1) is exactly zero → first test fires.
        let r = [1.0, 0.5, 0.25, 0.0];
        let mut a = [42.0; 4];
        let err = levinson_durbin(&r, &mut a, 3).unwrap_err();
        assert_eq!(err, LevinsonError::TrailingZero);
        // Predictor reset to the canonical all-pass shape.
        assert_eq!(a[0], 1.0);
        assert!(a[1..].iter().all(|&x| x == 0.0));
    }

    #[test]
    fn zero_signal_returns_error() {
        // R(1) = 0 with a non-zero trailing tap so the trailing-zero
        // guard doesn't fire first → second guard fires.
        let r = [0.0, 0.0, 0.0, 0.5];
        let mut a = [0.0; 4];
        let err = levinson_durbin(&r, &mut a, 3).unwrap_err();
        assert_eq!(err, LevinsonError::ZeroSignal);
    }

    #[test]
    fn ar_p2_recovers_known_predictor() {
        // Cross-check the recursion against an analytically-known
        // first-order AR(1) autocorrelation: r[k] = ρ^|k| with
        // ρ = 0.8. For an order-1 predictor the spec recursion
        // should set a[1] = -ρ = -0.8 and the residual energy
        // should be (1 - ρ²) · r[0]. This catches sign flips and
        // off-by-one errors in the autocorrelation index mapping
        // without consulting any external reference.
        let rho: f64 = 0.8;
        let r = [1.0, rho, rho * rho, rho.powi(3)];
        let mut a = [0.0; 4];
        levinson_durbin(&r, &mut a, 3).unwrap();
        assert_eq!(a[0], 1.0);
        assert!((a[1] - -rho).abs() < 1e-12);
        // For an AR(1) process the higher-order taps are zero (the
        // recursion drops to RC = 0 from order 2 onward).
        assert!(a[2].abs() < 1e-12);
        assert!(a[3].abs() < 1e-12);
    }

    #[test]
    fn ar_p2_with_two_poles_recovers_both_taps() {
        // Build an autocorrelation for a 2-pole AR process with
        // predictor (1, -a1, -a2) and check that the recursion
        // recovers (a1, a2). Autocorrelation for a stable AR(2) with
        // distinct real poles p1, p2 satisfies:
        //     r[0] = 1, r[k] = p1·r[k-1] + p2·r[k-2] - p1·p2·r[k-2]
        // We just synthesize r[] numerically from the predictor
        // (1, 0.6, -0.2) by solving the Yule-Walker system
        // directly — i.e. forward-substitute under the same
        // recurrence the Levinson is the inverse of.
        let a1_true = 0.6;
        let a2_true = -0.2;
        // Yule-Walker (forward): r[1] = -a1·r[0] - a2·r[1]
        // closed form for 2nd-order AR:
        let r0 = 1.0_f64;
        let r1 = -a1_true / (1.0 + a2_true);
        let r2 = -a1_true * r1 - a2_true * r0;
        let r = [r0, r1, r2];
        let mut a = [0.0; 3];
        levinson_durbin(&r, &mut a, 2).unwrap();
        assert!((a[1] - a1_true).abs() < 1e-10);
        assert!((a[2] - a2_true).abs() < 1e-10);
    }

    #[test]
    fn high_order_completes_for_well_conditioned_input() {
        // 50th-order recursion on a smoothly decaying autocorrelation
        // should complete without ill-conditioning. The test confirms
        // the recursion scales to the synthesis-filter order used by
        // G.728's block 50.
        let mut r = [0.0f64; 51];
        for k in 0..51 {
            // Mildly absorbing exponential — produces a well-defined
            // 50th-order AR autocorrelation.
            r[k] = (-0.05 * k as f64).exp();
        }
        let mut a = [0.0f64; 51];
        let result = levinson_durbin(&r, &mut a, 50);
        assert!(result.is_ok(), "recursion failed: {:?}", result);
        assert_eq!(a[0], 1.0);
    }
}
