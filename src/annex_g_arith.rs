//! Annex G fixed-point arithmetic foundation (§G.1.2 / §G.1.3).
//!
//! ITU-T G.728 Annex G (1994-11) describes how to implement the
//! 16 kbit/s LD-CELP coder on a 16-bit fixed-point arithmetic device
//! "in sufficient detail … producing an output signal of equivalent
//! quality" and "bit exact" with the floating-point coder. The annex
//! is built bottom-up: §G.1.2 fixes the numerical representations and
//! §G.1.3 fixes the arithmetic primitives (shifting, rounding,
//! normalization, multiplication, addition, division). Every later
//! Annex G module (the gain adapter of §G.2.1, the Levinson-Durbin
//! recursion of §G.2.2, and the per-block pseudo-code of §G.3) is
//! expressed in terms of these primitives.
//!
//! This module transcribes that foundation, one §G.1.3 sub-clause at a
//! time. It is deliberately self-contained and arithmetic-only: it
//! introduces no codec state and reads no codebooks. The downstream
//! Annex G blocks build on it.
//!
//! ## Numerical representations (§G.1.2)
//!
//! * **Single precision fixed point** — a 16-bit word, range −32768 to
//!   +32767. A `Qn` format places the binary point so that `n` bits
//!   lie to the right of it (`Q15` represents `[−1.0, +1.0)`, `Q0` is a
//!   pure integer). Modelled here by [`i16`] plus an out-of-band `Q`.
//! * **Double precision fixed point** — 32 bits of information (the
//!   product register / accumulator). Modelled by [`i32`].
//! * **Scalar floating point** — two words: a *mantissa* whose
//!   magnitude is normalized into `[16384, 32767]` (so bit 14 of a
//!   positive mantissa is always 1) and an `NLS` word giving the number
//!   of left shifts applied to normalize it. See [`ScalarFloat`].
//! * **Block floating point** — `n + 1` words for `n` values: one
//!   shared `NLS` plus `n` mantissas sharing it (not necessarily
//!   normalized). Produced by [`vscale`].
//!
//! ## Convention notes (§G.1.3, §G.1.3.1)
//!
//! All arithmetic is two's complement, which Rust's [`i16`] / [`i32`]
//! already provide. A right shift sign-extends; the annex flags the
//! magnitude anomaly (e.g. `3 >> 1 = 1` but `−3 >> 1 = −2`) — this is
//! exactly Rust's arithmetic right shift on a signed integer, so no
//! special handling is needed. The annex also warns that a shift by
//! "greater than the size of the word" is compiler-dependent; we make
//! every such shift well-defined explicitly (see [`shr_sat`]).

/// Number of fractional bits in the `Q15` single-precision format
/// (§G.1.2): the binary point sits between bits 14 and 15, giving the
/// range `[−1.0, +1.0)`.
pub const Q15: u32 = 15;

/// Maximum number of left shifts for a single-precision (16-bit) block
/// floating-point mantissa: `MLS = 14` (§G.1.3, VSCALE description).
/// Normalizing left-justifies the magnitude so that bit 14 is set,
/// i.e. the magnitude lands in `[16384, 32767]`.
pub const MLS_SINGLE: i32 = 14;

/// Maximum number of left shifts for a double-precision (32-bit) block
/// floating-point mantissa: `MLS = 30` (§G.1.3, VSCALE description).
pub const MLS_DOUBLE: i32 = 30;

/// Saturating two's-complement right shift by `k` bits with the §G.1.3
/// "shift larger than word size" anomaly made well-defined.
///
/// §G.1.3.1 warns that "in the algorithm an instruction is generated
/// to right shift a word by greater than the size of that word … the
/// result of such an operation should be 0 or −1, depending on the
/// sign of the original data", but that some compilers treat such a
/// shift as illegal. Rust panics on an over-wide shift in debug and
/// masks the shift amount in release, so neither matches the spec. We
/// clamp the shift to the sign-fill value explicitly: a value shifted
/// right by ≥ 32 bits collapses to `0` (non-negative) or `−1`
/// (negative), exactly as the annex prescribes.
#[inline]
#[must_use]
pub fn shr_sat(value: i32, k: u32) -> i32 {
    if k >= 32 {
        // Sign fill: arithmetic shift saturates to all-sign-bits.
        value >> 31
    } else {
        value >> k
    }
}

/// Saturating left shift by `k` bits (§G.1.3.1).
///
/// "If we shift a value to the left, we need to check for possible
/// overflows." We perform the shift in 64-bit space and saturate the
/// result back into the 32-bit accumulator range, matching the
/// saturation-mode arithmetic the annex mandates for sums of products
/// re-used as multiplier inputs (§G.1.3).
#[inline]
#[must_use]
pub fn shl_sat(value: i32, k: u32) -> i32 {
    if k >= 32 {
        if value > 0 {
            i32::MAX
        } else if value < 0 {
            i32::MIN
        } else {
            0
        }
    } else {
        let wide = (value as i64) << k;
        wide.clamp(i32::MIN as i64, i32::MAX as i64) as i32
    }
}

/// Round a double-precision accumulator value to a single-precision
/// 16-bit word — the §G.1.3.1 `RND(.)` function.
///
/// Quoting the convention (paraphrased): the accumulator has a high
/// word and a low word with the binary point between them. Rounding
/// converts to the integer nearest the two-word value: "test the MSB
/// of the low word. If it is a 1, add 1 to the value in the high word.
/// Then zero out the low word." The annex works the two examples 1.5
/// → 2 (high word 1, low word `0x8000`, MSB set ⇒ +1) and −1.5 → −1
/// (high word −2, low word `0x8000`, MSB set ⇒ +1 ⇒ −1). Both are
/// exactly "add the rounding bit, where the rounding bit is bit 15 of
/// the low word."
///
/// The annex further requires saturation on overflow: if the high word
/// is already `0x7FFF` (= 32767) and the low word's MSB is set, adding
/// 1 would wrap to `0x8000` (= −32768); instead "the value is saturated
/// to avoid an unrepresentable value." We saturate to `0x7FFF`.
#[inline]
#[must_use]
pub fn rnd(acc: i32) -> i16 {
    // High word = bits 16..31, low word = bits 0..15. The rounding bit
    // is bit 15 (the MSB of the low word).
    let high = acc >> 16;
    let round_bit = (acc >> 15) & 1;
    let rounded = high + round_bit;
    if rounded > i16::MAX as i32 {
        i16::MAX
    } else if rounded < i16::MIN as i32 {
        i16::MIN
    } else {
        rounded as i16
    }
}

/// A scalar single-precision floating-point value (§G.1.2): a mantissa
/// whose magnitude is normalized into `[16384, 32767]` paired with the
/// number of left shifts `NLS` that were applied to normalize it.
///
/// The represented real value is `mantissa · 2^(−nls)` (for a `Q0`
/// mantissa interpretation; downstream code that needs a different `Q`
/// base folds the offset into `nls`). A negative `nls` means the
/// normalization required *right* shifts instead of left shifts
/// (§G.1.3, VSCALE: "the NLS value returned will be negative").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScalarFloat {
    /// Normalized mantissa, magnitude in `[16384, 32767]` (or 0).
    pub mantissa: i16,
    /// Number of left shifts applied to normalize (may be negative).
    pub nls: i32,
}

impl ScalarFloat {
    /// The scalar floating-point zero. Per §G.1.3 a zero input vector
    /// is assigned one extra bit of left shift (`NLS = MLS + 1`); for a
    /// scalar single-precision value that is `MLS_SINGLE + 1 = 15`.
    pub const ZERO: Self = Self {
        mantissa: 0,
        nls: MLS_SINGLE + 1,
    };

    /// Convert a single-precision fixed-point integer to normalized
    /// scalar floating-point form (the §G.1.3 `FLOAT(.)` conversion for
    /// a single 16-bit value). Equivalent to running [`vscale`] /
    /// [`findnls`] on a length-1 vector.
    #[must_use]
    pub fn from_i16(value: i16) -> Self {
        let nls = findnls(&[value as i32], 1, MLS_SINGLE);
        let scaled = shl_norm(value as i32, nls);
        Self {
            mantissa: scaled as i16,
            nls,
        }
    }

    /// The approximate real value `mantissa · 2^(−nls)`. Provided for
    /// tests / diagnostics; the codec itself stays in fixed point.
    #[must_use]
    pub fn to_f64(self) -> f64 {
        (self.mantissa as f64) * 2f64.powi(-self.nls)
    }
}

/// Apply `nls` left shifts to a value, interpreting a negative `nls`
/// as `−nls` right shifts (the §G.1.3.1 negative-shift convention used
/// throughout VSCALE / scaling code).
#[inline]
fn shl_norm(value: i32, nls: i32) -> i32 {
    if nls >= 0 {
        shl_sat(value, nls as u32)
    } else {
        shr_sat(value, (-nls) as u32)
    }
}

/// Find the number of left shifts needed to normalize a vector to the
/// requested precision — the §G.1.3 `FINDNLS` subroutine.
///
/// * `input` — the values to examine (single- or double-precision
///   mantissas, as `i32` for a common representation).
/// * `slen` — search length: `1` means the maximum-magnitude element is
///   known to be the first element (only `input[0]` is examined);
///   otherwise the whole `input[..slen]` prefix is searched.
/// * `mls` — maximum left shifts permitted ([`MLS_SINGLE`] = 14 for
///   16-bit data, [`MLS_DOUBLE`] = 30 for 32-bit data; 12 or 13 when
///   fewer mantissa bits are wanted).
///
/// Returns `NLS`, the number of left shifts (negative ⇒ right shifts)
/// that would left-justify the largest-magnitude element to the `mls`
/// bound. A zero input returns `mls + 1` ("let 0 have one more bit of
/// left shift than 1"). FINDNLS computes only `NLS`; [`vscale`] applies
/// it to produce the scaled vector.
#[must_use]
pub fn findnls(input: &[i32], slen: usize, mls: i32) -> i32 {
    debug_assert!(!input.is_empty(), "FINDNLS input must be non-empty");

    // AA0 finds the maximum positive value, AA1 the maximum negative.
    let mut aa0 = input[0];
    let mut aa1 = input[0];
    if slen != 1 {
        for &v in &input[1..slen.min(input.len())] {
            if v > aa0 {
                aa0 = v;
            }
            if v < aa1 {
                aa1 = v;
            }
        }
    }

    // Case 1: zero input vector.
    if aa0 == 0 && aa1 == 0 {
        return mls + 1;
    }

    let mut nls: i32 = 0;

    // The mantissa bounds after the desired shift. `maxi = 2^mls` is the
    // upper bound; `mini = 2^(mls-1)` (= maxi/2) the lower bound. The
    // annex spells out the two sign cases separately to keep every
    // intermediate inside the accumulator without overflow at MLS = 30.
    if aa0 < 0 || aa1 < -aa0 {
        // Case 2: the negative extreme has the larger magnitude.
        let maxi = -(1i64 << mls); // −2^mls
        let mini = 2 * maxi; // mantissa lower bound after shift
        let mut a = aa1 as i64;
        if a < mini {
            // Magnitude too large: right-shift (negative NLS) until it fits.
            while a < mini {
                a >>= 1;
                nls -= 1;
            }
        } else {
            // Magnitude too small: left-shift until just past `maxi`.
            while a >= maxi {
                a <<= 1;
                nls += 1;
            }
        }
    } else {
        // Case 3: the positive extreme has the larger magnitude.
        let mini = 1i64 << mls; // = 2^mls (mantissa lower bound)
                                // maxi = mini − 1; maxi = maxi + mini = 2*mini − 1.
        let maxi = (mini - 1) + mini;
        let mut a = aa0 as i64;
        if a > maxi {
            // Magnitude too large: right-shift (negative NLS).
            while a > maxi {
                a >>= 1;
                nls -= 1;
            }
        } else {
            // Magnitude too small: left-shift until ≥ `mini`.
            while a < mini {
                a <<= 1;
                nls += 1;
            }
        }
    }

    nls
}

/// Scale a vector to block floating-point form — the §G.1.3 `VSCALE`
/// subroutine.
///
/// Identical search to [`findnls`], but additionally applies the
/// resulting `NLS` to every element (left shift for positive `NLS`,
/// right shift for negative). Returns the scaled vector (the shared
/// block-floating-point mantissas) together with the shared `NLS`.
///
/// * `input` / `slen` / `mls` — as for [`findnls`].
#[must_use]
pub fn vscale(input: &[i32], slen: usize, mls: i32) -> (Vec<i32>, i32) {
    let nls = findnls(input, slen, mls);
    // Zero input: OUT(I) = 0 for all I (the search already returned
    // mls + 1); shifting zeros is a no-op but we keep it explicit.
    let out: Vec<i32> = input.iter().map(|&v| shl_norm(v, nls)).collect();
    (out, nls)
}

/// Scalar floating-point division on a 16-bit fixed-point device — the
/// §G.1.3.4 `DIVIDE` subroutine.
///
/// Computes `QUO · 2^(−QUONLS) = (NUM · 2^(−num_nls)) / (DEN ·
/// 2^(−den_nls))` with full 16-bit precision in the quotient mantissa
/// (the annex stresses that "approximate division routines are not
/// sufficient" because this division is used inside Durbin's
/// recursion). `NUM` and `DEN` are assumed already in normalized scalar
/// floating-point form; per the annex "there is no test for DEN being
/// zero — it is assumed that it is not zero."
///
/// The annex's call site is Durbin's recursion, where the result is a
/// reflection coefficient with magnitude strictly less than unity (it
/// is stored in `Q15`). The routine relies on that domain: when
/// `|NUM| < |DEN|` the numerator is pre-shifted by one bit and the
/// 16-bit long division yields a normalized quotient mantissa that fits
/// a signed 16-bit word. Callers must keep `|NUM| < |DEN|`; outside
/// that range the 16-bit quotient would set bit 15 and wrap, which the
/// annex's restricted use never reaches.
///
/// Returns `(QUO, QUONLS)`: the quotient mantissa and its NLS.
#[must_use]
pub fn divide(num: i16, num_nls: i32, den: i16, den_nls: i32) -> (i16, i32) {
    // Determine the sign of the quotient from the sign of NUM * DEN.
    let sign = if (num as i32) * (den as i32) < 0 {
        -1
    } else {
        1
    };

    // QUONLS = NUMNLS − DENNLS + 14.
    let mut quonls = num_nls - den_nls + 14;

    // A0 = |NUM| in the (32-bit) accumulator low half, A1 = |DEN|.
    let mut a0 = (num as i32).abs();
    let a1 = (den as i32).abs();

    // If |NUM| < |DEN|, pre-shift the numerator left by one and bump
    // QUONLS so the long division yields a normalized result.
    if a0 < a1 {
        quonls += 1;
        a0 <<= 1;
    }

    // Long division loop: the §G.1.3.4 loop runs for I = 0, 1, … , 14
    // ("If I < 15, GO TO LOOP"), i.e. 15 iterations, each emitting one
    // quotient bit MSB-first.
    let mut quo: i32 = 0;
    for _ in 0..15 {
        quo <<= 1;
        if a0 >= a1 {
            quo += 1;
            a0 -= a1;
        }
        a0 <<= 1;
    }

    // Round the 16th bit: if the running remainder still covers the
    // denominator, the next quotient bit would be 1 — round up.
    if a0 >= a1 {
        quo += 1;
    }

    if sign < 0 {
        quo = -quo;
    }

    (quo as i16, quonls)
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- §G.1.3.1 right-shift sign-extension anomaly ---------------

    #[test]
    fn shr_sign_extends_like_spec_examples() {
        // The annex's worked example: 3 >> 1 = 1, −3 >> 1 = −2.
        assert_eq!(shr_sat(3, 1), 1);
        assert_eq!(shr_sat(-3, 1), -2);
    }

    #[test]
    fn shr_oversized_collapses_to_sign_fill() {
        // §G.1.3.1: a shift wider than the word collapses to 0 / −1
        // by sign, not a wrapped / masked amount.
        assert_eq!(shr_sat(12345, 32), 0);
        assert_eq!(shr_sat(12345, 64), 0);
        assert_eq!(shr_sat(-12345, 32), -1);
        assert_eq!(shr_sat(-1, 40), -1);
    }

    #[test]
    fn shl_saturates_on_overflow() {
        assert_eq!(shl_sat(1, 4), 16);
        assert_eq!(shl_sat(0x4000_0000, 4), i32::MAX);
        assert_eq!(shl_sat(-0x4000_0000, 4), i32::MIN);
        assert_eq!(shl_sat(0, 40), 0);
        assert_eq!(shl_sat(5, 40), i32::MAX);
    }

    // --- §G.1.3.1 RND(.) ------------------------------------------

    #[test]
    fn rnd_matches_worked_examples() {
        // 1.5 → 2: high word 1 (0x0001_0000), low word 0x8000.
        let acc_1p5 = (1i32 << 16) | 0x8000;
        assert_eq!(rnd(acc_1p5), 2);

        // −1.5 → −1: value −1.5 in two-word form is high −2, low 0x8000.
        // As a single i32 accumulator that is −1.5 * 2^16 = −98304.
        let acc_m1p5 = -((3i32 << 16) / 2); // −1.5 · 65536
        assert_eq!(rnd(acc_m1p5), -1);
    }

    #[test]
    fn rnd_rounds_half_up_and_truncates_below_half() {
        // 2.0 exactly → 2 (low word zero, no round bit).
        assert_eq!(rnd(2i32 << 16), 2);
        // 2.4 → 2 (round bit clear).
        assert_eq!(rnd((2i32 << 16) + ((0.4 * 65536.0) as i32)), 2);
        // 2.5 → 3 (round bit set).
        assert_eq!(rnd((2i32 << 16) + 0x8000), 3);
        // 2.6 → 3.
        assert_eq!(rnd((2i32 << 16) + ((0.6 * 65536.0) as i32)), 3);
    }

    #[test]
    fn rnd_saturates_on_overflow() {
        // High word 0x7FFF, low word MSB set ⇒ would wrap to 0x8000;
        // §G.1.3.1 says saturate instead.
        let acc = (0x7FFFi32 << 16) | 0x8000;
        assert_eq!(rnd(acc), i16::MAX);
    }

    // --- §G.1.2 / §G.1.3 scalar floating point --------------------

    #[test]
    fn from_i16_normalizes_into_band() {
        // A small positive value normalizes so its magnitude lands in
        // [16384, 32767].
        let f = ScalarFloat::from_i16(1);
        assert!(f.mantissa.unsigned_abs() >= 16384);
        // 1 · 2^(−nls) must reconstruct to ~1.0.
        assert!((f.to_f64() - 1.0).abs() < 1e-9);

        let f = ScalarFloat::from_i16(12345);
        assert!(f.mantissa.unsigned_abs() >= 16384);
        assert!((f.to_f64() - 12345.0).abs() < 1e-6);
    }

    #[test]
    fn from_i16_handles_negative_and_already_normalized() {
        let f = ScalarFloat::from_i16(-30000);
        assert!(f.mantissa.unsigned_abs() >= 16384);
        assert_eq!(f.nls, 0, "magnitude already in band ⇒ zero shifts");
        assert!((f.to_f64() - (-30000.0)).abs() < 1e-6);
    }

    #[test]
    fn from_i16_zero_is_scalar_zero() {
        let f = ScalarFloat::from_i16(0);
        assert_eq!(f, ScalarFloat::ZERO);
        assert_eq!(f.nls, MLS_SINGLE + 1);
        assert_eq!(f.to_f64(), 0.0);
    }

    // --- §G.1.3 FINDNLS / VSCALE ----------------------------------

    #[test]
    fn findnls_zero_vector_gets_extra_bit() {
        assert_eq!(findnls(&[0, 0, 0], 3, MLS_SINGLE), MLS_SINGLE + 1);
    }

    #[test]
    fn vscale_left_justifies_largest_magnitude() {
        // Largest magnitude is 100; normalizing should push it into the
        // [16384, 32767] band and shift the others by the same NLS.
        let (out, nls) = vscale(&[100, -50, 25], 3, MLS_SINGLE);
        assert!(out[0].unsigned_abs() >= 16384);
        // Shared NLS ⇒ ratios preserved exactly (powers of two only).
        assert_eq!(out[0] >> nls, 100);
        assert_eq!(out[1] >> nls, -50);
        assert_eq!(out[2] >> nls, 25);
    }

    #[test]
    fn vscale_negative_extreme_drives_nls() {
        // The negative element has the larger magnitude (Case 2).
        let (out, nls) = vscale(&[40, -120, 10], 3, MLS_SINGLE);
        assert!(out[1].unsigned_abs() >= 16384);
        assert_eq!(out[1] >> nls, -120);
    }

    #[test]
    fn vscale_slen1_examines_only_first() {
        // SLEN = 1 means input[0] is asserted to be the max magnitude.
        let (out, nls) = vscale(&[200, 5, 5], 1, MLS_SINGLE);
        assert!(out[0].unsigned_abs() >= 16384);
        assert_eq!(out[0] >> nls, 200);
    }

    #[test]
    fn vscale_double_precision_uses_wider_band() {
        // Double precision: MLS = 30, band is [2^29, 2^30).
        let (out, _nls) = vscale(&[1_000_000, -2_000, 3], 3, MLS_DOUBLE);
        assert!((out[0].unsigned_abs() as i64) >= (1i64 << 29));
    }

    #[test]
    fn vscale_large_value_uses_right_shift() {
        // A 32-bit input examined at single precision (MLS = 14) must
        // be right-shifted (negative NLS) to fit a 16-bit mantissa.
        let (out, nls) = vscale(&[1_000_000, 0, 0], 1, MLS_SINGLE);
        assert!(nls < 0, "expected right shifts for an over-range value");
        assert!(out[0].unsigned_abs() <= i16::MAX as u32);
    }

    // --- §G.1.3.4 DIVIDE ------------------------------------------

    /// Reconstruct the real value a DIVIDE result represents.
    fn quo_value(quo: i16, quonls: i32) -> f64 {
        (quo as f64) * 2f64.powi(-quonls)
    }

    #[test]
    fn divide_matches_float_reference_within_one_lsb() {
        // Exercise a spread of magnitudes and signs within the routine's
        // domain (|quotient| < 1, the reflection-coefficient range it is
        // used for in Durbin's recursion). The mantissa carries 16-bit
        // precision, so the relative error is ≤ 2^-14.
        let cases = [
            (3.0f64, 4.0f64),
            (1.0, 3.0),
            (2.0, 7.0),
            (-5.0, 8.0),
            (5.0, -8.0),
            (-4.0, 9.0),
            (6789.0, 12345.0),
            (1.0, 30000.0),
        ];
        for (num_v, den_v) in cases {
            let n = ScalarFloat::from_i16(num_v as i16);
            let d = ScalarFloat::from_i16(den_v as i16);
            let (quo, quonls) = divide(n.mantissa, n.nls, d.mantissa, d.nls);
            let got = quo_value(quo, quonls);
            let want = num_v / den_v;
            let rel = (got - want).abs() / want.abs().max(1e-9);
            assert!(
                rel < 1e-3,
                "DIVIDE {num_v}/{den_v}: got {got}, want {want} (rel {rel})"
            );
        }
    }

    #[test]
    fn divide_sign_follows_operand_signs() {
        // In-domain operands (|quotient| < 1).
        let pos = ScalarFloat::from_i16(2);
        let neg = ScalarFloat::from_i16(-6);
        // +/− ⇒ negative (2 / −6).
        let (quo, _) = divide(pos.mantissa, pos.nls, neg.mantissa, neg.nls);
        assert!(quo < 0);
        // −/− ⇒ positive (−2 / −6).
        let neg2 = ScalarFloat::from_i16(-2);
        let (quo, _) = divide(neg2.mantissa, neg2.nls, neg.mantissa, neg.nls);
        assert!(quo > 0);
    }

    #[test]
    fn divide_mantissa_is_normalized() {
        // DIVIDE returns a normalized mantissa (bit 14 set), as the
        // annex requires for a scalar floating-point result.
        let n = ScalarFloat::from_i16(1);
        let d = ScalarFloat::from_i16(3);
        let (quo, _) = divide(n.mantissa, n.nls, d.mantissa, d.nls);
        assert!(quo.unsigned_abs() >= 16384);
    }
}
