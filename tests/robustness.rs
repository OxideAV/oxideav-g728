//! Robustness / hardening gates for the public G.728 entry points.
//!
//! These tests carry no ITU reference material — they assert only
//! *self-consistency invariants* that must hold for **any** input:
//!
//! * the encoder never emits a channel index outside the 10-bit range
//!   the §5.11 packing layout can carry (`ICHAN < 1024`);
//! * neither decoder path ever panics, produces `NaN`/`Inf`, or lets the
//!   backward-adapted state run away, on the full index alphabet or on
//!   long adversarial streams;
//! * a fixed-point encode → decode round-trip over a long pseudo-random
//!   PCM stream stays bounded.
//!
//! They complement `tests/conformance.rs` (which pins bit-exactness
//! against the private ITU vectors) by guarding the code paths that the
//! vectors never exercise — every codeword value and a range of
//! degenerate inputs — so a future refactor that introduces a
//! divide-by-zero, an unchecked shift, or an unbounded feedback loop
//! fails here even when the private corpus is absent.

use oxideav_g728::{Decoder, DecoderFixed, Encoder, EncoderFixed, DEFAULT_MAX};

const IDIM: usize = 5;
/// Runaway detector for the fixed-point Q2 decoder output. Q2 speech is
/// bounded by roughly `±4095·2^3`; anything past `1 << 20` is a runaway.
const Q2_RUNAWAY: i32 = 1 << 20;

/// A tiny self-contained LCG so the "adversarial" streams are
/// deterministic and dependency-free.
struct Lcg(u32);
impl Lcg {
    fn new(seed: u32) -> Self {
        Self(seed)
    }
    fn next_u32(&mut self) -> u32 {
        self.0 = self.0.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        self.0
    }
    /// A 10-bit channel index.
    fn next_index(&mut self) -> u16 {
        (self.next_u32() >> 16) as u16 & 0x03FF
    }
    /// A pseudo-random i16 PCM sample across the full range.
    fn next_i16(&mut self) -> i16 {
        (self.next_u32() >> 16) as i16
    }
}

#[test]
fn fixed_encoder_never_emits_out_of_range_index() {
    // Adversarial Q2 inputs: extremes, alternating rails, DC steps and a
    // long pseudo-random tail. Every returned channel index must fit the
    // 10-bit §5.11 field regardless of input.
    let mut enc = EncoderFixed::new();
    let rails: [[i16; IDIM]; 6] = [
        [0; IDIM],
        [i16::MAX; IDIM],
        [i16::MIN; IDIM],
        [i16::MAX, i16::MIN, i16::MAX, i16::MIN, i16::MAX],
        [4095, -4096, 4095, -4096, 4095], // Q2 range rails
        [1, -1, 1, -1, 1],
    ];
    for r in &rails {
        for _ in 0..64 {
            let ichan = enc.encode_vector(r);
            assert!(
                ichan < 1024,
                "index {ichan} exceeds 10-bit range on rail {r:?}"
            );
        }
    }
    let mut rng = Lcg::new(0xC0FF_EE01);
    for _ in 0..4000 {
        let v: [i16; IDIM] = std::array::from_fn(|_| rng.next_i16() >> 1);
        let ichan = enc.encode_vector(&v);
        assert!(ichan < 1024, "index {ichan} exceeds 10-bit range");
    }
}

#[test]
fn fixed_decoder_full_alphabet_and_streams_do_not_run_away() {
    // Every one of the 1024 codewords in sequence, then reversed, then a
    // long pseudo-random stream. The decoder must never panic and its Q2
    // output + block-floating NLS must stay bounded.
    let mut dec = DecoderFixed::new();
    let forward = (0u16..1024).chain((0u16..1024).rev());
    let mut rng = Lcg::new(0x1357_9BDF);
    let random = std::iter::from_fn(|| Some(rng.next_index())).take(8000);
    for ichan in forward.chain(random) {
        let out = dec.decode_vector(ichan);
        for &v in &out {
            assert!(
                v.abs() < Q2_RUNAWAY,
                "fixed decoder Q2 runaway: {v} on index {ichan}"
            );
        }
        let (_st, nls) = dec.last_st();
        assert!(
            (-8..=40).contains(&nls),
            "block-floating NLS out of band: {nls}"
        );
    }
}

#[test]
fn float_decoder_outputs_stay_finite_and_bounded() {
    // Both float paths (postfilter off and on) over the full alphabet and
    // a long adversarial stream: every sample must be finite and within
    // the saturation envelope — this catches NaN/Inf leaking out of the
    // gain / pitch / postfilter adapters.
    for postfiltered in [false, true] {
        let mut dec = Decoder::new();
        // The raw synthesis path (block 32) is hard-clipped at
        // `DEFAULT_MAX`; the postfilter AGC (block 77) applies an output
        // gain that can legitimately overshoot that envelope, so the
        // postfiltered path is only checked for finiteness + no runaway.
        let bound = if postfiltered {
            4.0 * DEFAULT_MAX
        } else {
            DEFAULT_MAX
        };
        let full = 0u16..1024;
        let mut rng = Lcg::new(0x2468_ACE0);
        let random = std::iter::from_fn(|| Some(rng.next_index())).take(6000);
        for ichan in full.chain(random) {
            let out = if postfiltered {
                dec.decode_vector_postfiltered(ichan)
            } else {
                dec.decode_vector(ichan)
            };
            for &v in &out {
                assert!(
                    v.is_finite(),
                    "non-finite sample {v} (postfiltered={postfiltered})"
                );
                assert!(
                    v.abs() <= bound,
                    "sample {v} exceeds bound {bound} (postfiltered={postfiltered})"
                );
            }
        }
    }
}

#[test]
fn fixed_roundtrip_over_long_random_stream_is_stable() {
    // Fixed-point encode → decode over a long random PCM stream. The
    // reconstruction must stay bounded (LD-CELP is not lossless, so we
    // assert stability, not fidelity).
    let mut enc = EncoderFixed::new();
    let mut dec = DecoderFixed::new();
    let mut rng = Lcg::new(0x0BAD_F00D);
    for _ in 0..5000 {
        // Encoder input is Q2 (16-bit-linear >> 1 per §G.3.1).
        let s: [i16; IDIM] = std::array::from_fn(|_| rng.next_i16() >> 1);
        let ichan = enc.encode_vector(&s);
        assert!(ichan < 1024);
        let out = dec.decode_vector(ichan);
        for &v in &out {
            assert!(v.abs() < Q2_RUNAWAY, "roundtrip Q2 runaway: {v}");
        }
    }
}

#[test]
fn float_roundtrip_over_long_random_stream_is_stable() {
    // Float encode → decode over a long random PCM stream. Every stage
    // (encoder analysis-by-synthesis + decoder chain) must stay finite
    // and bounded.
    let mut enc = Encoder::new();
    let mut dec = Decoder::new();
    let mut rng = Lcg::new(0xFEED_BEEF);
    for _ in 0..5000 {
        // Float codec sample convention: (i16 >> 1) / 4.0 (see the
        // conformance harness). Feed the full amplitude range.
        let s: [f64; IDIM] = std::array::from_fn(|_| f64::from(rng.next_i16() >> 1) / 4.0);
        let ichan = enc
            .encode_vector(&s)
            .expect("encode never fails on valid vectors");
        assert!(ichan < 1024);
        let out = dec.decode_vector(ichan);
        for &v in &out {
            assert!(v.is_finite(), "non-finite reconstruction {v}");
            assert!(
                v.abs() <= DEFAULT_MAX,
                "reconstruction {v} exceeds envelope"
            );
        }
    }
}
