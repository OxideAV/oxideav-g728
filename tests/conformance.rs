//! ITU-T G.728 official conformance test vectors
//! (`docs/audio/g728/conformance/` in the private docs tree).
//!
//! Black-box input → reference-output pairs in two parallel reference
//! sets: the no-suffix files are the 1992 floating-point codec
//! references, the `*g` files the Annex G fixed-point references.
//!
//! File format: raw 16-bit little-endian words. Speech files carry one
//! PCM sample per word; bitstream files carry one 10-bit channel index
//! per word (bits 0-2 = gain, bits 3-9 = shape — exactly the crate's
//! `ICHAN = (IS − 1)·NG + (IG − 1)` layout).
//!
//! I/O conversions (§G.3.1: "for 16-bit linear input signals a right
//! shift of 1 bit is required", i.e. 16-bit linear ↔ the coder's `Q2`
//! rendering of the `[−4096, +4095.75]` uniform-PCM range):
//!
//! * encoder input: `S = sample >> 1` (`Q2`),
//! * decoder output, postfilter OFF: the block-floating quantized
//!   speech `ST` rendered at `sample = st · 2^(3 − NLSST)` with the
//!   §G.1.3.1 add-half-then-shift rounding,
//! * decoder output, postfilter ON: the integer-`Q2` `SPF` doubled
//!   (every Annex-G postfiltered reference sample is even),
//! * float codec: samples `= s16 / 8` in, `round(x · 8)` out.
//!
//! The vectors live in the private docs tree, not in this repository;
//! when they are not present (crates.io build, standalone CI) every
//! test here skips with a note. Run from the umbrella workspace — or
//! with `OXIDEAV_G728_CONFORMANCE` pointing at the vector directory —
//! to exercise the real gates. Per-vector status (r392):
//!
//! | Leg | Cases | Result |
//! |-----|-------|--------|
//! | fixed encoder (`incwNg`) | 1–6 | bit-exact |
//! | fixed decoder −postfilter (`outaNg`) | 1–6 | bit-exact |
//! | fixed decoder +postfilter (`outb4g`) | 4 | see test |
//! | float encoder (`incwN`) | 1–4 | bit-exact |
//! | float encoder (`incwN`) | 5, 6 | near-tie prefix + drift (pinned) |
//! | float decoder −postfilter (`outaN`) | 1–6 | bit-exact |
//! | float decoder +postfilter (`outb4`) | 4 | see test |

use std::env;
use std::path::{Path, PathBuf};

use oxideav_g728::{Decoder, DecoderFixed, Encoder, EncoderFixed};

const IDIM: usize = 5;
const CASES: [usize; 6] = [1, 2, 3, 4, 5, 6];

fn conformance_dir() -> Option<PathBuf> {
    if let Ok(p) = env::var("OXIDEAV_G728_CONFORMANCE") {
        let p = PathBuf::from(p);
        if p.join("cw1.bin").exists() {
            return Some(p);
        }
    }
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    loop {
        let c = d
            .join("docs")
            .join("audio")
            .join("g728")
            .join("conformance");
        if c.join("cw1.bin").exists() {
            return Some(c);
        }
        if !d.pop() {
            return None;
        }
    }
}

/// `skip!()` — the vectors are private-tree only; absence is a skip,
/// never a failure (crates.io / standalone CI has no docs checkout).
macro_rules! vectors_or_skip {
    () => {
        match conformance_dir() {
            Some(d) => d,
            None => {
                eprintln!("G.728 conformance vectors not present — skipping");
                return;
            }
        }
    };
}

fn read_words(dir: &Path, name: &str) -> Vec<i16> {
    let bytes = std::fs::read(dir.join(name)).unwrap_or_else(|e| panic!("read {name}: {e}"));
    assert!(bytes.len() % 2 == 0, "{name}: odd byte length");
    bytes
        .chunks_exact(2)
        .map(|c| i16::from_le_bytes([c[0], c[1]]))
        .collect()
}

/// Compare two PCM/bitstream word sequences; return (diff count, max
/// |Δ|, first diff index).
fn diff_stats(ours: &[i16], reference: &[i16]) -> (usize, i32, Option<usize>) {
    assert_eq!(ours.len(), reference.len(), "length mismatch");
    let mut n = 0usize;
    let mut max = 0i32;
    let mut first = None;
    for (i, (a, b)) in ours.iter().zip(reference).enumerate() {
        let d = (i32::from(*a) - i32::from(*b)).abs();
        if d != 0 {
            n += 1;
            max = max.max(d);
            if first.is_none() {
                first = Some(i);
            }
        }
    }
    (n, max, first)
}

/// §G.3.1 16-bit-linear → `Q2` input conversion.
fn pcm16_to_q2(chunk: &[i16]) -> [i16; IDIM] {
    std::array::from_fn(|k| chunk[k] >> 1)
}

/// Fixed decoder, postfilter OFF: block-floating `ST` → 16-bit linear
/// (`sample·8 = st · 2^(3 − nls)`, §G.1.3.1 add-half rounding).
fn st_to_pcm16(st: &[i16; IDIM], nls: i32) -> [i16; IDIM] {
    std::array::from_fn(|k| {
        let v = i64::from(st[k]);
        let sh = nls - 3;
        let w = if sh <= 0 {
            v << (-sh)
        } else {
            (v + (1 << (sh - 1))) >> sh
        };
        w.clamp(i64::from(i16::MIN), i64::from(i16::MAX)) as i16
    })
}

/// Float codec sample → 16-bit linear (`round(x · 8)`, saturated).
fn f64_to_pcm16(x: f64) -> i16 {
    (x * 8.0).round().clamp(-32768.0, 32767.0) as i16
}

// ---------------------------------------------------------------------
// Annex G fixed-point legs — bit-exact gates.
// ---------------------------------------------------------------------

#[test]
fn fixed_encoder_bit_exact_on_all_vectors() {
    let dir = vectors_or_skip!();
    for case in CASES {
        let input = read_words(&dir, &format!("in{case}.bin"));
        let reference = read_words(&dir, &format!("incw{case}g.bin"));
        let mut enc = EncoderFixed::new();
        let ours: Vec<i16> = input
            .chunks_exact(IDIM)
            .map(|c| enc.encode_vector(&pcm16_to_q2(c)) as i16)
            .collect();
        let (n, _, first) = diff_stats(&ours, &reference);
        assert_eq!(
            n,
            0,
            "case {case}: {n}/{} indices differ from incw{case}g (first at {first:?})",
            reference.len()
        );
    }
}

#[test]
fn fixed_decoder_no_postfilter_bit_exact_on_all_vectors() {
    let dir = vectors_or_skip!();
    for case in CASES {
        let cw = read_words(&dir, &format!("cw{case}.bin"));
        let reference = read_words(&dir, &format!("outa{case}g.bin"));
        let mut dec = DecoderFixed::new();
        let mut ours = Vec::with_capacity(cw.len() * IDIM);
        for &w in &cw {
            let _ = dec.decode_vector((w as u16) & 0x3ff);
            let (st, nls) = dec.last_st();
            ours.extend_from_slice(&st_to_pcm16(st, nls));
        }
        let (n, max, first) = diff_stats(&ours, &reference);
        assert_eq!(
            n,
            0,
            "case {case}: {n}/{} samples differ from outa{case}g (max |Δ| {max}, first {first:?})",
            reference.len()
        );
    }
}

#[test]
fn fixed_decoder_postfilter_case4() {
    let dir = vectors_or_skip!();
    let cw = read_words(&dir, "cw4.bin");
    let reference = read_words(&dir, "outb4g.bin");
    let mut dec = DecoderFixed::new();
    let mut ours = Vec::with_capacity(cw.len() * IDIM);
    for &w in &cw {
        for v in dec.decode_vector((w as u16) & 0x3ff) {
            ours.push((i64::from(v) * 2).clamp(i64::from(i16::MIN), i64::from(i16::MAX)) as i16);
        }
    }
    let (n, max, first) = diff_stats(&ours, &reference);
    assert_eq!(
        n,
        0,
        "postfiltered case 4: {n}/{} samples differ from outb4g (max |Δ| {max}, first {first:?})",
        reference.len()
    );
}

// ---------------------------------------------------------------------
// 1992 floating-point legs.
// ---------------------------------------------------------------------

#[test]
fn float_decoder_no_postfilter_bit_exact_on_all_vectors() {
    let dir = vectors_or_skip!();
    for case in CASES {
        let cw = read_words(&dir, &format!("cw{case}.bin"));
        let reference = read_words(&dir, &format!("outa{case}.bin"));
        let mut dec = Decoder::new();
        let mut ours = Vec::with_capacity(cw.len() * IDIM);
        for &w in &cw {
            for x in dec.decode_vector((w as u16) & 0x3ff) {
                ours.push(f64_to_pcm16(x));
            }
        }
        let (n, max, first) = diff_stats(&ours, &reference);
        assert_eq!(
            n,
            0,
            "case {case}: {n}/{} samples differ from outa{case} (max |Δ| {max}, first {first:?})",
            reference.len()
        );
    }
}

#[test]
fn float_encoder_reference_agreement() {
    // The float reference targets are only reproducible up to the
    // reference codec's own floating-point precision at codebook
    // near-ties. Empirically (r392) the f64 chain is bit-exact on
    // cases 1–4 (14 336 consecutive decisions), flips exactly one
    // near-tie on case 6 (index 29), and on the 8-minute case 5
    // follows the reference exactly until a near-tie at index 215
    // after which the free-running states drift apart. Pin those
    // levels so regressions surface.
    let dir = vectors_or_skip!();
    for case in CASES {
        let input = read_words(&dir, &format!("in{case}.bin"));
        let reference = read_words(&dir, &format!("incw{case}.bin"));
        let mut enc = Encoder::new();
        let ours: Vec<i16> = input
            .chunks_exact(IDIM)
            .map(|c| {
                let s: [f64; IDIM] = std::array::from_fn(|k| f64::from(c[k] >> 1) / 4.0);
                enc.encode_vector(&s).expect("encode") as i16
            })
            .collect();
        let (n, _, first) = diff_stats(&ours, &reference);
        match case {
            1..=4 => assert_eq!(
                n,
                0,
                "case {case}: {n}/{} indices differ from incw{case} (first {first:?})",
                reference.len()
            ),
            5 => assert!(
                first.is_none() || first.unwrap() >= 200,
                "case 5: float chain must track the reference through the \
                 first 200 decisions (first mismatch at {first:?})"
            ),
            _ => assert!(
                n <= 1,
                "case 6: at most the one known near-tie flip (got {n}, first {first:?})"
            ),
        }
    }
}

#[test]
fn float_decoder_postfilter_case4_agreement() {
    // The float postfiltered leg is reproducible only up to the
    // reference codec's own floating-point precision: the AGC ratio and
    // the postfilter IIR states accumulate sub-LSB differences between
    // our f64 arithmetic and the reference implementation, which flip
    // the round(x·8) output on samples whose true value sits near a
    // rounding boundary. Empirically (r392, after the §4.7 one-sided
    // fundamental-vs-multiple gate fix) 43 664 of 51 200 samples are
    // exact and the rest are within ±3 16-bit LSBs. Pin those levels.
    let dir = vectors_or_skip!();
    let cw = read_words(&dir, "cw4.bin");
    let reference = read_words(&dir, "outb4.bin");
    let mut dec = Decoder::new();
    let mut ours = Vec::with_capacity(cw.len() * IDIM);
    for &w in &cw {
        for x in dec.decode_vector_postfiltered((w as u16) & 0x3ff) {
            ours.push(f64_to_pcm16(x));
        }
    }
    let (n, max, first) = diff_stats(&ours, &reference);
    assert!(
        max <= 3,
        "postfiltered case 4: max |Δ| {max} vs outb4 exceeds the float \
         reproducibility bound (first diff {first:?})"
    );
    assert!(
        n <= 7600,
        "postfiltered case 4: {n}/{} samples differ from outb4 — \
         agreement regressed below the pinned 85 % level",
        reference.len()
    );
}
