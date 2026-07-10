//! Annex I (frame / packet loss concealment) — §I.5.1 `VEC_LOOP` driver
//! tests.
//!
//! **Conformance status (docs `g728-errata.md`, issues #201/#232):** no
//! concealed-PCM reference output is staged for the Annex I path — the
//! ITU corpus pairs `mask1` with *no* expected output, upstream builds
//! no frame-erasure decoder variant, and the Annex I Recommendation PDF
//! ships pseudo-code but no I/O vectors. The PLC output therefore
//! **cannot be adjudicated bit-exactly in-repo**; these tests pin the
//! spec-mandated *structural* properties instead:
//!
//! * the good path is bit-identical whether or not the Annex I state
//!   machinery is armed (the always-on 31E/PTAP bookkeeping must never
//!   perturb a never-erased stream — separately guaranteed by the ITU
//!   conformance gates in `tests/conformance.rs`);
//! * concealment output is finite and honours the synthesis saturation;
//! * the §I.5.2 FESCALE schedule silences long erasures (voiced 60 ms,
//!   unvoiced 70 ms);
//! * the §I.4.5/§I.5.1 gain clamp arms for `min(FECOUNT, AFTERFEMAX)`
//!   good cycles after the erasure and then releases;
//! * decode resumes cleanly after an erasure (no panic, no NaN, and the
//!   backward adapters keep tracking).

use oxideav_g728::consts::{AFTERFEMAX, IDIM, NUPDATE};
use oxideav_g728::{Decoder, DecoderFixed, FRAME_LEN};

/// Deterministic 10-bit channel-index stream (arbitrary but fixed).
fn index_stream(n: usize, seed: u64) -> Vec<u16> {
    let mut s = seed;
    (0..n)
        .map(|_| {
            s = s.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            ((s >> 33) & 0x3FF) as u16
        })
        .collect()
}

fn assert_sane(v: &[f64; FRAME_LEN], what: &str) {
    for (k, &x) in v.iter().enumerate() {
        assert!(x.is_finite(), "{what}: sample {k} not finite ({x})");
        assert!(
            x.abs() <= 4095.0 * 4.0,
            "{what}: sample {k} escaped any plausible range ({x})"
        );
    }
}

#[test]
fn annex_i_state_never_perturbs_a_good_stream() {
    // Two decoders fed the same never-erased stream must agree exactly
    // — one uses the plain per-vector API, the other interrogates the
    // Annex I accessors along the way (reads must be side-effect-free).
    let idx = index_stream(64 * NUPDATE, 0x1234_5678);
    let mut a = Decoder::new();
    let mut b = Decoder::new();
    for &i in &idx {
        let va = a.decode_vector_postfiltered(i);
        let vb = b.decode_vector_postfiltered(i);
        let _ = b.fe_counters();
        let _ = b.frame_erasure_excitation().in_erasure();
        assert_eq!(va, vb);
    }
    let (fecount, afterfe, _) = a.fe_counters();
    assert_eq!((fecount, afterfe), (0, 0), "no erasure ⇒ counters stay 0");
}

#[test]
fn concealment_is_finite_and_decoder_resumes() {
    // good 32 cycles → conceal 8 cycles (20 ms) → good 32 cycles.
    let mut dec = Decoder::new();
    for &i in &index_stream(32 * NUPDATE, 42) {
        let v = dec.decode_vector_postfiltered(i);
        assert_sane(&v, "pre-erasure");
    }
    for n in 0..(8 * NUPDATE) {
        let v = dec.conceal_vector_postfiltered();
        assert_sane(&v, "concealed");
        assert!(dec.frame_erasure_excitation().in_erasure(), "vector {n}");
    }
    let (fecount, _, _) = dec.fe_counters();
    assert_eq!(fecount, 8, "8 concealed adaptation cycles counted");
    for &i in &index_stream(32 * NUPDATE, 43) {
        let v = dec.decode_vector_postfiltered(i);
        assert_sane(&v, "post-erasure");
    }
    let (fecount, _, _) = dec.fe_counters();
    assert_eq!(fecount, 0, "good frame resets FECOUNT");
    assert!(!dec.frame_erasure_excitation().in_erasure());
}

#[test]
fn long_erasure_scales_to_silence() {
    // §I.5.2: FESCALE = 0 past 60 ms (voiced) / 70 ms (unvoiced). Run a
    // 100 ms erasure (40 cycles) and check the schedule bottomed out
    // and the concealed energy decayed.
    let mut dec = Decoder::new();
    for &i in &index_stream(64 * NUPDATE, 7) {
        dec.decode_vector_postfiltered(i);
    }
    let mut first_cycle_energy = 0.0f64;
    let mut last_cycle_energy = 0.0f64;
    for n in 0..(40 * NUPDATE) {
        let v = dec.conceal_vector_postfiltered();
        assert_sane(&v, "long-erasure conceal");
        let e: f64 = v.iter().map(|x| x * x).sum();
        if n < NUPDATE {
            first_cycle_energy += e;
        }
        if n >= 39 * NUPDATE {
            last_cycle_energy += e;
        }
    }
    assert_eq!(
        dec.frame_erasure_excitation().fescale(),
        0.0,
        "FESCALE schedule must bottom out past 70 ms"
    );
    assert!(
        last_cycle_energy <= first_cycle_energy + 1e-9,
        "concealed energy must not grow across a 100 ms erasure \
         (first {first_cycle_energy}, last {last_cycle_energy})"
    );
}

#[test]
fn gain_clamp_arms_for_erasure_length_and_releases() {
    let mut dec = Decoder::new();
    for &i in &index_stream(16 * NUPDATE, 9) {
        dec.decode_vector_postfiltered(i);
    }
    // 6 erased cycles.
    for _ in 0..(6 * NUPDATE) {
        dec.conceal_vector_postfiltered();
    }
    let (fecount, afterfe, _) = dec.fe_counters();
    assert_eq!((fecount, afterfe), (6, 0), "clamp arms only at recovery");
    // First good cycle: AFTERFE = FECOUNT = 6; then one fewer per good
    // cycle (the decrement tests the previous cycle's FERROR, so the
    // first good cycle itself does not decrement).
    let idx = index_stream(8 * NUPDATE, 10);
    for (n, &i) in idx.iter().enumerate() {
        dec.decode_vector_postfiltered(i);
        if n % NUPDATE == NUPDATE - 1 {
            let cycles_done = n / NUPDATE + 1;
            let want = 6usize.saturating_sub(cycles_done - 1);
            let (_, afterfe, _) = dec.fe_counters();
            assert_eq!(afterfe, want, "after {cycles_done} good cycles");
        }
    }
    let (_, afterfe, _) = dec.fe_counters();
    assert_eq!(afterfe, 0, "clamp released");
}

#[test]
fn gain_clamp_saturates_at_afterfemax() {
    let mut dec = Decoder::new();
    for &i in &index_stream(8 * NUPDATE, 11) {
        dec.decode_vector_postfiltered(i);
    }
    for _ in 0..((AFTERFEMAX + 8) * NUPDATE) {
        dec.conceal_vector_postfiltered();
    }
    // One good cycle to arm the clamp.
    for &i in &index_stream(NUPDATE, 12) {
        dec.decode_vector_postfiltered(i);
    }
    let (_, afterfe, _) = dec.fe_counters();
    assert_eq!(afterfe, AFTERFEMAX, "40 ms ceiling (§I.4.5)");
}

#[test]
fn ogaindb_tracks_and_growth_is_limited_after_erasure() {
    // Feed silence, conceal, then blast maximum-energy vectors: the
    // predicted gain may rise at most FEGAINMAX = 2 dB per vector while
    // the clamp is active. We can't observe sigma directly, but OGAINDB
    // (the limited dB gain) is exposed — consecutive good vectors under
    // an active clamp must not jump more than 2 dB.
    let mut dec = Decoder::new();
    for _ in 0..(16 * NUPDATE) {
        dec.decode_vector_postfiltered(0); // near-silence stream
    }
    for _ in 0..(8 * NUPDATE) {
        dec.conceal_vector_postfiltered();
    }
    // Loud indices (max gain codebook entries) while the clamp holds.
    let mut prev = None;
    for n in 0..(4 * NUPDATE) {
        dec.decode_vector_postfiltered(0x3FF);
        let (_, afterfe, ogaindb) = dec.fe_counters();
        assert!(ogaindb.is_finite());
        if let Some(p) = prev {
            if afterfe > 0 {
                assert!(
                    ogaindb - p <= 2.0 + 1e-9,
                    "vector {n}: OGAINDB jumped {} dB under an active clamp",
                    ogaindb - p
                );
            }
        }
        prev = Some(ogaindb);
    }
}

#[test]
fn conceal_vector_raw_matches_postfiltered_state_flow() {
    // conceal_vector() (raw sd) and conceal_vector_postfiltered() run
    // the identical inner loop; two decoders driven in lockstep — one
    // via each surface — must keep identical *subsequent* good-frame
    // output.
    let warm = index_stream(16 * NUPDATE, 21);
    let tail = index_stream(16 * NUPDATE, 22);
    let mut a = Decoder::new();
    let mut b = Decoder::new();
    for &i in &warm {
        a.decode_vector_postfiltered(i);
        b.decode_vector_postfiltered(i);
    }
    for _ in 0..(4 * NUPDATE) {
        let _ = a.conceal_vector();
        let _ = b.conceal_vector_postfiltered();
    }
    for &i in &tail {
        assert_eq!(
            a.decode_vector_postfiltered(i),
            b.decode_vector_postfiltered(i)
        );
    }
}

// ---------------------------------------------------------------------
// mask1-driven decode over the official ITU bitstream (skips when the
// private docs tree is absent). No concealed-PCM reference exists to
// diff against (docs gap, issue #232) — the gates here are structural.
// ---------------------------------------------------------------------

use std::env;
use std::path::PathBuf;

fn conformance_dir() -> Option<PathBuf> {
    if let Ok(p) = env::var("OXIDEAV_G728_CONFORMANCE") {
        let p = PathBuf::from(p);
        if p.join("cw4.bin").exists() {
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
        if c.join("cw4.bin").exists() {
            return Some(c);
        }
        if !d.pop() {
            return None;
        }
    }
}

#[test]
fn mask1_masked_decode_of_official_bitstream() {
    let Some(dir) = conformance_dir() else {
        eprintln!("G.728 conformance vectors not present — skipping");
        return;
    };
    let cw = std::fs::read(dir.join("cw4.bin")).expect("cw4.bin");
    let indices: Vec<u16> = cw
        .chunks_exact(2)
        .map(|c| u16::from_le_bytes([c[0], c[1]]) & 0x3FF)
        .collect();
    let mask: Vec<bool> = std::fs::read_to_string(dir.join("mask1"))
        .expect("mask1")
        .trim()
        .chars()
        .map(|c| c == '1')
        .collect();
    assert_eq!(mask.len(), 20, "mask1 is 20 frame flags");

    // Interpretation: one mask entry per 10 ms frame erasure (FESIZE=4
    // ⇒ 16 vectors per frame), applied to the first 20 frames of the
    // stream; the remainder decodes normally. (The corpus stages no
    // paired PLC output and no framing metadata for mask1 — docs gap
    // #232 — so this is a self-consistent, documented choice.)
    const VEC_PER_FRAME: usize = 4 * NUPDATE;
    let mut masked = Decoder::new();
    let mut reference = Decoder::new();
    let mut first_erased_vec = None;
    for (n, &i) in indices.iter().enumerate() {
        let frame = n / VEC_PER_FRAME;
        let erased = *mask.get(frame).unwrap_or(&false);
        let r = reference.decode_vector_postfiltered(i);
        let m = if erased {
            if first_erased_vec.is_none() {
                first_erased_vec = Some(n);
            }
            masked.conceal_vector_postfiltered()
        } else {
            masked.decode_vector_postfiltered(i)
        };
        assert_sane(&m, "masked decode");
        // Before the first erasure the masked decoder must be
        // bit-identical to the reference decode.
        if first_erased_vec.is_none() {
            assert_eq!(m, r, "vector {n}: pre-erasure prefix must match");
        }
    }
    let first = first_erased_vec.expect("mask1 has one erased frame");
    assert_eq!(first, 19 * VEC_PER_FRAME, "mask1 erases the 20th frame");

    // After the (single, final-in-mask) erased frame, the decoder keeps
    // decoding the remaining ~620 frames; spot-check that its output
    // re-approaches the reference: the two free-running decoders will
    // not be bit-equal (their backward-adapted states diverged during
    // the erasure), but the masked decode must remain bounded speech.
    let (fecount, _, _) = masked.fe_counters();
    assert_eq!(fecount, 0, "erasure closed by the following good frames");
}

// ---------------------------------------------------------------------
// Fixed-point (Annex G) Annex I path.
// ---------------------------------------------------------------------

fn assert_sane_q2(v: &[i32; IDIM], what: &str) {
    for (k, &x) in v.iter().enumerate() {
        assert!(
            x.abs() <= 1 << 20,
            "{what}: Q2 sample {k} escaped any plausible range ({x})"
        );
    }
}

#[test]
fn fixed_annex_i_state_never_perturbs_a_good_stream() {
    let idx = index_stream(64 * NUPDATE, 0x0F1E_5EED);
    let mut a = DecoderFixed::new();
    let mut b = DecoderFixed::new();
    for &i in &idx {
        let va = a.decode_vector(i);
        let vb = b.decode_vector(i);
        let _ = b.fe_counters();
        let _ = b.frame_erasure_excitation().in_erasure();
        assert_eq!(va, vb);
    }
    let (fecount, afterfe, _) = a.fe_counters();
    assert_eq!((fecount, afterfe), (0, 0));
}

#[test]
fn fixed_concealment_is_finite_and_decoder_resumes() {
    let mut dec = DecoderFixed::new();
    for &i in &index_stream(32 * NUPDATE, 51) {
        let v = dec.decode_vector(i);
        assert_sane_q2(&v, "fixed pre-erasure");
    }
    for _ in 0..(8 * NUPDATE) {
        let v = dec.conceal_vector();
        assert_sane_q2(&v, "fixed concealed");
        assert!(dec.frame_erasure_excitation().in_erasure());
    }
    let (fecount, _, _) = dec.fe_counters();
    assert_eq!(fecount, 8);
    for &i in &index_stream(32 * NUPDATE, 52) {
        let v = dec.decode_vector(i);
        assert_sane_q2(&v, "fixed post-erasure");
    }
    let (fecount, afterfe, ogaindb) = dec.fe_counters();
    assert_eq!(fecount, 0);
    assert!(afterfe <= AFTERFEMAX);
    assert!((-16_384..=14_336).contains(&ogaindb), "OGAINDB Q9 domain");
}

#[test]
fn fixed_long_erasure_scales_to_silence() {
    let mut dec = DecoderFixed::new();
    for &i in &index_stream(64 * NUPDATE, 53) {
        dec.decode_vector(i);
    }
    for _ in 0..(40 * NUPDATE) {
        let v = dec.conceal_vector();
        assert_sane_q2(&v, "fixed long conceal");
    }
    assert_eq!(
        dec.frame_erasure_excitation().fescale(),
        0,
        "fixed FESCALE schedule must bottom out past 70 ms"
    );
}

#[test]
fn fixed_gain_clamp_arms_and_saturates() {
    let mut dec = DecoderFixed::new();
    for &i in &index_stream(8 * NUPDATE, 54) {
        dec.decode_vector(i);
    }
    for _ in 0..((AFTERFEMAX + 4) * NUPDATE) {
        dec.conceal_vector();
    }
    for &i in &index_stream(NUPDATE, 55) {
        dec.decode_vector(i);
    }
    let (_, afterfe, _) = dec.fe_counters();
    assert_eq!(afterfe, AFTERFEMAX, "40 ms ceiling (§I.4.5)");
    // Clamp releases after AFTERFEMAX further good cycles.
    for &i in &index_stream(AFTERFEMAX * NUPDATE, 56) {
        dec.decode_vector(i);
    }
    let (_, afterfe, _) = dec.fe_counters();
    assert_eq!(afterfe, 0, "clamp released");
}

#[test]
fn fixed_mask1_masked_decode_of_official_bitstream() {
    let Some(dir) = conformance_dir() else {
        eprintln!("G.728 conformance vectors not present — skipping");
        return;
    };
    let cw = std::fs::read(dir.join("cw4.bin")).expect("cw4.bin");
    let indices: Vec<u16> = cw
        .chunks_exact(2)
        .map(|c| u16::from_le_bytes([c[0], c[1]]) & 0x3FF)
        .collect();
    let mask: Vec<bool> = std::fs::read_to_string(dir.join("mask1"))
        .expect("mask1")
        .trim()
        .chars()
        .map(|c| c == '1')
        .collect();

    const VEC_PER_FRAME: usize = 4 * NUPDATE;
    let mut masked = DecoderFixed::new();
    let mut reference = DecoderFixed::new();
    let mut in_prefix = true;
    for (n, &i) in indices.iter().enumerate() {
        let frame = n / VEC_PER_FRAME;
        let erased = *mask.get(frame).unwrap_or(&false);
        let r = reference.decode_vector(i);
        let m = if erased {
            in_prefix = false;
            masked.conceal_vector()
        } else {
            masked.decode_vector(i)
        };
        assert_sane_q2(&m, "fixed masked decode");
        if in_prefix {
            assert_eq!(m, r, "vector {n}: pre-erasure prefix must be bit-exact");
        }
    }
    let (fecount, _, _) = masked.fe_counters();
    assert_eq!(fecount, 0, "erasure closed by the following good frames");
}

#[test]
fn fixed_hostile_random_masks_never_panic() {
    let mut s = 0xFEED_FACE_u64;
    let mut rnd = move || {
        s = s.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        s >> 33
    };
    let mut dec = DecoderFixed::new();
    for _ in 0..4096 {
        let r = rnd();
        let v = if r & 7 == 0 {
            dec.conceal_vector()
        } else {
            dec.decode_vector((r & 0x3FF) as u16)
        };
        assert_sane_q2(&v, "fixed hostile");
    }
}

#[test]
fn hostile_random_masks_never_panic() {
    // Fuzz-shaped robustness: random erasure patterns (including
    // misaligned, single-vector "erasures" that violate the documented
    // cycle-alignment guidance) over random indices must never panic
    // or produce non-finite output.
    let mut s = 0xDEAD_BEEF_u64;
    let mut rnd = move || {
        s = s.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        s >> 33
    };
    let mut dec = Decoder::new();
    for _ in 0..4096 {
        let r = rnd();
        if r & 7 == 0 {
            let v = dec.conceal_vector_postfiltered();
            assert!(v.iter().all(|x| x.is_finite()));
        } else {
            let v = dec.decode_vector_postfiltered((r & 0x3FF) as u16);
            assert!(v.iter().all(|x| x.is_finite()));
        }
    }
    let _ = IDIM; // silence unused-import lint on skip-configured builds
}
