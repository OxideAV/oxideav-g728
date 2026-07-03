//! Annex G §G.7 — the fixed-point encoder and decoder main programs.
//!
//! [`EncoderFixed`] / [`DecoderFixed`] chain every landed §G.3
//! fixed-point block in exactly the §G.7 main-program execution order,
//! closing the full bit-exact coder loop:
//!
//! ```text
//! encoder, per vector (ICOUNT = 1..4):
//!   ICOUNT=3: commit block 51 (A) + block 38 (AWZ/AWP)
//!   blocks 46/98/99/48 (σ)  [45 at ICOUNT=2]
//!   blockzir → ZIR ; block 4 → SW
//!   blocks 11/16/12/14/15/13/17/18 → ICHAN ; 19/21 → ET
//!   blocks 9/10 memory update → ST
//!   blocks 93/94/96/97 → GSTATE(1)
//!   STTMP/STMP bookkeeping
//!   ICOUNT=4: blocks 49+50 (synthesis window + Levinson)
//!   ICOUNT=2: blocks 36+37 (weighting window + Levinson)
//!   ICOUNT=1: blocks 43+44 (log-gain window + Levinson)
//!
//! decoder, per vector:
//!   ICOUNT=3: commit block 51 (A)
//!   blocks 46/98/99/48 (σ)  [45 at ICOUNT=2]
//!   blocks 19/21 → ET ; block 32 → ST
//!   ICOUNT=1: block 85 (short-term postfilter coefficients)
//!   block 81 → SST/D ; ICOUNT=3: blocks 82/83/84
//!   blocks 71–77 → SPF
//!   blocks 93/94/96/97 → GSTATE(1) ; STTMP bookkeeping
//!   ICOUNT=4: block 49 + block 50 (order 1–10, save APF byproduct,
//!             then continue 11–50) + block 51 (pending)
//! ```
//!
//! The §G.7 note licenses any block ordering that stays bit-exact; the
//! two visible-commit points the drivers must preserve are block 51 /
//! block 38 becoming *live* at `ICOUNT = 3` and block 45 at
//! `ICOUNT = 2` (handled inside
//! [`crate::annex_g_gain_adapter::GainAdapterFixed`]). The encoder runs
//! block 50 to order 50 in one shot
//! ([`crate::annex_g_synth_adapter::SynthAdapterFixed`]); the decoder
//! interrupts it at order 10 to capture the postfilter predictor `APF`
//! (+ `RC1`), then resumes — §G.2.2 proves the interrupted recursion is
//! bit-exact with the one-shot, so encoder and decoder derive the same
//! `A`.

use crate::annex_g_codebook::{self, SearchResult};
use crate::annex_g_decoder::SynthesisFilterFixed;
use crate::annex_g_encoder::EncoderFiltersFixed;
use crate::annex_g_gain_adapter::GainAdapterFixed;
use crate::annex_g_hybrid::{BflSegment, HybridWindowFixed, HybridWindowFixedState, NLSATT50};
use crate::annex_g_levinson::{
    levinson_durbin_fixed, levinson_durbin_fixed_resume, LevinsonInput, RecursionResume,
};
use crate::annex_g_pf_coeff::short_term_coeff_fixed;
use crate::annex_g_pitch::{apf_to_q13, PitchAdapterFixed};
use crate::annex_g_postfilter::{
    LongTermCoeff, PostfilterFixed, ShortTermCoeff, PF_ORDER, SST_PAST,
};
use crate::annex_g_synth_adapter::{bandwidth_expand_q14, StSegment, SynthAdapterFixed, N5};
use crate::annex_g_weight_adapter::{WeightAdapterFixed, WeightCoeffFixed};
use crate::consts::{IDIM, LPC, LPCW, NFRSZ, NG, NUPDATE};
use crate::tables::FACV_Q14;

/// The full §G.7 fixed-point encoder: input `Q2` speech vectors, output
/// 10-bit channel indices.
#[derive(Debug, Clone)]
pub struct EncoderFixed {
    /// `ICOUNT` (0 before the first vector; 1..=4 inside the loop).
    icount: usize,
    /// §G.3.1 – §G.3.3 filter memories.
    filters: EncoderFiltersFixed,
    /// Blocks 49/50/51 — backward synthesis-filter adapter.
    synth: SynthAdapterFixed,
    /// Blocks 36/37/38 — backward weighting-filter adapter.
    weight: WeightAdapterFixed,
    /// Figure G.1 backward gain adapter.
    gain: GainAdapterFixed,
    /// Live `Q14` synthesis predictor `A` (block 51 visible-commit at
    /// `ICOUNT = 3`).
    a: [i16; LPC + 1],
    /// Live `Q14` weighting pair (block 38 visible-commit at
    /// `ICOUNT = 3`).
    wcoeff: WeightCoeffFixed,
    /// `STTMP` — the cycle's four `ST` segments (blocks 49/50 input).
    sttmp: [StSegment; N5],
    /// `STMP` — the block-36 input-speech window (`Q2`, the annex's
    /// `(ICOUNT − 3)·IDIM` layout: vectors 3, 4, 1, 2).
    stmp: [i16; NFRSZ],
}

impl Default for EncoderFixed {
    fn default() -> Self {
        Self::new()
    }
}

impl EncoderFixed {
    /// Fresh encoder with Table 2 / Table G.2 initial state everywhere.
    #[must_use]
    pub fn new() -> Self {
        let mut a = [0i16; LPC + 1];
        a[0] = 1 << 14;
        Self {
            icount: 0,
            filters: EncoderFiltersFixed::new(),
            synth: SynthAdapterFixed::new(),
            weight: WeightAdapterFixed::new(),
            gain: GainAdapterFixed::new(),
            a,
            wcoeff: WeightCoeffFixed::disabled(),
            sttmp: [StSegment::zero(16); N5],
            stmp: [0; NFRSZ],
        }
    }

    /// Encode one `Q2` input-speech vector (`S(1..=IDIM)`, the §G.3.1
    /// `[−16384, +16383]` range) into its 10-bit channel index.
    pub fn encode_vector(&mut self, s: &[i16; IDIM]) -> u16 {
        // | If ICOUNT = 4, set ICOUNT = 0; ICOUNT = ICOUNT + 1
        self.icount = self.icount % NUPDATE + 1;
        let icount = self.icount;

        // | If ICOUNT = 3: if !ILLCOND do block 51; if !ILLCONDW do
        // | block 38 (the adapters computed both at the end of the
        // | previous cycle; this is the visible commit), then blocks
        // | 12/14/15 (recomputed per call inside encode_vector below).
        if icount == 3 {
            self.a = *self.synth.predictor();
            self.wcoeff = *self.weight.coeff();
        }

        // | If ICOUNT = 2 and !ILLCONDG do block 45; then blocks
        // | 46/98/99/48.
        let (gain, nlsgain) = self.gain.predict(icount);

        // | do "blockzir" ; do block 4.
        let zir =
            self.filters
                .zero_input_response(&self.a, &self.wcoeff.q_gamma1, &self.wcoeff.q_gamma2);
        let sw = self
            .filters
            .weighting_filter(s, &self.wcoeff.q_gamma1, &self.wcoeff.q_gamma2);

        // | blocks 11 / 16 / 13 / 17 / 18 / 19 / 21.
        let a6: [i16; IDIM + 1] = std::array::from_fn(|i| self.a[i]);
        let awz6: [i16; IDIM + 1] = std::array::from_fn(|i| self.wcoeff.q_gamma1[i]);
        let awp6: [i16; IDIM + 1] = std::array::from_fn(|i| self.wcoeff.q_gamma2[i]);
        let out = annex_g_codebook::encode_vector(&sw, &zir, &a6, &awz6, &awp6, gain, nlsgain);

        // | blocks 9 and 10 during filter memory update → ST.
        let up = self.filters.memory_update(
            &out.et,
            out.nlset,
            &self.a,
            &self.wcoeff.q_gamma1,
            &self.wcoeff.q_gamma2,
        );

        // | blocks 93/94/96/97; GSTATE(1) = output of block 97 (and the
        // | ICOUNT = 1 GTMP + blocks 43/44 refresh).
        self.gain.update(icount, out.result.ig, out.result.is);

        // | I = (ICOUNT − 1)·IDIM; STTMP(I+1..I+5) = ST; NLSSTTMP(ICOUNT).
        self.sttmp[icount - 1] = StSegment {
            st: up.st,
            nls: up.nlsst,
        };
        // | I = (ICOUNT − 3)·IDIM; if ICOUNT < 3, I += 20; STMP(I+1..) = S.
        let i0 = if icount < 3 {
            (icount + 1) * IDIM
        } else {
            (icount - 3) * IDIM
        };
        self.stmp[i0..i0 + IDIM].copy_from_slice(s);

        // | If ICOUNT = 4: blocks 49 + 50 (block 51 is computed inside
        // | the adapter; its result becomes live at the next ICOUNT=3).
        if icount == 4 {
            let _ = self.synth.adapt(&self.sttmp);
        }
        // | If ICOUNT = 2: blocks 36 + 37 (+38 computed; live at 3).
        if icount == 2 {
            let stmp = self.stmp;
            let _ = self.weight.adapt(&stmp);
        }

        out.result.ichan
    }

    /// The most recent quantized-speech vector's mantissas + shift (for
    /// tests / decoder-lockstep verification).
    #[must_use]
    pub fn last_st(&self) -> (&[i16; IDIM], i32) {
        let seg = &self.sttmp[self.icount - 1];
        (&seg.st, seg.nls)
    }

    /// Encode a whole buffer of `Q2` speech vectors into the §5.11
    /// serial byte stream (MSB-first 10-bit-per-vector packing via
    /// [`crate::bitstream::pack_indices`]). The natural framing unit is
    /// one adaptation cycle = 4 indices = 40 bits = 5 bytes.
    pub fn encode_buffer(&mut self, vectors: &[[i16; IDIM]]) -> Vec<u8> {
        let indices: Vec<u16> = vectors.iter().map(|v| self.encode_vector(v)).collect();
        crate::bitstream::pack_indices(&indices)
    }
}

/// The full §G.7 fixed-point decoder: input 10-bit channel indices,
/// output postfiltered `Q2` speech vectors.
#[derive(Debug, Clone)]
pub struct DecoderFixed {
    icount: usize,
    /// Block 32 — the fixed-point synthesis filter.
    synth_filter: SynthesisFilterFixed,
    /// Block 49 — the decoder's own hybrid-window state (the block-50
    /// recursion is split at order 10 here, unlike the encoder's
    /// one-shot [`SynthAdapterFixed`]).
    hw49: HybridWindowFixedState,
    /// Figure G.1 backward gain adapter.
    gain: GainAdapterFixed,
    /// Blocks 81 – 84 — the long-term-postfilter adaptation chain.
    pitch: PitchAdapterFixed,
    /// Blocks 71 – 77 — the adaptive postfilter.
    postfilter: PostfilterFixed,
    /// Live `Q14` synthesis predictor `A` (commit at `ICOUNT = 3`).
    a: [i16; LPC + 1],
    /// Block 51 output pending its `ICOUNT = 3` commit.
    pending_a: Option<[i16; LPC + 1]>,
    /// `STTMP` — the cycle's four decoded `ST` segments.
    sttmp: [StSegment; N5],
    /// Order-10 Levinson by-product for block 85 (`ATMP(1..=11)` in the
    /// `nlsapf` format) + `RC1`.
    apf_atil: [i16; PF_ORDER + 1],
    /// `NLSAPF` — the by-product's Q signal.
    nlsapf: i32,
    /// `RC1` (Q15) captured with the by-product.
    rc1: i16,
    /// `APF` in flat `Q13` for block 81 (§G.3.28 LABEL conversion).
    apf_q13: [i16; PF_ORDER + 1],
    /// Live long-term postfilter coefficients (blocks 82 – 84 output).
    ltc: LongTermCoeff,
    /// Live short-term postfilter coefficients (block 85 output).
    stc: ShortTermCoeff,
}

impl Default for DecoderFixed {
    fn default() -> Self {
        Self::new()
    }
}

/// The §G.3.17 synthesis-filter window descriptor (block 49).
fn synth_window() -> HybridWindowFixed<'static> {
    HybridWindowFixed {
        order: LPC,
        l: NFRSZ,
        n: crate::consts::NONR,
        window: &crate::tables::WNR_Q15,
        nlsatt: NLSATT50,
    }
}

impl DecoderFixed {
    /// Fresh decoder with Table 2 / Table G.2 initial state: all-pass
    /// `A`, identity postfilter (`GL = 1`, `B = 0`, `AZ = AP = TILTZ =
    /// 0`), all-pass `APF`.
    #[must_use]
    pub fn new() -> Self {
        let mut a = [0i16; LPC + 1];
        a[0] = 1 << 14;
        let mut apf_atil = [0i16; PF_ORDER + 1];
        apf_atil[0] = i16::MAX; // Q15 unity head (NLSAPF = 15)
        let mut apf_q13 = [0i16; PF_ORDER + 1];
        apf_q13[0] = 1 << 13;
        let mut az_q14 = [0i32; PF_ORDER + 1];
        az_q14[0] = 1 << 14;
        Self {
            icount: 0,
            synth_filter: SynthesisFilterFixed::new(),
            hw49: HybridWindowFixedState::new(&synth_window()),
            gain: GainAdapterFixed::new(),
            pitch: PitchAdapterFixed::new(),
            postfilter: PostfilterFixed::new(),
            a,
            pending_a: None,
            sttmp: [StSegment::zero(16); N5],
            apf_atil,
            nlsapf: 15,
            rc1: 0,
            apf_q13,
            ltc: LongTermCoeff {
                kp: crate::annex_g_pitch::KP1_INIT,
                gl_q14: 1 << 14,
                glb_q16: 0,
            },
            stc: ShortTermCoeff {
                az_q14,
                ap_q14: az_q14,
                tiltz_q14: 0,
            },
        }
    }

    /// Decode one 10-bit channel index into the postfiltered `Q2`
    /// speech vector `SPF(1..=IDIM)`.
    pub fn decode_vector(&mut self, ichan: u16) -> [i32; IDIM] {
        // | If ICOUNT = 4, set ICOUNT = 0; ICOUNT = ICOUNT + 1
        self.icount = self.icount % NUPDATE + 1;
        let icount = self.icount;

        // | Obtain IS and IG from ICHAN (1-based).
        let is = (ichan as usize) / NG + 1;
        let ig = (ichan as usize) % NG + 1;

        // | If ICOUNT = 3: If ILLCOND = .FALSE., do block 51 (commit the
        // | pending expansion from the previous cycle's ICOUNT = 4).
        if icount == 3 {
            if let Some(a) = self.pending_a.take() {
                self.a = a;
            }
        }

        // | do blocks 46, 98, 99, and 48 (45 at ICOUNT = 2 inside).
        let (gain, nlsgain) = self.gain.predict(icount);

        // | do blocks 19 and 21 → ET.
        let result = SearchResult { is, ig, ichan };
        let (et, nlset) = annex_g_codebook::excitation(&result, gain, nlsgain);

        // | do block 32 → ST (BFL, NLSST).
        let et_i32: [i32; IDIM] = std::array::from_fn(|k| i32::from(et[k]));
        let a_i32: [i32; LPC + 1] = std::array::from_fn(|i| i32::from(self.a[i]));
        let st = self.synth_filter.synthesize(&et_i32, nlset, &a_i32);
        let nlsst = self.synth_filter.nlsst();

        // | If ICOUNT = 1, do block 85 (short-term postfilter coeffs
        // | from the order-10 by-product; keep the old set on decline).
        if icount == 1 {
            if let Some(sc) = short_term_coeff_fixed(&self.apf_atil, self.rc1, self.nlsapf) {
                self.stc = sc;
            }
        }

        // | do block 81 → SST (Q2) + residual D.
        let sst = self.pitch.inverse_filter(&st, nlsst, &self.apf_q13);

        // | If ICOUNT = 3: blocks 82 (pitch), 83 (tap), 84 (GL/GLB).
        if icount == 3 {
            self.pitch.pitch_extract();
            let ptap = self.pitch.pitch_tap(self.postfilter.sst(), SST_PAST - 1);
            self.ltc = self.pitch.long_term_coeff(ptap);
        }

        // | blocks 71 – 77 → SPF.
        let spf = self.postfilter.filter_vector(&sst, &self.ltc, &self.stc);

        // | blocks 93/94/96/97 (+ ICOUNT = 1 GTMP + 43/44).
        self.gain.update(icount, ig, is);

        // | STTMP(I+1..I+5) = ST; NLSSTTMP(ICOUNT) = NLSST.
        let st_i16: [i16; IDIM] = std::array::from_fn(|k| st[k] as i16);
        self.sttmp[icount - 1] = StSegment {
            st: st_i16,
            nls: nlsst,
        };

        // | If ICOUNT = 4: block 49, block 50 order 1–10 (save APF/RC1),
        // | continue block 50 order 11–50, compute block 51 (pending).
        if icount == 4 {
            self.adapt_synthesis();
        }

        spf
    }

    /// The decoder's `ICOUNT = 4` once-per-frame adaptation: block 49 →
    /// block 50 interrupted at order 10 (APF by-product) → resumed to
    /// order 50 → block 51 (pending commit at `ICOUNT = 3`).
    fn adapt_synthesis(&mut self) {
        let segs: [BflSegment; N5] = std::array::from_fn(|j| BflSegment {
            samples: self.sttmp[j].st,
            nls: self.sttmp[j].nls,
        });
        let core = self.hw49.run(&synth_window(), &segs);

        // Block 50, order 1 – 10.
        let mut atmp10 = [0i16; PF_ORDER + 1];
        let st10 = levinson_durbin_fixed(
            &LevinsonInput {
                rtmp: &core.rtmp,
                illcond: core.illcond,
                lpc: LPCW,
            },
            &mut atmp10,
        );
        // | NLSAPF = NLSATMP; copy ATMP(2:11) to APF — gated on
        // | ILLCONDP (§G.3.28: "ATMP will not be copied to APF" when
        // | the 10th-order results are invalid).
        if !st10.illcondp {
            self.apf_atil = atmp10;
            self.nlsapf = st10.nlsatmp;
            self.rc1 = st10.rc1;
            // §G.3.28 LABEL — the flat Q13 rendering block 81 consumes.
            self.apf_q13 = apf_to_q13(&atmp10, st10.nlsatmp);
        }
        if st10.illcond {
            // The 10th-order recursion already failed ⇒ the order-50
            // continuation cannot proceed either; keep the old A.
            self.pending_a = None;
            return;
        }

        // | continue block 50, order 11 – 50.
        let resume = RecursionResume {
            atmp: atmp10.to_vec(),
            nrs: st10.nrs,
            alphatmp: st10.alphatmp,
        };
        let mut atmp = [0i16; LPC + 1];
        let st50 = levinson_durbin_fixed_resume(
            &LevinsonInput {
                rtmp: &core.rtmp,
                illcond: core.illcond,
                lpc: LPC,
            },
            &mut atmp,
            &resume,
        );

        // Block 51 (compute now; visible at the next ICOUNT = 3).
        self.pending_a = if st50.illcond {
            None
        } else {
            bandwidth_expand_q14(&atmp[..=LPC], &FACV_Q14, st50.nlsatmp).map(|v| {
                let mut a = [0i16; LPC + 1];
                a.copy_from_slice(&v);
                a
            })
        };
    }

    /// The most recent decoded (pre-postfilter) quantized-speech vector
    /// (mantissas + shift) — for encoder-lockstep verification.
    #[must_use]
    pub fn last_st(&self) -> (&[i16; IDIM], i32) {
        let seg = &self.sttmp[self.icount - 1];
        (&seg.st, seg.nls)
    }

    /// Decode a §5.11 serial byte stream (whole 10-bit indices only —
    /// see [`crate::bitstream::unpack_indices`]) into postfiltered `Q2`
    /// speech vectors.
    pub fn decode_bytes(&mut self, bytes: &[u8]) -> crate::Result<Vec<[i32; IDIM]>> {
        let indices = crate::bitstream::unpack_indices(bytes)?;
        Ok(indices.iter().map(|&i| self.decode_vector(i)).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Synthetic `Q2` speech: a two-tone swept signal with realistic
    /// amplitude (|s| ≤ ~2000 in sample units ⇒ ≤ 8000 in `Q2`).
    fn speech_vector(n: usize) -> [i16; IDIM] {
        std::array::from_fn(|k| {
            let t = (n * IDIM + k) as f64;
            let v =
                (t * 0.089).sin() * 1400.0 + (t * 0.031).sin() * 500.0 + (t * 0.53).sin() * 140.0;
            (v * 4.0).round() as i16
        })
    }

    #[test]
    fn encoder_decoder_quantized_speech_lockstep() {
        // The heart of Annex G: the decoder, fed only the channel
        // indices, must reproduce the encoder's quantized speech ST
        // BIT-EXACTLY (same mantissas at the same shift). 200 vectors
        // = 50 adaptation cycles, exercising every once-per-frame
        // block many times.
        let mut enc = EncoderFixed::new();
        let mut dec = DecoderFixed::new();

        for n in 0..200 {
            let s = speech_vector(n);
            let ichan = enc.encode_vector(&s);
            let _spf = dec.decode_vector(ichan);

            let (st_e, nls_e) = enc.last_st();
            let (st_d, nls_d) = dec.last_st();
            assert_eq!(
                nls_e, nls_d,
                "vector {n}: NLSST diverged (enc {nls_e} vs dec {nls_d})"
            );
            assert_eq!(
                st_e, st_d,
                "vector {n}: quantized speech diverged at NLS {nls_e}"
            );
        }
    }

    #[test]
    fn decoder_output_tracks_input_speech() {
        // End-to-end sanity: after convergence the postfiltered output
        // should correlate strongly with the input (LD-CELP at
        // 16 kbit/s is near-transparent for a stationary two-tone).
        let mut enc = EncoderFixed::new();
        let mut dec = DecoderFixed::new();

        let mut num = 0.0f64;
        let mut den_s = 0.0f64;
        let mut den_o = 0.0f64;
        for n in 0..240 {
            let s = speech_vector(n);
            let ichan = enc.encode_vector(&s);
            let spf = dec.decode_vector(ichan);
            if n >= 40 {
                for k in 0..IDIM {
                    let x = f64::from(s[k]);
                    let y = spf[k] as f64;
                    num += x * y;
                    den_s += x * x;
                    den_o += y * y;
                }
            }
        }
        let corr = num / (den_s.sqrt() * den_o.sqrt());
        assert!(
            corr > 0.9,
            "decoded speech should track the input (corr {corr})"
        );
        // And the energy should be in the same ballpark (postfilter AGC).
        let gain_db = 10.0 * (den_o / den_s).log10();
        assert!(
            gain_db.abs() < 3.0,
            "output energy within 3 dB of input ({gain_db} dB)"
        );
    }

    #[test]
    fn silence_stays_bounded_and_settles() {
        // All-zero input exercises every ILLCOND path each cycle.
        // LD-CELP has no silence escape — the encoder keeps transmitting
        // the least-distortion codevector — and the cold-start backward
        // adapters ride a transient before the loop settles (the ±4095
        // quantized-speech clip bounds it by construction; the encoder
        // and decoder stay in bit-exact lockstep throughout, which the
        // dedicated lockstep test proves on speech and this one re-checks
        // on the degenerate input). Assert (a) hard boundedness, (b) the
        // limit cycle settles to near-silence, and (c) ST lockstep.
        let mut enc = EncoderFixed::new();
        let mut dec = DecoderFixed::new();
        let mut tail_peak = 0i32;
        for n in 0..400 {
            let ichan = enc.encode_vector(&[0i16; IDIM]);
            let spf = dec.decode_vector(ichan);
            let (st_e, nls_e) = enc.last_st();
            let (st_d, nls_d) = dec.last_st();
            assert_eq!((st_e, nls_e), (st_d, nls_d), "lockstep at vector {n}");
            for &v in spf.iter() {
                assert!(v.abs() <= i32::from(i16::MAX), "SPF overflow at {n}: {v}");
            }
            if n >= 200 {
                tail_peak = tail_peak.max(spf.iter().map(|v| v.abs()).max().unwrap());
            }
        }
        assert!(
            tail_peak <= 1024,
            "idle limit cycle must settle to near-silence (tail peak {tail_peak} Q2)"
        );
    }

    #[test]
    fn byte_stream_roundtrip_preserves_lockstep() {
        // §5.11 framing: encode_buffer → bytes → decode_bytes must be
        // exactly the per-vector path (2 cycles = 10 bytes per chunk).
        let mut enc_a = EncoderFixed::new();
        let mut enc_b = EncoderFixed::new();
        let mut dec = DecoderFixed::new();
        let mut dec_ref = DecoderFixed::new();

        let vectors: Vec<[i16; IDIM]> = (0..40).map(speech_vector).collect();
        let bytes = enc_a.encode_buffer(&vectors);
        assert_eq!(bytes.len(), 40 * 10 / 8, "40 vectors = 50 bytes");
        let decoded = dec.decode_bytes(&bytes).expect("whole indices");
        assert_eq!(decoded.len(), 40);

        for (n, v) in vectors.iter().enumerate() {
            let ichan = enc_b.encode_vector(v);
            let spf = dec_ref.decode_vector(ichan);
            assert_eq!(decoded[n], spf, "vector {n} through the byte path");
        }
    }

    #[test]
    fn decoder_survives_arbitrary_index_streams() {
        // Robustness: a decoder fed arbitrary (even adversarial) 10-bit
        // channel indices must never panic and never overflow its Q2
        // output — every ICHAN value is a legal codeword by
        // construction, but the resulting gain trajectory can ride the
        // +28 dB limiter for long stretches.
        let mut dec = DecoderFixed::new();
        let mut seed = 0x0dd5_eed5_0000_1111u64;
        for _ in 0..400 {
            seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            let ichan = ((seed >> 33) % 1024) as u16;
            let spf = dec.decode_vector(ichan);
            for &v in spf.iter() {
                assert!(v.abs() <= i32::from(i16::MAX), "SPF escaped 16 bits: {v}");
            }
        }
        // And every single codeword at least once, in order.
        let mut dec = DecoderFixed::new();
        for ichan in 0..1024u16 {
            let _ = dec.decode_vector(ichan);
        }
    }

    #[test]
    fn channel_index_roundtrip_layout() {
        // ICHAN = (IS − 1)·NG + (IG − 1): the decoder must recover the
        // same (IS, IG) the encoder's SearchResult packed.
        for is in [1usize, 2, 64, 128] {
            for ig in 1..=NG {
                let ichan = ((is - 1) * NG + (ig - 1)) as u16;
                assert_eq!((ichan as usize) / NG + 1, is);
                assert_eq!((ichan as usize) % NG + 1, ig);
            }
        }
    }
}
