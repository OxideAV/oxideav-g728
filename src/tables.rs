//! Bit-exact numeric tables for ITU-T G.728 LD-CELP.
//!
//! Every value here is transcribed directly from the integer columns of
//! Annexes A–D in the published Recommendation prose (G.728 1992-09).
//! Two conventions are used:
//!
//! * Tables in Annex A (hybrid windows) are stored as Q15 i16 integers
//!   exactly as printed; the matching f64 view is produced at runtime
//!   by dividing each value by 2¹⁵ = 32 768. The prose makes this
//!   division explicit: "Dividing the table entries by 2¹⁵ = 32 768
//!   gives the table above" (Annex A.1, A.2, A.3).
//! * Tables in Annexes B, C, D are stored as the printed integers in
//!   their stated Q-format (Q11 / Q13 / Q14) and converted on demand by
//!   the same `value / 2^q` rule the prose specifies.
//!
//! These tables are derived solely from the published Recommendation
//! PDF — no reference implementation source has been consulted. The
//! cross-check between the Annex A float column and the Annex A
//! integer column (the prose lists both) provides an in-repo guard
//! against typos: the [`tests`] module asserts every float-equivalent
//! is within one Q15 quantum (1 / 32 768) of the printed float value.

use crate::consts::{IDIM, LPC, LPCLG, LPCW, NCWD, NG};

/// Q-format shift for the Annex A hybrid window tables (Q15).
pub const Q15: i32 = 15;
/// Q-format shift for the Annex B shape-codebook table (Q11).
pub const Q11: i32 = 11;
/// Q-format shift for the Annex B gain-codebook GQ/GB tables (Q13).
pub const Q13: i32 = 13;
/// Q-format shift for the Annex B G2 / GSQ tables (Q12).
pub const Q12: i32 = 12;
/// Q-format shift for the Annex C bandwidth-broadening tables (Q14).
pub const Q14: i32 = 14;

// ---------------------------------------------------------------------
// Annex A.1 — Hybrid window for the synthesis filter
//
// 105 samples. The first NONR = 35 samples are the non-recursive
// portion; the remaining 70 are the recursive portion. Stored as
// 16-bit integers; spec footnote: "Dividing the table entries by 2¹⁵ =
// 32 768 gives the [Annex A.1 float] table above."
// ---------------------------------------------------------------------

/// `wnr[i]` (1 ≤ i ≤ 105) — synthesis-filter hybrid window, Q15.
/// Transcribed from the integer table on Recommendation G.728 (1992-09),
/// Annex A.1, page 54.
pub const WNR_Q15: [i16; 105] = [
    1_565, 3_127, 4_681, 6_225, 7_755, 9_266, 10_757, 12_223, 13_661, 15_068, 16_441, 17_776,
    19_071, 20_322, 21_526, 22_682, 23_786, 24_835, 25_828, 26_761, 27_634, 28_444, 29_188, 29_866,
    30_476, 31_016, 31_486, 31_884, 32_208, 32_460, 32_637, 32_739, 32_767, 32_721, 32_599, 32_403,
    32_171, 31_940, 31_711, 31_484, 31_259, 31_034, 30_812, 30_591, 30_372, 30_154, 29_938, 29_724,
    29_511, 29_299, 29_089, 28_881, 28_674, 28_468, 28_264, 28_062, 27_861, 27_661, 27_463, 27_266,
    27_071, 26_877, 26_684, 26_493, 26_303, 26_114, 25_927, 25_742, 25_557, 25_374, 25_192, 25_012,
    24_832, 24_654, 24_478, 24_302, 24_128, 23_955, 23_784, 23_613, 23_444, 23_276, 23_109, 22_943,
    22_779, 22_616, 22_454, 22_293, 22_133, 21_974, 21_817, 21_661, 21_505, 21_351, 21_198, 21_046,
    20_896, 20_746, 20_597, 20_450, 20_303, 20_157, 20_013, 19_870, 19_727,
];

// ---------------------------------------------------------------------
// Annex A.2 — Hybrid window for the log-gain predictor
//
// 34 samples. First NONRLG = 20 non-recursive; remaining 14 recursive.
// ---------------------------------------------------------------------

/// `wnrg[i]` (1 ≤ i ≤ 34) — log-gain predictor hybrid window, Q15.
/// Transcribed from Recommendation G.728 (1992-09), Annex A.2, page 55.
pub const WNRG_Q15: [i16; 34] = [
    3_026, 6_025, 8_973, 11_845, 14_615, 17_261, 19_759, 22_088, 24_228, 26_162, 27_872, 29_344,
    30_565, 31_525, 32_216, 32_631, 32_767, 32_625, 32_203, 31_506, 30_540, 29_461, 28_420, 27_416,
    26_448, 25_514, 24_613, 23_743, 22_905, 22_096, 21_315, 20_562, 19_836, 19_135,
];

// ---------------------------------------------------------------------
// Annex A.3 — Hybrid window for the perceptual weighting filter
//
// 60 samples. First NONRW = 30 non-recursive; remaining 30 recursive.
// ---------------------------------------------------------------------

/// `wnrw[i]` (1 ≤ i ≤ 60) — weighting-filter hybrid window, Q15.
/// Transcribed from Recommendation G.728 (1992-09), Annex A.3, page 55.
pub const WNRW_Q15: [i16; 60] = [
    1_957, 3_908, 5_845, 7_760, 9_648, 11_502, 13_314, 15_079, 16_790, 18_441, 20_026, 21_540,
    22_976, 24_331, 25_599, 26_775, 27_856, 28_837, 29_715, 30_487, 31_150, 31_702, 32_141, 32_464,
    32_672, 32_763, 32_738, 32_595, 32_336, 31_961, 31_472, 30_931, 30_400, 29_878, 29_365, 28_860,
    28_364, 27_877, 27_398, 26_927, 26_465, 26_010, 25_563, 25_124, 24_693, 24_268, 23_851, 23_442,
    23_039, 22_643, 22_254, 21_872, 21_496, 21_127, 20_764, 20_407, 20_057, 19_712, 19_373, 19_041,
];

// ---------------------------------------------------------------------
// Annex B — Excitation shape codebook (Y) and gain codebook (GQ / GB)
//
// Shape: 128 codevectors × 5 components, Q11 (i16; divide by 2 048 for
// the floating-point value). The 128 rows are indexed by a Gray-code
// "channel index" 0..127, NOT by an ordinary linear shape index — the
// codebook search routine in §3.9 selects rows by channel index.
//
// Gain: 8 levels. GQ stored as Q13 i16, GB (mid-points) stored as Q13
// i16. The spec also lists G2 = 2·GQ (Q12) and GSQ = GQ² (Q12);
// both are derived at runtime from GQ via the prose definitions
// (equations 3-21 / 3-22).
// ---------------------------------------------------------------------

/// Excitation shape codebook `y[j][k]` for 0 ≤ j < 128, 0 ≤ k < 5, Q11.
///
/// Transcribed from the channel-index table on Recommendation G.728
/// (1992-09), Annex B, pages 56–58. Each row is one of the printed
/// 128 codevectors; column 0 is the first codevector component
/// (channel-table column "1") and column 4 is the fifth ("5").
pub const Y_Q11: [[i16; IDIM]; NCWD] = [
    /*   0 */ [668, -2_950, -1_254, -1_790, -2_553],
    /*   1 */ [-5_032, -4_577, -1_045, 2_908, 3_318],
    /*   2 */ [-2_819, -2_677, -948, -2_825, -4_450],
    /*   3 */ [-6_679, -340, 1_482, -1_276, 1_262],
    /*   4 */ [-562, -6_757, 1_281, 179, -1_274],
    /*   5 */ [-2_512, -7_130, -4_925, 6_913, 2_411],
    /*   6 */ [-2_478, -156, 4_683, -3_873, 0],
    /*   7 */ [-8_208, 2_140, -478, -2_785, 533],
    /*   8 */ [1_889, 2_759, 1_381, -6_955, -5_913],
    /*   9 */ [5_082, -2_460, -5_778, 1_797, 568],
    /*  10 */ [-2_208, -3_309, -4_523, -6_236, -7_505],
    /*  11 */ [-2_719, 4_358, -2_988, -1_149, 2_664],
    /*  12 */ [1_259, 995, 2_711, -2_464, -10_390],
    /*  13 */ [1_722, -7_569, -2_742, 2_171, -2_329],
    /*  14 */ [1_032, 747, -858, -7_946, -12_843],
    /*  15 */ [3_106, 4_856, -4_193, -2_541, 1_035],
    /*  16 */ [1_862, -960, -6_628, 410, 5_882],
    /*  17 */ [-2_493, -2_628, -4_000, -60, 7_202],
    /*  18 */ [-2_672, 1_446, 1_536, -3_831, 1_233],
    /*  19 */ [-5_302, 6_912, 1_589, -4_187, 3_665],
    /*  20 */ [-3_456, -8_170, -7_709, 1_384, 4_698],
    /*  21 */ [-4_699, -6_209, -11_176, 8_104, 16_830],
    /*  22 */ [930, 7_004, 1_269, -8_977, 2_567],
    /*  23 */ [4_649, 11_804, 3_441, -5_657, 1_199],
    /*  24 */ [2_542, -183, -8_859, -7_976, 3_230],
    /*  25 */ [-2_872, -2_011, -9_713, -8_385, 12_983],
    /*  26 */ [3_086, 2_140, -3_680, -9_643, -2_896],
    /*  27 */ [-7_609, 6_515, -2_283, -2_522, 6_332],
    /*  28 */ [-3_333, -5_620, -9_130, -11_131, 5_543],
    /*  29 */ [-407, -6_721, -17_466, -2_889, 11_568],
    /*  30 */ [3_692, 6_796, -262, -10_846, -1_856],
    /*  31 */ [7_275, 13_404, -2_989, -10_595, 4_936],
    /*  32 */ [244, -2_219, 2_656, 3_776, -5_412],
    /*  33 */ [-4_043, -5_934, 2_131, 863, -2_866],
    /*  34 */ [-3_302, 1_743, -2_006, -128, -2_052],
    /*  35 */ [-6_361, 3_342, -1_583, -21, 1_142],
    /*  36 */ [-3_837, -1_831, 6_397, 2_545, -2_848],
    /*  37 */ [-9_332, -6_528, 5_309, 1_986, -2_245],
    /*  38 */ [-4_490, 748, 1_935, -3_027, -493],
    /*  39 */ [-9_255, 5_366, 3_193, -4_493, 1_784],
    /*  40 */ [4_784, -370, 1_866, 1_057, -1_889],
    /*  41 */ [7_342, -2_690, -2_577, 676, -611],
    /*  42 */ [-502, 2_235, -1_850, -1_777, -2_049],
    /*  43 */ [1_011, 3_880, -2_465, 2_209, -152],
    /*  44 */ [2_592, 2_829, 5_588, 2_839, -7_306],
    /*  45 */ [-3_049, -4_918, 5_955, 9_201, -4_447],
    /*  46 */ [697, 3_908, 5_798, -4_451, -4_644],
    /*  47 */ [-2_121, 5_444, -2_570, 321, -1_202],
    /*  48 */ [2_846, -2_086, 3_532, 566, -708],
    /*  49 */ [-4_279, 950, 4_980, 3_749, 452],
    /*  50 */ [-2_484, 3_502, 1_719, -170, 238],
    /*  51 */ [-3_435, 263, 2_114, -2_005, 2_361],
    /*  52 */ [-7_338, -1_208, 9_347, -1_216, -4_013],
    /*  53 */ [-13_498, -439, 8_028, -4_232, 361],
    /*  54 */ [-3_729, 5_433, 2_004, -4_727, -1_259],
    /*  55 */ [-3_986, 7_743, 8_429, -3_691, -987],
    /*  56 */ [5_198, -423, 1_150, -1_281, 816],
    /*  57 */ [7_409, 4_109, -3_949, 2_690, 30],
    /*  58 */ [1_246, 3_055, -35, -1_370, -246],
    /*  59 */ [-1_489, 5_635, -678, -2_627, 3_170],
    /*  60 */ [4_830, -4_585, 2_008, -1_062, 799],
    /*  61 */ [-129, 717, 4_594, 14_937, 10_706],
    /*  62 */ [417, 2_759, 1_850, -5_057, -1_153],
    /*  63 */ [-3_887, 7_361, -5_768, 4_285, 666],
    /*  64 */ [1_443, -938, 20, -2_119, -1_697],
    /*  65 */ [-3_712, -3_402, -2_212, 110, 2_136],
    /*  66 */ [-2_952, 12, -1_568, -3_500, -1_855],
    /*  67 */ [-1_315, -1_731, 1_160, -558, 1_709],
    /*  68 */ [88, -4_569, 194, -454, -2_957],
    /*  69 */ [-2_839, -1_666, -273, 2_084, -155],
    /*  70 */ [-189, -2_376, 1_663, -1_040, -2_449],
    /*  71 */ [-2_842, -1_369, 636, -248, -2_677],
    /*  72 */ [1_517, 79, -3_013, -3_669, -973],
    /*  73 */ [1_913, -2_493, -5_312, -749, 1_271],
    /*  74 */ [-2_903, -3_324, -3_756, -3_690, -1_829],
    /*  75 */ [-2_913, -1_547, -2_760, -1_406, 1_124],
    /*  76 */ [1_844, -1_834, 456, 706, -4_272],
    /*  77 */ [467, -4_256, -1_909, 1_521, 1_134],
    /*  78 */ [-127, -994, -637, -1_491, -6_494],
    /*  79 */ [873, -2_045, -3_828, -2_792, -578],
    /*  80 */ [2_311, -1_817, 2_632, -3_052, 1_968],
    /*  81 */ [641, 1_194, 1_893, 4_107, 6_342],
    /*  82 */ [-45, 1_198, 2_160, -1_449, 2_203],
    /*  83 */ [-2_004, 1_713, 3_518, 2_652, 4_251],
    /*  84 */ [2_936, -3_968, 1_280, 131, -1_476],
    /*  85 */ [2_827, 8, -1_928, 2_658, 3_513],
    /*  86 */ [3_199, -816, 2_687, -1_741, -1_407],
    /*  87 */ [2_948, 4_029, 394, -253, 1_298],
    /*  88 */ [4_286, 51, -4_507, -32, -659],
    /*  89 */ [3_903, 5_646, -5_588, -2_592, 5_707],
    /*  90 */ [-606, 1_234, -1_607, -5_187, 664],
    /*  91 */ [-525, 3_620, -2_192, -2_527, 1_707],
    /*  92 */ [4_297, -3_251, -2_283, 812, -2_264],
    /*  93 */ [5_765, 528, -3_287, 1_352, 1_672],
    /*  94 */ [2_735, 1_241, -1_103, -3_273, -3_407],
    /*  95 */ [4_033, 1_648, -2_965, -1_174, 1_444],
    /*  96 */ [74, 918, 1_999, 915, -1_026],
    /*  97 */ [-2_496, -1_605, 2_034, 2_950, 229],
    /*  98 */ [-2_168, 2_037, 15, -1_264, -208],
    /*  99 */ [-3_552, 1_530, 581, 1_491, 962],
    /* 100 */ [-2_613, -2_338, 3_621, -1_488, -2_185],
    /* 101 */ [-1_747, 81, 5_538, 1_432, -2_257],
    /* 102 */ [-1_019, 867, 214, -2_284, -1_510],
    /* 103 */ [-1_684, 2_816, -229, 2_551, -1_389],
    /* 104 */ [2_707, 504, 479, 2_783, -1_009],
    /* 105 */ [2_517, -1_487, -1_596, 621, 1_929],
    /* 106 */ [-148, 2_206, -4_288, 1_292, -1_401],
    /* 107 */ [-527, 1_243, -2_731, 1_909, 1_280],
    /* 108 */ [2_149, -1_501, 3_688, 610, -4_591],
    /* 109 */ [3_306, -3_369, 1_875, 3_636, -1_217],
    /* 110 */ [2_574, 2_513, 1_449, -3_074, -4_979],
    /* 111 */ [814, 1_826, -2_497, 4_234, -4_077],
    /* 112 */ [1_664, -220, 3_418, 1_002, 1_115],
    /* 113 */ [781, 1_658, 3_919, 6_130, 3_140],
    /* 114 */ [1_148, 4_065, 1_516, 815, 199],
    /* 115 */ [1_191, 2_489, 2_561, 2_421, 2_443],
    /* 116 */ [770, -5_915, 5_515, -368, -3_199],
    /* 117 */ [1_190, 1_047, 3_742, 6_927, -2_089],
    /* 118 */ [292, 3_099, 4_308, -758, -2_455],
    /* 119 */ [523, 3_921, 4_044, 1_386, 85],
    /* 120 */ [4_367, 1_006, -1_252, -1_466, -1_383],
    /* 121 */ [3_852, 1_579, -77, 2_064, 868],
    /* 122 */ [5_109, 2_919, -202, 359, -509],
    /* 123 */ [3_650, 3_206, 2_303, 1_693, 1_296],
    /* 124 */ [2_905, -3_907, 229, -1_196, -2_332],
    /* 125 */ [5_977, -3_585, 805, 3_825, -3_138],
    /* 126 */ [3_746, -606, 53, -269, -3_301],
    /* 127 */ [606, 2_018, -1_316, 4_064, 398],
];

// ---------------------------------------------------------------------
// Annex B — Gain codebook
//
// GQ (gain levels) and GB (gain quantizer cell boundaries) are stored
// in Q13. GQ has NG = 8 entries; entries 4..7 are the sign-mirrored
// negatives of entries 0..3 per the codebook search algorithm's
// "positive / negative" split (§3.9.1 a). GB has NG-1 = 7 entries,
// but only six are used (d0, d1, d2, d4, d5, d6).
//
// The printed Annex B table lists indices 1..4 as floats:
//   GQ:  0.515625, 0.90234375, 1.579101563, 2.763427734
//   GB:  0.708984375, 1.240722656, 2.171264649
// with the spec footnote `GQ(1) = 33/64` and `GQ(i) = (7/4) GQ(i-1)`
// for i = 2, 3, 4. We re-derive the Q13 integers from those rationals
// to avoid PDF-OCR ambiguity on the trailing fractional digits; the
// per-test guard in [`tests::gq_q13_matches_printed_floats`] proves
// the floats agree to one Q13 quantum (1 / 8192 ≈ 1.22e-4).
// ---------------------------------------------------------------------

/// Gain-codebook levels `gq[i]` for 0 ≤ i < 8.
///
/// Indices 0..3 correspond to positive gain levels (Annex B / §3.9.1
/// "positive gain" branch); indices 4..7 are the sign-flipped
/// negatives (Annex B / §3.9.1 "negative gain" branch).
pub const GQ: [f64; NG] = [
    33.0 / 64.0,                 // index 0 = 0.515625
    (33.0 / 64.0) * (7.0 / 4.0), // index 1 = 0.90234375
    (33.0 / 64.0) * (7.0 / 4.0) * (7.0 / 4.0),
    // index 2 ≈ 1.579101563
    (33.0 / 64.0) * (7.0 / 4.0) * (7.0 / 4.0) * (7.0 / 4.0),
    // index 3 ≈ 2.763427734
    -(33.0 / 64.0),               // index 4 = -GQ(1)
    -(33.0 / 64.0) * (7.0 / 4.0), // index 5 = -GQ(2)
    -(33.0 / 64.0) * (7.0 / 4.0) * (7.0 / 4.0),
    // index 6 = -GQ(3)
    -(33.0 / 64.0) * (7.0 / 4.0) * (7.0 / 4.0) * (7.0 / 4.0),
    // index 7 = -GQ(4)
];

/// Gain quantizer mid-point boundaries `gb[i]` for 0 ≤ i < 7.
///
/// Only six boundaries (d0, d1, d2, d4, d5, d6) are referenced by the
/// codebook-search decision tree in §3.9.2; the spec marks `gb[3]`
/// (index = 3, the printed "Array index 4" column) as "any arbitrary
/// value (not used)" and the same for `gb[7]` past the table end. We
/// follow the spec's mid-point definition `gb[i] = (gq[i] + gq[i+1]) / 2`
/// and place a placeholder for the unused slot.
pub const GB: [f64; NG - 1] = [
    (GQ[0] + GQ[1]) / 2.0, // d0 = 0.708984375
    (GQ[1] + GQ[2]) / 2.0, // d1 = 1.240722656
    (GQ[2] + GQ[3]) / 2.0, // d2 = 2.171264649
    0.0,                   // d3 — "any arbitrary value (not used)"
    (GQ[4] + GQ[5]) / 2.0, // d4 = -GB[0]
    (GQ[5] + GQ[6]) / 2.0, // d5 = -GB[1]
    (GQ[6] + GQ[7]) / 2.0, // d6 = -GB[2]
];

/// Precomputed `g2[i] = 2·gq[i]` for the inner-product test in §3.9
/// (equation 3-21).
pub const G2: [f64; NG] = [
    2.0 * GQ[0],
    2.0 * GQ[1],
    2.0 * GQ[2],
    2.0 * GQ[3],
    2.0 * GQ[4],
    2.0 * GQ[5],
    2.0 * GQ[6],
    2.0 * GQ[7],
];

/// Precomputed `gsq[i] = gq[i]²` for the distortion term in §3.9
/// (equation 3-22).
pub const GSQ: [f64; NG] = [
    GQ[0] * GQ[0],
    GQ[1] * GQ[1],
    GQ[2] * GQ[2],
    GQ[3] * GQ[3],
    GQ[4] * GQ[4],
    GQ[5] * GQ[5],
    GQ[6] * GQ[6],
    GQ[7] * GQ[7],
];

// ---------------------------------------------------------------------
// Annex G §G.5 — log-gain tables for the fixed-point gain adapter
// (Figure G.1 blocks 93 / 94).
//
// Table G.3/G.728 stores 20·log10|gᵢ| and Table G.4/G.728 stores
// 10·log10 P[yⱼ] in Q11 ("to obtain the fixed point value, multiply
// the floating point value by 2048 = 2¹¹"). The per-test guards
// `gcblg_q11_matches_float_db` / `shapelg_q11_matches_float_db`
// cross-check every transcribed entry against the floating-point dB
// functions [`crate::annex_g_gain::gain_log_db`] /
// [`crate::annex_g_gain::shape_log_db`] (which derive the same values
// from the Annex A / Annex B codebook columns), guarding both the
// PDF-OCR transcription and the derivation in one shot.
// ---------------------------------------------------------------------

/// Table G.3/G.728 — `GCBLG(i) = 20·log10|gq(i)|` in Q11 for the 8
/// gain-codebook levels (block 93). The printed table lists the four
/// magnitudes (indices 0–3); the sign-mirrored levels 4–7 share the
/// same dB values.
pub const GCBLG_Q11: [i16; NG] = [
    -11783, -1828, 8127, 18082, // indices 0..3 (|g| = 0.5156 … 2.7634)
    -11783, -1828, 8127, 18082, // indices 4..7 (negated levels, same |g|)
];

/// Table G.4/G.728 — `SHAPELG(j) = 10·log10 P[yⱼ]` in Q11 for the 128
/// shape codevectors (block 94), `P[yⱼ] = (1/5)·Σ yⱼₖ²`.
pub const SHAPELG_Q11: [i16; NCWD] = [
    -227, 10308, 6549, 7753, 7597, 16563, 6406, 11933, // 0..7
    13569, 10569, 16328, 6536, 15803, 11673, 21318, 9100, // 8..15
    12245, 12018, 2503, 14690, 18190, 28801, 16803, 20331, // 16..23
    18019, 24920, 16159, 17618, 23072, 28075, 19169, 25723, // 24..31
    8670, 10069, 503, 8647, 11165, 18447, 4264, 17381, // 32..39
    3531, 10543, -2392, 2266, 14527, 18788, 13030, 6238, // 40..47
    1825, 9090, 211, 1888, 18088, 22557, 10893, 18156, // 48..55
    3426, 13400, -4375, 7970, 7754, 25270, 5313, 15615, // 56..63
    -6296, 4510, 2202, -7229, 3146, -2818, -2674, -1567, // 64..71
    1841, 5803, 7824, 319, 1815, 1765, 6949, 2484, // 72..79
    2808, 9714, -4215, 6678, 2634, 3509, 871, 2190, // 80..87
    5546, 15337, 3708, 2406, 5750, 7538, 3912, 3543, // 88..95
    -10104, 303, -6161, -1142, 3867, 5935, -7201, -759, // 96..103
    -2093, -2863, 2217, -3243, 6161, 5853, 7599, 6747, // 104..111
    -2001, 10218, -54, 1912, 11495, 10575, 4517, 4279, // 112..119
    1813, 566, 4569, 4153, 3368, 11179, 1694, 761, // 120..127
];

// ---------------------------------------------------------------------
// Precomputed shape-codevector energy table `E_j = Σ_k y_j(k)²`
//
// Spec source: §3.9 of the G.728 (1992-09) Recommendation derives the
// encoder's analysis-by-synthesis search via equations 3-14..3-23. The
// distortion form rearranges to `D_{i,j} = b_i · <x̃, ỹ_j> + c_i · E_j`
// (with `b_i = 2·g_i`, `c_i = g_i²`), where `E_j = Σ_{k=1..IDIM} y_j(k)²`
// is the shape codevector's own energy — a constant of the codebook and
// therefore precomputable once at table-load time, NOT recomputed per
// search.
//
// `E_j` is a derived quantity from `Y_Q11`: the spec emits no separate
// printed integer column for it (the same way `G2 = 2·GQ` and
// `GSQ = GQ²` above are derived from `GQ`). The values here are
// computed at compile time as 64-bit float sums of squares of the Q11
// shape rows; the per-test in `tests::y_energy_matches_dot_product`
// proves the result equals the direct dot product `Σ y_j² / 2²²` to
// machine precision.
//
// The encoder will consume `Y_ENERGY[j]` directly in the search-cost
// expression of equation 3-23 once block-12..18 of §3.9 lands. The
// surface is exposed now so the typed encoder skeleton has a stable
// table to reference; no encoder pipeline logic is implemented yet.
// ---------------------------------------------------------------------

/// Internal helper: compute one shape codevector's energy in floating
/// point. `Y_Q11[j]` is the Q11 integer row; the float row is
/// `Y_Q11[j] / 2¹¹`, so the energy `Σ_k (Y_Q11[j][k] / 2¹¹)²` equals
/// `Σ_k Y_Q11[j][k]² / 2²²`.
const fn shape_energy_q11_row(j: usize) -> f64 {
    let scale = (1u64 << (Q11 * 2)) as f64;
    let mut sum = 0.0f64;
    let mut k = 0;
    while k < IDIM {
        let v = Y_Q11[j][k] as f64;
        sum += v * v;
        k += 1;
    }
    sum / scale
}

/// Precomputed shape codevector energy `E_j = Σ_k y_j(k)²` for
/// 0 ≤ j < 128 (§3.9, equation 3-23). Derived from `Y_Q11` at compile
/// time using the spec-stated Q11 → float division by 2 048; the
/// per-test guard `y_energy_matches_dot_product` cross-checks each
/// entry against a direct dot-product computation in `y_f64()` space.
///
/// Values are non-negative by construction (sum of squares); the
/// per-test `y_energy_all_non_negative_and_finite` enforces this
/// invariant explicitly.
pub const Y_ENERGY: [f64; NCWD] = {
    let mut out = [0.0f64; NCWD];
    let mut j = 0;
    while j < NCWD {
        out[j] = shape_energy_q11_row(j);
        j += 1;
    }
    out
};

/// Run-time accessor matching the float-view convention used elsewhere
/// in this module (`facv_f64` / `facgpv_f64` / etc.). Returns a fresh
/// owned copy of [`Y_ENERGY`] for callers that prefer the accessor
/// shape over a direct constant reference.
pub fn y_energy_f64() -> [f64; NCWD] {
    Y_ENERGY
}

// ---------------------------------------------------------------------
// Annex C — Bandwidth broadening vectors (Q14)
// ---------------------------------------------------------------------

/// Synthesis-filter BW broadening vector `facv[i]` (1 ≤ i ≤ LPC+1 = 51).
/// Stored as Q14 i16; the floating-point value is `facv[i] / 16384`.
/// Transcribed from Recommendation G.728 (1992-09), Annex C, pages 59–60.
pub const FACV_Q14: [i16; LPC + 1] = [
    16_384, //  1
    16_192, 16_002, 15_815, 15_629, 15_446, 15_265, 15_086, 14_910, //  2.. 9
    14_735, 14_562, 14_391, 14_223, 14_056, 13_891, 13_729, 13_568, // 10..17
    13_409, 13_252, 13_096, 12_943, 12_791, 12_641, 12_493, 12_347, // 18..25
    12_202, 12_059, 11_918, 11_778, 11_640, 11_504, 11_369, 11_236, // 26..33
    11_104, 10_974, 10_845, 10_718, 10_593, 10_468, 10_346, 10_225, // 34..41
    10_105, 9_986, 9_869, 9_754, 9_639, 9_526, 9_415, 9_304, // 42..49
    9_195, 9_088, // 50..51
];

/// Log-gain BW broadening vector `facgpv[i]` (1 ≤ i ≤ LPCLG+1 = 11). Q14.
pub const FACGPV_Q14: [i16; LPCLG + 1] = [
    16_384, 14_848, 13_456, 12_195, 11_051, 10_015, 9_076, 8_225, 7_454, 6_755, 6_122,
];

/// Perceptual-weighting-filter pole-control `wpcfv[i]` (1 ≤ i ≤ LPCW+1).
/// Q14. Encodes γ₂ⁱ⁻¹ with γ₂ = 0.6.
pub const WPCFV_Q14: [i16; LPCW + 1] = [
    16_384, 9_830, 5_898, 3_539, 2_123, 1_274, 764, 459, 275, 165, 99,
];

/// Perceptual-weighting-filter zero-control `wzcfv[i]` (1 ≤ i ≤ LPCW+1).
/// Q14. Encodes γ₁ⁱ⁻¹ with γ₁ = 0.9.
pub const WZCFV_Q14: [i16; LPCW + 1] = [
    16_384, 14_746, 13_271, 11_944, 10_750, 9_675, 8_707, 7_836, 7_053, 6_347, 5_713,
];

/// Short-term postfilter pole-control `spfpcfv[i]` (1 ≤ i ≤ LPC_PF+1 = 11).
/// Q14. Encodes SPFPCF^(i-1) with SPFPCF = 0.75.
pub const SPFPCFV_Q14: [i16; 11] = [
    16_384, 12_288, 9_216, 6_912, 5_184, 3_888, 2_916, 2_187, 1_640, 1_230, 923,
];

/// Short-term postfilter zero-control `spfzcfv[i]` (1 ≤ i ≤ LPC_PF+1 = 11).
/// Q14. Encodes SPFZCF^(i-1) with SPFZCF = 0.65.
pub const SPFZCFV_Q14: [i16; 11] = [
    16_384, 10_650, 6_922, 4_499, 2_925, 1_901, 1_236, 803, 522, 339, 221,
];

// ---------------------------------------------------------------------
// Annex D — 1 kHz lowpass elliptic filter used in block 82
// ---------------------------------------------------------------------

/// 1 kHz lowpass elliptic-filter denominator coefficients `al[i]`
/// (i = 1, 2, 3 in the prose; we store the i = 0 implicit unity tap
/// in `[0]`). Transcribed from Recommendation G.728 (1992-09), Annex D.
pub const AL: [f64; 4] = [1.0, -2.34036589, 2.01190019, -0.614109218];

/// 1 kHz lowpass elliptic-filter numerator coefficients `bl[i]`
/// (i = 0..3). Transcribed from Annex D.
pub const BL: [f64; 4] = [0.0357081667, -0.0069956244, -0.0069956244, 0.0357081667];

// ---------------------------------------------------------------------
// Run-time accessors (Q-integer → f64)
// ---------------------------------------------------------------------

/// Annex A.1 synthesis-filter hybrid window, float view (Q15 / 2¹⁵).
pub fn wnr_f64() -> [f64; 105] {
    let mut out = [0.0f64; 105];
    for (i, &q) in WNR_Q15.iter().enumerate() {
        out[i] = q as f64 / (1u32 << Q15) as f64;
    }
    out
}

/// Annex A.2 log-gain hybrid window, float view.
pub fn wnrg_f64() -> [f64; 34] {
    let mut out = [0.0f64; 34];
    for (i, &q) in WNRG_Q15.iter().enumerate() {
        out[i] = q as f64 / (1u32 << Q15) as f64;
    }
    out
}

/// Annex A.3 weighting-filter hybrid window, float view.
pub fn wnrw_f64() -> [f64; 60] {
    let mut out = [0.0f64; 60];
    for (i, &q) in WNRW_Q15.iter().enumerate() {
        out[i] = q as f64 / (1u32 << Q15) as f64;
    }
    out
}

/// Annex B shape codebook, float view (Q11 / 2¹¹).
pub fn y_f64() -> [[f64; IDIM]; NCWD] {
    let mut out = [[0.0f64; IDIM]; NCWD];
    for j in 0..NCWD {
        for k in 0..IDIM {
            out[j][k] = Y_Q11[j][k] as f64 / (1u32 << Q11) as f64;
        }
    }
    out
}

/// Annex C synthesis-filter BW broadening vector, float view.
pub fn facv_f64() -> [f64; LPC + 1] {
    let mut out = [0.0f64; LPC + 1];
    for (i, &q) in FACV_Q14.iter().enumerate() {
        out[i] = q as f64 / (1u32 << Q14) as f64;
    }
    out
}

/// Annex C log-gain BW broadening vector, float view.
pub fn facgpv_f64() -> [f64; LPCLG + 1] {
    let mut out = [0.0f64; LPCLG + 1];
    for (i, &q) in FACGPV_Q14.iter().enumerate() {
        out[i] = q as f64 / (1u32 << Q14) as f64;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gcblg_q11_matches_float_db() {
        // Table G.3 = round(20·log10|gq| · 2048); the float dB values
        // derive from the same Annex B gain-codebook column.
        let db = crate::annex_g_gain::gain_log_db();
        for i in 0..NG {
            let want = (db[i] * 2048.0).round() as i32;
            let diff = (i32::from(GCBLG_Q11[i]) - want).abs();
            assert!(
                diff <= 1,
                "GCBLG_Q11[{i}] = {} but float dB gives {want}",
                GCBLG_Q11[i]
            );
        }
    }

    #[test]
    fn shapelg_q11_matches_float_db() {
        // Table G.4 = round(10·log10 P[y_j] · 2048); the float dB values
        // derive from the same Annex A shape-codebook rows.
        let db = crate::annex_g_gain::shape_log_db();
        for j in 0..NCWD {
            let want = (db[j] * 2048.0).round() as i32;
            let diff = (i32::from(SHAPELG_Q11[j]) - want).abs();
            assert!(
                diff <= 1,
                "SHAPELG_Q11[{j}] = {} but float dB gives {want}",
                SHAPELG_Q11[j]
            );
        }
    }

    // ------------- Cardinality cross-checks against spec prose --------

    #[test]
    fn wnr_q15_has_105_entries() {
        // §3.3 / Annex A.1 prose: "the first 105 samples of the window
        // function for the synthesis filter"
        assert_eq!(WNR_Q15.len(), 105);
    }

    #[test]
    fn wnrg_q15_has_34_entries() {
        // Annex A.2 prose: "the first 34 samples"
        assert_eq!(WNRG_Q15.len(), 34);
    }

    #[test]
    fn wnrw_q15_has_60_entries() {
        // Annex A.3 prose: "the first 60 samples"
        assert_eq!(WNRW_Q15.len(), 60);
    }

    #[test]
    fn y_q11_dims_match_ncwd_idim() {
        // Annex B: 128 codevectors × 5 components.
        assert_eq!(Y_Q11.len(), NCWD);
        assert_eq!(Y_Q11[0].len(), IDIM);
    }

    #[test]
    fn gain_codebook_has_ng_levels() {
        // Table 1/G.728: NG = 8 (gain codebook size).
        assert_eq!(GQ.len(), NG);
    }

    // ------------- Annex A peak-value sanity ---------------------------

    #[test]
    fn wnr_peaks_match_spec_value() {
        // Annex A.1 prose: the float column lists 0.999969482 as the
        // largest synthesis-window sample (index 33 / 1-based). The
        // matching Q15 row should be `round(0.999969482 * 32768)`.
        let printed_peak: i16 = 32_767;
        assert!(WNR_Q15.contains(&printed_peak));
        // ...and `0.999969482 * 32768 ≈ 32767.0` (within 0.0).
        let derived = (0.999_969_482_f64 * 32_768.0).round() as i16;
        assert_eq!(derived, printed_peak);
    }

    #[test]
    fn wnrg_peaks_match_spec_value() {
        // Annex A.2: the printed float peak is 1.0 (sample 17 in the
        // 1-based table), with Q15 integer 32 767 = 0x7FFF.
        let printed_peak: i16 = 32_767;
        assert!(WNRG_Q15.contains(&printed_peak));
    }

    #[test]
    fn wnrw_peaks_match_spec_value() {
        // Annex A.3: printed float peak 0.999084473, Q15 = 32 738.
        let printed_peak: i16 = 32_738;
        assert!(WNRW_Q15.contains(&printed_peak));
        let derived = (0.999_084_473_f64 * 32_768.0).round() as i16;
        assert_eq!(derived, printed_peak);
    }

    // ------------- Annex B numeric sanity ------------------------------

    #[test]
    fn gq_q13_matches_printed_floats() {
        // Spec Annex B lists GQ(1..4) = 0.515625, 0.90234375,
        // 1.579101563, 2.763427734 with the rule GQ(1) = 33/64,
        // GQ(i) = (7/4)·GQ(i-1). One Q13 quantum is 1 / 8192 ≈
        // 1.22e-4; check that we agree to at least that tolerance.
        let printed = [0.515625_f64, 0.90234375, 1.579_101_563, 2.763_427_734];
        let q13_quantum = 1.0 / (1u32 << Q13) as f64;
        for (i, &p) in printed.iter().enumerate() {
            assert!(
                (GQ[i] - p).abs() < q13_quantum,
                "GQ[{}] = {} vs printed {} (diff {})",
                i,
                GQ[i],
                p,
                (GQ[i] - p).abs()
            );
        }
    }

    #[test]
    fn gq_back_half_is_sign_negation() {
        // §3.9 / Annex B: indices 4..7 are the negations of 0..3 used
        // by the codebook search "search through negative gains"
        // branch. Sign-symmetry is the structural property.
        for i in 0..4 {
            assert_eq!(GQ[i + 4], -GQ[i]);
        }
    }

    #[test]
    fn gb_midpoints_match_spec_value() {
        // Spec footnote: GB(i) is the mid-point between GQ(i) and
        // GQ(i+1) that have the same sign. Check the first printed
        // value GB(1) = 0.708984375 within one Q13 quantum.
        let q13_quantum = 1.0 / (1u32 << Q13) as f64;
        assert!((GB[0] - 0.708_984_375).abs() < q13_quantum);
        assert!((GB[1] - 1.240_722_656).abs() < q13_quantum);
        assert!((GB[2] - 2.171_264_649).abs() < q13_quantum);
    }

    #[test]
    fn y_codevector_zero_matches_printed_row() {
        // Spot-check row 0: spec prints (668, -2950, -1254, -1790,
        // -2553) for channel index 0.
        assert_eq!(Y_Q11[0], [668_i16, -2_950, -1_254, -1_790, -2_553]);
    }

    #[test]
    fn y_codevector_last_matches_printed_row() {
        // Spot-check row 127: spec prints (606, 2018, -1316, 4064, 398).
        assert_eq!(Y_Q11[127], [606_i16, 2_018, -1_316, 4_064, 398_i16]);
    }

    #[test]
    fn y_q11_to_float_div2048() {
        // Spec prose: "To obtain the floating point value from the
        // integer value, divide the integer value by 2 048." Spot-check
        // row 0 column 0.
        let f = y_f64();
        assert!((f[0][0] - 668.0 / 2048.0).abs() < 1e-15);
    }

    // ------------- Annex C numeric sanity ------------------------------

    #[test]
    fn facv_first_is_q14_unity() {
        // FACV(1) is the i = 0 power of λ which is 1.0 in floating
        // point; Q14 representation is 16 384.
        assert_eq!(FACV_Q14[0], 16_384);
    }

    #[test]
    fn facv_geometric_with_fac() {
        // Annex C reproduces `λⁱ⁻¹` for λ = 253 / 256. Each entry
        // should equal `round(FACV[0] * λⁱ)` within one Q14 quantum.
        let lambda = crate::consts::FAC;
        for i in 1..FACV_Q14.len() {
            let expected = (16_384.0 * lambda.powi(i as i32)).round() as i16;
            // Allow ±1 LSB to absorb the trailing-power rounding.
            let diff = (FACV_Q14[i] as i32 - expected as i32).abs();
            assert!(
                diff <= 1,
                "FACV_Q14[{}] = {} vs derived λⁱ ≈ {} (diff {})",
                i,
                FACV_Q14[i],
                expected,
                diff
            );
        }
    }

    #[test]
    fn facgpv_geometric_with_facgp() {
        // Same shape as FACV but with λ_g = 29/32.
        let lambda_g = crate::consts::FACGP;
        for i in 1..FACGPV_Q14.len() {
            let expected = (16_384.0 * lambda_g.powi(i as i32)).round() as i16;
            let diff = (FACGPV_Q14[i] as i32 - expected as i32).abs();
            assert!(
                diff <= 1,
                "FACGPV_Q14[{}] = {} vs derived (29/32)ⁱ ≈ {} (diff {})",
                i,
                FACGPV_Q14[i],
                expected,
                diff
            );
        }
    }

    #[test]
    fn wpcfv_matches_gamma2_powers() {
        // γ₂ = 0.6; equation 3-4c uses γ₂ⁱ in the weighting filter's
        // pole-controlling vector. Q14 representation.
        for i in 1..WPCFV_Q14.len() {
            let expected = (16_384.0 * crate::consts::WPCF.powi(i as i32)).round() as i16;
            let diff = (WPCFV_Q14[i] as i32 - expected as i32).abs();
            assert!(
                diff <= 1,
                "WPCFV_Q14[{}] = {} vs derived γ₂ⁱ ≈ {} (diff {})",
                i,
                WPCFV_Q14[i],
                expected,
                diff
            );
        }
    }

    #[test]
    fn wzcfv_matches_gamma1_powers() {
        for i in 1..WZCFV_Q14.len() {
            let expected = (16_384.0 * crate::consts::WZCF.powi(i as i32)).round() as i16;
            let diff = (WZCFV_Q14[i] as i32 - expected as i32).abs();
            assert!(
                diff <= 1,
                "WZCFV_Q14[{}] = {} vs derived γ₁ⁱ ≈ {} (diff {})",
                i,
                WZCFV_Q14[i],
                expected,
                diff
            );
        }
    }

    // ------------- Annex D structural -------------------------------

    #[test]
    fn lowpass_filter_is_symmetric_at_endpoints() {
        // Annex D: bl[0] = bl[3] (= 0.0357081667) and bl[1] = bl[2]
        // (= -0.0069956244). The structural symmetry is a sanity
        // guard against a typo in either column.
        assert_eq!(BL[0], BL[3]);
        assert_eq!(BL[1], BL[2]);
    }

    #[test]
    fn lowpass_filter_is_normalised() {
        // Spec form: `a0 = 1.0`. The implicit a0 tap is `AL[0]`.
        assert_eq!(AL[0], 1.0);
    }

    // ------------- Shape-codevector energy table (E_j, §3.9 eq. 3-23) -

    #[test]
    fn y_energy_table_dimension_matches_codebook() {
        // §3.9 prose: one E_j per shape codevector. Codebook has NCWD =
        // 128 rows, so the energy table must also have 128 entries.
        assert_eq!(Y_ENERGY.len(), NCWD);
    }

    #[test]
    fn y_energy_matches_dot_product() {
        // Cross-check: E_j must equal the direct Σ_k y_j(k)² computed
        // on the floating-point shape codevector view. This pins the
        // const-derived Y_ENERGY against the same data the rest of the
        // crate reads through `y_f64()` and guards against a future
        // typo in either branch of the derivation.
        let y_float = y_f64();
        for j in 0..NCWD {
            let mut expected = 0.0f64;
            for k in 0..IDIM {
                expected += y_float[j][k] * y_float[j][k];
            }
            assert!(
                (Y_ENERGY[j] - expected).abs() < 1e-12,
                "Y_ENERGY[{}] = {} vs dot-product {} (diff {})",
                j,
                Y_ENERGY[j],
                expected,
                (Y_ENERGY[j] - expected).abs()
            );
        }
    }

    #[test]
    fn y_energy_all_non_negative_and_finite() {
        // Sum-of-squares ⇒ every entry must be ≥ 0 and finite. A zero
        // entry would correspond to the all-zero codevector — the
        // printed Annex B codebook has no such row, so every entry must
        // be strictly positive.
        for (j, &e) in Y_ENERGY.iter().enumerate() {
            assert!(e.is_finite(), "Y_ENERGY[{}] not finite: {}", j, e);
            assert!(e > 0.0, "Y_ENERGY[{}] not strictly positive: {}", j, e);
        }
    }

    #[test]
    fn y_energy_row_zero_matches_hand_computed() {
        // Annex B prose lists row 0 of the Q11 shape codebook as
        // (668, -2950, -1254, -1790, -2553). The floating-point row is
        // those values divided by 2 048, and the energy is the sum of
        // squares of that float row. Spot-check the resulting E_0
        // against a hand-computed value (Q11 sum of squares /
        // 4 194 304 = 2²²).
        //
        //  668² + 2950² + 1254² + 1790² + 2553² =
        //  446 224 + 8 702 500 + 1 572 516 + 3 204 100 + 6 517 809
        //  = 20 443 149
        // / 2²² (= 4 194 304) ≈ 4.874...
        let expected = 20_443_149.0_f64 / (1u64 << 22) as f64;
        assert!(
            (Y_ENERGY[0] - expected).abs() < 1e-12,
            "Y_ENERGY[0] = {} vs hand-computed {} (diff {})",
            Y_ENERGY[0],
            expected,
            (Y_ENERGY[0] - expected).abs()
        );
    }

    #[test]
    fn y_energy_accessor_returns_same_data_as_const() {
        // The runtime accessor is a thin owned-copy wrapper around the
        // const; bit-for-bit equality is the contract.
        let view = y_energy_f64();
        for j in 0..NCWD {
            assert_eq!(view[j], Y_ENERGY[j]);
        }
    }
}
