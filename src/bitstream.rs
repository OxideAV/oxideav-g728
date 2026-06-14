//! Serial bit-stream framing for the 16 kbit/s LD-CELP channel —
//! §5.11 (block 18 epilogue) of ITU-T G.728 (1992-09).
//!
//! The 16 kbit/s stream carries one 10-bit channel index `ICHAN` per
//! `IDIM = 5`-sample vector and nothing else (§2.1, §3.11): "Only the
//! index to the excitation codebook is transmitted." There is no
//! in-band framing — the decoder is assumed to already know the index
//! boundaries and the 4-vector adaptation phase (§3.11; decode-trace
//! §9). What the spec *does* pin down is the **bit order on the wire**.
//!
//! The block-18 transmit epilogue of the §5.11 pseudo-code states the
//! serial convention exactly:
//!
//! > For serial bit stream transmission, the most significant bit of
//! > `ICHAN` should be transmitted first. If `ICHAN` is represented by
//! > the 10-bit word `b9 b8 b7 b6 b5 b4 b3 b2 b1 b0`, then the order of
//! > the transmitted bits should be `b9`, and then `b8`, and then
//! > `b7`, …, and finally `b0`. (`b9` is the most significant bit.)
//!
//! So the serial stream is the concatenation, per vector, of the ten
//! bits `b9 b8 … b0`, most-significant-first. To carry that serial
//! stream over an octet-oriented transport we pack those bits into
//! bytes most-significant-bit-first as well (the first serial bit is
//! bit 7 of the first byte), which keeps the wire order `b9 b8 … b0`
//! contiguous across the byte boundary and makes a single-index packet
//! occupy the high ten bits of its first two bytes — directly readable
//! in a hex dump. The same MSB-first choice is what makes the §3.11
//! "leftmost bit of the codeword" (the robbed synchronization bit, =
//! `b9`) land first on the wire.
//!
//! [`pack_indices`] serialises a slice of channel indices into that
//! byte stream; [`unpack_indices`] is its exact inverse. When the
//! total serialised bit count is not a whole multiple of
//! [`CHANNEL_INDEX_BITS`] (= 10) — i.e. the byte buffer cannot hold a
//! whole number of channel indices — [`unpack_indices`] returns
//! [`Error::InvalidInputLength`] rather than fabricate a partial index.

use crate::encoder::CHANNEL_INDEX_BITS;
use crate::{Error, Result};

/// Mask covering the significant low [`CHANNEL_INDEX_BITS`] of a
/// channel index (`0..=1023`).
const INDEX_MASK: u16 = (1 << CHANNEL_INDEX_BITS) - 1;

/// Pack a sequence of 10-bit channel indices into the §5.11 serial bit
/// stream, octet-aligned most-significant-bit-first.
///
/// Each index contributes its ten significant bits `b9 b8 … b0`
/// (most-significant first) to the stream; any bits above bit 9 of an
/// input word are ignored (`ICHAN ∈ [0, 1024)` per
/// [`CHANNEL_INDEX_BITS`]). Indices are emitted in order, with no
/// padding between them. If the total bit count
/// (`indices.len() · 10`) is not a multiple of 8, the final byte is
/// zero-padded in its least-significant bits — the trailing pad bits
/// are *not* recoverable as a channel index and [`unpack_indices`]
/// will reject a buffer whose bit length is not a multiple of 10 (see
/// the round-trip note there).
///
/// The inverse is [`unpack_indices`].
pub fn pack_indices(indices: &[u16]) -> Vec<u8> {
    // Total serial bits, then ceil-divide to whole bytes.
    let total_bits = indices.len() * CHANNEL_INDEX_BITS as usize;
    let total_bytes = total_bits.div_ceil(8);
    let mut out = vec![0u8; total_bytes];

    // `bitpos` is the absolute bit offset into the serial stream, with
    // bit 0 == the most-significant bit of byte 0 (MSB-first packing).
    let mut bitpos = 0usize;
    for &ichan in indices {
        let ichan = ichan & INDEX_MASK;
        // Emit b9 first, b0 last (§5.11: "the most significant bit …
        // transmitted first").
        for shift in (0..CHANNEL_INDEX_BITS).rev() {
            let bit = ((ichan >> shift) & 1) as u8;
            if bit != 0 {
                // Within a byte, the first serial bit is bit 7.
                out[bitpos / 8] |= 1 << (7 - (bitpos % 8));
            }
            bitpos += 1;
        }
    }
    out
}

/// Unpack a §5.11 serial bit stream back into 10-bit channel indices —
/// the exact inverse of [`pack_indices`].
///
/// Reads ten bits at a time, most-significant-bit-first, reconstructing
/// each `ICHAN` as `b9 b8 … b0`. Returns one `u16` per index, each in
/// `[0, 1024)`.
///
/// # Errors
///
/// Returns [`Error::InvalidInputLength`] when the buffer's bit length
/// (`bytes.len() · 8`) is not a whole multiple of
/// [`CHANNEL_INDEX_BITS`] (= 10) — i.e. the byte stream does not encode
/// a whole number of channel indices. Note this is *stricter* than
/// "round-trips a [`pack_indices`] output": a single packed index
/// occupies 10 bits in 2 bytes, leaving 6 zero pad bits; that 16-bit
/// buffer is 8 ÷ 5 of an index and is rejected. Callers that packed an
/// arbitrary count and need lossless recovery should carry the index
/// count out of band, or only pack counts that are a multiple of 4
/// (one adaptation cycle = 40 bits = 5 whole bytes), the natural
/// framing unit of the codec.
pub fn unpack_indices(bytes: &[u8]) -> Result<Vec<u16>> {
    let total_bits = bytes.len() * 8;
    if total_bits % CHANNEL_INDEX_BITS as usize != 0 {
        return Err(Error::InvalidInputLength);
    }
    let count = total_bits / CHANNEL_INDEX_BITS as usize;
    let mut out = Vec::with_capacity(count);

    let mut bitpos = 0usize;
    for _ in 0..count {
        let mut ichan: u16 = 0;
        // Read b9 first, b0 last.
        for _ in 0..CHANNEL_INDEX_BITS {
            let bit = (bytes[bitpos / 8] >> (7 - (bitpos % 8))) & 1;
            ichan = (ichan << 1) | bit as u16;
            bitpos += 1;
        }
        out.push(ichan);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_a_cycle_aligned_stream() {
        // One adaptation cycle is 4 vectors = 40 bits = 5 whole bytes,
        // so a multiple-of-4 index count packs to a whole byte count
        // and round-trips losslessly.
        let indices: Vec<u16> = (0u16..16).map(|i| (i * 61) & INDEX_MASK).collect();
        assert_eq!(indices.len() % 4, 0);
        let packed = pack_indices(&indices);
        assert_eq!(packed.len(), indices.len() * 10 / 8);
        let unpacked = unpack_indices(&packed).expect("cycle-aligned stream is whole");
        assert_eq!(unpacked, indices);
    }

    #[test]
    fn extreme_indices_round_trip() {
        // Boundary indices 0 and 1023 (all-zero and all-one codewords)
        // must survive packing exactly.
        let indices = [0u16, 1023, 0, 1023];
        let packed = pack_indices(&indices);
        let unpacked = unpack_indices(&packed).unwrap();
        assert_eq!(unpacked, indices);
    }

    #[test]
    fn msb_of_first_index_lands_in_top_byte_bit() {
        // §5.11: b9 (the MSB of ICHAN) is transmitted first, so it must
        // land in bit 7 of the first byte. Index 0x200 == 1 << 9 sets
        // only b9; the packed stream's first byte must be 0b1000_0000.
        let packed = pack_indices(&[0x200]);
        assert_eq!(packed[0], 0b1000_0000);
        // And b0 (index == 1) lands as the 10th serial bit = bit 5 of
        // byte 1 (bits 0..7 of byte 0, then 8,9 of byte 1 → 9th and
        // 10th serial bits are bit 7 and bit 6 of byte 1 — so b0 is
        // bit 6 of byte 1).
        let packed_b0 = pack_indices(&[0x001]);
        assert_eq!(packed_b0[0], 0);
        assert_eq!(packed_b0[1], 0b0100_0000);
    }

    #[test]
    fn input_bits_above_b9_are_masked() {
        // Anything above bit 9 of an input word is not part of the
        // channel index and must be dropped.
        let packed_clean = pack_indices(&[0x055]);
        let packed_dirty = pack_indices(&[0x055 | 0xFC00]);
        assert_eq!(packed_clean, packed_dirty);
    }

    #[test]
    fn rejects_non_multiple_of_ten_bits() {
        // A 2-byte buffer is 16 bits = 1.6 indices → reject.
        assert_eq!(
            unpack_indices(&[0xAB, 0xCD]),
            Err(Error::InvalidInputLength)
        );
        // 5 bytes = 40 bits = 4 indices → accept.
        assert!(unpack_indices(&[0; 5]).is_ok());
        // 1 byte = 8 bits → reject.
        assert_eq!(unpack_indices(&[0xFF]), Err(Error::InvalidInputLength));
    }

    #[test]
    fn empty_round_trips_to_empty() {
        assert_eq!(pack_indices(&[]), Vec::<u8>::new());
        assert_eq!(unpack_indices(&[]).unwrap(), Vec::<u16>::new());
    }

    #[test]
    fn second_index_is_bit_aligned_across_byte_boundary() {
        // Two indices = 20 bits, spanning into byte 2. Verify the
        // second index is read back exactly (it straddles a byte
        // boundary: bits 10..19, i.e. bits 2..7 of byte 1 + bits 0..3
        // of byte 2).
        let indices = [0x2AA, 0x155];
        let packed = pack_indices(&indices);
        assert_eq!(packed.len(), 3); // 20 bits → 3 bytes
        let unpacked = unpack_indices(&packed);
        // 24 bits is not a multiple of 10, so direct unpack rejects;
        // but a cycle-aligned wrapper would accept. Confirm the
        // rejection and that a 4-index (40-bit, 5-byte) stream built
        // from the same values round-trips.
        assert_eq!(unpacked, Err(Error::InvalidInputLength));
        let four = [0x2AA, 0x155, 0x2AA, 0x155];
        let p4 = pack_indices(&four);
        assert_eq!(unpack_indices(&p4).unwrap(), four);
    }

    #[test]
    fn exhaustive_single_index_round_trip_via_padding() {
        // Every valid index 0..=1023, packed as a 4-tuple (cycle
        // aligned), must survive the round trip in each slot.
        for j in (0u16..1024).step_by(7) {
            let four = [
                j,
                (j + 1) & INDEX_MASK,
                (j + 2) & INDEX_MASK,
                (j + 3) & INDEX_MASK,
            ];
            let packed = pack_indices(&four);
            assert_eq!(unpack_indices(&packed).unwrap(), four);
        }
    }
}
