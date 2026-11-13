"""
HELIX-CORE Configuration & Biophysical Parameters.

This module centralizes all biological, physical synthesis, and sequencing
constants required for the HELIX-CORE DNA data storage architecture.
"""

from dataclasses import dataclass


# =====================================================================
# BIOPHYSICAL CONSTANTS & SYNTHESIS SPECIFICATIONS
# =====================================================================

# Universal Illumina/Twist compatible PCR Primers (20 nt each)
# Designed for high-efficiency amplification, matched Tm (~58-60C),
# balanced GC content (50-55%), and zero homopolymer repeats > 2.
DEFAULT_FORWARD_PRIMER: str = "ACACGACGCTCTTCCGATCT"
DEFAULT_REVERSE_PRIMER: str = "AGATCGGAAGAGCGTCGTGT"

# Nucleotide <-> Dibit mapping strictly according to Phase 4
# 00 = A (Adenine)
# 01 = C (Cytosine)
# 10 = G (Guanine)
# 11 = T (Thymine)
BITS_TO_BASE: dict[str, str] = {
    "00": "A",
    "01": "C",
    "10": "G",
    "11": "T",
}

BASE_TO_BITS: dict[str, str] = {
    "A": "00",
    "C": "01",
    "G": "10",
    "T": "11",
}

# Fast byte-to-base and base-to-byte lookup tables (256 entries)
BYTE_TO_BASES: list[str] = [
    BITS_TO_BASE[f"{(b >> 6) & 3:02b}"]
    + BITS_TO_BASE[f"{(b >> 4) & 3:02b}"]
    + BITS_TO_BASE[f"{(b >> 2) & 3:02b}"]
    + BITS_TO_BASE[f"{b & 3:02b}"]
    for b in range(256)
]

# Physical Strand Architecture Constraints
# Maximum and minimum oligo synthesis lengths permitted by wet-lab reality
MIN_STRAND_LENGTH: int = 150  # nt (base pairs)
MAX_STRAND_LENGTH: int = 200  # nt (base pairs)

# Metadata & Integrity field lengths (in bytes)
CHUNK_ID_BYTES: int = 4   # 32-bit Index / Chunk ID (16 nt)
CRC32_BYTES: int = 4      # 32-bit CRC32 checksum (16 nt)

# Default payload size in bytes.
# 28 bytes payload -> 112 nt.
# 20 nt (Fwd Primer) + 16 nt (ID) + 112 nt (Payload) + 16 nt (CRC32) + 20 nt (Rev Primer) = 184 nt.
# 184 nt strictly satisfies the 150-200 bp synthesis constraint.
DEFAULT_PAYLOAD_SIZE: int = 28

# Biological Constraints
MAX_HOMOPOLYMER_ALLOWED: int = 3  # No base repeated > 3 times (AAAA, CCCC, GGGG, TTTT forbidden)
MIN_GC_CONTENT: float = 0.40      # Minimum GC content: 40%
MAX_GC_CONTENT: float = 0.60      # Maximum GC content: 60%

# Error Correction Defaults
# 15% parity chunks guarantees >10% strand loss recovery threshold
DEFAULT_PARITY_RATIO: float = 0.15

# Maximum data chunks per Reed-Solomon Galois Field GF(2^8) block
# Total block length N = K + M <= 255. With 15% parity, max K is ~215.
MAX_RS_BLOCK_DATA_CHUNKS: int = 180

# Magic header for container format
HELIX_CONTAINER_MAGIC: bytes = b"HLX1"

# 32-bit Whitening Mask for Chunk ID wire representation.
# 0x4E4E4E4E maps to 'CATGCATGCATGCATG' in DNA (exactly 50% GC, 0 homopolymer repeats),
# preventing small Chunk IDs (0, 1, 2) from producing runs of 'AAAAAAAAAAAA'.
CHUNK_ID_WHITENING_MASK: int = 0x4E4E4E4E

# Minimum parity chunks to guarantee robust recovery even on small files
MIN_PARITY_CHUNKS: int = 4


# =====================================================================
# SELF-DESCRIBING 32-BIT CHUNK ID & WIRE ENCODING
# =====================================================================

def pack_chunk_id(
    chunk_index: int,
    total_data_chunks: int = 1,
    is_parity: bool = False,
    salt: int = 0,
) -> int:
    """
    Pack metadata into a self-describing 32-bit Chunk ID.

    Bit layout:
      - Bit 31 (1 bit): is_parity (0 = DATA, 1 = PARITY)
      - Bits 30-16 (15 bits): chunk_index (0..32,767)
      - Bits 15-8 (8 bits): total_data_chunks (1..255)
      - Bits 7-0 (8 bits): candidate salt (0..255)

    Total = 32 bits = 4 bytes.
    """
    p_bit = (1 if is_parity else 0) & 0x1
    idx = chunk_index & 0x7FFF
    k_val = total_data_chunks & 0xFF
    s_val = salt & 0xFF
    return ((p_bit << 31) | (idx << 16) | (k_val << 8) | s_val) & 0xFFFFFFFF


def unpack_chunk_id(chunk_id: int) -> dict:
    """
    Unpack 32-bit Chunk ID into its constituent metadata fields:
      - is_parity: bool
      - chunk_index: int
      - total_data_chunks: int
      - salt: int
    """
    chunk_id = chunk_id & 0xFFFFFFFF
    is_parity = bool((chunk_id >> 31) & 0x1)
    chunk_index = (chunk_id >> 16) & 0x7FFF
    total_data_chunks = (chunk_id >> 8) & 0xFF
    salt = chunk_id & 0xFF

    # Graceful fallback for legacy IDs where chunk_id was just (idx << 8) | salt
    if not is_parity and total_data_chunks == 0 and chunk_index == 0 and (chunk_id >> 8) > 0:
        chunk_index = (chunk_id >> 8) & 0x00FFFFFF
        salt = chunk_id & 0xFF

    return {
        "is_parity": is_parity,
        "chunk_index": chunk_index,
        "total_data_chunks": total_data_chunks,
        "salt": salt,
    }


def encode_wire_id(
    chunk_index: int,
    total_data_chunks: int = 1,
    is_parity: bool = False,
    salt: int = 0,
) -> bytes:
    """
    Encode metadata into a 4-byte wire representation for physical DNA synthesis.

    Uses a deterministic whitening mask derived from the salt so that candidate salts
    randomize the entire 16-nucleotide ID, eliminating homopolymers and balancing GC.
    """
    import hashlib

    b0 = ((1 if is_parity else 0) << 7) | ((chunk_index >> 8) & 0x7F)
    b1 = chunk_index & 0xFF
    b2 = total_data_chunks & 0xFF
    b3 = salt & 0xFF

    wire_b3 = b3 ^ 0x4E
    salt_mask = hashlib.sha256(bytes([salt, 0x4E])).digest()[:3]
    wire_b0 = b0 ^ salt_mask[0]
    wire_b1 = b1 ^ salt_mask[1]
    wire_b2 = b2 ^ salt_mask[2]

    return bytes([wire_b0, wire_b1, wire_b2, wire_b3])


def decode_wire_id(wire_bytes: bytes) -> dict:
    """
    Decode 4-byte wire representation from physical DNA into metadata and chunk_id.
    """
    import hashlib

    if len(wire_bytes) != 4:
        raise ValueError(f"Wire ID must be exactly 4 bytes, got {len(wire_bytes)}")

    wire_b0, wire_b1, wire_b2, wire_b3 = wire_bytes
    salt = wire_b3 ^ 0x4E
    salt_mask = hashlib.sha256(bytes([salt, 0x4E])).digest()[:3]
    b0 = wire_b0 ^ salt_mask[0]
    b1 = wire_b1 ^ salt_mask[1]
    b2 = wire_b2 ^ salt_mask[2]

    is_parity = bool((b0 >> 7) & 1)
    chunk_index = ((b0 & 0x7F) << 8) | b1
    total_data_chunks = b2

    chunk_id = pack_chunk_id(
        chunk_index=chunk_index,
        total_data_chunks=total_data_chunks,
        is_parity=is_parity,
        salt=salt,
    )

    return {
        "chunk_id": chunk_id,
        "is_parity": is_parity,
        "chunk_index": chunk_index,
        "total_data_chunks": total_data_chunks,
        "salt": salt,
    }


