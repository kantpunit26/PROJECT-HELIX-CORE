"""
HELIX-CORE Domain Models & Data Transfer Objects.

Defines all core dataclasses representing biological constraints, data fragments,
masked chunks, physical DNA oligo strands, and validation reports.
"""

from dataclasses import dataclass, field
from typing import Optional

from .config import (
    DEFAULT_FORWARD_PRIMER,
    DEFAULT_REVERSE_PRIMER,
    MAX_HOMOPOLYMER_ALLOWED,
    MAX_GC_CONTENT,
    MAX_STRAND_LENGTH,
    MIN_GC_CONTENT,
    MIN_STRAND_LENGTH,
)


@dataclass(frozen=True)
class BioConstraints:
    """Biophysical synthesis and sequencing constraints for DNA oligos."""

    max_homopolymer: int = MAX_HOMOPOLYMER_ALLOWED
    min_gc: float = MIN_GC_CONTENT
    max_gc: float = MAX_GC_CONTENT
    min_length: int = MIN_STRAND_LENGTH
    max_length: int = MAX_STRAND_LENGTH

    def validate_values(self) -> None:
        """Validate logical consistency of constraints."""
        if not (0.0 <= self.min_gc <= self.max_gc <= 1.0):
            raise ValueError(f"Invalid GC range: [{self.min_gc}, {self.max_gc}]")
        if self.max_homopolymer < 1:
            raise ValueError("max_homopolymer must be at least 1")
        if self.min_length > self.max_length:
            raise ValueError("min_length cannot exceed max_length")


@dataclass
class ValidationResult:
    """Results of biophysical validation on a DNA sequence."""

    is_valid: bool
    sequence: str
    length: int
    gc_content: float
    max_homopolymer: int
    homopolymer_violations: list[tuple[str, int, int]] = field(default_factory=list)
    gc_violation: bool = False
    length_violation: bool = False
    error_messages: list[str] = field(default_factory=list)

    def summary(self) -> str:
        """Human-readable validation summary."""
        status = "PASSED" if self.is_valid else "FAILED"
        return (
            f"Validation [{status}] Length={self.length} bp, "
            f"GC={self.gc_content:.2%}, MaxHomopolymer={self.max_homopolymer}. "
            f"Errors: {', '.join(self.error_messages) if self.error_messages else 'None'}"
        )


@dataclass
class Chunk:
    """Represents a fragmented payload before obfuscation."""

    chunk_id: int          # 32-bit unique chunk ID / seed
    chunk_index: int       # Logical sequential index (0..K-1 or K..K+M-1)
    payload: bytes         # Raw byte payload (fixed size)
    is_parity: bool = False
    salt: int = 0          # Candidate salt for biological constraint compliance
    block_index: int = 0   # RS block index for multi-block payloads

    @property
    def payload_size(self) -> int:
        return len(self.payload)


@dataclass
class MaskedChunk:
    """Represents a chunk after seeded XOR bitmasking and CRC32 injection."""

    chunk_id: int          # 32-bit unique identifier / seed
    chunk_index: int       # Logical sequential index
    masked_payload: bytes  # Obfuscated payload bytes
    crc32: int             # 32-bit unsigned CRC checksum
    is_parity: bool = False
    salt: int = 0
    block_index: int = 0

    @property
    def payload_size(self) -> int:
        return len(self.masked_payload)


@dataclass
class DnaStrand:
    """Represents a synthesized DNA oligo ready for physical wet-lab handling."""

    chunk_id: int
    chunk_index: int
    forward_primer: str
    id_sequence: str
    payload_sequence: str
    crc_sequence: str
    reverse_primer: str
    full_sequence: str
    gc_content: float
    max_homopolymer: int
    is_parity: bool = False

    @property
    def length(self) -> int:
        return len(self.full_sequence)


@dataclass
class PipelineStats:
    """End-to-end metrics for the DNA storage pipeline execution."""

    original_filename: str
    original_size_bytes: int
    original_sha256: str
    total_binary_bits: int
    data_chunk_count: int
    parity_chunk_count: int
    total_strand_count: int
    payload_size_bytes: int
    strand_length_bp: int
    mean_gc_content: float
    min_gc_content: float
    max_gc_content: float
    max_homopolymer_observed: int
    recovered_sha256: Optional[str] = None
    is_bit_perfect: bool = False

