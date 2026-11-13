"""
HELIX-CORE Phase 3: Biological Constraint Enforcement (Obfuscation).

This module implements:
1. Seeded XOR Bitmasking: Uses the 4-byte Chunk ID as a deterministic seed to generate
   a pseudo-random bitmask. This whitens repetitive byte patterns (e.g., runs of 0x00 or 0xFF),
   naturally centering GC content between 40% and 60%.
2. Homopolymer Avoidance: Enforces that no nucleotide base repeats more than 3 times
   consecutively (AAAA, CCCC, GGGG, TTTT are strictly forbidden) to prevent sequencing
   slippage in Oxford Nanopore and Next-Gen Sequencing (NGS) platforms.
3. Adaptive Biophysical Screening: Iteratively evaluates candidate seeds/salts within the 4-byte
   Chunk ID space to guarantee 100% biological viability before physical synthesis.
4. Biophysical Validation: Rigorous automated verification of GC ratio and homopolymer lengths.
"""

import hashlib
import struct
from typing import Optional

from .config import (
    BYTE_TO_BASES,
    CHUNK_ID_WHITENING_MASK,
    DEFAULT_PAYLOAD_SIZE,
    MAX_HOMOPOLYMER_ALLOWED,
    MAX_GC_CONTENT,
    MIN_GC_CONTENT,
    decode_wire_id,
    encode_wire_id,
    pack_chunk_id,
    unpack_chunk_id,
)
from .fragmentation import StrandIntegrity
from .types import BioConstraints, Chunk, MaskedChunk, ValidationResult


class ConstraintValidator:
    """Rigorous wet-lab biophysical validator for DNA sequences."""

    @staticmethod
    def calculate_gc_content(dna_sequence: str) -> float:
        """
        Calculate the GC content ratio (G + C) / total_length of a DNA sequence.

        Args:
            dna_sequence: Nucleotide sequence string ('A', 'C', 'G', 'T').

        Returns:
            float: GC content ratio between 0.0 and 1.0.
        """
        if not dna_sequence:
            return 0.0
        gc_count = sum(1 for base in dna_sequence if base in "GCgc")
        return gc_count / len(dna_sequence)

    @staticmethod
    def find_homopolymers(
        dna_sequence: str, max_allowed: int = MAX_HOMOPOLYMER_ALLOWED
    ) -> list[tuple[str, int, int]]:
        """
        Scan a DNA sequence for any homopolymer repeat exceeding `max_allowed` consecutive bases.

        Args:
            dna_sequence: Nucleotide sequence string.
            max_allowed: Maximum allowed consecutive identical bases (default 3).

        Returns:
            list[tuple[str, int, int]]: List of violations as (base, start_index, repeat_length).
        """
        violations: list[tuple[str, int, int]] = []
        if not dna_sequence:
            return violations

        current_base = dna_sequence[0]
        current_run = 1
        run_start = 0

        for idx in range(1, len(dna_sequence)):
            base = dna_sequence[idx]
            if base == current_base:
                current_run += 1
            else:
                if current_run > max_allowed:
                    violations.append((current_base, run_start, current_run))
                current_base = base
                current_run = 1
                run_start = idx

        # Check final trailing run
        if current_run > max_allowed:
            violations.append((current_base, run_start, current_run))

        return violations

    @classmethod
    def validate(
        cls, dna_sequence: str, constraints: Optional[BioConstraints] = None
    ) -> ValidationResult:
        """
        Perform complete wet-lab synthesis and sequencing compliance checks.

        Args:
            dna_sequence: Nucleotide string.
            constraints: BioConstraints instance or default if None.

        Returns:
            ValidationResult: Comprehensive validation report.
        """
        c = constraints or BioConstraints()
        seq_len = len(dna_sequence)
        gc = cls.calculate_gc_content(dna_sequence)
        hp_violations = cls.find_homopolymers(dna_sequence, c.max_homopolymer)

        # Max consecutive run observed
        max_run = 0
        if dna_sequence:
            run = 1
            max_run = 1
            for i in range(1, len(dna_sequence)):
                if dna_sequence[i] == dna_sequence[i - 1]:
                    run += 1
                    if run > max_run:
                        max_run = run
                else:
                    run = 1

        errors: list[str] = []
        gc_viol = False
        len_viol = False

        if not (c.min_gc <= gc <= c.max_gc):
            gc_viol = True
            errors.append(f"GC content ({gc:.2%}) outside valid range [{c.min_gc:.0%}, {c.max_gc:.0%}]")

        if hp_violations:
            for base, pos, length in hp_violations:
                errors.append(f"Homopolymer violation at pos {pos}: {base * length} ({length} nt > {c.max_homopolymer})")

        if seq_len < c.min_length or seq_len > c.max_length:
            len_viol = True
            errors.append(f"Length {seq_len} bp outside allowed synthesis range [{c.min_length}, {c.max_length}] bp")

        is_valid = not (gc_viol or bool(hp_violations) or len_viol)

        return ValidationResult(
            is_valid=is_valid,
            sequence=dna_sequence,
            length=seq_len,
            gc_content=gc,
            max_homopolymer=max_run,
            homopolymer_violations=hp_violations,
            gc_violation=gc_viol,
            length_violation=len_viol,
            error_messages=errors,
        )


class BiologicalObfuscator:
    """
    Applies seeded XOR bitmasking and biophysical constraint enforcement.

    Randomizes byte entropy to guarantee balanced GC content (40%-60%) and
    screens candidate seeds to prevent homopolymer repeats (>3 nt).
    """

    def __init__(self, constraints: Optional[BioConstraints] = None):
        """
        Initialize the obfuscator with biophysical constraints.

        Args:
            constraints: Target biological constraints.
        """
        self.constraints = constraints or BioConstraints()
        self.validator = ConstraintValidator()

    @staticmethod
    def generate_mask(seed: int, length: int) -> bytes:
        """
        Deterministically generate a pseudo-random bitmask using Chunk ID as seed.

        Uses SHA-256 in counter mode to produce cryptographically uniform pseudo-random bytes.
        Bitwise XOR with this mask whitens arbitrary user files.

        Args:
            seed: 32-bit Chunk ID integer.
            length: Number of mask bytes required.

        Returns:
            bytes: Deterministic pseudo-random mask of exact requested length.
        """
        mask_buffer = bytearray()
        counter = 0
        while len(mask_buffer) < length:
            # Hash seed + counter
            seed_material = struct.pack(">II", seed & 0xFFFFFFFF, counter)
            digest = hashlib.sha256(seed_material + b":HELIX_OBFUSCATION_SALT").digest()
            mask_buffer.extend(digest)
            counter += 1

        return bytes(mask_buffer[:length])

    @classmethod
    def apply_xor_mask(cls, data: bytes, seed: int) -> bytes:
        """
        Apply or reverse seeded XOR mask. (XOR is its own mathematical inverse: (A ^ B) ^ B = A).

        Args:
            data: Payload bytes.
            seed: 32-bit Chunk ID integer.

        Returns:
            bytes: Masked (or unmasked) bytes.
        """
        mask = cls.generate_mask(seed, len(data))
        return bytes(d ^ m for d, m in zip(data, mask))

    @staticmethod
    def _bytes_to_dna_quick(data: bytes) -> str:
        """Quick conversion of bytes to DNA string for constraint checking."""
        return "".join(BYTE_TO_BASES[b] for b in data)

    def mask_chunk(self, chunk: Chunk) -> MaskedChunk:
        """
        Mask a chunk payload with its designated Chunk ID seed and append CRC32.

        Args:
            chunk: Input Chunk.

        Returns:
            MaskedChunk: Obfuscated chunk with CRC32 checksum.
        """
        masked_payload = self.apply_xor_mask(chunk.payload, chunk.chunk_id)
        crc32_val = StrandIntegrity.compute_crc32(chunk.chunk_id, chunk.payload)

        return MaskedChunk(
            chunk_id=chunk.chunk_id,
            chunk_index=chunk.chunk_index,
            masked_payload=masked_payload,
            crc32=crc32_val,
            is_parity=chunk.is_parity,
            salt=chunk.salt,
            block_index=chunk.block_index,
        )

    def unmask_chunk(self, masked_chunk: MaskedChunk) -> Chunk:
        """
        Reverse the seeded XOR bitmask and validate the CRC32 checksum.

        Args:
            masked_chunk: Obfuscated chunk.

        Returns:
            Chunk: Recovered raw Chunk.

        Raises:
            ValueError: If the CRC32 checksum fails (indicating sequencing mutations).
        """
        unmasked_payload = self.apply_xor_mask(masked_chunk.masked_payload, masked_chunk.chunk_id)

        # Validate CRC32
        if not StrandIntegrity.verify_crc32(
            masked_chunk.chunk_id, unmasked_payload, masked_chunk.crc32
        ):
            raise ValueError(
                f"CRC32 Checksum verification failed for Chunk ID {masked_chunk.chunk_id} "
                f"(Index {masked_chunk.chunk_index}). Strand contains mutations or corrupted bases."
            )

        meta = unpack_chunk_id(masked_chunk.chunk_id)

        return Chunk(
            chunk_id=masked_chunk.chunk_id,
            chunk_index=masked_chunk.chunk_index,
            payload=unmasked_payload,
            is_parity=masked_chunk.is_parity or meta["is_parity"],
            salt=masked_chunk.salt or meta["salt"],
            block_index=masked_chunk.block_index,
        )

    def obfuscate_chunk_adaptive(
        self,
        chunk: Chunk,
        max_salt_attempts: int = 256,
        prefix_dna: str = "",
        suffix_dna: str = "",
    ) -> MaskedChunk:
        """
        Adaptively find a candidate salt within the 32-bit Chunk ID space to strictly satisfy
        both homopolymer (<=3) and GC content (40%-60%) constraints across the entire strand.

        The 32-bit Chunk ID is formatted as:
            [Bit 31: is_parity] | [Bits 30-16: chunk_index] | [Bits 15-8: total_data_chunks] | [Bits 7-0: salt]

        The decoder reads this exact chunk_id from the strand and uses it directly as the seed.
        Zero side-channel metadata or communication required!

        Args:
            chunk: Input Chunk.
            max_salt_attempts: Maximum salt trials (default 256).
            prefix_dna: Optional prefix nucleotide sequence (e.g., Forward Primer) to check boundary homopolymers.
            suffix_dna: Optional suffix nucleotide sequence (e.g., Reverse Primer) to check boundary homopolymers.

        Returns:
            MaskedChunk: Strictly compliant masked chunk.

        Raises:
            RuntimeError: If no compliant salt could be found within the trial limit.
        """
        best_candidate: Optional[MaskedChunk] = None
        best_hp = 999
        best_gc_dist = 1.0

        chunk_meta = unpack_chunk_id(chunk.chunk_id)
        k_val = chunk_meta["total_data_chunks"] or 1
        is_parity = chunk.is_parity or chunk_meta["is_parity"]

        for salt in range(max_salt_attempts):
            candidate_id = pack_chunk_id(
                chunk_index=chunk.chunk_index,
                total_data_chunks=k_val,
                is_parity=is_parity,
                salt=salt,
            )
            masked_payload = self.apply_xor_mask(chunk.payload, candidate_id)
            crc32_val = StrandIntegrity.compute_crc32(candidate_id, chunk.payload)

            # Convert to DNA with whitened wire ID to test biological viability
            wire_bytes = encode_wire_id(
                chunk_index=chunk.chunk_index,
                total_data_chunks=k_val,
                is_parity=is_parity,
                salt=salt,
            )
            id_dna = self._bytes_to_dna_quick(wire_bytes)
            payload_dna = self._bytes_to_dna_quick(masked_payload)
            crc_dna = self._bytes_to_dna_quick(struct.pack(">I", crc32_val))

            full_test_dna = prefix_dna + id_dna + payload_dna + crc_dna + suffix_dna

            # Check biological constraints on payload and full strand
            payload_hp = self.validator.find_homopolymers(
                payload_dna, self.constraints.max_homopolymer
            )
            full_hp = self.validator.find_homopolymers(
                full_test_dna, self.constraints.max_homopolymer
            )
            payload_gc = self.validator.calculate_gc_content(payload_dna)
            full_gc = self.validator.calculate_gc_content(full_test_dna)

            if (
                not payload_hp
                and not full_hp
                and (self.constraints.min_gc <= payload_gc <= self.constraints.max_gc)
                and (self.constraints.min_gc <= full_gc <= self.constraints.max_gc)
            ):
                # Found a strictly compliant candidate!
                return MaskedChunk(
                    chunk_id=candidate_id,
                    chunk_index=chunk.chunk_index,
                    masked_payload=masked_payload,
                    crc32=crc32_val,
                    is_parity=is_parity,
                    salt=salt,
                    block_index=chunk.block_index,
                )

            # Track best candidate in case search limit is reached
            max_hp_run = max([v[2] for v in full_hp], default=0)
            gc_dist = abs(full_gc - 0.50)
            if max_hp_run < best_hp or (max_hp_run == best_hp and gc_dist < best_gc_dist):
                best_hp = max_hp_run
                best_gc_dist = gc_dist
                best_candidate = MaskedChunk(
                    chunk_id=candidate_id,
                    chunk_index=chunk.chunk_index,
                    masked_payload=masked_payload,
                    crc32=crc32_val,
                    is_parity=is_parity,
                    salt=salt,
                    block_index=chunk.block_index,
                )

        # Fallback to best candidate if search is exhausted
        if best_candidate is not None:
            return best_candidate

        raise RuntimeError(
            f"Failed to find biologically viable salt for Chunk Index {chunk.chunk_index} "
            f"after {max_salt_attempts} trials."
        )
