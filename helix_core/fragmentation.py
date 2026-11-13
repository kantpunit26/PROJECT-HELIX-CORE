"""
HELIX-CORE Phase 2: Fragmentation & Reed-Solomon Error Correction.

This module provides:
1. Fragmentation of binary byte streams into fixed-size payloads respecting synthesis limits (150-200 bp).
2. Metadata injection: prepending 4-byte Index/Chunk IDs to allow reassembly from a free-floating liquid pool.
3. Reed-Solomon (RS) forward error correction: generates parity chunks so that if >=10% of synthesized
   strands are lost, degraded, or destroyed, 100% of the original data is recovered.
4. CRC32 strand-level checksumming to immediately detect mutations (SNPs, indels, substitutions).
"""

import math
import struct
from typing import Optional, Union
import zlib

import reedsolo

from .config import (
    CHUNK_ID_BYTES,
    CRC32_BYTES,
    DEFAULT_PARITY_RATIO,
    DEFAULT_PAYLOAD_SIZE,
    MAX_RS_BLOCK_DATA_CHUNKS,
    MIN_PARITY_CHUNKS,
    pack_chunk_id,
    unpack_chunk_id,
)
from .types import Chunk


class StrandIntegrity:
    """Computes and validates 32-bit CRC32 checksums for physical strand error detection."""

    @staticmethod
    def compute_crc32(chunk_id: int, payload: bytes) -> int:
        """
        Calculate a 32-bit CRC32 checksum over the Chunk ID and payload bytes.

        Args:
            chunk_id: 32-bit chunk identifier.
            payload: Chunk payload bytes.

        Returns:
            int: Unsigned 32-bit integer checksum.
        """
        id_bytes = struct.pack(">I", chunk_id)
        return zlib.crc32(id_bytes + payload) & 0xFFFFFFFF

    @classmethod
    def verify_crc32(cls, chunk_id: int, payload: bytes, expected_crc32: int) -> bool:
        """
        Verify whether the received payload matches the expected checksum.

        If sequencing introduces SNPs, insertions, deletions, or bit-flips,
        the CRC32 will fail and the strand must be discarded as an erasure.
        """
        actual_crc32 = cls.compute_crc32(chunk_id, payload)
        return actual_crc32 == expected_crc32


class FragmentationEngine:
    """Handles fragmentation of byte streams into fixed payloads and reassembly."""

    def __init__(self, payload_size: int = DEFAULT_PAYLOAD_SIZE):
        """
        Initialize the fragmentation engine.

        Args:
            payload_size: Payload length in bytes per fragment.
                          Default 28 bytes = 112 nt (bases), which yields a final
                          strand length of 184 nt (within 150-200 bp synthesis limit).
        """
        if payload_size < 4:
            raise ValueError(f"Payload size ({payload_size}) must be at least 4 bytes.")
        self.payload_size = payload_size

    def fragment(self, data: bytes) -> list[Chunk]:
        """
        Fragment raw bytes into fixed-size Chunk objects.

        If the final fragment is smaller than `payload_size`, it is zero-padded
        deterministically to preserve uniform column dimensions for Reed-Solomon coding.

        Args:
            data: Raw continuous byte stream.

        Returns:
            list[Chunk]: Ordered list of data chunks with sequential chunk indices.
        """
        if not data:
            raise ValueError("Cannot fragment empty data stream.")

        chunks: list[Chunk] = []
        total_len = len(data)
        num_chunks = math.ceil(total_len / self.payload_size)

        for i in range(num_chunks):
            start = i * self.payload_size
            end = min(start + self.payload_size, total_len)
            payload = data[start:end]

            # Zero-pad final chunk if necessary to maintain fixed column dimensions
            if len(payload) < self.payload_size:
                payload = payload.ljust(self.payload_size, b"\x00")

            # Pack self-describing 32-bit chunk_id (records index, total chunks, parity flag, and salt=0)
            chunk_id = pack_chunk_id(
                chunk_index=i,
                total_data_chunks=num_chunks,
                is_parity=False,
                salt=0,
            )
            chunk = Chunk(
                chunk_id=chunk_id,
                chunk_index=i,
                payload=payload,
                is_parity=False,
                block_index=i // MAX_RS_BLOCK_DATA_CHUNKS,
            )
            chunks.append(chunk)

        return chunks

    @staticmethod
    def reassemble(chunks: list[Chunk]) -> bytes:
        """
        Reassemble data chunks in order of their Index ID into a continuous byte array.

        Args:
            chunks: List of Chunk objects (must be data chunks, not parity chunks).

        Returns:
            bytes: Continuous reconstructed byte stream.
        """
        # Filter and sort strictly by chunk_index
        data_chunks = [c for c in chunks if not c.is_parity]
        data_chunks.sort(key=lambda c: c.chunk_index)

        # Verify contiguous indices
        for expected_idx, chunk in enumerate(data_chunks):
            if chunk.chunk_index != expected_idx:
                raise ValueError(
                    f"Missing chunk in reassembly: expected index {expected_idx}, got {chunk.chunk_index}"
                )

        return b"".join(c.payload for c in data_chunks)


class ReedSolomonCodec:
    """
    Inter-strand Reed-Solomon forward error correction codec.

    Operates across chunk payloads column-wise using Galois Field GF(2^8) arithmetic.
    Treats lost or corrupted strands (identified by missing IDs or failed CRC32) as
    erasures, enabling 100% mathematical recovery of missing data as long as the
    number of lost strands does not exceed the number of parity strands.
    """

    def __init__(self, parity_ratio: float = DEFAULT_PARITY_RATIO):
        """
        Initialize the Reed-Solomon codec.

        Args:
            parity_ratio: Ratio of parity chunks to generate (e.g. 0.15 = 15%).
                          Must be >= 0.10 to satisfy wet-lab 10% recovery constraint.
        """
        if parity_ratio < 0.10:
            raise ValueError(
                f"Parity ratio ({parity_ratio:.2%}) must be >= 10% to guarantee wet-lab recovery."
            )
        self.parity_ratio = parity_ratio

    def encode_block(self, data_chunks: list[Chunk], block_index: int = 0) -> list[Chunk]:
        """
        Generate Reed-Solomon parity chunks for a single block of data chunks.

        Uses column-wise encoding: for byte position j in [0..payload_size-1],
        the vector [C_0[j], C_1[j], ..., C_{K-1}[j]] is encoded to generate M parity bytes.

        Args:
            data_chunks: List of K data chunks in this block (K <= 200).
            block_index: Numerical index of this block.

        Returns:
            list[Chunk]: List of M parity chunks with indices starting at K.
        """
        k = len(data_chunks)
        if k == 0:
            return []

        payload_size = data_chunks[0].payload_size
        # Verify uniform payload sizes
        for c in data_chunks:
            if len(c.payload) != payload_size:
                raise ValueError("All chunks within an RS block must have identical payload sizes.")

        # Compute number of parity chunks M (at least MIN_PARITY_CHUNKS to ensure >=10% loss recovery)
        m = max(MIN_PARITY_CHUNKS, math.ceil(k * self.parity_ratio))

        # Check Galois Field GF(2^8) constraint: N = K + M <= 255
        if k + m > 255:
            raise ValueError(
                f"RS block size K={k} + M={m} = {k+m} exceeds GF(2^8) maximum block length 255. "
                "Reduce block size using multi-block partitioning."
            )

        rsc = reedsolo.RSCodec(m)
        parity_payloads = [bytearray(payload_size) for _ in range(m)]

        # Column-wise encoding across all K chunks
        for col in range(payload_size):
            message_col = bytes(data_chunks[row].payload[col] for row in range(k))
            encoded_col = rsc.encode(message_col)
            # Extracted parity bytes are located at encoded_col[k : k + m]
            for p in range(m):
                parity_payloads[p][col] = encoded_col[k + p]

        # Construct parity Chunk objects with self-describing 32-bit chunk_id
        base_index = data_chunks[-1].chunk_index + 1
        parity_chunks: list[Chunk] = []
        for p in range(m):
            parity_idx = base_index + p
            chunk_id = pack_chunk_id(
                chunk_index=parity_idx,
                total_data_chunks=k,
                is_parity=True,
                salt=0,
            )
            parity_chunk = Chunk(
                chunk_id=chunk_id,
                chunk_index=parity_idx,
                payload=bytes(parity_payloads[p]),
                is_parity=True,
                block_index=block_index,
            )
            parity_chunks.append(parity_chunk)

        return parity_chunks

    def decode_block(
        self,
        received_chunks: dict[int, Chunk],
        k_data_count: int,
        m_parity_count: int,
        start_index: int = 0,
        payload_size: int = DEFAULT_PAYLOAD_SIZE,
        block_index: int = 0,
    ) -> list[Chunk]:
        """
        Reconstruct lost or erased data chunks in an RS block.

        Args:
            received_chunks: Dict mapping chunk_index -> Chunk for all surviving, valid chunks.
            k_data_count: Number of original data chunks in this block.
            m_parity_count: Number of parity chunks generated for this block.
            start_index: First chunk_index of this block (usually 0 for single block).
            payload_size: Byte length of each chunk's payload.
            block_index: Numerical index of this block.

        Returns:
            list[Chunk]: Reconstructed list of all K data chunks for this block.

        Raises:
            reedsolo.ReedSolomonError: If too many chunks were lost (> M erasures).
        """
        total_chunks = k_data_count + m_parity_count
        rsc = reedsolo.RSCodec(m_parity_count)

        # Check which chunks are missing (erasures)
        missing_indices: list[int] = []
        for local_idx in range(total_chunks):
            global_idx = start_index + local_idx
            if global_idx not in received_chunks:
                missing_indices.append(local_idx)

        # If no data chunks are missing, return existing data chunks directly
        if all((start_index + i) in received_chunks for i in range(k_data_count)):
            return [received_chunks[start_index + i] for i in range(k_data_count)]

        # Check erasure capacity limit
        if len(missing_indices) > m_parity_count:
            raise reedsolo.ReedSolomonError(
                f"Unrecoverable strand loss: {len(missing_indices)} chunks missing, "
                f"exceeding parity capacity M={m_parity_count}."
            )

        # Prepare reconstructed column buffers
        reconstructed_data = [bytearray(payload_size) for _ in range(k_data_count)]

        # Column-wise decoding
        for col in range(payload_size):
            full_col = bytearray(total_chunks)
            for local_idx in range(total_chunks):
                global_idx = start_index + local_idx
                if global_idx in received_chunks:
                    full_col[local_idx] = received_chunks[global_idx].payload[col]
                else:
                    full_col[local_idx] = 0  # Erasure placeholder

            # reedsolo decode with explicit erasure positions
            decoded_col, _, _ = rsc.decode(full_col, erase_pos=missing_indices)

            for row in range(k_data_count):
                reconstructed_data[row][col] = decoded_col[row]

        # Build output Chunk objects
        result_chunks: list[Chunk] = []
        for i in range(k_data_count):
            global_idx = start_index + i
            if global_idx in received_chunks:
                result_chunks.append(received_chunks[global_idx])
            else:
                recovered_chunk = Chunk(
                    chunk_id=global_idx,
                    chunk_index=global_idx,
                    payload=bytes(reconstructed_data[i]),
                    is_parity=False,
                    block_index=block_index,
                )
                result_chunks.append(recovered_chunk)

        return result_chunks
