"""
HELIX-CORE Phase 1: Raw File Ingestion & Binary Conversion.

This module is responsible for:
1. Ingesting arbitrary digital files (.png, .jpg, .pdf, .mp4, etc.) via standard binary I/O.
2. Converting the raw byte stream into an unbroken binary string.
3. Cryptographically fingerprinting files via SHA-256 hashes to guarantee bit-perfect fidelity.
4. Packaging payloads into self-describing HELIX containers to preserve exact file lengths
   across fixed-size biological chunk boundaries.
"""

import hashlib
import os
from pathlib import Path
import struct
from typing import Optional, Union

from .config import HELIX_CONTAINER_MAGIC


class BinaryIngestor:
    """Production-grade engine for file ingestion, binary conversion, and integrity hashing."""

    @staticmethod
    def read_file(filepath: Union[str, Path]) -> bytes:
        """
        Ingest a file from the physical disk as a raw byte array.

        Strict constraint: Forbidden from using mock data or random generators.
        Must use standard Python I/O open(file, 'rb').

        Args:
            filepath: Path to the target file on disk.

        Returns:
            bytes: Raw binary content of the file.

        Raises:
            FileNotFoundError: If the file does not exist on disk.
            PermissionError: If read permissions are lacking.
            ValueError: If the file is completely empty (0 bytes).
        """
        path = Path(filepath).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Source file not found on disk: {path}")

        file_size = path.stat().st_size
        if file_size == 0:
            raise ValueError(f"Target file is empty (0 bytes): {path}")

        with open(path, "rb") as fh:
            data = fh.read()

        return data

    @staticmethod
    def write_file(filepath: Union[str, Path], data: bytes) -> Path:
        """
        Write a raw byte array back to physical disk.

        Args:
            filepath: Destination file path.
            data: Raw byte array to write.

        Returns:
            Path: Resolved path of the written file.
        """
        dest_path = Path(filepath).resolve()
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        with open(dest_path, "wb") as fh:
            fh.write(data)
        return dest_path

    @staticmethod
    def compute_sha256(data: Union[bytes, str, Path]) -> str:
        """
        Calculate the cryptographic SHA-256 checksum.

        Supports in-memory byte arrays or direct file path streaming to handle
        large video/audio files without excessive RAM consumption.

        Args:
            data: Raw bytes OR a file path.

        Returns:
            str: 64-character lowercase hexadecimal SHA-256 hash.
        """
        hasher = hashlib.sha256()

        if isinstance(data, (str, Path)):
            path = Path(data).resolve()
            if not path.is_file():
                raise FileNotFoundError(f"File not found for hash calculation: {path}")
            with open(path, "rb") as fh:
                while chunk := fh.read(65536):
                    hasher.update(chunk)
        elif isinstance(data, (bytes, bytearray)):
            hasher.update(data)
        else:
            raise TypeError(f"Unsupported data type for SHA-256 computation: {type(data)}")

        return hasher.hexdigest()

    @staticmethod
    def bytes_to_binary_string(data: bytes) -> str:
        """
        Convert an entire byte stream into an unbroken binary string.

        Each byte is formatted as an 8-bit binary representation ('0' and '1'),
        preserving leading zeros and order.

        Args:
            data: Raw byte array.

        Returns:
            str: Continuous unbroken bit string (e.g., '0110100001100101...').
        """
        # Fast formatting via list comprehension
        return "".join(f"{b:08b}" for b in data)

    @staticmethod
    def binary_string_to_bytes(bit_string: str) -> bytes:
        """
        Convert an unbroken binary string back into a raw byte array.

        Args:
            bit_string: Binary string containing only '0' and '1' characters.

        Returns:
            bytes: Reconstructed raw byte array.

        Raises:
            ValueError: If bit_string length is not a multiple of 8 or contains invalid characters.
        """
        if len(bit_string) % 8 != 0:
            raise ValueError(
                f"Binary string length ({len(bit_string)}) must be an exact multiple of 8."
            )

        # Slice in 8-bit steps and convert to integer
        byte_list = [
            int(bit_string[i : i + 8], 2)
            for i in range(0, len(bit_string), 8)
        ]
        return bytes(byte_list)

    @classmethod
    def pack_container(cls, raw_data: bytes, filename: str = "") -> bytes:
        """
        Package raw file bytes into a self-describing HELIX container stream.

        Structure:
          - Magic Header (4 bytes): b'HLX1'
          - Original File Size (8 bytes, uint64 big-endian)
          - Original SHA-256 Digest (32 bytes raw binary)
          - Filename Length (2 bytes, uint16 big-endian)
          - Filename (UTF-8 encoded bytes)
          - Payload (Raw file bytes)

        This ensures that when data is fragmented into fixed-size DNA chunks (which
        may require zero-padding on the final chunk), the exact byte boundary and
        hash are preserved with 100% mathematical fidelity.
        """
        magic = HELIX_CONTAINER_MAGIC
        file_size = len(raw_data)
        sha256_raw = hashlib.sha256(raw_data).digest()
        fname_bytes = filename.encode("utf-8")
        fname_len = len(fname_bytes)

        header = struct.pack(
            ">4sQ32sH",
            magic,
            file_size,
            sha256_raw,
            fname_len,
        ) + fname_bytes

        return header + raw_data

    @classmethod
    def unpack_container(cls, container_data: bytes) -> tuple[bytes, dict]:
        """
        Unpack a HELIX container stream and extract the exact original file bytes.

        Validates the magic header and checks the cryptographic SHA-256 digest
        against the unpacked data.

        Returns:
            tuple[bytes, dict]: (original_file_bytes, metadata_dict)

        Raises:
            ValueError: If container magic is invalid, corrupted, or SHA-256 does not match.
        """
        min_header_size = 4 + 8 + 32 + 2
        if len(container_data) < min_header_size:
            raise ValueError("Corrupt container: byte stream shorter than minimum header size.")

        magic, file_size, sha256_expected, fname_len = struct.unpack(
            ">4sQ32sH", container_data[:min_header_size]
        )

        if magic != HELIX_CONTAINER_MAGIC:
            raise ValueError(f"Invalid container signature: {magic!r} != {HELIX_CONTAINER_MAGIC!r}")

        offset = min_header_size
        filename = container_data[offset : offset + fname_len].decode("utf-8", errors="replace")
        offset += fname_len

        payload = container_data[offset : offset + file_size]
        if len(payload) != file_size:
            raise ValueError(
                f"Container payload truncation: expected {file_size} bytes, found {len(payload)} bytes."
            )

        sha256_actual = hashlib.sha256(payload).digest()
        if sha256_actual != sha256_expected:
            raise ValueError("Container cryptographic mismatch: SHA-256 verification failed.")

        metadata = {
            "filename": filename,
            "original_size": file_size,
            "sha256": sha256_actual.hex(),
        }
        return payload, metadata

