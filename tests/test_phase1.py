"""
Unit tests for Phase 1: Raw File Ingestion & Binary Conversion.
"""

import hashlib
from pathlib import Path
import unittest

from helix_core.ingestion import BinaryIngestor


class TestPhase1Ingestion(unittest.TestCase):
    """Test suite for Phase 1 functionality."""

    def setUp(self):
        self.ingestor = BinaryIngestor()
        self.test_png_path = Path("test_sample.png").resolve()

    def test_read_real_file(self):
        """Verify standard Python I/O reads actual disk file without mock data."""
        self.assertTrue(self.test_png_path.is_file(), "Test sample PNG must exist on disk")
        raw_bytes = self.ingestor.read_file(self.test_png_path)
        self.assertGreater(len(raw_bytes), 0)
        # Check standard PNG header signature: \x89PNG\r\n\x1a\n
        self.assertEqual(raw_bytes[:8], b"\x89PNG\r\n\x1a\n")

    def test_sha256_computation(self):
        """Verify cryptographic SHA-256 calculation matches standard hashlib."""
        raw_bytes = self.ingestor.read_file(self.test_png_path)
        expected_hash = hashlib.sha256(raw_bytes).hexdigest()
        actual_hash = self.ingestor.compute_sha256(self.test_png_path)
        self.assertEqual(actual_hash, expected_hash)

        # Test hashing in-memory bytes
        byte_hash = self.ingestor.compute_sha256(raw_bytes)
        self.assertEqual(byte_hash, expected_hash)

    def test_bytes_to_binary_string_and_reverse(self):
        """Verify unbroken binary string conversion is 100% bit-perfect and reversible."""
        raw_bytes = self.ingestor.read_file(self.test_png_path)
        bit_string = self.ingestor.bytes_to_binary_string(raw_bytes)

        # Length must be exactly len(raw_bytes) * 8
        self.assertEqual(len(bit_string), len(raw_bytes) * 8)
        self.assertTrue(set(bit_string).issubset({"0", "1"}))

        # Reverse back to bytes
        reconstructed_bytes = self.ingestor.binary_string_to_bytes(bit_string)
        self.assertEqual(reconstructed_bytes, raw_bytes)
        self.assertEqual(
            hashlib.sha256(reconstructed_bytes).hexdigest(),
            hashlib.sha256(raw_bytes).hexdigest(),
        )

    def test_container_packaging_and_unpacking(self):
        """Verify self-describing container preserves exact length and validates SHA-256."""
        raw_bytes = self.ingestor.read_file(self.test_png_path)
        filename = "test_sample.png"

        container_bytes = self.ingestor.pack_container(raw_bytes, filename=filename)
        self.assertTrue(container_bytes.startswith(b"HLX1"))

        unpacked_bytes, metadata = self.ingestor.unpack_container(container_bytes)
        self.assertEqual(unpacked_bytes, raw_bytes)
        self.assertEqual(metadata["filename"], filename)
        self.assertEqual(metadata["original_size"], len(raw_bytes))
        self.assertEqual(metadata["sha256"], self.ingestor.compute_sha256(raw_bytes))


if __name__ == "__main__":
    unittest.main()

