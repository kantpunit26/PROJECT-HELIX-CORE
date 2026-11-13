"""
Unit tests for Phase 2: Fragmentation & Reed-Solomon Error Correction.
"""

from pathlib import Path
import unittest

from helix_core.fragmentation import FragmentationEngine, ReedSolomonCodec, StrandIntegrity
from helix_core.ingestion import BinaryIngestor


class TestPhase2FragmentationAndRS(unittest.TestCase):
    """Test suite for Phase 2 functionality."""

    def setUp(self):
        self.ingestor = BinaryIngestor()
        self.test_png_path = Path("test_sample.png").resolve()
        self.raw_bytes = self.ingestor.read_file(self.test_png_path)
        self.payload_size = 28
        self.engine = FragmentationEngine(payload_size=self.payload_size)
        self.rs_codec = ReedSolomonCodec(parity_ratio=0.15)

    def test_fragmentation_fixed_payload(self):
        """Verify fragments have uniform payload lengths matching synthesis limits."""
        chunks = self.engine.fragment(self.raw_bytes)
        self.assertGreater(len(chunks), 0)

        for c in chunks:
            # Each fragment payload must exactly match payload_size (28 bytes = 112 nt)
            self.assertEqual(len(c.payload), self.payload_size)
            # Chunk ID must be 32-bit integer (0 <= id < 2^32)
            self.assertGreaterEqual(c.chunk_id, 0)
            self.assertLess(c.chunk_id, 2**32)

    def test_fragmentation_and_reassembly(self):
        """Verify reassembly reconstructs exact bytes in order of Index ID."""
        chunks = self.engine.fragment(self.raw_bytes)
        reassembled = self.engine.reassemble(chunks)

        # Truncated to original length (due to zero-padding on last chunk)
        self.assertEqual(reassembled[: len(self.raw_bytes)], self.raw_bytes)

    def test_strand_integrity_crc32(self):
        """Verify CRC32 correctly detects mutations and bit flips."""
        chunk_id = 42
        payload = b"SamplePayloadForCRCTesting12"
        crc = StrandIntegrity.compute_crc32(chunk_id, payload)

        # Match verification
        self.assertTrue(StrandIntegrity.verify_crc32(chunk_id, payload, crc))

        # Mutation in payload (SNP / bit-flip)
        mutated_payload = b"SamplePayloadForCRCTesting13"
        self.assertFalse(StrandIntegrity.verify_crc32(chunk_id, mutated_payload, crc))

        # Mutation in chunk_id
        self.assertFalse(StrandIntegrity.verify_crc32(chunk_id + 1, payload, crc))

    def test_reed_solomon_parity_generation_and_recovery(self):
        """
        Verify Reed-Solomon generates parity chunks capable of recovering 100% of data
        when >=10% of strands are destroyed or lost in a simulated sequencing pool.
        """
        chunks = self.engine.fragment(self.raw_bytes)
        k_data = len(chunks)

        # Generate parity chunks (15% redundancy)
        parity_chunks = self.rs_codec.encode_block(chunks, block_index=0)
        m_parity = len(parity_chunks)
        self.assertGreaterEqual(m_parity, 1)
        self.assertGreaterEqual(m_parity / k_data, 0.10, "Must satisfy >= 10% parity threshold")

        all_chunks = {c.chunk_index: c for c in (chunks + parity_chunks)}

        # Simulate 10% loss by dropping data chunks
        drop_count = max(1, int(k_data * 0.10))
        dropped_indices = list(range(drop_count))  # Drop first drop_count chunks

        received_chunks = {idx: c for idx, c in all_chunks.items() if idx not in dropped_indices}

        # Reconstruct missing data chunks using RS parity
        reconstructed = self.rs_codec.decode_block(
            received_chunks=received_chunks,
            k_data_count=k_data,
            m_parity_count=m_parity,
            start_index=0,
            payload_size=self.payload_size,
        )

        self.assertEqual(len(reconstructed), k_data)
        for i in range(k_data):
            self.assertEqual(
                reconstructed[i].payload,
                chunks[i].payload,
                f"Chunk {i} must be bit-perfectly reconstructed",
            )


if __name__ == "__main__":
    unittest.main()

