"""
Unit tests for Phase 3: Biological Constraint Enforcement (Obfuscation).
"""

from pathlib import Path
import unittest

from helix_core.fragmentation import FragmentationEngine
from helix_core.ingestion import BinaryIngestor
from helix_core.obfuscation import BiologicalObfuscator, ConstraintValidator
from helix_core.types import BioConstraints, Chunk


class TestPhase3Obfuscation(unittest.TestCase):
    """Test suite for Phase 3 functionality."""

    def setUp(self):
        self.validator = ConstraintValidator()
        self.constraints = BioConstraints(max_homopolymer=3, min_gc=0.40, max_gc=0.60)
        self.obfuscator = BiologicalObfuscator(constraints=self.constraints)

    def test_homopolymer_detection(self):
        """Verify detector flags runs of > 3 consecutive identical bases."""
        # AAAA (4 nt) should be flagged
        violations = self.validator.find_homopolymers("ACGTACGTAAAATGC", max_allowed=3)
        self.assertEqual(len(violations), 1)
        self.assertEqual(violations[0], ("A", 8, 4))

        # GGGGG (5 nt) should be flagged
        violations = self.validator.find_homopolymers("ACGTGGGGGACT", max_allowed=3)
        self.assertEqual(len(violations), 1)
        self.assertEqual(violations[0], ("G", 4, 5))

        # AAA (3 nt) is allowed
        violations = self.validator.find_homopolymers("ACGTAAATGCA", max_allowed=3)
        self.assertEqual(len(violations), 0)

    def test_gc_content_calculation(self):
        """Verify exact calculation of GC ratio."""
        self.assertEqual(self.validator.calculate_gc_content("ACGT"), 0.50)
        self.assertEqual(self.validator.calculate_gc_content("GGCC"), 1.00)
        self.assertEqual(self.validator.calculate_gc_content("AATT"), 0.00)
        self.assertEqual(self.validator.calculate_gc_content("AGCTAGCTAG"), 0.50)

    def test_seeded_xor_mask_determinism_and_reversibility(self):
        """Verify XOR mask is deterministic and mathematically self-inverting."""
        seed = 12345
        payload = b"HelloDnaWorldWithBioPayload!"
        mask = self.obfuscator.generate_mask(seed, len(payload))
        self.assertEqual(len(mask), len(payload))

        # Same seed must yield identical mask
        mask2 = self.obfuscator.generate_mask(seed, len(payload))
        self.assertEqual(mask, mask2)

        # Reversibility
        masked = self.obfuscator.apply_xor_mask(payload, seed)
        self.assertNotEqual(masked, payload)
        unmasked = self.obfuscator.apply_xor_mask(masked, seed)
        self.assertEqual(unmasked, payload)

    def test_whitening_repetitive_file_bytes(self):
        """
        Verify seeded XOR mask whitens repeating zero bytes (empty space in files)
        to prevent long runs of AAAAA (since 00=A).
        """
        seed = 999
        empty_space_bytes = b"\x00" * 28  # Unmasked would be 112 consecutive 'A' bases!
        raw_dna = BiologicalObfuscator._bytes_to_dna_quick(empty_space_bytes)
        self.assertEqual(raw_dna, "A" * 112)
        self.assertEqual(self.validator.calculate_gc_content(raw_dna), 0.0)

        # Apply seeded mask
        masked = self.obfuscator.apply_xor_mask(empty_space_bytes, seed)
        masked_dna = BiologicalObfuscator._bytes_to_dna_quick(masked)

        # After masking, homopolymers are broken up and GC is balanced!
        gc = self.validator.calculate_gc_content(masked_dna)
        self.assertGreaterEqual(gc, 0.35)
        self.assertLessEqual(gc, 0.65)
        self.assertNotEqual(masked_dna, "A" * 112)

    def test_adaptive_biological_screening(self):
        """
        Verify adaptive screening finds a candidate salt within 4-byte Chunk ID space
        that guarantees zero homopolymers > 3 and 40%-60% GC content.
        """
        # Test on 10 arbitrary chunks
        for i in range(10):
            chunk = Chunk(
                chunk_id=i,
                chunk_index=i,
                payload=bytes([(i * 13 + j * 7) % 256 for j in range(28)]),
            )
            masked_chunk = self.obfuscator.obfuscate_chunk_adaptive(chunk)

            # Check unmasking integrity
            recovered_chunk = self.obfuscator.unmask_chunk(masked_chunk)
            self.assertEqual(recovered_chunk.payload, chunk.payload)

            # Check biological constraints on payload DNA
            dna_seq = BiologicalObfuscator._bytes_to_dna_quick(masked_chunk.masked_payload)
            violations = self.validator.find_homopolymers(dna_seq, max_allowed=3)
            gc = self.validator.calculate_gc_content(dna_seq)

            self.assertEqual(len(violations), 0, f"Found homopolymer violations in chunk {i}")
            self.assertGreaterEqual(gc, 0.40, f"GC below 40% in chunk {i}: {gc:.2%}")
            self.assertLessEqual(gc, 0.60, f"GC above 60% in chunk {i}: {gc:.2%}")


if __name__ == "__main__":
    unittest.main()

