"""
HELIX-CORE Production Pipeline: Bi-directional DNA Data Storage Architecture.

Orchestrates Phase 1 (Ingestion & Binary Conversion), Phase 2 (Fragmentation & RS Error Correction),
Phase 3 (Biological Obfuscation & Homopolymer/GC Balancing), Phase 4 (Translation & Assembly),
Phase 5 (Disassembly, CRC Verification & RS Recovery), and Phase 6 (Reassembly & Cryptographic Validation).
"""

import math
from pathlib import Path
import random
from typing import Optional, Union

from .config import (
    DEFAULT_FORWARD_PRIMER,
    DEFAULT_PARITY_RATIO,
    DEFAULT_PAYLOAD_SIZE,
    DEFAULT_REVERSE_PRIMER,
    MIN_PARITY_CHUNKS,
    unpack_chunk_id,
)
from .fragmentation import FragmentationEngine, ReedSolomonCodec, StrandIntegrity
from .ingestion import BinaryIngestor
from .obfuscation import BiologicalObfuscator, ConstraintValidator
from .translation import DnaTranslator
from .types import BioConstraints, Chunk, DnaStrand, MaskedChunk, PipelineStats


class HelixPipeline:
    """Master orchestrator for the HELIX-CORE DNA data storage pipeline."""

    def __init__(
        self,
        payload_size: int = DEFAULT_PAYLOAD_SIZE,
        parity_ratio: float = DEFAULT_PARITY_RATIO,
        constraints: Optional[BioConstraints] = None,
        forward_primer: str = DEFAULT_FORWARD_PRIMER,
        reverse_primer: str = DEFAULT_REVERSE_PRIMER,
        adaptive_screening: bool = True,
    ):
        """
        Initialize the HELIX-CORE pipeline.

        Args:
            payload_size: Fragment payload size in bytes (default 28 -> 184 nt strand).
            parity_ratio: Parity chunk ratio (default 0.15 = 15% redundancy).
            constraints: Target biophysical constraints (homopolymer <= 3, GC 40-60%).
            forward_primer: 5' universal PCR primer.
            reverse_primer: 3' universal PCR primer.
            adaptive_screening: If True, uses salt-based biological screening to guarantee
                                zero homopolymers > 3 and 40-60% GC on 100% of strands.
        """
        self.payload_size = payload_size
        self.parity_ratio = parity_ratio
        self.constraints = constraints or BioConstraints()
        self.adaptive_screening = adaptive_screening

        # Instantiate pipeline engines
        self.ingestor = BinaryIngestor()
        self.fragmenter = FragmentationEngine(payload_size=payload_size)
        self.rs_codec = ReedSolomonCodec(parity_ratio=parity_ratio)
        self.obfuscator = BiologicalObfuscator(constraints=self.constraints)
        self.translator = DnaTranslator(
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            payload_size=payload_size,
        )
        self.validator = ConstraintValidator()

    def encode_file(self, input_filepath: Union[str, Path]) -> tuple[list[DnaStrand], dict]:
        """
        Execute Phase 1, 2, 3, and 4: Convert any raw file into synthesis-ready DNA strands.

        Workflow:
          1. Phase 1: Ingest raw bytes from disk and calculate SHA-256.
          2. Phase 1: Pack into container stream preserving length and hash.
          3. Phase 1: Convert entire stream into unbroken binary string.
          4. Phase 2: Fragment binary into fixed payloads (150-200 bp).
          5. Phase 2: Generate Reed-Solomon parity chunks (15% redundancy).
          6. Phase 3: Apply seeded XOR bitmask with Chunk ID, avoid homopolymers (<=3), balance GC (40-60%).
          7. Phase 4: Map binary dibits to nucleotides (00=A, 01=C, 10=G, 11=T) and attach universal primers.

        Args:
            input_filepath: Source file path on disk (NO mock data).

        Returns:
            tuple[list[DnaStrand], dict]: Synthesized strands and encoding metadata.
        """
        input_path = Path(input_filepath).resolve()
        raw_bytes = self.ingestor.read_file(input_path)
        original_sha256 = self.ingestor.compute_sha256(raw_bytes)

        # Containerize to preserve file boundary across fixed chunks
        container_bytes = self.ingestor.pack_container(raw_bytes, filename=input_path.name)

        # Convert to unbroken binary bit string
        bit_string = self.ingestor.bytes_to_binary_string(container_bytes)

        # Phase 2: Fragmentation
        data_chunks = self.fragmenter.fragment(container_bytes)
        k_data_chunks = len(data_chunks)

        # Phase 2: Reed-Solomon Error Correction (Generate Parity Chunks)
        parity_chunks = self.rs_codec.encode_block(data_chunks, block_index=0)
        m_parity_chunks = len(parity_chunks)

        all_raw_chunks = data_chunks + parity_chunks

        # Phase 3 & 4: Obfuscation and Physical Strand Assembly
        strands: list[DnaStrand] = []
        for chunk in all_raw_chunks:
            if self.adaptive_screening:
                masked_chunk = self.obfuscator.obfuscate_chunk_adaptive(
                    chunk,
                    max_salt_attempts=256,
                    prefix_dna=self.translator.forward_primer,
                    suffix_dna=self.translator.reverse_primer,
                )
            else:
                masked_chunk = self.obfuscator.mask_chunk(chunk)

            strand = self.translator.assemble_strand(masked_chunk)
            strands.append(strand)

        metadata = {
            "filename": input_path.name,
            "original_size": len(raw_bytes),
            "original_sha256": original_sha256,
            "total_bits": len(bit_string),
            "data_chunks": k_data_chunks,
            "parity_chunks": m_parity_chunks,
            "total_strands": len(strands),
            "strand_length_bp": strands[0].length if strands else 0,
        }

        return strands, metadata

    def decode_strands_to_bytes(
        self,
        received_strands: list[DnaStrand | str],
        total_data_chunks: Optional[int] = None,
        total_parity_chunks: Optional[int] = None,
    ) -> tuple[bytes, dict]:
        """
        Execute Phase 5 and 6 in-memory: Recover exact raw file bytes and metadata from sequenced DNA strands.

        Uses self-describing 32-bit strand metadata to automatically determine K (data count) and M (parity count)
        if not explicitly provided.
        """
        valid_chunks: dict[int, Chunk] = {}
        discarded_mutated_count = 0

        # Phase 5: Disassembly, Unmasking, and CRC32 Validation
        for item in received_strands:
            raw_seq = item.full_sequence if isinstance(item, DnaStrand) else item
            try:
                masked_chunk = self.translator.disassemble_strand(raw_seq)
            except Exception:
                discarded_mutated_count += 1
                continue

            try:
                # Reverse XOR mask with Chunk ID seed & verify CRC32
                unmasked_chunk = self.obfuscator.unmask_chunk(masked_chunk)
                valid_chunks[unmasked_chunk.chunk_index] = unmasked_chunk
            except ValueError:
                # CRC32 failed: Strand suffered mutations (SNPs, indels); safely discard
                discarded_mutated_count += 1
                continue

        if not valid_chunks:
            raise ValueError("No valid DNA strands could be decoded (primers or CRC32 failed).")

        # Automatically determine K (data chunks) from self-describing strand IDs
        if total_data_chunks is None or total_data_chunks <= 0:
            k_candidates = [
                unpack_chunk_id(c.chunk_id)["total_data_chunks"]
                for c in valid_chunks.values()
                if unpack_chunk_id(c.chunk_id)["total_data_chunks"] > 0
            ]
            if k_candidates:
                total_data_chunks = max(k_candidates)
            else:
                # Fallback: non-parity chunks count
                non_parity = [c.chunk_index for c in valid_chunks.values() if not c.is_parity]
                total_data_chunks = max(non_parity) + 1 if non_parity else len(valid_chunks)

        # Automatically determine M (parity chunks)
        if total_parity_chunks is None or total_parity_chunks <= 0:
            parity_count = sum(1 for c in valid_chunks.values() if c.is_parity)
            total_parity_chunks = max(
                parity_count,
                MIN_PARITY_CHUNKS,
                math.ceil(total_data_chunks * self.parity_ratio),
            )

        # Check if all data chunks are already present without needing RS
        all_data_present = all(i in valid_chunks for i in range(total_data_chunks))
        if all_data_present:
            reconstructed_data_chunks = [valid_chunks[i] for i in range(total_data_chunks)]
        else:
            # Phase 5: Reed-Solomon Erasure Reconstruction
            reconstructed_data_chunks = self.rs_codec.decode_block(
                received_chunks=valid_chunks,
                k_data_count=total_data_chunks,
                m_parity_count=total_parity_chunks,
                start_index=0,
                payload_size=self.payload_size,
                block_index=0,
            )

        # Phase 6: Reassemble Continuous Stream
        container_bytes = self.fragmenter.reassemble(reconstructed_data_chunks)

        # Unpack Container to recover original raw bytes and metadata
        recovered_bytes, metadata = self.ingestor.unpack_container(container_bytes)

        metadata["total_data_chunks"] = total_data_chunks
        metadata["total_parity_chunks"] = total_parity_chunks
        metadata["discarded_mutations"] = discarded_mutated_count
        metadata["valid_strands_received"] = len(valid_chunks)

        return recovered_bytes, metadata

    def decode_strands(
        self,
        received_strands: list[DnaStrand | str],
        output_filepath: Union[str, Path],
        total_data_chunks: Optional[int] = None,
        total_parity_chunks: Optional[int] = None,
    ) -> Path:
        """
        Execute Phase 5 and 6: Recover exact raw file from sequenced DNA strands and write to disk.

        Workflow:
          1. Phase 5: Strip primers, translate nucleotides back to XOR-masked binary.
          2. Phase 5: Reverse XOR mask using Chunk ID seed.
          3. Phase 5: Validate CRC32 checksum. Discard mutated strands.
          4. Phase 5: Reconstruct missing/erased data chunks using Reed-Solomon parity chunks.
          5. Phase 6: Order data chunks by Index ID and reassemble continuous byte stream.
          6. Phase 6: Unpack container and write recovered file to disk.
          7. Phase 6: Run cryptographic SHA-256 hash comparison.

        Args:
            received_strands: List of DnaStrand objects or raw sequence strings.
            output_filepath: Destination path for the recovered file.
            total_data_chunks: Optional expected number of data chunks (K). Inferred from strands if omitted.
            total_parity_chunks: Optional number of parity chunks (M). Inferred if omitted.

        Returns:
            Path: Path of the recovered file on disk.

        Raises:
            ValueError: If file reconstruction fails or hash does not match.
        """
        recovered_bytes, metadata = self.decode_strands_to_bytes(
            received_strands=received_strands,
            total_data_chunks=total_data_chunks,
            total_parity_chunks=total_parity_chunks,
        )

        # Write to physical disk
        dest_path = self.ingestor.write_file(output_filepath, recovered_bytes)

        return dest_path

    @staticmethod
    def simulate_sequencing_loss(
        strands: list[DnaStrand],
        drop_rate: float = 0.10,
        mutation_rate: float = 0.00,
        random_seed: int = 42,
    ) -> list[DnaStrand]:
        """
        Simulate wet-lab Next-Gen Sequencing degradation:
        1. Pool Dropout: Drops `drop_rate` fraction of strands from the liquid pool.
        2. Point Mutations (SNPs): Mutates nucleotides on random strands to test CRC32 rejection.

        Args:
            strands: List of synthesized DNA strands.
            drop_rate: Fraction of strands destroyed/lost (e.g. 0.10 = 10%).
            mutation_rate: Fraction of surviving strands that suffer a point mutation.
            random_seed: Seed for reproducible wet-lab simulation.

        Returns:
            list[DnaStrand]: Surviving sequencing output pool.
        """
        rng = random.Random(random_seed)
        n = len(strands)
        num_drops = int(round(n * drop_rate))
        drop_indices = set(rng.sample(range(n), min(num_drops, n)))

        # Identify candidate indices for point mutations from surviving pool
        surviving_pool_indices = [i for i in range(n) if i not in drop_indices]
        num_mutations = int(round(len(surviving_pool_indices) * mutation_rate))
        mutation_indices = set(
            rng.sample(surviving_pool_indices, min(num_mutations, len(surviving_pool_indices)))
        )

        surviving: list[DnaStrand] = []
        for idx in surviving_pool_indices:
            strand = strands[idx]
            seq_chars = list(strand.full_sequence)

            # Check if this strand suffers a point mutation (SNP)
            if idx in mutation_indices:
                mut_pos = rng.randint(20, len(seq_chars) - 21)  # mutate payload region
                orig_base = seq_chars[mut_pos]
                bases = [b for b in "ACGT" if b != orig_base]
                seq_chars[mut_pos] = rng.choice(bases)

            mutated_seq = "".join(seq_chars)
            surviving.append(
                DnaStrand(
                    chunk_id=strand.chunk_id,
                    chunk_index=strand.chunk_index,
                    forward_primer=strand.forward_primer,
                    id_sequence=strand.id_sequence,
                    payload_sequence=strand.payload_sequence,
                    crc_sequence=strand.crc_sequence,
                    reverse_primer=strand.reverse_primer,
                    full_sequence=mutated_seq,
                    gc_content=strand.gc_content,
                    max_homopolymer=strand.max_homopolymer,
                    is_parity=strand.is_parity,
                )
            )

        return surviving

    def run_full_cycle(
        self,
        input_filepath: Union[str, Path],
        output_filepath: Union[str, Path],
        fasta_export_path: Optional[Union[str, Path]] = None,
        drop_rate: float = 0.10,
        mutation_rate: float = 0.02,
    ) -> PipelineStats:
        """
        Execute the complete 6-Phase bi-directional DNA data storage pipeline on a real file.

        Verifies bit-perfect reconstruction via cryptographic SHA-256 hash comparison.
        """
        input_path = Path(input_filepath).resolve()
        output_path = Path(output_filepath).resolve()

        # Phase 1, 2, 3, 4: Encoding
        strands, meta = self.encode_file(input_path)

        # Optional FASTA Export
        if fasta_export_path:
            self.translator.export_fasta(strands, fasta_export_path)

        # Wet-Lab Simulation: 10% Strand Loss + 2% Point Mutations
        sequenced_strands = self.simulate_sequencing_loss(
            strands, drop_rate=drop_rate, mutation_rate=mutation_rate
        )

        # Phase 5 & 6: Decoding & Bit-Perfect Reassembly
        self.decode_strands(
            received_strands=sequenced_strands,
            output_filepath=output_path,
            total_data_chunks=meta["data_chunks"],
            total_parity_chunks=meta["parity_chunks"],
        )

        recovered_sha256 = self.ingestor.compute_sha256(output_path)
        is_bit_perfect = recovered_sha256 == meta["original_sha256"]

        # Biophysical metrics across synthesized library
        gc_list = [s.gc_content for s in strands]
        max_hp_all = max(s.max_homopolymer for s in strands) if strands else 0

        return PipelineStats(
            original_filename=meta["filename"],
            original_size_bytes=meta["original_size"],
            original_sha256=meta["original_sha256"],
            total_binary_bits=meta["total_bits"],
            data_chunk_count=meta["data_chunks"],
            parity_chunk_count=meta["parity_chunks"],
            total_strand_count=meta["total_strands"],
            payload_size_bytes=self.payload_size,
            strand_length_bp=meta["strand_length_bp"],
            mean_gc_content=sum(gc_list) / len(gc_list) if gc_list else 0.0,
            min_gc_content=min(gc_list) if gc_list else 0.0,
            max_gc_content=max(gc_list) if gc_list else 0.0,
            max_homopolymer_observed=max_hp_all,
            recovered_sha256=recovered_sha256,
            is_bit_perfect=is_bit_perfect,
        )
