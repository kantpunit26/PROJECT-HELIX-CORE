"""
HELIX-CORE Phase 4 & 5: Digital-to-Biological & Biological-to-Digital Translation.

This module is responsible for:
1. Precise mapping between binary dibits and nucleotide bases (00=A, 01=C, 10=G, 11=T).
2. Assembling full physical DNA oligo strands with Universal Forward and Reverse Primers.
3. Disassembling raw sequenced DNA strands: stripping primers, parsing 4-byte Chunk IDs,
   masked payloads, and CRC32 checksums.
4. Exporting synthesis-ready .fasta and .csv files formatted for Twist Bioscience and IDT.
5. Ingesting simulated or physical sequencing FASTA files.
"""

from pathlib import Path
import struct
from typing import Optional, Union

from .config import (
    BASE_TO_BITS,
    BITS_TO_BASE,
    BYTE_TO_BASES,
    CHUNK_ID_BYTES,
    CHUNK_ID_WHITENING_MASK,
    CRC32_BYTES,
    DEFAULT_FORWARD_PRIMER,
    DEFAULT_REVERSE_PRIMER,
    DEFAULT_PAYLOAD_SIZE,
    decode_wire_id,
    encode_wire_id,
    pack_chunk_id,
    unpack_chunk_id,
)
from .obfuscation import ConstraintValidator
from .types import DnaStrand, MaskedChunk


class DnaTranslator:
    """Production translation engine bridging digital bits and biological DNA strands."""

    def __init__(
        self,
        forward_primer: str = DEFAULT_FORWARD_PRIMER,
        reverse_primer: str = DEFAULT_REVERSE_PRIMER,
        payload_size: int = DEFAULT_PAYLOAD_SIZE,
    ):
        """
        Initialize the DNA translation engine.

        Args:
            forward_primer: 5' universal PCR amplification primer (default 20 nt).
            reverse_primer: 3' universal PCR amplification primer (default 20 nt).
            payload_size: Payload length in bytes (default 28 bytes = 112 nt).
        """
        self.forward_primer = forward_primer.upper()
        self.reverse_primer = reverse_primer.upper()
        self.payload_size = payload_size
        self.validator = ConstraintValidator()

        # Build reverse base lookup for decoding 4 bases -> 1 byte
        self._base4_to_byte = {bases: b for b, bases in enumerate(BYTE_TO_BASES)}

    @staticmethod
    def bytes_to_dna(data: bytes) -> str:
        """
        Translate raw bytes into nucleotide sequence using strict dibit mapping:
        00 = A, 01 = C, 10 = G, 11 = T.

        Args:
            data: Raw byte array.

        Returns:
            str: Nucleotide string of length len(data) * 4.
        """
        return "".join(BYTE_TO_BASES[b] for b in data)

    def dna_to_bytes(self, dna_sequence: str) -> bytes:
        """
        Translate a nucleotide sequence back into raw bytes.

        Args:
            dna_sequence: Nucleotide string (A, C, G, T) of length multiple of 4.

        Returns:
            bytes: Reconstructed byte array.

        Raises:
            ValueError: If sequence length is not a multiple of 4 or contains invalid characters.
        """
        clean_seq = dna_sequence.strip().upper()
        if len(clean_seq) % 4 != 0:
            raise ValueError(f"DNA sequence length ({len(clean_seq)}) must be a multiple of 4.")

        byte_list: list[int] = []
        for i in range(0, len(clean_seq), 4):
            kmer = clean_seq[i : i + 4]
            val = self._base4_to_byte.get(kmer)
            if val is None:
                # Handle any unexpected non-ACGT characters
                bits = "".join(BASE_TO_BITS[b] for b in kmer)
                val = int(bits, 2)
            byte_list.append(val)

        return bytes(byte_list)

    def assemble_strand(self, masked_chunk: MaskedChunk) -> DnaStrand:
        """
        Assemble a complete physical DNA strand from a MaskedChunk.

        Strand Architecture (5' to 3'):
          [Forward Primer (20 nt)]
          + [4-byte Chunk ID (16 nt)]
          + [Masked Payload (112 nt)]
          + [4-byte CRC32 (16 nt)]
          + [Reverse Primer (20 nt)]
          = 184 nt total (strictly within 150-200 bp synthesis limit).

        Args:
            masked_chunk: Obfuscated chunk with CRC32.

        Returns:
            DnaStrand: Complete biological strand.
        """
        # Self-describing wire representation for physical synthesis
        meta = unpack_chunk_id(masked_chunk.chunk_id)
        k_val = meta["total_data_chunks"] or 1
        is_p = masked_chunk.is_parity or meta["is_parity"]
        salt_val = masked_chunk.salt or meta["salt"]

        wire_bytes = encode_wire_id(
            chunk_index=masked_chunk.chunk_index,
            total_data_chunks=k_val,
            is_parity=is_p,
            salt=salt_val,
        )
        crc_bytes = struct.pack(">I", masked_chunk.crc32 & 0xFFFFFFFF)

        id_dna = self.bytes_to_dna(wire_bytes)
        payload_dna = self.bytes_to_dna(masked_chunk.masked_payload)
        crc_dna = self.bytes_to_dna(crc_bytes)

        full_seq = self.forward_primer + id_dna + payload_dna + crc_dna + self.reverse_primer

        gc = self.validator.calculate_gc_content(full_seq)

        # Max consecutive repeat
        max_run = 1
        run = 1
        for i in range(1, len(full_seq)):
            if full_seq[i] == full_seq[i - 1]:
                run += 1
                if run > max_run:
                    max_run = run
            else:
                run = 1

        return DnaStrand(
            chunk_id=masked_chunk.chunk_id,
            chunk_index=masked_chunk.chunk_index,
            forward_primer=self.forward_primer,
            id_sequence=id_dna,
            payload_sequence=payload_dna,
            crc_sequence=crc_dna,
            reverse_primer=self.reverse_primer,
            full_sequence=full_seq,
            gc_content=gc,
            max_homopolymer=max_run,
            is_parity=masked_chunk.is_parity,
        )

    def disassemble_strand(self, raw_sequence: str) -> MaskedChunk:
        """
        Disassemble a raw sequenced DNA strand into a MaskedChunk.

        Strips universal Forward and Reverse primers, extracts the 4-byte Chunk ID,
        the masked payload bytes, and the 4-byte CRC32 checksum.

        Args:
            raw_sequence: Nucleotide sequence string from sequencing output.

        Returns:
            MaskedChunk: Extracted masked chunk.

        Raises:
            ValueError: If primers cannot be detected or strand length is corrupt.
        """
        seq = raw_sequence.strip().upper()
        fwd_len = len(self.forward_primer)
        rev_len = len(self.reverse_primer)

        # Verify and strip primers
        if not seq.startswith(self.forward_primer):
            raise ValueError("Strand rejected: 5' Forward Primer mismatch or missing.")
        if not seq.endswith(self.reverse_primer):
            raise ValueError("Strand rejected: 3' Reverse Primer mismatch or missing.")

        core_seq = seq[fwd_len : len(seq) - rev_len]

        expected_core_len = (CHUNK_ID_BYTES + self.payload_size + CRC32_BYTES) * 4
        if len(core_seq) != expected_core_len:
            raise ValueError(
                f"Strand core length ({len(core_seq)} nt) != expected ({expected_core_len} nt)."
            )

        # Slice sub-sequences
        id_dna_len = CHUNK_ID_BYTES * 4      # 16 nt
        crc_dna_len = CRC32_BYTES * 4        # 16 nt
        payload_dna_len = self.payload_size * 4

        id_dna = core_seq[:id_dna_len]
        payload_dna = core_seq[id_dna_len : id_dna_len + payload_dna_len]
        crc_dna = core_seq[id_dna_len + payload_dna_len :]

        # Decode to bytes
        id_bytes = self.dna_to_bytes(id_dna)
        masked_payload = self.dna_to_bytes(payload_dna)
        crc_bytes = self.dna_to_bytes(crc_dna)

        # Unpack self-describing wire ID
        wire_meta = decode_wire_id(id_bytes)
        chunk_id = wire_meta["chunk_id"]
        crc32_val = struct.unpack(">I", crc_bytes)[0]

        return MaskedChunk(
            chunk_id=chunk_id,
            chunk_index=wire_meta["chunk_index"],
            masked_payload=masked_payload,
            crc32=crc32_val,
            is_parity=wire_meta["is_parity"],
            salt=wire_meta["salt"],
        )

    @staticmethod
    def export_fasta(strands: list[DnaStrand], filepath: Union[str, Path]) -> Path:
        """
        Export synthesized sequences into a standard FASTA file for synthesis providers.

        Format compatible with Twist Bioscience, IDT, and CustomArray.
        """
        out_path = Path(filepath).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)

        with open(out_path, "w", encoding="utf-8") as fh:
            for s in strands:
                kind = "parity" if s.is_parity else "data"
                header = (
                    f">HELIX_{s.chunk_id:08X}|index={s.chunk_index}|type={kind}|"
                    f"len={s.length}|gc={s.gc_content:.1%}|max_hp={s.max_homopolymer}"
                )
                fh.write(f"{header}\n{s.full_sequence}\n")

        return out_path

    @staticmethod
    def export_csv(strands: list[DnaStrand], filepath: Union[str, Path]) -> Path:
        """Export synthesis order CSV with standard Twist/IDT columns."""
        out_path = Path(filepath).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)

        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write("Oligo_Name,Sequence,Length_bp,GC_Percent,Type,Chunk_ID,Chunk_Index\n")
            for s in strands:
                kind = "parity" if s.is_parity else "data"
                name = f"HELIX_{s.chunk_id:08X}"
                fh.write(
                    f"{name},{s.full_sequence},{s.length},{s.gc_content:.2%},{kind},{s.chunk_id},{s.chunk_index}\n"
                )

        return out_path

    @staticmethod
    def import_fasta(filepath: Union[str, Path]) -> list[str]:
        """Read sequences from a FASTA file."""
        path = Path(filepath).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"FASTA file not found: {path}")

        sequences: list[str] = []
        current_seq: list[str] = []

        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_seq:
                        sequences.append("".join(current_seq))
                        current_seq = []
                else:
                    current_seq.append(line)

            if current_seq:
                sequences.append("".join(current_seq))

        return sequences
