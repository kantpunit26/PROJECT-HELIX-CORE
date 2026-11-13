"""
HELIX-CORE On-The-Spot DNA Sequence Verification Tool.

Allows users to instantly verify:
1. An individual DNA strand (passed via command line or selected from FASTA).
2. The exact biophysical properties (length, GC content, homopolymers, primers).
3. The underlying binary bits, hexadecimal bytes, and decoded data within that specific strand.
4. A complete .fasta oligo pool, proving that the physical DNA sequences decode
   back into the exact bit-perfect original file.
"""

import argparse
import hashlib
from pathlib import Path
import struct
import sys

from helix_core.config import (
    BASE_TO_BITS,
    BITS_TO_BASE,
    DEFAULT_FORWARD_PRIMER,
    DEFAULT_PAYLOAD_SIZE,
    DEFAULT_REVERSE_PRIMER,
)
from helix_core.fragmentation import FragmentationEngine, ReedSolomonCodec, StrandIntegrity
from helix_core.ingestion import BinaryIngestor
from helix_core.obfuscation import BiologicalObfuscator, ConstraintValidator
from helix_core.pipeline import HelixPipeline
from helix_core.translation import DnaTranslator


def format_hex_dump(data: bytes, bytes_per_line: int = 16) -> str:
    """Format bytes as a classic side-by-side hex dump and ASCII preview."""
    lines = []
    for i in range(0, len(data), bytes_per_line):
        chunk = data[i : i + bytes_per_line]
        hex_part = " ".join(f"{b:02x}" for b in chunk)
        ascii_part = "".join(chr(b) if 32 <= b < 127 else "." for b in chunk)
        lines.append(f"  {i:04x}: {hex_part:<48} |{ascii_part}|")
    return "\n".join(lines)


def verify_single_strand(seq: str, payload_size: int = DEFAULT_PAYLOAD_SIZE):
    """Deep inspection and decoding of a single 184 nt DNA oligo sequence."""
    translator = DnaTranslator(payload_size=payload_size)
    validator = ConstraintValidator()
    obfuscator = BiologicalObfuscator()

    clean_seq = seq.strip().upper()
    print("=" * 80)
    print("                ON-THE-SPOT DNA STRAND VERIFICATION REPORT                 ")
    print("=" * 80)
    print(f"[*] Input DNA Sequence ({len(clean_seq)} nt):")
    print(f"    {clean_seq}")

    # 1. Biophysical Analysis
    gc = validator.calculate_gc_content(clean_seq)
    hp_violations = validator.find_homopolymers(clean_seq, max_allowed=3)
    max_run = 1
    run = 1
    for i in range(1, len(clean_seq)):
        if clean_seq[i] == clean_seq[i - 1]:
            run += 1
            if run > max_run:
                max_run = run
        else:
            run = 1

    print("\n--- 1. Biophysical Synthesis Validation ---")
    print(f"[+] Total Strand Length: {len(clean_seq)} bp (Budget: 150-200 bp) -> {'VALID' if 150 <= len(clean_seq) <= 200 else 'INVALID'}")
    print(f"[+] GC Content:          {gc:.2%} (Wet-Lab Limit: 40.0% - 60.0%) -> {'BALANCED' if 0.40 <= gc <= 0.60 else 'OUT OF RANGE'}")
    print(f"[+] Max Homopolymer Run: {max_run} consecutive bases -> {'PASSED (<=3)' if max_run <= 3 else 'FAILED (>3)'}")
    if hp_violations:
        for base, pos, length in hp_violations:
            print(f"    [!] Violation at pos {pos}: {base * length} ({length} bases)")

    # 2. Primer Verification
    print("\n--- 2. Primer & Sub-Region Disassembly ---")
    fwd_len = len(DEFAULT_FORWARD_PRIMER)
    rev_len = len(DEFAULT_REVERSE_PRIMER)
    has_fwd = clean_seq.startswith(DEFAULT_FORWARD_PRIMER)
    has_rev = clean_seq.endswith(DEFAULT_REVERSE_PRIMER)

    print(f"[+] 5' Forward Primer:   {clean_seq[:fwd_len]} [{'MATCH' if has_fwd else 'MISMATCH'}]")
    id_dna = clean_seq[fwd_len : fwd_len + 16]
    payload_dna = clean_seq[fwd_len + 16 : fwd_len + 16 + (payload_size * 4)]
    crc_dna = clean_seq[fwd_len + 16 + (payload_size * 4) : len(clean_seq) - rev_len]
    print(f"[+] 4-Byte Chunk ID:     {id_dna} (16 nt)")
    print(f"[+] Masked Payload:      {payload_dna} ({len(payload_dna)} nt)")
    print(f"[+] 4-Byte CRC32:        {crc_dna} (16 nt)")
    print(f"[+] 3' Reverse Primer:   {clean_seq[-rev_len:]} [{'MATCH' if has_rev else 'MISMATCH'}]")

    # 3. Disassembly & CRC32 Verification
    masked_chunk = translator.disassemble_strand(clean_seq)
    print("\n--- 3. Chunk Header & Integrity ---")
    print(f"[+] Chunk ID (Hex):      0x{masked_chunk.chunk_id:08X} (Index: {masked_chunk.chunk_index}, Salt: {masked_chunk.salt})")
    print(f"[+] Strand CRC32:        0x{masked_chunk.crc32:08X}")

    unmasked_payload = obfuscator.apply_xor_mask(masked_chunk.masked_payload, masked_chunk.chunk_id)
    crc_valid = StrandIntegrity.verify_crc32(masked_chunk.chunk_id, unmasked_payload, masked_chunk.crc32)
    print(f"[+] CRC32 Verification:  {'PASSED (Zero Mutations Detected)' if crc_valid else 'FAILED (Mutated Strand)'}")

    # 4. Decoded Data Payload
    print("\n--- 4. Decoded Binary & Raw Data Payload ---")
    bit_string = BinaryIngestor.bytes_to_binary_string(unmasked_payload)
    print(f"[+] Binary Bitstream ({len(bit_string)} bits):")
    print(f"    {bit_string}")
    print(f"\n[+] Hex & ASCII Dump ({len(unmasked_payload)} bytes):")
    print(format_hex_dump(unmasked_payload))
    print("=" * 80)


def verify_fasta_pool(fasta_file: Path, original_file: Path = None, output_file: Path = None):
    """Verify an entire FASTA file by decoding it back into a file on the spot."""
    if not fasta_file.is_file():
        print(f"[ERROR] FASTA file not found: {fasta_file}", file=sys.stderr)
        sys.exit(1)

    print("=" * 80)
    print("             ON-THE-SPOT FASTA OLIGO POOL VERIFICATION             ")
    print("=" * 80)
    print(f"[*] Reading FASTA File: {fasta_file.name}")

    sequences = DnaTranslator.import_fasta(fasta_file)
    print(f"[+] Ingested {len(sequences)} DNA sequences from FASTA.")

    pipeline = HelixPipeline(payload_size=28, parity_ratio=0.15)

    # Disassemble and count data vs parity chunks
    valid_chunks = {}
    total_data = 0
    total_parity = 0

    for seq in sequences:
        try:
            mc = pipeline.translator.disassemble_strand(seq)
            chunk = pipeline.obfuscator.unmask_chunk(mc)
            valid_chunks[chunk.chunk_index] = chunk
            # Check type from index
        except Exception as e:
            print(f"[!] Warning: sequence rejected: {e}")

    # Count data and parity chunks directly from FASTA headers or sequence count
    k_data = 0
    m_parity = 0
    with open(fasta_file, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if "type=data" in line:
                    k_data += 1
                elif "type=parity" in line:
                    m_parity += 1

    if k_data == 0:
        # Fallback if no type metadata in header
        k_data = len(sequences)
        m_parity = 0

    print(f"[+] Identified {k_data} Data Chunks + {m_parity} Parity Chunks from Pool Headers")

    temp_output = output_file or (fasta_file.parent / f"decoded_from_{fasta_file.stem}.bin")

    pipeline.decode_strands(
        received_strands=sequences,
        output_filepath=temp_output,
        total_data_chunks=k_data,
        total_parity_chunks=m_parity,
    )

    rec_hash = BinaryIngestor.compute_sha256(temp_output)
    print(f"\n[+] Reconstructed File: {temp_output}")
    print(f"[+] Size:               {temp_output.stat().st_size:,} bytes")
    print(f"[+] SHA-256 Digest:     {rec_hash}")

    if original_file and original_file.is_file():
        orig_hash = BinaryIngestor.compute_sha256(original_file)
        print(f"\n[*] Comparing with Original File: {original_file.name}")
        print(f"    Original SHA-256:  {orig_hash}")
        print(f"    Recovered SHA-256: {rec_hash}")
        if orig_hash == rec_hash:
            print("\n[SUCCESS] 100.000% BIT-PERFECT MATCH CONFIRMED ON THE SPOT!")
        else:
            print("\n[FAILURE] CHECKSUM MISMATCH!", file=sys.stderr)
            sys.exit(1)
    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description="HELIX-CORE: On-the-spot DNA sequence verification tool"
    )
    parser.add_argument(
        "--seq",
        type=str,
        help="A raw 184 nt DNA sequence string to decode and inspect immediately",
    )
    parser.add_argument(
        "--fasta",
        type=str,
        help="Path to a .fasta file containing synthesized DNA sequences",
    )
    parser.add_argument(
        "--strand-index",
        type=int,
        default=0,
        help="When --fasta is provided without --seq, inspect this specific strand index (default: 0)",
    )
    parser.add_argument(
        "--original",
        type=str,
        help="Path to the original file to compare SHA-256 hashes against",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Path to save the decoded output file",
    )

    args = parser.parse_args()

    if args.seq:
        verify_single_strand(args.seq)
    elif args.fasta:
        fasta_path = Path(args.fasta).resolve()
        if args.original:
            orig_path = Path(args.original).resolve()
            out_path = Path(args.output).resolve() if args.output else None
            verify_fasta_pool(fasta_path, original_file=orig_path, output_file=out_path)
        else:
            # Inspect specific strand from the FASTA file
            sequences = DnaTranslator.import_fasta(fasta_path)
            idx = args.strand_index
            if 0 <= idx < len(sequences):
                print(f"[*] Inspecting Strand #{idx} from {fasta_path.name}:")
                verify_single_strand(sequences[idx])
            else:
                print(f"[ERROR] Strand index {idx} out of range (0..{len(sequences)-1})", file=sys.stderr)
    else:
        # Default: inspect strand 0 from synthesis_pool.fasta if it exists
        default_fasta = Path("synthesis_pool.fasta").resolve()
        if default_fasta.is_file():
            sequences = DnaTranslator.import_fasta(default_fasta)
            if sequences:
                print("[*] No arguments provided. Verifying Strand #0 from synthesis_pool.fasta:")
                verify_single_strand(sequences[0])
                return
        parser.print_help()


if __name__ == "__main__":
    main()
