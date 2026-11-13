# HELIX-CORE: Production-Grade Bi-Directional DNA Data Storage Pipeline

**HELIX-CORE** is an industrial-strength, bi-directional DNA data storage pipeline engineered to translate arbitrary digital files (`.jpg`, `.png`, `.pdf`, `.mp4`, etc.) into biologically viable, physically synthesizable DNA oligo sequences and reconstruct the exact original file with **100% zero data loss (bit-perfect SHA-256 validation)**.

---

## 🧬 Biophysical Architecture & Wet-Lab Reality

DNA data storage differs fundamentally from traditional magnetic or optical media:
1. **Free-Floating Pool**: Synthesized oligos exist in liquid suspension with no physical tracks or sectors. Every fragment requires an immutable address identifier.
2. **Physical Synthesis Limits**: Chemical and enzymatic synthesis by providers like Twist Bioscience or Integrated DNA Technologies (IDT) degrades rapidly above 200 base pairs. HELIX-CORE enforces a target strand length of **184 nt**, strictly within the 150–200 bp synthesis limit.
3. **Sequencing Chemistry Failure Modes**:
   - **Homopolymers**: Repetitive single-base runs (`AAAA`, `CCCC`, `GGGG`, `TTTT`) induce polymerase slippage and optical/electrical saturation in Oxford Nanopore and Illumina flow cells. HELIX-CORE guarantees **no base repeats more than 3 times consecutively**.
   - **GC Content Skew**: Sequences with `<40%` or `>60%` GC content suffer from severe PCR amplification bias and secondary hairpin formation. HELIX-CORE guarantees GC content between **40.0% and 60.0%**.
4. **Stochastic Strand Dropout & Mutations**: Enzymatic degradation, PCR stochasticity, and pipetting loss routinely cause 5%–10% of oligos to vanish. HELIX-CORE integrates inter-strand **Reed-Solomon erasure coding** combined with **strand-level CRC32 checksums** to recover 100% of the original file even if >10% of strands are completely lost or corrupted.

---

## 🔬 Pipeline Workflow

```
[Raw Digital File] (.png, .pdf, .mp4)
       │
       ▼ (Phase 1)
[Binary Ingestion & Containerization] ──► Cryptographic SHA-256 Hash
       │
       ▼ (Phase 2)
[Fixed-Length Fragmentation] ──────────► 28-byte payloads (112 nt)
       │
       ▼ (Phase 2)
[Reed-Solomon Erasure Coding] ─────────► Column-wise GF(2^8) parity chunks (>=15%)
       │
       ▼ (Phase 3)
[Seeded XOR Obfuscation] ──────────────► Chunk ID seed, GC balancing (40-60%), Homopolymers (<=3)
       │
       ▼ (Phase 4)
[DNA Translation & Assembly] ──────────► 00=A, 01=C, 10=G, 11=T + Universal Primers (184 bp)
       │
       ▼ (Export)
[.fasta / .csv Synthesis Order] ───────► Compatible with Twist Bioscience / IDT
       │
  ═══════════════════════════════════════════════════════════════════════════
       WET-LAB REALITY SIMULATION: 10% Strand Loss + Point Mutations (SNPs)
  ═══════════════════════════════════════════════════════════════════════════
       │
       ▼ (Phase 5)
[Sequencing Ingestion & Primer Strip] ─► Extract Chunk ID, Payload, CRC32
       │
       ▼ (Phase 5)
[CRC32 Mutation Filtering] ────────────► Corrupted strands discarded as erasures
       │
       ▼ (Phase 5)
[Reed-Solomon Erasure Reconstruction] ──► Rebuilds missing data from surviving parity
       │
       ▼ (Phase 6)
[Reassembly & Container Unpack] ───────► Reorder by Index ID & trim zero-padding
       │
       ▼ (Validation)
[Bit-Perfect Output File] ─────────────► SHA-256 input == SHA-256 output (100% Match)
```

---

## 📐 Strand Anatomy (184 bp Oligo Budget)

Each synthesized oligo is architected as follows (5' to 3'):

| Region | Length (nt) | Length (bytes) | Description |
| :--- | :---: | :---: | :--- |
| **Forward Primer** | 20 nt | — | `5'-ACACGACGCTCTTCCGATCT-3'` (Twist/Illumina universal PCR) |
| **Chunk ID / Index** | 16 nt | 4 bytes | 32-bit big-endian index + biological screening salt |
| **Masked Payload** | 112 nt | 28 bytes | Whitened payload (XOR masked with Chunk ID seed) |
| **CRC32 Checksum** | 16 nt | 4 bytes | 32-bit checksum protecting Chunk ID + Payload |
| **Reverse Primer** | 20 nt | — | `5'-AGATCGGAAGAGCGTCGTGT-3'` (Universal PCR reverse) |
| **Total Strand** | **184 nt** | **36 bytes** | **Strictly within 150–200 bp synthesis limit** |

---

## 🚀 Quickstart & Usage

### 1. Installation

Ensure Python 3.10+ is installed along with dependencies:
```bash
pip install reedsolo
```

### 2. Launch the Web Dashboard (Recommended)

Start the interactive bi-directional web dashboard:
```bash
python dashboard.py
```
Open **`http://localhost:8080`** in any browser.

The dashboard features two dedicated tabs:
1. **Encode (File ➔ DNA)**: Upload any digital file to synthesize 184 nt DNA oligos and export FASTA/CSV.
2. **Decode (DNA ➔ File)**: Upload a `.fasta` file or paste DNA sequences to reverse the biological translation, unmask payloads, verify CRC32, and reconstruct the exact original digital file for download!

### 3. Verify DNA Sequences On-The-Spot (CLI Verification)

Inspect a single 184 nt DNA strand or verify an entire `.fasta` synthesis pool on the spot:
```bash
# Verify an individual 184 nt sequence (biophysical check, unmasking, bitstream & hex dump)
python verify_dna.py --seq ACACGACGCTCTTCCGATCTGAGGCGTATTACCATCTCGTTGGACGGACGTAATCTTAGACTAGCTGGGAAATGCGATACCTTAGTACACGGATAGATCACTAACATGCAGACATGCAGTGGATCTGAAATGTGGAATCCTGGGATAGGATGAACGGGTCTGAGAGATCGGAAGAGCGTCGTGT

# Verify an entire synthesis pool against an original file
python verify_dna.py --fasta sample_english.fasta --orig sample_english.txt
```

### 4. Run the CLI Demonstration

Run a full round-trip on any file with simulated wet-lab degradation (10% strand dropout + 2% point mutations):
```bash
python run_helix_demo.py --input test_sample.png --output recovered_sample.png --drop-rate 0.10 --mutation-rate 0.02
```

### 5. Run Test Suite

```bash
python -m unittest discover -s tests -p "test_*.py" -v
```

---

## 🧪 Object-Oriented Component Reference

- **`helix_core.BinaryIngestor`**: Ingests files using standard Python I/O (`open(file, 'rb')`), computes SHA-256 digests, converts bytes to unbroken binary strings, and packs/unpacks container streams.
- **`helix_core.FragmentationEngine`**: Fragments continuous data into uniform payloads (28 bytes) and reassembles ordered chunks.
- **`helix_core.ReedSolomonCodec`**: Column-wise Galois Field $GF(2^8)$ forward error correction. Generates parity chunks and recovers erased chunks.
- **`helix_core.StrandIntegrity`**: Computes and verifies 32-bit CRC32 checksums for mutation detection.
- **`helix_core.BiologicalObfuscator`**: Applies seeded PRNG XOR bitmasking and adaptive candidate screening to eliminate homopolymers ($>3$) and maintain 40%–60% GC content.
- **`helix_core.ConstraintValidator`**: Biophysical synthesis validation engine checking GC content, homopolymers, and oligo length constraints.
- **`helix_core.DnaTranslator`**: Maps dibits (`00=A`, `01=C`, `10=G`, `11=T`) to nucleotides, encodes/decodes self-describing 32-bit wire IDs, and attaches universal PCR primers.
- **`helix_core.HelixPipeline`**: High-level coordinator executing end-to-end encoding, simulated wet-lab degradation, and bit-perfect reconstruction.


