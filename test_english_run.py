from pathlib import Path
from helix_core.pipeline import HelixPipeline
from helix_core.ingestion import BinaryIngestor

input_file = Path("sample_english.txt")
output_file = Path("recovered_sample_english.txt")
fasta_file = Path("sample_english.fasta")

pipeline = HelixPipeline(payload_size=28, parity_ratio=0.15)
strands, meta = pipeline.encode_file(input_file)
print(f"Encoded into {len(strands)} strands ({meta['data_chunks']} data + {meta['parity_chunks']} parity)")

pipeline.translator.export_fasta(strands, fasta_file)

# Decode from strands
pipeline.decode_strands(
    received_strands=strands,
    output_filepath=output_file,
    total_data_chunks=meta["data_chunks"],
    total_parity_chunks=meta["parity_chunks"],
)

with open(output_file, "r", encoding="utf-8") as f:
    text = f.read()

print("--- Recovered File Content ---")
print(text)
print("--- End Recovered File Content ---")

orig_hash = BinaryIngestor.compute_sha256(input_file)
rec_hash = BinaryIngestor.compute_sha256(output_file)
print("Orig SHA-256:           ", orig_hash)
print("Rec  SHA-256:           ", rec_hash)
print("Hashes Match:           ", orig_hash == rec_hash)
assert orig_hash == rec_hash

# Autonomous Decode from FASTA (Zero external metadata passed!)
fasta_seqs = pipeline.translator.import_fasta(fasta_file)
recovered_bytes, auto_meta = pipeline.decode_strands_to_bytes(fasta_seqs)

print("\n--- Autonomous Reconstructed Metadata ---")
print(f"Filename from DNA:      {auto_meta['filename']}")
print(f"Data Chunks Inferred:   {auto_meta['total_data_chunks']}")
print(f"Parity Chunks Inferred: {auto_meta['total_parity_chunks']}")
print(f"Bytes Decoded:          {len(recovered_bytes)}")

auto_output_file = Path(f"autonomous_{auto_meta['filename']}")
with open(auto_output_file, "wb") as f:
    f.write(recovered_bytes)

auto_hash = BinaryIngestor.compute_sha256(auto_output_file)
print("Auto Rec SHA-256:       ", auto_hash)
print("Autonomous Match:       ", orig_hash == auto_hash)
assert orig_hash == auto_hash
print("\n[SUCCESS] 100% Autonomous DNA-to-File Restoration Confirmed!")
