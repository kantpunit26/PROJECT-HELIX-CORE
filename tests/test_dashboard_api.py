"""
Integration test for HELIX-CORE Dashboard Server API.
"""

import json
from pathlib import Path
import subprocess
import time
import unittest
import urllib.request
import uuid


class TestDashboardAPI(unittest.TestCase):
    """Test dashboard HTTP server endpoints."""

    @classmethod
    def setUpClass(cls):
        cls.server_proc = subprocess.Popen(
            ["C:\\Users\\Punit\\AppData\\Local\\Python\\bin\\python.exe", "dashboard.py"]
        )
        time.sleep(1.5)

    @classmethod
    def tearDownClass(cls):
        cls.server_proc.terminate()
        cls.server_proc.wait()

    def test_get_dashboard_html(self):
        """Verify GET / serves the dashboard web page."""
        with urllib.request.urlopen("http://localhost:8080/") as res:
            self.assertEqual(res.status, 200)
            html = res.read().decode("utf-8")
            self.assertIn("HELIX-CORE", html)
            self.assertIn("Select or Drop Your Digital File", html)

    def test_post_upload_file_api(self):
        """Verify POST /api/upload ingests a real file and returns bit-perfect pipeline metrics."""
        boundary = "----TestBoundary" + uuid.uuid4().hex
        with open("test_sample.png", "rb") as fh:
            file_bytes = fh.read()

        parts = [
            f"--{boundary}\r\nContent-Disposition: form-data; name=\"file\"; filename=\"test_sample.png\"\r\nContent-Type: image/png\r\n\r\n".encode("utf-8")
            + file_bytes
            + b"\r\n",
            f"--{boundary}\r\nContent-Disposition: form-data; name=\"parity_ratio\"\r\n\r\n0.15\r\n".encode("utf-8"),
            f"--{boundary}\r\nContent-Disposition: form-data; name=\"drop_rate\"\r\n\r\n0.10\r\n".encode("utf-8"),
            f"--{boundary}\r\nContent-Disposition: form-data; name=\"mutation_rate\"\r\n\r\n0.02\r\n".encode("utf-8"),
            f"--{boundary}--\r\n".encode("utf-8"),
        ]
        body = b"".join(parts)

        req = urllib.request.Request(
            "http://localhost:8080/api/upload",
            data=body,
            headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
        )

        with urllib.request.urlopen(req) as res:
            self.assertEqual(res.status, 200)
            data = json.loads(res.read().decode("utf-8"))
            self.assertTrue(data["success"])
            self.assertTrue(data["is_bit_perfect"])
            self.assertEqual(data["original_sha256"], data["recovered_sha256"])
            self.assertEqual(data["strand_length_bp"], 184)
            self.assertGreater(len(data["sample_strands"]), 0)

    def test_post_decode_dna_api(self):
        """Verify POST /api/decode ingests DNA sequence strings and restores the digital file."""
        fasta_path = Path("dashboard_storage/exports/synthesis_pool.fasta")
        if not fasta_path.is_file():
            # Run upload first to generate fasta
            self.test_post_upload_file_api()

        fasta_text = fasta_path.read_text(encoding="utf-8")

        payload = json.dumps({"dna_text": fasta_text}).encode("utf-8")
        req = urllib.request.Request(
            "http://localhost:8080/api/decode",
            data=payload,
            headers={"Content-Type": "application/json"},
        )

        with urllib.request.urlopen(req) as res:
            self.assertEqual(res.status, 200)
            data = json.loads(res.read().decode("utf-8"))
            self.assertTrue(data["success"])
            self.assertGreater(data["file_size"], 0)
            self.assertGreater(data["total_sequences_ingested"], 0)
            self.assertIn("filename", data)
            self.assertIn("sha256", data)


if __name__ == "__main__":
    unittest.main()

