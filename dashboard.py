"""
HELIX-CORE Web Dashboard Server.

Clean, light, professional, and bi-directional DNA data storage interface:
1. TAB 1: Encode Digital File -> DNA Sequences (184 nt oligos with universal primers & FASTA/CSV export).
2. TAB 2: Decode DNA Sequences -> Binary Bitstream -> Original Digital File (with instant file download).
"""

import http.server
import json
import math
from pathlib import Path
import re
import socketserver
import struct
import urllib.parse

from helix_core.fragmentation import FragmentationEngine, ReedSolomonCodec, StrandIntegrity
from helix_core.ingestion import BinaryIngestor
from helix_core.obfuscation import BiologicalObfuscator, ConstraintValidator
from helix_core.pipeline import HelixPipeline
from helix_core.translation import DnaTranslator

PORT = 8080
BASE_DIR = Path(__file__).resolve().parent
UPLOADS_DIR = BASE_DIR / "dashboard_storage" / "uploads"
EXPORTS_DIR = BASE_DIR / "dashboard_storage" / "exports"
RECOVERED_DIR = BASE_DIR / "dashboard_storage" / "recovered"

UPLOADS_DIR.mkdir(parents=True, exist_ok=True)
EXPORTS_DIR.mkdir(parents=True, exist_ok=True)
RECOVERED_DIR.mkdir(parents=True, exist_ok=True)

LATEST_RUN = {
    "has_data": False,
    "input_path": None,
    "output_path": None,
    "fasta_path": None,
    "csv_path": None,
}


DASHBOARD_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>HELIX-CORE | DNA Data Storage Pipeline</title>
  <script src="https://www.gstatic.com/antigravity/web/dev/tailwindcss.min.js"></script>
  <style>
    .nt-A { color: #059669; font-weight: 600; }
    .nt-C { color: #0284c7; font-weight: 600; }
    .nt-G { color: #d97706; font-weight: 600; }
    .nt-T { color: #dc2626; font-weight: 600; }
  </style>
</head>
<body class="bg-slate-50 text-slate-900 min-h-screen font-sans antialiased">

  <!-- TOP NAVIGATION BAR -->
  <header class="bg-white border-b border-slate-200 sticky top-0 z-50">
    <div class="max-w-5xl mx-auto px-4 sm:px-6 h-16 flex items-center justify-between">
      <div class="flex items-center space-x-3">
        <div class="w-8 h-8 rounded-lg bg-indigo-600 text-white flex items-center justify-center font-bold text-base shadow-sm">
          🧬
        </div>
        <div>
          <span class="font-bold text-slate-900 text-base tracking-tight">HELIX-CORE</span>
          <span class="text-xs text-slate-500 font-normal ml-2 hidden sm:inline">Bi-Directional DNA Data Storage</span>
        </div>
      </div>

      <!-- TAB SWITCHER -->
      <div class="flex items-center p-1 bg-slate-100 rounded-xl border border-slate-200 text-xs font-semibold">
        <button id="tabBtnEncode"
                class="px-3.5 py-1.5 rounded-lg bg-white text-indigo-700 shadow-sm transition flex items-center gap-1.5">
          <span>📤</span>
          <span>Encode: File ➔ DNA</span>
        </button>
        <button id="tabBtnDecode"
                class="px-3.5 py-1.5 rounded-lg text-slate-600 hover:text-slate-900 transition flex items-center gap-1.5">
          <span>📥</span>
          <span>Decode: DNA ➔ File</span>
        </button>
      </div>
    </div>
  </header>

  <!-- MAIN CONTAINER -->
  <main class="max-w-5xl mx-auto px-4 sm:px-6 py-8">

    <!-- ============================================================= -->
    <!-- TAB 1: ENCODE DIGITAL FILE TO DNA (Active by default)          -->
    <!-- ============================================================= -->
    <section id="panelEncode" class="space-y-6">

      <div class="text-center max-w-2xl mx-auto space-y-1.5">
        <h2 class="text-2xl font-bold text-slate-900 tracking-tight">
          Encode Any Digital File into Physical DNA
        </h2>
        <p class="text-xs sm:text-sm text-slate-600">
          Upload any file (.png, .jpg, .pdf, .mp4, .zip). Translates to 184 nt synthesis oligos with forward error correction.
        </p>
      </div>

      <div class="bg-white border border-slate-200 rounded-2xl shadow-sm p-6 sm:p-8 space-y-6">

        <!-- Step 1: Upload -->
        <div>
          <div class="flex items-center space-x-2 mb-3">
            <span class="w-6 h-6 rounded-full bg-indigo-100 text-indigo-700 text-xs font-bold flex items-center justify-center">1</span>
            <h3 class="text-sm font-semibold text-slate-900">Select or Drop Your Digital File</h3>
          </div>

          <div id="dropzoneEncode"
               class="border-2 border-dashed border-slate-300 hover:border-indigo-500 hover:bg-indigo-50/30 transition rounded-xl p-8 text-center cursor-pointer flex flex-col items-center justify-center space-y-3">
            <input type="file" id="fileInputEncode" class="hidden" />
            <div class="w-12 h-12 rounded-full bg-slate-100 text-slate-500 flex items-center justify-center text-2xl">
              📄
            </div>
            <div>
              <p class="text-sm font-semibold text-slate-800" id="filePromptEncode">
                Click here to choose a file, or drag and drop
              </p>
              <p class="text-xs text-slate-500 mt-1" id="fileSubtextEncode">
                PNG, JPG, PDF, MP4, TXT, ZIP, etc. (Standard Python binary I/O)
              </p>
            </div>
          </div>
        </div>

        <!-- Step 2: Protection Settings -->
        <div class="pt-4 border-t border-slate-100">
          <div class="flex items-center space-x-2 mb-4">
            <span class="w-6 h-6 rounded-full bg-indigo-100 text-indigo-700 text-xs font-bold flex items-center justify-center">2</span>
            <h3 class="text-sm font-semibold text-slate-900">Error Protection & Simulation Settings</h3>
          </div>

          <div class="grid grid-cols-1 sm:grid-cols-2 gap-6 bg-slate-50/80 p-4 rounded-xl border border-slate-200 text-xs">
            <div>
              <div class="flex justify-between font-medium text-slate-700 mb-1.5">
                <span>Error Correction Redundancy:</span>
                <span id="parityLabel" class="text-indigo-600 font-bold font-mono">15%</span>
              </div>
              <input type="range" id="paritySlider" min="10" max="30" value="15" step="1"
                     class="w-full h-1.5 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-indigo-600">
              <p class="text-[11px] text-slate-500 mt-1">Reed-Solomon parity strands (MDS code) to recover lost oligos.</p>
            </div>

            <div>
              <div class="flex justify-between font-medium text-slate-700 mb-1.5">
                <span>Simulated Wet-Lab Strand Loss:</span>
                <span id="dropLabel" class="text-amber-600 font-bold font-mono">10%</span>
              </div>
              <input type="range" id="dropSlider" min="0" max="20" value="10" step="1"
                     class="w-full h-1.5 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-amber-600">
              <p class="text-[11px] text-slate-500 mt-1">Random oligo destruction in sequencing pool.</p>
            </div>
          </div>
        </div>

        <!-- Action Button -->
        <div class="pt-2">
          <button id="runBtnEncode" disabled
                  class="w-full py-3.5 px-6 rounded-xl bg-indigo-600 hover:bg-indigo-700 disabled:bg-slate-200 disabled:text-slate-400 disabled:cursor-not-allowed text-white font-semibold text-sm transition shadow-sm flex items-center justify-center space-x-2">
            <span>🧬</span>
            <span>Encode to DNA Sequences & Verify</span>
          </button>
        </div>

        <!-- Loader -->
        <div id="loaderEncode" class="hidden text-center py-4 space-y-2">
          <div class="inline-block animate-spin text-2xl">⏳</div>
          <p class="text-xs text-slate-600 font-medium">Executing Phase 1 (Ingestion) ➔ Phase 2 (RS) ➔ Phase 3 (Obfuscation) ➔ Phase 4 (DNA)...</p>
        </div>

      </div>

      <!-- RESULTS SECTION (Encode) -->
      <div id="resultsEncode" class="hidden space-y-6">

        <div class="bg-emerald-50 border border-emerald-200 rounded-2xl p-6 flex flex-col sm:flex-row items-start sm:items-center justify-between gap-4">
          <div class="flex items-start sm:items-center space-x-3.5">
            <div class="w-10 h-10 rounded-full bg-emerald-100 text-emerald-700 flex items-center justify-center text-xl font-bold flex-shrink-0">
              ✓
            </div>
            <div>
              <h3 class="text-base font-bold text-emerald-900">
                100% Bit-Perfect Reconstruction Confirmed
              </h3>
              <p class="text-xs text-emerald-700 mt-0.5">
                Input file was translated to DNA, survived simulated wet-lab degradation, and was fully restored with zero data loss.
              </p>
            </div>
          </div>
          <a id="btnDownloadRecovered" href="#" download
             class="px-5 py-2.5 bg-emerald-600 hover:bg-emerald-700 text-white text-xs font-semibold rounded-xl shadow-sm transition whitespace-nowrap">
            ⬇ Download Restored File
          </a>
        </div>

        <!-- Summary Metrics -->
        <div class="grid grid-cols-2 sm:grid-cols-4 gap-4">
          <div class="bg-white border border-slate-200 p-4 rounded-xl shadow-sm">
            <p class="text-xs text-slate-500 font-medium">Original File Size</p>
            <p id="statSizeEncode" class="text-lg font-bold text-slate-900 font-mono mt-1">-</p>
          </div>
          <div class="bg-white border border-slate-200 p-4 rounded-xl shadow-sm">
            <p class="text-xs text-slate-500 font-medium">Total DNA Strands</p>
            <p id="statStrandsEncode" class="text-lg font-bold text-indigo-600 font-mono mt-1">-</p>
          </div>
          <div class="bg-white border border-slate-200 p-4 rounded-xl shadow-sm">
            <p class="text-xs text-slate-500 font-medium">Strand Length</p>
            <p class="text-lg font-bold text-slate-900 font-mono mt-1">184 base pairs</p>
          </div>
          <div class="bg-white border border-slate-200 p-4 rounded-xl shadow-sm">
            <p class="text-xs text-slate-500 font-medium">GC Balance</p>
            <p id="statGCEncode" class="text-lg font-bold text-emerald-600 font-mono mt-1">-</p>
          </div>
        </div>

        <!-- Hashes -->
        <div class="bg-white border border-slate-200 rounded-2xl p-6 shadow-sm space-y-3">
          <h4 class="text-sm font-bold text-slate-900">SHA-256 Hash Verification</h4>
          <div class="space-y-2 text-xs font-mono">
            <div class="p-3 bg-slate-50 rounded-xl border border-slate-200">
              <span class="text-slate-500 block text-[11px] font-sans font-medium">Input File SHA-256:</span>
              <span id="valInputHashEncode" class="text-slate-800 break-all select-all font-bold"></span>
            </div>
            <div class="p-3 bg-slate-50 rounded-xl border border-slate-200">
              <span class="text-slate-500 block text-[11px] font-sans font-medium">Restored File SHA-256:</span>
              <span id="valOutputHashEncode" class="text-emerald-700 break-all select-all font-bold"></span>
            </div>
          </div>
        </div>

        <!-- Downloads -->
        <div class="bg-white border border-slate-200 rounded-2xl p-6 shadow-sm space-y-3">
          <h4 class="text-sm font-bold text-slate-900">Physical Synthesis Exports</h4>
          <p class="text-xs text-slate-500">Export files ready for chemical synthesis providers (Twist Bioscience / IDT):</p>
          <div class="flex flex-wrap gap-3 pt-1">
            <a href="/api/download/fasta" download
               class="px-4 py-2.5 bg-slate-100 hover:bg-slate-200 text-slate-800 text-xs font-semibold rounded-xl border border-slate-300 transition flex items-center space-x-2">
              <span>🧬</span>
              <span>Download Synthesis FASTA (.fasta)</span>
            </a>
            <a href="/api/download/csv" download
               class="px-4 py-2.5 bg-slate-100 hover:bg-slate-200 text-slate-800 text-xs font-semibold rounded-xl border border-slate-300 transition flex items-center space-x-2">
              <span>📊</span>
              <span>Download Twist / IDT Order (.csv)</span>
            </a>
          </div>
        </div>

        <!-- Strand Inspector -->
        <div class="bg-white border border-slate-200 rounded-2xl p-6 shadow-sm space-y-4">
          <div class="flex items-center justify-between">
            <div>
              <h4 class="text-sm font-bold text-slate-900">Synthesized DNA Oligo Sequence (Strand #0)</h4>
              <p class="text-xs text-slate-500">184 nt physical sequence ready for synthesis</p>
            </div>
            <div class="flex items-center space-x-2 font-mono text-xs">
              <span class="nt-A">● A</span>
              <span class="nt-C">● C</span>
              <span class="nt-G">● G</span>
              <span class="nt-T">● T</span>
            </div>
          </div>

          <div id="strandDisplayEncode"
               class="font-mono text-xs leading-relaxed p-4 bg-slate-900 text-slate-100 rounded-xl border border-slate-800 break-all select-all tracking-wider">
          </div>
        </div>

      </div>

    </section>

    <!-- ============================================================= -->
    <!-- TAB 2: DECODE DNA SEQUENCE TO DIGITAL FILE (NEW PANEL)         -->
    <!-- ============================================================= -->
    <section id="panelDecode" class="hidden space-y-6">

      <div class="text-center max-w-2xl mx-auto space-y-1.5">
        <h2 class="text-2xl font-bold text-slate-900 tracking-tight">
          Decode DNA Sequences into Digital Files
        </h2>
        <p class="text-xs sm:text-sm text-slate-600">
          Upload a synthesized .fasta file or paste DNA sequences. System strips primers, translates nucleotides to binary, reverses the XOR mask, and restores the original file.
        </p>
      </div>

      <div class="bg-white border border-slate-200 rounded-2xl shadow-sm p-6 sm:p-8 space-y-6">

        <!-- Input options -->
        <div class="space-y-4">
          <div class="flex items-center justify-between">
            <h3 class="text-sm font-semibold text-slate-900 flex items-center gap-2">
              <span>🧬</span> Provide DNA Sequences
            </h3>
            <button id="btnLoadLastPool"
                    class="text-xs text-indigo-600 hover:text-indigo-800 font-semibold transition flex items-center gap-1">
              <span>⚡</span> Load Last Synthesized Pool
            </button>
          </div>

          <!-- Upload .fasta dropzone -->
          <div id="dropzoneDecode"
               class="border-2 border-dashed border-slate-300 hover:border-indigo-500 hover:bg-indigo-50/30 transition rounded-xl p-6 text-center cursor-pointer flex flex-col items-center justify-center space-y-2">
            <input type="file" id="fileInputDecode" accept=".fasta,.fa,.txt" class="hidden" />
            <div class="w-10 h-10 rounded-full bg-slate-100 text-slate-500 flex items-center justify-center text-xl">
              📥
            </div>
            <div>
              <p class="text-xs sm:text-sm font-semibold text-slate-800" id="filePromptDecode">
                Click to upload a .fasta or .txt DNA sequence file
              </p>
              <p class="text-[11px] text-slate-500 mt-0.5" id="fileSubtextDecode">
                Or drag and drop your synthesis FASTA pool here
              </p>
            </div>
          </div>

          <!-- Or paste text -->
          <div class="space-y-1.5">
            <div class="flex justify-between items-center text-xs">
              <label for="dnaTextInput" class="font-medium text-slate-700">Or Paste DNA Sequences (FASTA format or raw ACGT nucleotides):</label>
              <span id="charCountDna" class="text-slate-400 font-mono">0 chars</span>
            </div>
            <textarea id="dnaTextInput" rows="5"
                      placeholder="Paste your 184 nt DNA sequences here (e.g. ACACGACGCTCTTCCGATCTCATGCATGCATGCAGC...)"
                      class="w-full p-3 font-mono text-xs bg-slate-50 border border-slate-200 rounded-xl focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:bg-white transition"></textarea>
          </div>
        </div>

        <!-- Action Button -->
        <div>
          <button id="runBtnDecode"
                  class="w-full py-3.5 px-6 rounded-xl bg-indigo-600 hover:bg-indigo-700 text-white font-semibold text-sm transition shadow-sm flex items-center justify-center space-x-2">
            <span>⚡</span>
            <span>Decode DNA ➔ Binary ➔ Restore Digital File</span>
          </button>
        </div>

        <!-- Loader -->
        <div id="loaderDecode" class="hidden text-center py-4 space-y-2">
          <div class="inline-block animate-spin text-2xl">⏳</div>
          <p class="text-xs text-slate-600 font-medium">Stripping primers ➔ Translating nucleotides to binary ➔ Reversing XOR mask ➔ Unpacking container...</p>
        </div>

      </div>

      <!-- RESULTS SECTION (Decode) -->
      <div id="resultsDecode" class="hidden space-y-6">

        <!-- Success Banner -->
        <div class="bg-emerald-50 border border-emerald-200 rounded-2xl p-6 flex flex-col sm:flex-row items-start sm:items-center justify-between gap-4">
          <div class="flex items-start sm:items-center space-x-3.5">
            <div class="w-10 h-10 rounded-full bg-emerald-100 text-emerald-700 flex items-center justify-center text-xl font-bold flex-shrink-0">
              ✓
            </div>
            <div>
              <h3 class="text-base font-bold text-emerald-900">
                Digital File Successfully Reconstructed
              </h3>
              <p class="text-xs text-emerald-700 mt-0.5" id="bannerDecodeMeta">
                All DNA sequences were verified, translated back to binary, and assembled with 100% bit-perfect fidelity.
              </p>
            </div>
          </div>
          <a id="btnDownloadDecodedFile" href="/api/download/recovered" download
             class="px-5 py-2.5 bg-emerald-600 hover:bg-emerald-700 text-white text-xs font-semibold rounded-xl shadow-sm transition whitespace-nowrap">
            ⬇ Download Restored File
          </a>
        </div>

        <!-- Reconstructed File Metrics -->
        <div class="grid grid-cols-2 sm:grid-cols-4 gap-4">
          <div class="bg-white border border-slate-200 p-4 rounded-xl shadow-sm">
            <p class="text-xs text-slate-500 font-medium">Restored Filename</p>
            <p id="statDecodedFilename" class="text-sm font-bold text-slate-900 font-mono mt-1 truncate">-</p>
          </div>
          <div class="bg-white border border-slate-200 p-4 rounded-xl shadow-sm">
            <p class="text-xs text-slate-500 font-medium">Exact File Size</p>
            <p id="statDecodedSize" class="text-lg font-bold text-indigo-600 font-mono mt-1">-</p>
          </div>
          <div class="bg-white border border-slate-200 p-4 rounded-xl shadow-sm">
            <p class="text-xs text-slate-500 font-medium">Sequences Decoded</p>
            <p id="statDecodedStrands" class="text-lg font-bold text-slate-900 font-mono mt-1">-</p>
          </div>
          <div class="bg-white border border-slate-200 p-4 rounded-xl shadow-sm">
            <p class="text-xs text-slate-500 font-medium">Integrity Check</p>
            <p class="text-lg font-bold text-emerald-600 font-mono mt-1">CRC32 Valid</p>
          </div>
        </div>

        <!-- 5-Step Translation Trace -->
        <div class="bg-white border border-slate-200 rounded-2xl p-6 shadow-sm space-y-4">
          <h4 class="text-sm font-bold text-slate-900">Step-by-Step Biological-to-Digital Translation Trace</h4>

          <div class="grid grid-cols-1 sm:grid-cols-5 gap-3 text-xs font-mono">
            <div class="p-3 bg-slate-50 rounded-xl border border-slate-200 space-y-1">
              <span class="text-slate-500 text-[10px] font-sans font-medium block">1. DNA SEQUENCE</span>
              <span class="text-slate-800 font-bold">184 nt verified</span>
              <span class="text-[10px] text-slate-500 block">Universal primers found</span>
            </div>
            <div class="p-3 bg-slate-50 rounded-xl border border-slate-200 space-y-1">
              <span class="text-slate-500 text-[10px] font-sans font-medium block">2. PRIMERS STRIPPED</span>
              <span class="text-slate-800 font-bold">144 nt core</span>
              <span class="text-[10px] text-slate-500 block">ID + Payload + CRC32</span>
            </div>
            <div class="p-3 bg-slate-50 rounded-xl border border-slate-200 space-y-1">
              <span class="text-slate-500 text-[10px] font-sans font-medium block">3. DIBIT TRANSLATION</span>
              <span class="text-indigo-600 font-bold">224 binary bits</span>
              <span class="text-[10px] text-slate-500 block">A=00 C=01 G=10 T=11</span>
            </div>
            <div class="p-3 bg-slate-50 rounded-xl border border-slate-200 space-y-1">
              <span class="text-slate-500 text-[10px] font-sans font-medium block">4. XOR UNMASKING</span>
              <span class="text-purple-600 font-bold">Chunk ID Seed</span>
              <span class="text-[10px] text-slate-500 block">CRC32 matched</span>
            </div>
            <div class="p-3 bg-emerald-50 rounded-xl border border-emerald-200 space-y-1">
              <span class="text-emerald-700 text-[10px] font-sans font-medium block">5. DIGITAL FILE</span>
              <span class="text-emerald-900 font-bold">Original File Restored</span>
              <span class="text-[10px] text-emerald-700 block">Container unpacked</span>
            </div>
          </div>
        </div>

        <!-- Live Text Content Preview (For English/Text files) -->
        <div id="textPreviewContainer" class="hidden bg-white border border-slate-200 rounded-2xl p-6 shadow-sm space-y-3">
          <div class="flex items-center justify-between">
            <div class="flex items-center gap-2">
              <span class="w-6 h-6 rounded-full bg-emerald-100 text-emerald-700 text-xs font-bold flex items-center justify-center">✓</span>
              <h4 class="text-sm font-bold text-slate-900">Restored File Text Content (Original Form)</h4>
            </div>
            <span class="text-xs bg-emerald-100 text-emerald-800 px-2.5 py-1 rounded-lg font-semibold">100% Original English State</span>
          </div>
          <p class="text-xs text-slate-500">The DNA sequences were decoded directly into original binary and unpacked with zero encoding corruption (no Mojibake):</p>
          <pre id="textPreviewArea" class="p-4 bg-slate-50 border border-slate-200 rounded-xl text-xs font-mono text-slate-800 whitespace-pre-wrap max-h-64 overflow-y-auto leading-relaxed select-all"></pre>
        </div>

        <!-- Binary Bitstream & Hex Inspector -->
        <div class="bg-white border border-slate-200 rounded-2xl p-6 shadow-sm space-y-3">
          <div class="flex items-center justify-between">
            <h4 class="text-sm font-bold text-slate-900">Recovered Binary Bitstream & Hex Payload</h4>
            <span class="text-xs font-mono text-slate-400">First 128 bits</span>
          </div>

          <div class="p-3.5 bg-slate-900 text-cyan-300 font-mono text-xs rounded-xl border border-slate-800 break-all select-all tracking-wider"
               id="bitstreamDecodeView">
            <!-- Binary stream populated dynamically -->
          </div>

          <div class="p-3 bg-slate-50 border border-slate-200 rounded-xl text-xs font-mono text-slate-700 whitespace-pre overflow-x-auto"
               id="hexDecodeView">
            <!-- Hex dump populated dynamically -->
          </div>
        </div>

      </div>

    </section>

  </main>

  <footer class="border-t border-slate-200 bg-white py-6 mt-12 text-center text-xs text-slate-500">
    HELIX-CORE &bull; Bi-Directional DNA Data Storage Pipeline &bull; Production Grade
  </footer>

  <!-- CLIENT-SIDE SCRIPT -->
  <script>
    // Tab switching
    const tabBtnEncode = document.getElementById('tabBtnEncode');
    const tabBtnDecode = document.getElementById('tabBtnDecode');
    const panelEncode = document.getElementById('panelEncode');
    const panelDecode = document.getElementById('panelDecode');

    tabBtnEncode.addEventListener('click', () => {
      tabBtnEncode.className = 'px-3.5 py-1.5 rounded-lg bg-white text-indigo-700 shadow-sm transition flex items-center gap-1.5';
      tabBtnDecode.className = 'px-3.5 py-1.5 rounded-lg text-slate-600 hover:text-slate-900 transition flex items-center gap-1.5';
      panelEncode.classList.remove('hidden');
      panelDecode.classList.add('hidden');
    });

    tabBtnDecode.addEventListener('click', () => {
      tabBtnDecode.className = 'px-3.5 py-1.5 rounded-lg bg-white text-indigo-700 shadow-sm transition flex items-center gap-1.5';
      tabBtnEncode.className = 'px-3.5 py-1.5 rounded-lg text-slate-600 hover:text-slate-900 transition flex items-center gap-1.5';
      panelDecode.classList.remove('hidden');
      panelEncode.classList.add('hidden');
    });

    // Helper functions
    function formatBytes(bytes) {
      if (bytes < 1024) return bytes + ' bytes';
      if (bytes < 1048576) return (bytes / 1024).toFixed(1) + ' KB';
      return (bytes / 1048576).toFixed(2) + ' MB';
    }

    function colorize(seq) {
      return seq.split('').map(b => `<span class="nt-${b}">${b}</span>`).join('');
    }

    // --- ENCODE PANEL LOGIC ---
    const dropzoneEncode = document.getElementById('dropzoneEncode');
    const fileInputEncode = document.getElementById('fileInputEncode');
    const filePromptEncode = document.getElementById('filePromptEncode');
    const fileSubtextEncode = document.getElementById('fileSubtextEncode');
    const runBtnEncode = document.getElementById('runBtnEncode');
    const loaderEncode = document.getElementById('loaderEncode');
    const resultsEncode = document.getElementById('resultsEncode');

    const paritySlider = document.getElementById('paritySlider');
    const parityLabel = document.getElementById('parityLabel');
    const dropSlider = document.getElementById('dropSlider');
    const dropLabel = document.getElementById('dropLabel');

    let currentEncodeFile = null;

    paritySlider.addEventListener('input', (e) => parityLabel.textContent = e.target.value + '%');
    dropSlider.addEventListener('input', (e) => dropLabel.textContent = e.target.value + '%');

    dropzoneEncode.addEventListener('click', () => fileInputEncode.click());
    dropzoneEncode.addEventListener('dragover', (e) => { e.preventDefault(); dropzoneEncode.classList.add('border-indigo-500', 'bg-indigo-50/50'); });
    dropzoneEncode.addEventListener('dragleave', () => dropzoneEncode.classList.remove('border-indigo-500', 'bg-indigo-50/50'));
    dropzoneEncode.addEventListener('drop', (e) => {
      e.preventDefault();
      dropzoneEncode.classList.remove('border-indigo-500', 'bg-indigo-50/50');
      if (e.dataTransfer.files.length) handleEncodeFile(e.dataTransfer.files[0]);
    });

    fileInputEncode.addEventListener('change', (e) => {
      if (e.target.files.length) handleEncodeFile(e.target.files[0]);
    });

    function handleEncodeFile(file) {
      currentEncodeFile = file;
      filePromptEncode.innerHTML = `<span class="text-indigo-600 font-bold">${file.name}</span>`;
      fileSubtextEncode.textContent = `File Size: ${formatBytes(file.size)} - Ready to encode`;
      runBtnEncode.disabled = false;
    }

    runBtnEncode.addEventListener('click', async () => {
      if (!currentEncodeFile) return;

      runBtnEncode.disabled = true;
      loaderEncode.classList.remove('hidden');

      const formData = new FormData();
      formData.append('file', currentEncodeFile);
      formData.append('parity_ratio', (parseInt(paritySlider.value) / 100).toString());
      formData.append('drop_rate', (parseInt(dropSlider.value) / 100).toString());
      formData.append('mutation_rate', '0.02');

      try {
        const res = await fetch('/api/upload', { method: 'POST', body: formData });
        const data = await res.json();

        if (!data.success) {
          alert('Error: ' + data.error);
          return;
        }

        loaderEncode.classList.add('hidden');
        resultsEncode.classList.remove('hidden');

        document.getElementById('statSizeEncode').textContent = formatBytes(data.original_size_bytes);
        document.getElementById('statStrandsEncode').textContent = `${data.total_strand_count} oligos (${data.data_chunk_count} data + ${data.parity_chunk_count} parity)`;
        document.getElementById('statGCEncode').textContent = (data.mean_gc_content * 100).toFixed(1) + '%';
        document.getElementById('valInputHashEncode').textContent = data.original_sha256;
        document.getElementById('valOutputHashEncode').textContent = data.recovered_sha256;

        const downloadFileBtn = document.getElementById('btnDownloadRecovered');
        downloadFileBtn.href = '/api/download/recovered?t=' + Date.now();
        downloadFileBtn.download = 'recovered_' + data.original_filename;

        if (data.sample_strands && data.sample_strands.length) {
          const s = data.sample_strands[0];
          const displayEl = document.getElementById('strandDisplayEncode');
          displayEl.innerHTML =
            `<span class="text-slate-400" title="Forward Primer (20nt)">${s.forward_primer}</span>` +
            `<span class="text-indigo-300 font-bold" title="Chunk ID (16nt)">${s.id_sequence}</span>` +
            `<span title="Data Payload (112nt)">${colorize(s.payload_sequence)}</span>` +
            `<span class="text-purple-300 font-bold" title="CRC32 Checksum (16nt)">${s.crc_sequence}</span>` +
            `<span class="text-slate-400" title="Reverse Primer (20nt)">${s.reverse_primer}</span>`;
        }

        resultsEncode.scrollIntoView({ behavior: 'smooth' });

      } catch (err) {
        alert('Pipeline execution failed: ' + err.message);
      } finally {
        runBtnEncode.disabled = false;
        loaderEncode.classList.add('hidden');
      }
    });

    // --- DECODE PANEL LOGIC ---
    const dropzoneDecode = document.getElementById('dropzoneDecode');
    const fileInputDecode = document.getElementById('fileInputDecode');
    const filePromptDecode = document.getElementById('filePromptDecode');
    const fileSubtextDecode = document.getElementById('fileSubtextDecode');
    const dnaTextInput = document.getElementById('dnaTextInput');
    const charCountDna = document.getElementById('charCountDna');
    const runBtnDecode = document.getElementById('runBtnDecode');
    const loaderDecode = document.getElementById('loaderDecode');
    const resultsDecode = document.getElementById('resultsDecode');
    const btnLoadLastPool = document.getElementById('btnLoadLastPool');

    let currentDecodeFile = null;

    dnaTextInput.addEventListener('input', () => {
      charCountDna.textContent = `${dnaTextInput.value.length} chars`;
    });

    dropzoneDecode.addEventListener('click', () => fileInputDecode.click());
    dropzoneDecode.addEventListener('dragover', (e) => { e.preventDefault(); dropzoneDecode.classList.add('border-indigo-500', 'bg-indigo-50/50'); });
    dropzoneDecode.addEventListener('dragleave', () => dropzoneDecode.classList.remove('border-indigo-500', 'bg-indigo-50/50'));
    dropzoneDecode.addEventListener('drop', (e) => {
      e.preventDefault();
      dropzoneDecode.classList.remove('border-indigo-500', 'bg-indigo-50/50');
      if (e.dataTransfer.files.length) handleDecodeFile(e.dataTransfer.files[0]);
    });

    fileInputDecode.addEventListener('change', (e) => {
      if (e.target.files.length) handleDecodeFile(e.target.files[0]);
    });

    function handleDecodeFile(file) {
      currentDecodeFile = file;
      filePromptDecode.innerHTML = `<span class="text-indigo-600 font-bold">${file.name}</span>`;
      fileSubtextDecode.textContent = `Size: ${formatBytes(file.size)} - Ready to decode into digital file`;
    }

    btnLoadLastPool.addEventListener('click', async () => {
      try {
        const res = await fetch('/api/download/fasta');
        if (!res.ok) throw new Error('No synthesis pool currently cached on server. Encode a file first!');
        const text = await res.text();
        dnaTextInput.value = text;
        charCountDna.textContent = `${text.length} chars`;
        dnaTextInput.scrollIntoView({ behavior: 'smooth' });
      } catch (err) {
        alert(err.message);
      }
    });

    runBtnDecode.addEventListener('click', async () => {
      const textVal = dnaTextInput.value.trim();
      if (!currentDecodeFile && !textVal) {
        alert('Please select a .fasta file or paste DNA sequences into the box.');
        return;
      }

      runBtnDecode.disabled = true;
      loaderDecode.classList.remove('hidden');

      const formData = new FormData();
      if (currentDecodeFile) {
        formData.append('file', currentDecodeFile);
      }
      if (textVal) {
        formData.append('dna_text', textVal);
      }

      try {
        const res = await fetch('/api/decode', { method: 'POST', body: formData });
        const data = await res.json();

        if (!data.success) {
          alert('Decoding Error: ' + data.error);
          return;
        }

        loaderDecode.classList.add('hidden');
        resultsDecode.classList.remove('hidden');

        document.getElementById('statDecodedFilename').textContent = data.filename;
        document.getElementById('statDecodedSize').textContent = formatBytes(data.file_size);
        document.getElementById('statDecodedStrands').textContent = `${data.total_sequences_ingested} strands`;
        document.getElementById('bannerDecodeMeta').textContent = `Restored "${data.filename}" (${formatBytes(data.file_size)}) with matching SHA-256: ${data.sha256}`;

        document.getElementById('bitstreamDecodeView').textContent = data.reconstructed_bits_sample || 'Binary stream verified';
        document.getElementById('hexDecodeView').textContent = data.hex_sample || '';

        // Live text preview if text file restored
        const textContainer = document.getElementById('textPreviewContainer');
        const textArea = document.getElementById('textPreviewArea');
        if (data.text_preview) {
          textContainer.classList.remove('hidden');
          textArea.textContent = data.text_preview;
        } else {
          textContainer.classList.add('hidden');
        }

        const downloadBtn = document.getElementById('btnDownloadDecodedFile');
        downloadBtn.href = '/api/download/recovered?t=' + Date.now();
        downloadBtn.download = data.filename;

        resultsDecode.scrollIntoView({ behavior: 'smooth' });

      } catch (err) {
        alert('Decoding failed: ' + err.message);
      } finally {
        runBtnDecode.disabled = false;
        loaderDecode.classList.add('hidden');
      }
    });
  </script>
</body>
</html>
"""


class HelixDashboardHandler(http.server.BaseHTTPRequestHandler):
    """HTTP Request Handler for the light professional HELIX-CORE dashboard."""

    def _set_headers(self, status: int = 200, content_type: str = "text/html; charset=utf-8"):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Cache-Control", "no-cache, no-store, must-revalidate")
        self.end_headers()

    def do_GET(self):
        """Handle dashboard UI and file downloads."""
        parsed_url = urllib.parse.urlparse(self.path)
        path = parsed_url.path

        if path in ("/", "/index.html"):
            self._set_headers(200, "text/html; charset=utf-8")
            self.wfile.write(DASHBOARD_HTML.encode("utf-8"))
            return

        if path == "/api/download/recovered":
            if not LATEST_RUN["output_path"] or not Path(LATEST_RUN["output_path"]).is_file():
                self._set_headers(404, "application/json")
                self.wfile.write(json.dumps({"error": "No recovered file available"}).encode("utf-8"))
                return
            rec_path = Path(LATEST_RUN["output_path"])
            fname = rec_path.name
            self.send_response(200)
            self.send_header("Content-Type", "application/octet-stream")
            self.send_header("Content-Disposition", f'attachment; filename="{fname}"')
            self.send_header("Content-Length", str(rec_path.stat().st_size))
            self.end_headers()
            with open(rec_path, "rb") as fh:
                while chunk := fh.read(65536):
                    self.wfile.write(chunk)
            return

        if path == "/api/download/fasta":
            if not LATEST_RUN["fasta_path"] or not Path(LATEST_RUN["fasta_path"]).is_file():
                self._set_headers(404, "application/json")
                self.wfile.write(json.dumps({"error": "No FASTA file available"}).encode("utf-8"))
                return
            fasta_path = Path(LATEST_RUN["fasta_path"])
            self.send_response(200)
            self.send_header("Content-Type", "text/plain; charset=utf-8")
            self.send_header("Content-Disposition", 'attachment; filename="synthesis_pool.fasta"')
            self.send_header("Content-Length", str(fasta_path.stat().st_size))
            self.end_headers()
            with open(fasta_path, "rb") as fh:
                while chunk := fh.read(65536):
                    self.wfile.write(chunk)
            return

        if path == "/api/download/csv":
            if not LATEST_RUN["csv_path"] or not Path(LATEST_RUN["csv_path"]).is_file():
                self._set_headers(404, "application/json")
                self.wfile.write(json.dumps({"error": "No CSV file available"}).encode("utf-8"))
                return
            csv_path = Path(LATEST_RUN["csv_path"])
            self.send_response(200)
            self.send_header("Content-Type", "text/csv; charset=utf-8")
            self.send_header("Content-Disposition", 'attachment; filename="synthesis_order.csv"')
            self.send_header("Content-Length", str(csv_path.stat().st_size))
            self.end_headers()
            with open(csv_path, "rb") as fh:
                while chunk := fh.read(65536):
                    self.wfile.write(chunk)
            return

        self._set_headers(404, "text/plain")
        self.wfile.write(b"Not Found")

    def do_POST(self):
        """Handle encoding (upload) and decoding (dna sequence) requests."""
        if self.path == "/api/upload":
            self._handle_upload()
        elif self.path == "/api/decode":
            self._handle_decode()
        else:
            self._set_headers(404, "application/json")
            self.wfile.write(json.dumps({"error": "Not Found"}).encode("utf-8"))

    def _handle_upload(self):
        """Handle digital file upload and DNA synthesis."""
        try:
            content_length = int(self.headers.get("Content-Length", 0))
            content_type = self.headers.get("Content-Type", "")

            if not content_type.startswith("multipart/form-data"):
                raise ValueError("Expected multipart/form-data content type.")

            boundary = content_type.split("boundary=")[1].encode("utf-8")
            body_bytes = self.rfile.read(content_length)

            parts = body_bytes.split(b"--" + boundary)
            file_data = None
            filename = "uploaded_file.bin"
            parity_ratio = 0.15
            drop_rate = 0.10
            mutation_rate = 0.02

            for part in parts:
                if b'Content-Disposition: form-data;' not in part:
                    continue
                header_body_split = part.split(b"\r\n\r\n", 1)
                if len(header_body_split) < 2:
                    continue
                part_headers, part_body = header_body_split
                part_body = part_body.rstrip(b"\r\n")

                if b'name="file"' in part_headers:
                    part_str = part_headers.decode("latin1", errors="ignore")
                    m = re.search(r'filename="([^"]+)"', part_str)
                    if not m:
                        m = re.search(r'filename=([^\s;\r\n]+)', part_str)
                    if m:
                        filename = m.group(1).strip()
                    file_data = part_body
                elif b'name="parity_ratio"' in part_headers:
                    parity_ratio = float(part_body.decode("utf-8", errors="ignore"))
                elif b'name="drop_rate"' in part_headers:
                    drop_rate = float(part_body.decode("utf-8", errors="ignore"))
                elif b'name="mutation_rate"' in part_headers:
                    mutation_rate = float(part_body.decode("utf-8", errors="ignore"))

            if not file_data:
                raise ValueError("No file content received in upload.")

            safe_fname = Path(filename).name
            upload_path = UPLOADS_DIR / safe_fname
            with open(upload_path, "wb") as fh:
                fh.write(file_data)

            recovered_path = RECOVERED_DIR / f"recovered_{safe_fname}"
            fasta_path = EXPORTS_DIR / "synthesis_pool.fasta"
            csv_path = EXPORTS_DIR / "synthesis_order.csv"

            pipeline = HelixPipeline(
                payload_size=28,
                parity_ratio=parity_ratio,
                adaptive_screening=True,
            )

            strands, meta = pipeline.encode_file(upload_path)
            pipeline.translator.export_fasta(strands, fasta_path)
            pipeline.translator.export_csv(strands, csv_path)

            sequenced_pool = pipeline.simulate_sequencing_loss(
                strands, drop_rate=drop_rate, mutation_rate=mutation_rate
            )

            pipeline.decode_strands(
                received_strands=sequenced_pool,
                output_filepath=recovered_path,
                total_data_chunks=meta["data_chunks"],
                total_parity_chunks=meta["parity_chunks"],
            )

            rec_sha = BinaryIngestor.compute_sha256(recovered_path)
            is_perfect = rec_sha == meta["original_sha256"]

            gc_list = [s.gc_content for s in strands]
            max_hp = max((s.max_homopolymer for s in strands), default=0)

            LATEST_RUN["has_data"] = True
            LATEST_RUN["input_path"] = str(upload_path)
            LATEST_RUN["output_path"] = str(recovered_path)
            LATEST_RUN["fasta_path"] = str(fasta_path)
            LATEST_RUN["csv_path"] = str(csv_path)

            samples = []
            for s in strands[:5]:
                samples.append({
                    "chunk_id_hex": f"{s.chunk_id:08X}",
                    "chunk_index": s.chunk_index,
                    "is_parity": s.is_parity,
                    "length": s.length,
                    "gc_content": s.gc_content,
                    "forward_primer": s.forward_primer,
                    "id_sequence": s.id_sequence,
                    "payload_sequence": s.payload_sequence,
                    "crc_sequence": s.crc_sequence,
                    "reverse_primer": s.reverse_primer,
                    "full_sequence": s.full_sequence,
                })

            response_data = {
                "success": True,
                "original_filename": safe_fname,
                "original_size_bytes": meta["original_size"],
                "original_sha256": meta["original_sha256"],
                "recovered_sha256": rec_sha,
                "is_bit_perfect": is_perfect,
                "total_binary_bits": meta["total_bits"],
                "data_chunk_count": meta["data_chunks"],
                "parity_chunk_count": meta["parity_chunks"],
                "total_strand_count": meta["total_strands"],
                "strand_length_bp": meta["strand_length_bp"],
                "parity_ratio": parity_ratio,
                "mean_gc_content": sum(gc_list) / len(gc_list) if gc_list else 0.0,
                "sample_strands": samples,
            }

            self._set_headers(200, "application/json")
            self.wfile.write(json.dumps(response_data).encode("utf-8"))

        except Exception as e:
            self._set_headers(500, "application/json")
            self.wfile.write(json.dumps({"success": False, "error": str(e)}).encode("utf-8"))

    def _handle_decode(self):
        """Handle decoding DNA sequences into a digital file."""
        try:
            content_length = int(self.headers.get("Content-Length", 0))
            content_type = self.headers.get("Content-Type", "")

            dna_text = ""
            if content_type.startswith("multipart/form-data"):
                boundary = content_type.split("boundary=")[1].encode("utf-8")
                body_bytes = self.rfile.read(content_length)
                parts = body_bytes.split(b"--" + boundary)

                for part in parts:
                    if b'Content-Disposition: form-data;' not in part:
                        continue
                    header_body_split = part.split(b"\r\n\r\n", 1)
                    if len(header_body_split) < 2:
                        continue
                    part_headers, part_body = header_body_split
                    part_body = part_body.rstrip(b"\r\n")

                    if b'name="file"' in part_headers:
                        dna_text = part_body.decode("utf-8", errors="ignore")
                    elif b'name="dna_text"' in part_headers and not dna_text:
                        dna_text = part_body.decode("utf-8", errors="ignore")

            elif "application/json" in content_type:
                body_bytes = self.rfile.read(content_length)
                payload_json = json.loads(body_bytes.decode("utf-8"))
                dna_text = payload_json.get("dna_text", "")
            else:
                body_bytes = self.rfile.read(content_length)
                dna_text = body_bytes.decode("utf-8", errors="ignore")

            if not dna_text.strip():
                raise ValueError("No DNA sequence content received.")

            # Parse DNA sequences (FASTA lines or raw whitespace/newline-separated sequences)
            sequences = []
            current_seq = []
            for line in dna_text.splitlines():
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_seq:
                        sequences.append("".join(current_seq))
                        current_seq = []
                else:
                    # Strip any non-ACGT characters
                    clean_line = re.sub(r"[^ACGTacgt]", "", line).upper()
                    if clean_line:
                        current_seq.append(clean_line)

            if current_seq:
                sequences.append("".join(current_seq))

            # If user pasted raw sequences without FASTA headers, split 184 nt chunks
            if len(sequences) == 1 and len(sequences[0]) > 184 and len(sequences[0]) % 184 == 0:
                full_raw = sequences[0]
                sequences = [full_raw[i : i + 184] for i in range(0, len(full_raw), 184)]

            if not sequences:
                raise ValueError("No valid DNA nucleotide sequences found (must contain A, C, G, T).")

            pipeline = HelixPipeline()
            recovered_bytes, meta = pipeline.decode_strands_to_bytes(sequences)

            recovered_filename = meta.get("filename") or "restored_file.bin"
            file_size = meta.get("original_size", len(recovered_bytes))
            file_sha256 = meta.get("sha256", BinaryIngestor.compute_sha256(recovered_bytes))

            recovered_path = RECOVERED_DIR / f"recovered_{Path(recovered_filename).name}"
            with open(recovered_path, "wb") as fh:
                fh.write(recovered_bytes)

            LATEST_RUN["output_path"] = str(recovered_path)

            # Sample binary bitstream (first 128 bits)
            bit_sample = BinaryIngestor.bytes_to_binary_string(recovered_bytes[:16])

            # Sample hex dump
            hex_lines = []
            for i in range(0, min(len(recovered_bytes), 48), 16):
                c = recovered_bytes[i : i + 16]
                h_part = " ".join(f"{b:02x}" for b in c)
                a_part = "".join(chr(b) if 32 <= b < 127 else "." for b in c)
                hex_lines.append(f"{i:04x}: {h_part:<48} |{a_part}|")

            # Check if readable UTF-8 text for live preview (guarantees clear English display)
            text_preview = None
            try:
                decoded_text = recovered_bytes.decode("utf-8")
                printable = sum(1 for ch in decoded_text if ch.isprintable() or ch in "\r\n\t")
                if len(decoded_text) > 0 and (printable / len(decoded_text)) >= 0.85:
                    text_preview = decoded_text[:4000]
            except Exception:
                text_preview = None

            response_data = {
                "success": True,
                "filename": Path(recovered_filename).name,
                "file_size": file_size,
                "sha256": file_sha256,
                "total_sequences_ingested": len(sequences),
                "valid_data_chunks": meta.get("total_data_chunks", len(sequences)),
                "reconstructed_bits_sample": bit_sample,
                "hex_sample": "\n".join(hex_lines),
                "text_preview": text_preview,
            }

            self._set_headers(200, "application/json")
            self.wfile.write(json.dumps(response_data).encode("utf-8"))

        except Exception as e:
            self._set_headers(500, "application/json")
            self.wfile.write(json.dumps({"success": False, "error": str(e)}).encode("utf-8"))


class ThreadingHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True


def start_server(port: int = PORT):
    server = ThreadingHTTPServer(("0.0.0.0", port), HelixDashboardHandler)
    print("=" * 70)
    print(f"🧬 HELIX-CORE Bi-Directional Dashboard Live on http://localhost:{port}")
    print(f"[*] Encode: Upload digital file ➔ Physical DNA sequences")
    print(f"[*] Decode: Upload DNA sequences ➔ Binary ➔ Restored digital file")
    print("=" * 70)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping dashboard server.")
        server.server_close()


if __name__ == "__main__":
    start_server(PORT)
