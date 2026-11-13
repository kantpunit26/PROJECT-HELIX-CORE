"""
HELIX-CORE: Production-Grade Bi-Directional DNA Data Storage Pipeline.

Developed by Master Bio-Computational Engineering.
"""

from .config import (
    BITS_TO_BASE,
    BASE_TO_BITS,
    DEFAULT_FORWARD_PRIMER,
    DEFAULT_PARITY_RATIO,
    DEFAULT_PAYLOAD_SIZE,
    DEFAULT_REVERSE_PRIMER,
    MAX_HOMOPOLYMER_ALLOWED,
    MAX_GC_CONTENT,
    MAX_STRAND_LENGTH,
    MIN_GC_CONTENT,
    MIN_STRAND_LENGTH,
)
from .fragmentation import FragmentationEngine, ReedSolomonCodec, StrandIntegrity
from .ingestion import BinaryIngestor
from .obfuscation import BiologicalObfuscator, ConstraintValidator
from .pipeline import HelixPipeline
from .translation import DnaTranslator
from .types import (
    BioConstraints,
    Chunk,
    DnaStrand,
    MaskedChunk,
    PipelineStats,
    ValidationResult,
)

__version__ = "1.0.0"
__all__ = [
    "BinaryIngestor",
    "FragmentationEngine",
    "ReedSolomonCodec",
    "StrandIntegrity",
    "BiologicalObfuscator",
    "ConstraintValidator",
    "DnaTranslator",
    "HelixPipeline",
    "BioConstraints",
    "Chunk",
    "MaskedChunk",
    "DnaStrand",
    "PipelineStats",
    "ValidationResult",
    "DEFAULT_FORWARD_PRIMER",
    "DEFAULT_REVERSE_PRIMER",
    "DEFAULT_PAYLOAD_SIZE",
    "DEFAULT_PARITY_RATIO",
    "MAX_HOMOPOLYMER_ALLOWED",
    "MIN_GC_CONTENT",
    "MAX_GC_CONTENT",
    "MIN_STRAND_LENGTH",
    "MAX_STRAND_LENGTH",
]

