"""
Top-level package for CaverDock-style lower-bound reproduction via Vina.

Exposes key abstractions:
- ChargeProvider: strategy interface for per-disc charge preparation
- run_pipeline: orchestrate discs extraction, charge prep and Vina docking
"""

from .pipeline import run_pipeline
from .charge_providers import ChargeProvider, MGLToolsGasteigerChargeProvider

__all__ = [
    "run_pipeline",
    "ChargeProvider",
    "MGLToolsGasteigerChargeProvider",
]
