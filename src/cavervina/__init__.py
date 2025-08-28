"""
cavervina: Recreate CaverDock upper-bound energies using CAVER3 + AutoDock Vina.

This package orchestrates:
- Tunnel/disc parsing (from CAVER3 outputs)
- Receptor/ligand preparation to PDBQT with Gasteiger charges
- Per-disc Vina docking (parallel)
- Pose filtering by proximity to disc
- Energy profile extraction and comparison to CaverDock outputs

CLI entry point: python -m cavervina.cli
"""

__all__ = [
    "discs",
    "prepare",
    "vina",
    "caver",
    "profile",
]

