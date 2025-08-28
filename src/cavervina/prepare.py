from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Optional


def _run(cmd: list[str]) -> None:
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stdout}")


def prepare_receptor_pdbqt(
    receptor_input: str | Path,
    out_pdbqt: str | Path,
    adt_bin: Optional[Path] = None,
    obabel: Optional[Path] = None,
) -> None:
    """
    Prepare receptor PDBQT with Gasteiger charges.
    Preference order: AutoDockTools (prepare_receptor4.py) > Open Babel (obabel -xr -xp)
    """
    out_pdbqt = Path(out_pdbqt)
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    if adt_bin is not None:
        prep = Path(adt_bin) / "prepare_receptor4.py"
        if not prep.exists():
            raise FileNotFoundError(f"Missing ADT script: {prep}")
        _run([
            "python",
            str(prep),
            "-r",
            str(receptor_input),
            "-o",
            str(out_pdbqt),
            "-A",
            "checkhydrogens",  # ADT will add Gasteiger by default for receptor
        ])
        return

    if obabel is not None:
        # Use Open Babel as a fallback; try to assign Gasteiger charges
        _run([
            str(obabel),
            "-i",
            Path(receptor_input).suffix.lstrip("."),
            str(receptor_input),
            "-xr",
            "-xp",
            "-o",
            "pdbqt",
            "-O",
            str(out_pdbqt),
        ])
        return

    raise RuntimeError("Provide either --adt-bin to MGLTools scripts or --obabel path")


def prepare_ligand_pdbqt(
    ligand_input: str | Path,
    out_pdbqt: str | Path,
    adt_bin: Optional[Path] = None,
    obabel: Optional[Path] = None,
) -> None:
    """
    Prepare ligand PDBQT with Gasteiger charges and rotatable bonds.
    Preference order: AutoDockTools (prepare_ligand4.py) > Open Babel.
    """
    out_pdbqt = Path(out_pdbqt)
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    if adt_bin is not None:
        prep = Path(adt_bin) / "prepare_ligand4.py"
        if not prep.exists():
            raise FileNotFoundError(f"Missing ADT script: {prep}")
        _run([
            "python",
            str(prep),
            "-l",
            str(ligand_input),
            "-o",
            str(out_pdbqt),
            "-A",
            "gasteiger",
        ])
        return

    if obabel is not None:
        _run([
            str(obabel),
            "-i",
            Path(ligand_input).suffix.lstrip("."),
            str(ligand_input),
            "--partialcharge",
            "gasteiger",
            "-o",
            "pdbqt",
            "-O",
            str(out_pdbqt),
        ])
        return

    raise RuntimeError("Provide either --adt-bin to MGLTools scripts or --obabel path")

