from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class DiscContext:
    index: int
    x: float
    y: float
    z: float
    radius: float
    distance: float


class ChargeProvider:
    """
    Strategy interface to provide receptor/ligand PDBQT per disc context.
    Implementations can recompute charges per sphere (XTB, AMOEBA, E-field, etc.).
    """

    def get_receptor_pdbqt(self, context: DiscContext) -> str:
        raise NotImplementedError

    def get_ligand_pdbqt(self, context: DiscContext) -> str:
        raise NotImplementedError


class MGLToolsGasteigerChargeProvider(ChargeProvider):
    """
    Default provider using MGLTools prepare_{receptor,ligand}4.py with Gasteiger charges.

    It supports optional per-sphere recomputation if requested; otherwise caches
    a single preparation for all discs.
    """

    def __init__(
        self,
        receptor_input: str,
        ligand_input: str,
        out_dir: str,
        prepare_receptor_script: str = "prepare_receptor4.py",
        prepare_ligand_script: str = "prepare_ligand4.py",
        mgltools_pythonsh: Optional[str] = None,
        per_sphere_recompute_receptor: bool = False,
        per_sphere_recompute_ligand: bool = False,
        receptor_pdbqt_basename: str = "receptor.pdbqt",
        ligand_pdbqt_basename: str = "ligand.pdbqt",
        extra_receptor_args: Optional[list[str]] = None,
        extra_ligand_args: Optional[list[str]] = None,
    ) -> None:
        self.receptor_input = str(Path(receptor_input).resolve())
        self.ligand_input = str(Path(ligand_input).resolve())
        self.out_dir = Path(out_dir).resolve()
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.prepare_receptor_script = prepare_receptor_script
        self.prepare_ligand_script = prepare_ligand_script
        self.mgltools_pythonsh = mgltools_pythonsh
        self.per_sphere_recompute_receptor = per_sphere_recompute_receptor
        self.per_sphere_recompute_ligand = per_sphere_recompute_ligand
        self.receptor_pdbqt_basename = receptor_pdbqt_basename
        self.ligand_pdbqt_basename = ligand_pdbqt_basename
        self.extra_receptor_args = extra_receptor_args or []
        self.extra_ligand_args = extra_ligand_args or []

        self._receptor_cached: Optional[str] = None
        self._ligand_cached: Optional[str] = None

    def _prepare_receptor(self, out_path: Path) -> str:
        rec_in = self.receptor_input
        out_path = out_path.resolve()
        rec_in = Path(self.receptor_input)
        rec_dir, rec_name = rec_in.parent, rec_in.name
        out_path = out_path.resolve()
        if self.mgltools_pythonsh:
            cmd = [
                self.mgltools_pythonsh,
                self.prepare_receptor_script,
                "-r",
                rec_name,
                "-o",
                str(out_path),
            ] + self.extra_receptor_args
        else:
            cmd = [
                self.prepare_receptor_script,
                "-r",
                rec_name,
                "-o",
                str(out_path),
            ] + self.extra_receptor_args
        subprocess.run(cmd, check=True, cwd=str(rec_dir))
        return str(out_path)

    def _prepare_ligand(self, out_path: Path) -> str:
        lig_in = Path(self.ligand_input)
        lig_dir, lig_name = lig_in.parent, lig_in.name
        out_path = out_path.resolve()
        if self.mgltools_pythonsh:
            cmd = [
                self.mgltools_pythonsh,
                self.prepare_ligand_script,
                "-l",
                lig_name,
                "-o",
                str(out_path),
            ] + self.extra_ligand_args
        else:
            cmd = [
                self.prepare_ligand_script,
                "-l",
                lig_name,
                "-o",
                str(out_path),
            ] + self.extra_ligand_args
        subprocess.run(cmd, check=True, cwd=str(lig_dir))
        return str(out_path)

    def get_receptor_pdbqt(self, context: DiscContext) -> str:
        if not self.per_sphere_recompute_receptor and self._receptor_cached:
            return self._receptor_cached
        out_name = (
            f"disc_{context.index}_" + self.receptor_pdbqt_basename
            if self.per_sphere_recompute_receptor
            else self.receptor_pdbqt_basename
        )
        out_path = self.out_dir / out_name
        result = self._prepare_receptor(out_path)
        if not self.per_sphere_recompute_receptor:
            self._receptor_cached = result
        return result

    def get_ligand_pdbqt(self, context: DiscContext) -> str:
        if not self.per_sphere_recompute_ligand and self._ligand_cached:
            return self._ligand_cached
        out_name = (
            f"disc_{context.index}_" + self.ligand_pdbqt_basename
            if self.per_sphere_recompute_ligand
            else self.ligand_pdbqt_basename
        )
        out_path = self.out_dir / out_name
        result = self._prepare_ligand(out_path)
        if not self.per_sphere_recompute_ligand:
            self._ligand_cached = result
        return result


class ExternalCommandChargeProvider(ChargeProvider):
    """
    Generic provider that shells out to user-specified commands to produce
    per-sphere receptor/ligand PDBQT files.

    Use format strings with {index}, {x}, {y}, {z}, {radius}, {distance}, {out}
    placeholders to build commands.
    """

    def __init__(
        self,
        out_dir: str,
        receptor_cmd_template: str,
        ligand_cmd_template: str,
    ) -> None:
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.receptor_cmd_template = receptor_cmd_template
        self.ligand_cmd_template = ligand_cmd_template

    def _run_template(self, template: str, context: DiscContext, out_name: str) -> str:
        out_path = self.out_dir / out_name
        mapping = {
            "index": context.index,
            "x": context.x,
            "y": context.y,
            "z": context.z,
            "radius": context.radius,
            "distance": context.distance,
            "out": str(out_path),
        }
        cmd = template.format(**mapping)
        subprocess.run(cmd, shell=True, check=True)
        return str(out_path)

    def get_receptor_pdbqt(self, context: DiscContext) -> str:
        return self._run_template(self.receptor_cmd_template, context, f"disc_{context.index}_receptor.pdbqt")

    def get_ligand_pdbqt(self, context: DiscContext) -> str:
        return self._run_template(self.ligand_cmd_template, context, f"disc_{context.index}_ligand.pdbqt")


class StaticPDBQTProvider(ChargeProvider):
    """
    Returns precomputed receptor/ligand PDBQT paths for all discs.
    Useful when preparation is done externally once.
    """

    def __init__(self, receptor_pdbqt: str, ligand_pdbqt: str) -> None:
        self.receptor_pdbqt = str(Path(receptor_pdbqt).resolve())
        self.ligand_pdbqt = str(Path(ligand_pdbqt).resolve())

    def get_receptor_pdbqt(self, context: DiscContext) -> str:
        return self.receptor_pdbqt

    def get_ligand_pdbqt(self, context: DiscContext) -> str:
        return self.ligand_pdbqt
