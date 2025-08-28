from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass
class VinaBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float


def box_from_sphere(x: float, y: float, z: float, radius: float, margin: float = 2.0) -> VinaBox:
    size = 2.0 * radius + margin
    return VinaBox(center_x=x, center_y=y, center_z=z, size_x=size, size_y=size, size_z=size)


def run_vina(
    vina_bin: str,
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    box: VinaBox,
    out_pdbqt: str,
    log_path: Optional[str] = None,
    exhaustiveness: int = 8,
    num_modes: int = 1,
    seed: Optional[int] = None,
    cpu: Optional[int] = 1,
    mode: str = "dock",  # "dock", "local_only", "score_only"
    scoring: str = "vina",  # "vina", "ad4", "vinardo"
) -> Tuple[Optional[float], str]:
    cmd = [
        vina_bin,
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        str(box.center_x),
        "--center_y",
        str(box.center_y),
        "--center_z",
        str(box.center_z),
        "--size_x",
        str(box.size_x),
        "--size_y",
        str(box.size_y),
        "--size_z",
        str(box.size_z),
    ]
    if mode == "dock":
        cmd += [
            "--exhaustiveness",
            str(exhaustiveness),
            "--num_modes",
            str(num_modes),
        ]
    elif mode == "local_only":
        cmd += ["--local_only"]
    elif mode == "score_only":
        cmd += ["--score_only"]
    else:
        cmd += [
            "--exhaustiveness",
            str(exhaustiveness),
            "--num_modes",
            str(num_modes),
        ]
    cmd += ["--out", str(out_pdbqt)]
    if scoring and scoring != "vina":
        cmd += ["--scoring", scoring]
    if seed is not None and mode == "dock":
        cmd += ["--seed", str(seed)]
    if cpu is not None:
        cmd += ["--cpu", str(cpu)]

    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    stdout = proc.stdout
    if log_path:
        with open(log_path, "w") as f:
            f.write(stdout)

    score = None
    m = re.search(r"^\s*1\s+(-?\d+\.\d+)", stdout, flags=re.M)
    if m:
        score = float(m.group(1))
    # Fallback 1: parse from score_only stdout
    if score is None:
        m2 = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.\d+)", stdout)
        if m2:
            score = float(m2.group(1))
    # Fallback 2: parse from output PDBQT REMARK if needed
    if score is None or abs(score) > 50:
        try:
            with open(out_pdbqt, "r") as fp:
                for line in fp:
                    mm = re.search(r"^REMARK\s+VINA RESULT:\s+(-?\d+\.\d+)", line)
                    if mm:
                        score = float(mm.group(1))
                        break
        except Exception:
            pass
    return score, stdout
