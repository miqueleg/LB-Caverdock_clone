from __future__ import annotations

import csv
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List


def _run(cmd: list[str], cwd: Path | None = None) -> str:
    proc = subprocess.run(
        cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stdout}")
    return proc.stdout


def generate_disc_boxes(discs: List[dict], box_margin: float = 3.0) -> List[dict]:
    # passthrough if already boxes
    from .discs import generate_disc_boxes as _g, Disc

    if discs and isinstance(discs[0], dict) and "size_x" in discs[0]:
        return discs  # already in box dict format
    return _g(discs, box_margin=box_margin)


def run_vina_for_discs(
    boxes_csv: str | Path,
    receptor_pdbqt: str | Path,
    ligand_pdbqt: str | Path,
    vina_bin: str = "vina",
    exhaustiveness: int = 16,
    num_poses: int = 20,
    threads: int = os.cpu_count() or 4,
    outdir: str | Path = "vina_runs",
) -> Dict[int, float]:
    """
    Run vina for each disc (box) listed in boxes_csv. Results are saved under
    {outdir}/disc_{index}/ (out.pdbqt, log.txt). Returns a mapping disc_index->best_energy.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(boxes_csv) as fh:
        boxes = list(csv.DictReader(fh))

    def _job(box: dict) -> tuple[int, float]:
        idx = int(box["index"])
        ddir = outdir / f"disc_{idx:04d}"
        ddir.mkdir(exist_ok=True)
        out_pdbqt = ddir / "out.pdbqt"
        log = ddir / "log.txt"
        cmd = [
            vina_bin,
            "--receptor",
            str(receptor_pdbqt),
            "--ligand",
            str(ligand_pdbqt),
            "--center_x",
            str(box["center_x"]),
            "--center_y",
            str(box["center_y"]),
            "--center_z",
            str(box["center_z"]),
            "--size_x",
            str(box["size_x"]),
            "--size_y",
            str(box["size_y"]),
            "--size_z",
            str(box["size_z"]),
            "--exhaustiveness",
            str(exhaustiveness),
            "--num_modes",
            str(num_poses),
            "--out",
            str(out_pdbqt),
            "--log",
            str(log),
        ]
        _run(cmd)
        best = parse_best_energy_from_log(log)
        # Save metadata for later filtering
        with (ddir / "box.json").open("w") as fh:
            json.dump(box, fh)
        return idx, best

    best_map: Dict[int, float] = {}
    with ThreadPoolExecutor(max_workers=max(1, int(threads))) as ex:
        futs = [ex.submit(_job, box) for box in boxes]
        for fut in as_completed(futs):
            idx, val = fut.result()
            best_map[idx] = val
    return best_map


def parse_best_energy_from_log(log_path: str | Path) -> float:
    """
    Parse Vina log for the best (lowest) energy. Fallback to +inf if not found.
    """
    best = float("inf")
    with open(log_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("   1 "):
                # Expected Vina table line: rank, energy, rmsd1, rmsd2
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        best = float(parts[1])
                        break
                    except ValueError:
                        continue
    return best


def parse_energy_profile(vina_outdir: str | Path) -> Dict[int, float]:
    """Collect best energies from disc_* folders."""
    outdir = Path(vina_outdir)
    res: Dict[int, float] = {}
    for d in sorted(outdir.glob("disc_*")):
        try:
            idx = int(d.name.split("_")[-1])
        except Exception:
            continue
        log = d / "log.txt"
        if log.exists():
            res[idx] = parse_best_energy_from_log(log)
    return res

