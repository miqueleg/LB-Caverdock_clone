from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, List, Tuple


def _unit(v: Tuple[float, float, float]) -> Tuple[float, float, float]:
    x, y, z = v
    n = math.sqrt(x * x + y * y + z * z) or 1.0
    return (x / n, y / n, z / n)


def _dot(a, b) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _parse_boxes_csv(boxes_csv: str | Path) -> List[dict]:
    with open(boxes_csv) as fh:
        rows = list(csv.DictReader(fh))
        # cast
        for r in rows:
            for k in [
                "index",
                "center_x",
                "center_y",
                "center_z",
                "size_x",
                "size_y",
                "size_z",
                "radius",
                "nx",
                "ny",
                "nz",
            ]:
                if k in r and r[k] != "":
                    r[k] = float(r[k]) if k != "index" else int(float(r[k]))
        return rows


def _iter_pdbqt_models(pdbqt_path: Path):
    """
    Yield tuples (energy, coords) for each MODEL in a Vina PDBQT. coords is a list of (x,y,z).
    Energy is parsed from "REMARK VINA RESULT:" lines preceding coordinates for the model.
    """
    energy = None
    coords: List[Tuple[float, float, float]] = []
    in_model = False
    with pdbqt_path.open() as fh:
        for line in fh:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                try:
                    energy = float(parts[3])
                except Exception:
                    energy = None
            elif line.startswith("MODEL"):
                in_model = True
                coords = []
            elif line.startswith("ENDMDL"):
                if in_model:
                    yield energy, coords
                in_model = False
                energy = None
                coords = []
            elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except Exception:
                    continue


def _geom_center(coords: List[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    if not coords:
        return (0.0, 0.0, 0.0)
    sx = sum(c[0] for c in coords)
    sy = sum(c[1] for c in coords)
    sz = sum(c[2] for c in coords)
    n = len(coords)
    return (sx / n, sy / n, sz / n)


def filter_poses_near_disc(
    disc_row: dict,
    poses_pdbqt: Path,
    disc_half_thickness: float = 2.0,
) -> List[Tuple[float, float]]:
    """
    Return list of (energy, radial_distance) for poses whose centroid lies within:
    - plane distance <= disc_half_thickness from disc center, and
    - radial distance <= disc radius.
    """
    center = (float(disc_row["center_x"]), float(disc_row["center_y"]), float(disc_row["center_z"]))
    normal = _unit((float(disc_row["nx"]), float(disc_row["ny"]), float(disc_row["nz"])) )
    radius = float(disc_row["radius"])
    accepted: List[Tuple[float, float]] = []
    for energy, coords in _iter_pdbqt_models(poses_pdbqt):
        if energy is None:
            continue
        c = _geom_center(coords)
        v = (c[0] - center[0], c[1] - center[1], c[2] - center[2])
        # Signed distance to plane along normal
        plane_dist = abs(_dot(v, normal))
        # Radial component after subtracting the normal projection
        proj = _dot(v, normal)
        rad_vec = (v[0] - proj * normal[0], v[1] - proj * normal[1], v[2] - proj * normal[2])
        radial = math.sqrt(rad_vec[0] ** 2 + rad_vec[1] ** 2 + rad_vec[2] ** 2)
        if plane_dist <= disc_half_thickness and radial <= radius:
            accepted.append((energy, radial))
    return accepted


def export_profile_csv(
    boxes_csv: str | Path,
    vina_outdir: str | Path,
    disc_half_thickness: float,
    out_csv: str | Path,
) -> None:
    """For each disc_X folder, filter poses near disc and export the best energy per disc."""
    boxes = _parse_boxes_csv(boxes_csv)
    by_idx = {int(r["index"]): r for r in boxes}
    outdir = Path(vina_outdir)
    rows: List[dict] = []
    for d in sorted(outdir.glob("disc_*")):
        try:
            idx = int(d.name.split("_")[-1])
        except Exception:
            continue
        if idx not in by_idx:
            continue
        poses = d / "out.pdbqt"
        if not poses.exists():
            continue
        accepted = filter_poses_near_disc(by_idx[idx], poses, disc_half_thickness)
        if accepted:
            best_e = min(e for e, _ in accepted)
        else:
            # Fallback to best overall energy from log
            from .vina import parse_best_energy_from_log

            best_e = parse_best_energy_from_log(d / "log.txt")
        rows.append({"index": idx, "energy": best_e})

    with Path(out_csv).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["index", "energy"])
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["index"]):
            w.writerow(r)


def compare_profiles(ours_csv: str | Path, cdock_csv: str | Path) -> Dict[str, float]:
    def _read(path: str | Path) -> Dict[int, float]:
        with open(path) as fh:
            rd = csv.DictReader(fh)
            out: Dict[int, float] = {}
            for r in rd:
                idx = int(r.get("index") or r.get("disc") or r.get("step"))
                e = float(r.get("energy") or r.get("upper_bound") or r.get("Eub"))
                out[idx] = e
            return out

    a = _read(ours_csv)
    b = _read(cdock_csv)
    keys = sorted(set(a) & set(b))
    if not keys:
        return {"n": 0}
    diffs = [a[k] - b[k] for k in keys]
    mae = sum(abs(d) for d in diffs) / len(diffs)
    rmse = math.sqrt(sum(d * d for d in diffs) / len(diffs))
    return {"n": len(keys), "mae": mae, "rmse": rmse}

