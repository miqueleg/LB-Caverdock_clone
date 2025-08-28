from __future__ import annotations

from pathlib import Path
from typing import List, Tuple


def read_pdbqt_coords(path: str | Path) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    lines = []
    coords = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    x = y = z = 0.0
                coords.append((x, y, z))
            lines.append(line)
    return lines, coords


def write_pdbqt_coords(lines: List[str], coords: List[Tuple[float, float, float]], out_path: str | Path) -> None:
    out_lines = []
    ci = 0
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x, y, z = coords[ci]
            new = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            out_lines.append(new)
            ci += 1
        else:
            out_lines.append(line)
    with open(out_path, "w") as f:
        f.writelines(out_lines)


def compute_centroid(coords: List[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    if not coords:
        return (0.0, 0.0, 0.0)
    sx = sum(c[0] for c in coords)
    sy = sum(c[1] for c in coords)
    sz = sum(c[2] for c in coords)
    n = float(len(coords))
    return (sx / n, sy / n, sz / n)


def translate_coords(coords: List[Tuple[float, float, float]], dx: float, dy: float, dz: float) -> List[Tuple[float, float, float]]:
    return [(x + dx, y + dy, z + dz) for (x, y, z) in coords]


def recenter_ligand_to(coords: List[Tuple[float, float, float]], target: Tuple[float, float, float]) -> List[Tuple[float, float, float]]:
    cx, cy, cz = compute_centroid(coords)
    dx = target[0] - cx
    dy = target[1] - cy
    dz = target[2] - cz
    return translate_coords(coords, dx, dy, dz)

