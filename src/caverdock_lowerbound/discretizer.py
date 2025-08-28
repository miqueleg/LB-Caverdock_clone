from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from pathlib import Path
from typing import List, Tuple

import numpy as np


@dataclass
class Disc:
    index: int
    x: float
    y: float
    z: float
    radius: float
    distance: float
    nx: float | None = None
    ny: float | None = None
    nz: float | None = None


def discretize_from_profile(
    distances: np.ndarray,
    coords: Tuple[np.ndarray, np.ndarray, np.ndarray],
    radii: np.ndarray,
    spacing: float = 0.5,
) -> List[Disc]:
    """
    Discretize a tunnel profile (X,Y,Z,R vs cumulative distance) into evenly spaced discs.

    distances: shape (N,), monotonically increasing cumulative distances
    coords: (X, Y, Z) arrays shape (N,)
    radii: shape (N,)
    spacing: target spacing between discs in Angstrom
    """
    s = np.asarray(distances, dtype=float)
    x, y, z = (np.asarray(c, dtype=float) for c in coords)
    r = np.asarray(radii, dtype=float)
    if s.ndim != 1 or x.shape != s.shape or y.shape != s.shape or z.shape != s.shape or r.shape != s.shape:
        raise ValueError("Input arrays must be 1D and of equal length")
    if len(s) < 2:
        raise ValueError("Need at least 2 points to discretize")
    s_total = s[-1]
    m = max(2, int(np.floor(s_total / spacing)) + 1)
    s_grid = np.linspace(0.0, s_total, m)

    xi = np.interp(s_grid, s, x)
    yi = np.interp(s_grid, s, y)
    zi = np.interp(s_grid, s, z)
    ri = np.interp(s_grid, s, r)

    # Tangents by finite differences on interpolated curve
    tx = np.gradient(xi, s_grid)
    ty = np.gradient(yi, s_grid)
    tz = np.gradient(zi, s_grid)
    norms = np.sqrt(tx**2 + ty**2 + tz**2) + 1e-12
    nx, ny, nz = tx / norms, ty / norms, tz / norms

    discs: List[Disc] = []
    for i, (sg, cx, cy, cz, rr, vx, vy, vz) in enumerate(zip(s_grid, xi, yi, zi, ri, nx, ny, nz), start=1):
        discs.append(Disc(index=i, x=cx, y=cy, z=cz, radius=float(rr), distance=float(sg), nx=float(vx), ny=float(vy), nz=float(vz)))
    return discs


def write_dsd(discs: List[Disc], out_path: str | Path) -> None:
    """
    Write discs to CaverDock-like .dsd file: x y z nx ny nz radius
    """
    p = Path(out_path)
    with p.open("w") as f:
        for d in discs:
            vx = d.nx if d.nx is not None else 0.0
            vy = d.ny if d.ny is not None else 0.0
            vz = d.nz if d.nz is not None else 1.0
            f.write(f"{d.x} {d.y} {d.z} {vx} {vy} {vz} {d.radius}\n")


def parse_dsd(path: str | Path) -> List[Disc]:
    discs: List[Disc] = []
    prev = None
    dist = 0.0
    with open(path, "r") as f:
        for i, line in enumerate(f, start=1):
            parts = line.split()
            if len(parts) < 7:
                continue
            x, y, z = map(float, parts[:3])
            vx, vy, vz = map(float, parts[3:6])
            r = float(parts[6])
            if prev is not None:
                dist += sqrt((x - prev[0]) ** 2 + (y - prev[1]) ** 2 + (z - prev[2]) ** 2)
            discs.append(Disc(i, x, y, z, r, dist, vx, vy, vz))
            prev = (x, y, z)
    return discs

