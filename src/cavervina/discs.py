from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple


@dataclass
class Disc:
    index: int
    x: float
    y: float
    z: float
    radius: float
    # Optional: direction vector (unit tangent). If not provided, computed later.
    nx: float | None = None
    ny: float | None = None
    nz: float | None = None


def load_discs(path: str | Path) -> List[Disc]:
    """
    Load discs from common CAVER/CaverDock formats.
    Supported rows (whitespace or CSV):
      x y z r [nx ny nz]
    Headers are ignored.
    """
    p = Path(path)
    discs: List[Disc] = []
    with p.open() as fh:
        reader = csv.reader(fh)
        for raw in reader:
            if not raw:
                continue
            # Try split on whitespace if CSV parser returned a single token line
            if len(raw) == 1:
                raw = raw[0].strip().split()
            # Skip header-like lines
            if raw[0].startswith("#") or raw[0].lower().startswith("x"):
                continue
            try:
                vals = list(map(float, raw[:7]))
            except ValueError:
                # Non-numeric line
                continue
            x, y, z, r = vals[:4]
            nx = ny = nz = None
            if len(vals) >= 7:
                nx, ny, nz = vals[4:7]
            discs.append(Disc(len(discs), x, y, z, r, nx, ny, nz))

    if len(discs) < 2:
        raise ValueError("Expected at least 2 discs to define a tunnel")

    # Fill missing tangents as finite differences along the path
    for i, d in enumerate(discs):
        if d.nx is None or d.ny is None or d.nz is None:
            if i == 0:
                dx, dy, dz = (
                    discs[i + 1].x - d.x,
                    discs[i + 1].y - d.y,
                    discs[i + 1].z - d.z,
                )
            elif i == len(discs) - 1:
                dx, dy, dz = (
                    d.x - discs[i - 1].x,
                    d.y - discs[i - 1].y,
                    d.z - discs[i - 1].z,
                )
            else:
                dx, dy, dz = (
                    discs[i + 1].x - discs[i - 1].x,
                    discs[i + 1].y - discs[i - 1].y,
                    discs[i + 1].z - discs[i - 1].z,
                )
            norm = math.sqrt(dx * dx + dy * dy + dz * dz) or 1.0
            d.nx, d.ny, d.nz = dx / norm, dy / norm, dz / norm
    return discs


def generate_disc_boxes(discs: List[Disc], box_margin: float = 3.0) -> List[dict]:
    """
    Generate Vina box parameters per disc. Vina uses an axis-aligned grid box, so we
    approximate disc coverage by using a cubic box centered at the disc center with
    size = 2 * (radius + margin).
    Returns dicts: {index, center_x, center_y, center_z, size_x, size_y, size_z, radius}
    """
    boxes: List[dict] = []
    for d in discs:
        size = max(2.0 * (d.radius + box_margin), 8.0)  # enforce a reasonable minimum
        boxes.append(
            {
                "index": d.index,
                "center_x": d.x,
                "center_y": d.y,
                "center_z": d.z,
                "size_x": size,
                "size_y": size,
                "size_z": size,
                "radius": d.radius,
                "nx": d.nx,
                "ny": d.ny,
                "nz": d.nz,
            }
        )
    return boxes


def write_boxes_csv(boxes: Iterable[dict], out_csv: str | Path) -> None:
    fieldnames = [
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
    ]
    p = Path(out_csv)
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for row in boxes:
            w.writerow(row)


