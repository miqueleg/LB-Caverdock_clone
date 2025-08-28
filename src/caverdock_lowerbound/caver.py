from __future__ import annotations

import csv
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Tuple


@dataclass
class Disc:
    index: int
    x: float
    y: float
    z: float
    radius: float
    distance: float


def run_caver(
    pdb_path: str | os.PathLike,
    output_dir: str | os.PathLike,
    caver_cmd: Optional[List[str]] = None,
    config_path: Optional[str | os.PathLike] = None,
    cwd: Optional[str | os.PathLike] = None,
) -> None:
    """
    Run CAVER 3 to compute tunnels.

    Parameters
    ----------
    pdb_path : PDB structure path
    output_dir : CAVER output directory
    caver_cmd : command to invoke CAVER (default: ["caver"]). Examples:
        - ["caver"] if available on PATH
        - ["java", "-jar", "/path/to/Caver3.jar"]
    config_path : optional CAVER config file
    cwd : working directory
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    if caver_cmd is None:
        caver_cmd = ["caver"]

    cmd = list(caver_cmd)
    # Heuristic: common CLIs accept "--pdb" and "--out" or read from config.
    # We support both: if a config is supplied, pass it; otherwise pass pdb and out.
    if config_path:
        cmd += ["--config", str(config_path)]
    cmd += ["--pdb", str(pdb_path), "--out", str(output_dir)]

    subprocess.run(cmd, check=True, cwd=cwd)


def find_tunnel_dir(root: str | os.PathLike, tunnel_index: int) -> Optional[Path]:
    root_p = Path(root)
    if not root_p.exists():
        return None
    # Common CAVER 3 patterns
    candidates = [
        root_p / f"tunnel_{tunnel_index}",
        root_p / f"Tunnel_{tunnel_index}",
        root_p / f"tunnel{tunnel_index}",
    ]
    for c in candidates:
        if c.exists():
            return c
    # Try auto-detect a single tunnel directory
    dirs = [d for d in root_p.iterdir() if d.is_dir() and re.search(r"tunnel", d.name, re.I)]
    return dirs[0] if len(dirs) == 1 else None


def parse_profile_csv(path: Path) -> List[Disc]:
    discs: List[Disc] = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        # Guess column names
        cols = {k.lower().strip(): k for k in reader.fieldnames or []}
        def first(*names: str) -> Optional[str]:
            for n in names:
                if n in cols:
                    return cols[n]
            return None

        xk = first("x", "coordx", "center_x")
        yk = first("y", "coordy", "center_y")
        zk = first("z", "coordz", "center_z")
        rk = first("r", "radius")
        dk = first("distance", "dist", "s")
        ik = first("index", "i")
        idx = 0
        for row in reader:
            idx += 1
            discs.append(
                Disc(
                    index=int(row.get(ik) or idx),
                    x=float(row[xk]),
                    y=float(row[yk]),
                    z=float(row[zk]),
                    radius=float(row[rk]),
                    distance=float(row.get(dk) or (idx - 1)),
                )
            )
    return discs


def parse_centerline_pdb(path: Path) -> List[Disc]:
    """
    Parse a CAVER centerline/spheres PDB. Heuristics:
    - ATOM/HETATM coords are disc centers
    - Radius may be stored in B-factor or occupancy; try both and fallback to 1.0
    - Distance is cumulative along sequence index
    """
    discs: List[Disc] = []
    prev = None
    from math import sqrt

    def dist(a, b):
        return sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)

    with open(path, "r") as f:
        idx = 0
        cumulative = 0.0
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                idx += 1
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                occ = float(line[54:60]) if len(line) >= 60 else 1.0
                bfac = float(line[60:66]) if len(line) >= 66 else 1.0
                r = bfac if bfac > 0 else (occ if occ > 0 else 1.0)
                if prev is not None:
                    cumulative += dist(prev, (x, y, z))
                discs.append(Disc(index=idx, x=x, y=y, z=z, radius=r, distance=cumulative))
                prev = (x, y, z)
    return discs


def load_discs(
    caver_output_dir: str | os.PathLike,
    tunnel_index: int = 1,
) -> List[Disc]:
    """
    Load discs from CAVER outputs by trying common file names:
    - profile.csv / profile.txt within tunnel dir
    - centerline.pdb / spheres.pdb within tunnel dir
    - profiles.csv at root
    """
    root = Path(caver_output_dir)
    tdir = find_tunnel_dir(root, tunnel_index) or root

    # First, try CAVER 3 analysis summary which encodes X/Y/Z/R/distance as rows
    analysis_profiles = tdir / "analysis" / "tunnel_profiles.csv"
    if analysis_profiles.exists():
        return _parse_analysis_tunnel_profiles(analysis_profiles, tunnel_index)

    candidates = [
        tdir / "profile.csv",
        tdir / "profile.txt",
        tdir / "profiles.csv",
        root / "profiles.csv",
    ]
    for p in candidates:
        if p.exists():
            return parse_profile_csv(p)

    pdb_candidates = [
        tdir / "centerline.pdb",
        tdir / "spheres.pdb",
        tdir / f"tunnel_{tunnel_index}.pdb",
    ]
    for p in pdb_candidates:
        if p.exists():
            return parse_centerline_pdb(p)

    raise FileNotFoundError(
        f"Could not find CAVER discs in {caver_output_dir} (tunnel {tunnel_index})."
    )


def _parse_analysis_tunnel_profiles(path: Path, tunnel_index: int) -> List[Disc]:
    """
    Parse CAVER 3 analysis/tunnel_profiles.csv which stores rows for axes and radius.
    We select the first snapshot and requested tunnel index, then zip X,Y,Z,R,distance.
    """
    import csv

    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        for r in reader:
            # Expect: Snapshot, Tunnel cluster, Tunnel, ..., Axis, v1, v2, ...
            if len(r) < 6:
                continue
            snap = r[0].strip()
            clus = r[1].strip()
            tun = r[2].strip()
            axis = r[12].strip() if len(r) > 12 else ""
            values = [v.strip() for v in r[13:] if v.strip() != "-"]
            try:
                tun_i = int(tun)
            except ValueError:
                continue
            if tun_i != tunnel_index:
                continue
            rows.append((axis, values))

    def to_floats(vals):
        out = []
        for v in vals:
            try:
                out.append(float(v))
            except ValueError:
                pass
        return out

    # Collect series by axis name
    series = {axis: to_floats(vals) for axis, vals in rows}
    required = ["X", "Y", "Z", "R", "distance"]
    if not all(k in series for k in required):
        raise FileNotFoundError(
            f"Missing required series {required} in {path}. Found: {list(series.keys())}"
        )
    n = min(len(series["X"]), len(series["Y"]), len(series["Z"]), len(series["R"]), len(series["distance"]))
    discs: List[Disc] = []
    for i in range(n):
        discs.append(
            Disc(
                index=i + 1,
                x=series["X"][i],
                y=series["Y"][i],
                z=series["Z"][i],
                radius=series["R"][i],
                distance=series["distance"][i],
            )
        )
    return discs


def write_centerline_pdb(discs: List[Disc], out_path: str | os.PathLike) -> None:
    """
    Write a simple PDB with one pseudo-atom per disc. Stores radius in B-factor and
    cumulative distance in occupancy.
    """
    p = Path(out_path)
    lines = []
    for d in discs:
        # PDB columns: HETATM, serial, name, resName, chain, resSeq, x, y, z, occupancy, tempFactor, element
        # Use OCC (occupancy) to store distance and B-factor to store radius to be consistent with our parser.
        lines.append(
            f"HETATM{d.index:5d}  S{d.index%100:2d} TUN  A{d.index:4d}    "
            f"{d.x:8.3f}{d.y:8.3f}{d.z:8.3f}{d.distance:6.2f}{d.radius:6.2f}          S\n"
        )
    lines.append("END\n")
    p.write_text("".join(lines))


def load_profile_series(caver_output_dir: str | os.PathLike, tunnel_index: int = 1) -> Tuple[List[float], Tuple[List[float], List[float], List[float]], List[float]]:
    """
    Load (distance, (X,Y,Z), R) series for a tunnel suitable for discretization.
    """
    root = Path(caver_output_dir)
    tdir = find_tunnel_dir(root, tunnel_index) or root
    analysis_profiles = tdir / "analysis" / "tunnel_profiles.csv"
    if not analysis_profiles.exists():
        raise FileNotFoundError(f"Missing analysis/tunnel_profiles.csv in {tdir}")
    discs = _parse_analysis_tunnel_profiles(analysis_profiles, tunnel_index)
    s = [d.distance for d in discs]
    x = [d.x for d in discs]
    y = [d.y for d in discs]
    z = [d.z for d in discs]
    r = [d.radius for d in discs]
    return s, (x, y, z), r
