#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

# Use non-interactive backend and single-threaded math to avoid SHM issues
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import matplotlib.pyplot as plt
import pandas as pd


def parse_cd_lb_pdbqt(path: str) -> pd.DataFrame:
    energies = []
    distances = []
    radii = []
    indices = []

    model_idx = 0
    re_result = re.compile(r"^REMARK\s+CAVERDOCK RESULT:\s+(-?\d+\.\d+)\s+")
    re_tunnel = re.compile(r"^REMARK\s+CAVERDOCK TUNNEL:\s+(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")

    with open(path, "r") as f:
        for line in f:
            if line.startswith("MODEL"):
                model_idx += 1
            m2 = re_tunnel.match(line)
            if m2:
                idx = int(m2.group(1)) + 1  # make 1-based
                e = float(m2.group(2))
                r = float(m2.group(3))
                d = float(m2.group(4))
                indices.append(idx)
                energies.append(e)
                radii.append(r)
                distances.append(d)

    df = pd.DataFrame({
        "index": indices,
        "distance": distances,
        "radius": radii,
        "cd_lb_energy": energies,
    })
    return df.sort_values("index").reset_index(drop=True)


def parse_dsd_local(path: str) -> pd.DataFrame:
    """Self-contained DSD parser: returns DataFrame with index and cumulative distance.

    DSD format: x y z nx ny nz radius. We compute cumulative distance along (x,y,z).
    """
    xs, ys, zs = [], [], []
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 7:
                continue
            x, y, z = map(float, parts[:3])
            xs.append(x); ys.append(y); zs.append(z)
    dist = [0.0]
    import math
    for i in range(1, len(xs)):
        dx = xs[i] - xs[i-1]
        dy = ys[i] - ys[i-1]
        dz = zs[i] - zs[i-1]
        dist.append(dist[-1] + math.sqrt(dx*dx + dy*dy + dz*dz))
    return pd.DataFrame({"index": list(range(1, len(xs)+1)), "distance": dist})


def main():
    ap = argparse.ArgumentParser(description="Compare Vina per-disc LB vs CaverDock LB (PDBQT)")
    ap.add_argument("--vina-profile", required=True, help="runs/.../vina_profile.csv")
    ap.add_argument("--cd-pdbqt", required=True, help="CaverDock analysis-lb.pdbqt path")
    ap.add_argument("--out", required=True, help="Output PNG path")
    ap.add_argument("--cd-dsd", help="Optional CaverDock tunnel.dsd to enforce distance x-axis")
    ap.add_argument("--x-axis", choices=["distance","index"], default="distance")
    args = ap.parse_args()

    vina_df = pd.read_csv(args.vina_profile)
    cd_df = parse_cd_lb_pdbqt(args.cd_pdbqt)
    # If provided, use the DSD distances (cumulative) instead of any distances inferred from PDBQT
    if args.cd_dsd and Path(args.cd_dsd).exists():
        dsd_df = parse_dsd_local(args.cd_dsd)
        cd_df = cd_df.drop(columns=[c for c in ["distance"] if c in cd_df.columns])
        cd_df = pd.merge(cd_df, dsd_df, on="index", how="inner")

    # Plot using chosen x-axis to avoid misleading overlays
    fig, ax = plt.subplots(figsize=(8, 4))
    if args.x_axis == "distance":
        vina_df = vina_df.sort_values("distance")
        cd_df = cd_df.sort_values("distance")
        ax.plot(vina_df["distance"], vina_df["vina_score"], label="Lower-bound (Vina per-disc)", color="tab:blue")
        ax.plot(cd_df["distance"], cd_df["cd_lb_energy"], label="CaverDock LB (PDBQT)", color="tab:green")
        ax.set_xlabel("Distance along tunnel (Å)")
    else:
        ax.plot(vina_df["index"], vina_df["vina_score"], label="Lower-bound (Vina per-disc)", color="tab:blue")
        ax.plot(cd_df["index"], cd_df["cd_lb_energy"], label="CaverDock LB (PDBQT)", color="tab:green")
        ax.set_xlabel("Disc index")

    ax.set_ylabel("Binding energy (kcal/mol)")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
