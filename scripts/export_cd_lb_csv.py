#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import re


def parse_cd_lb_pdbqt(path: str) -> pd.DataFrame:
    indices = []
    energies = []
    distances = []
    radii = []
    re_tunnel = re.compile(r"^REMARK\s+CAVERDOCK TUNNEL:\s+(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")
    with open(path, "r") as f:
        for line in f:
            m = re_tunnel.match(line)
            if m:
                idx = int(m.group(1)) + 1
                e = float(m.group(2))
                r = float(m.group(3))
                d = float(m.group(4))
                indices.append(idx)
                energies.append(e)
                distances.append(d)
                radii.append(r)
    return pd.DataFrame({
        "index": indices,
        "distance": distances,
        "radius": radii,
        "cd_lb_energy": energies,
    }).sort_values("index").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description="Export CaverDock LB PDBQT energies to CSV")
    ap.add_argument("--pdbqt", required=True)
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()

    df = parse_cd_lb_pdbqt(args.pdbqt)
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv} with {len(df)} rows")


if __name__ == "__main__":
    main()

