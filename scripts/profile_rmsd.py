#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pandas as pd
from pathlib import Path


def rmsd(a, b):
    import numpy as np
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n = min(len(a), len(b))
    if n == 0:
        return float("nan")
    a = a[:n]
    b = b[:n]
    return float(((a - b) ** 2).mean() ** 0.5)


def main():
    ap = argparse.ArgumentParser(description="Compute RMSD between Vina LB and CaverDock LB profiles (aligned by index)")
    ap.add_argument("--vina-csv", required=True)
    ap.add_argument("--cd-csv", required=True)
    args = ap.parse_args()

    v = pd.read_csv(args.vina_csv).sort_values("index")
    c = pd.read_csv(args.cd_csv).sort_values("index")

    merged = pd.merge(v[["index","vina_score"]], c[["index","cd_lb_energy"]], on="index", how="inner")
    r = rmsd(merged["vina_score"].values, merged["cd_lb_energy"].values)
    print(f"RMSD: {r:.4f} kcal/mol over {len(merged)} points")


if __name__ == "__main__":
    main()

