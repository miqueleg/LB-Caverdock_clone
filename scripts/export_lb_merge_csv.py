#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="Merge Vina LB profile with CaverDock LB (by nearest distance)")
    ap.add_argument("--vina-profile", required=True)
    ap.add_argument("--cd-lb-csv", required=True)
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()

    vina = pd.read_csv(args.vina_profile).sort_values("distance")
    cd = pd.read_csv(args.cd_lb_csv).sort_values("distance")
    merged = pd.merge_asof(vina, cd, on="distance", direction="nearest", suffixes=("_vina","_cd"))
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv} with {len(merged)} rows")


if __name__ == "__main__":
    main()

