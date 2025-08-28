#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from caverdock_lowerbound.vina import VinaBox, run_vina


def main():
    ap = argparse.ArgumentParser(description="Score existing per-disc poses with Vina (score_only)")
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--vina-bin", default="vina")
    ap.add_argument("--out-csv", default=None)
    args = ap.parse_args()

    run = Path(args.run_dir)
    discs_csv = run / "discs.csv"
    out_csv = Path(args.out_csv) if args.out_csv else (run / "vina_profile.csv")

    rows_out = []
    with open(discs_csv, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            idx = int(row["index"]) if "index" in row else int(float(row["index"]))
            lig = run / f"dock_disc_{idx}.pdbqt"
            if not lig.exists():
                continue
            cx = float(row["x"]) ; cy = float(row["y"]) ; cz = float(row["z"]) ; rad = float(row["radius"]) ; dist = float(row["distance"]) ;
            size = max(6.0, min(12.0, 2.0*rad + 0.5))
            box = VinaBox(cx, cy, cz, size, size, size)
            score, _ = run_vina(
                vina_bin=args.vina_bin,
                receptor_pdbqt=str(args.receptor),
                ligand_pdbqt=str(lig),
                box=box,
                out_pdbqt=str(lig),
                log_path=None,
                cpu=1,
                mode="score_only",
                scoring="vina",
            )
            rows_out.append({
                "index": idx,
                "distance": dist,
                "radius": rad,
                "x": cx,
                "y": cy,
                "z": cz,
                "vina_score": score,
            })

    rows_out.sort(key=lambda d: d["index"]) 
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["index","distance","radius","x","y","z","vina_score"]) 
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"Wrote {out_csv} with {len(rows_out)} rows")


if __name__ == "__main__":
    main()

