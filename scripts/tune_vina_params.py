#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import subprocess as sp
import os
import pandas as pd


def run(cmd: list[str]):
    env = os.environ.copy()
    # ensure local package is importable
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1] / "src") + os.pathsep + env.get("PYTHONPATH", "")
    cp = sp.run(cmd, capture_output=True, text=True, env=env)
    if cp.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{cp.stdout}\n{cp.stderr}")
    return cp.stdout


def compute_rmsd(vina_csv: Path, cd_csv: Path) -> float:
    import numpy as np
    v = pd.read_csv(vina_csv).sort_values("index")
    c = pd.read_csv(cd_csv).sort_values("index")
    m = pd.merge(v[["index","vina_score"]], c[["index","cd_lb_energy"]], on="index", how="inner")
    if len(m) == 0:
        return float("inf")
    a = m["vina_score"].to_numpy()
    b = m["cd_lb_energy"].to_numpy()
    return float(np.sqrt(((a-b)**2).mean()))


def main():
    ap = argparse.ArgumentParser(description="Parameter sweep to match CaverDock LB profile")
    ap.add_argument("--dsd", required=True)
    ap.add_argument("--receptor-pdbqt", required=True)
    ap.add_argument("--ligand-pdbqt", required=True)
    ap.add_argument("--cd-pdbqt", required=True)
    ap.add_argument("--vina-bin", default="vina")
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--scoring", nargs="*", default=["vina","vinardo"]) 
    ap.add_argument("--exhaustiveness", nargs="*", type=int, default=[8])
    ap.add_argument("--margin", nargs="*", type=float, default=[1.5])
    ap.add_argument("--min-box", nargs="*", type=float, default=[8.0])
    ap.add_argument("--mode", choices=["independent","centered_local"], default="centered_local")
    # centered_local parameters
    ap.add_argument("--local-cycles", nargs="*", type=int, default=[1,2,3])
    ap.add_argument("--local-margin", nargs="*", type=float, default=[0.25,0.5,0.75])
    ap.add_argument("--local-min-box", nargs="*", type=float, default=[8.0,10.0])
    ap.add_argument("--local-max-box", nargs="*", type=float, default=[10.0,12.0])
    ap.add_argument("--cpu", type=int, default=4)
    args = ap.parse_args()

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    # Export CD LB profile once
    cd_csv = out_root / "cd_lb_profile.csv"
    if not cd_csv.exists():
        run(["python3", "scripts/export_cd_lb_csv.py", "--pdbqt", args.cd_pdbqt, "--out-csv", str(cd_csv)])

    results = []
    if args.mode == "independent":
        grid = itertools.product(args.scoring, args.exhaustiveness, args.margin, args.min_box)
    else:
        grid = itertools.product(args.scoring, args.local_cycles, args.local_margin, args.local_min_box, args.local_max_box)

    for combo in grid:
        if args.mode == "independent":
            scoring, ex, margin, minbox = combo
            label = f"ind_sc_{scoring}_ex{ex}_m{margin}_min{minbox}"
            cmd = [
                "python3", "scripts/run_lb_on_cd_discs.py",
                "--mode", "independent",
                "--dsd", args.dsd,
                "--receptor-pdbqt", args.receptor_pdbqt,
                "--ligand-pdbqt", args.ligand_pdbqt,
                "--out-dir", str(out_root / label),
                "--vina-bin", args.vina_bin,
                "--exhaustiveness", str(ex),
                "--cpu", str(args.cpu),
                "--margin", str(margin),
                "--min-box-size", str(minbox),
                "--scoring", scoring,
            ]
        else:
            scoring, cycles, lmargin, lmin, lmax = combo
            label = f"cent_sc_{scoring}_cyc{cycles}_lm{lmargin}_min{lmin}_max{lmax}"
            cmd = [
                "python3", "scripts/run_lb_on_cd_discs.py",
                "--mode", "centered_local",
                "--dsd", args.dsd,
                "--receptor-pdbqt", args.receptor_pdbqt,
                "--ligand-pdbqt", args.ligand_pdbqt,
                "--out-dir", str(out_root / label),
                "--vina-bin", args.vina_bin,
                "--cpu", str(args.cpu),
                "--local-cycles", str(cycles),
                "--local-margin", str(lmargin),
                "--local-min-box", str(lmin),
                "--local-max-box", str(lmax),
                "--scoring", scoring,
            ]
        out_dir = out_root / label
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            run(cmd)
            rmsd = compute_rmsd(out_dir / "vina_profile.csv", cd_csv)
            results.append({
                "label": label,
                **({"scoring": scoring} if True else {}),
                **({"exhaustiveness": ex, "margin": margin, "min_box": minbox} if args.mode=="independent" else {}),
                **({"local_cycles": cycles, "local_margin": lmargin, "local_min_box": lmin, "local_max_box": lmax} if args.mode=="centered_local" else {}),
                "rmsd": rmsd,
            })
            print(f"{label}: RMSD={rmsd:.3f}")
        except Exception as e:
            print(f"{label}: failed: {e}")
            continue

    res_df = pd.DataFrame(results).sort_values("rmsd")
    res_csv = out_root / "tuning_results.csv"
    res_df.to_csv(res_csv, index=False)
    print(f"Wrote {res_csv}")

    if len(res_df):
        best = res_df.iloc[0]
        best_dir = out_root / best["label"]
        # Also write a plot for the best
        run(["python3", "scripts/plot_lb_svg.py", "--vina-csv", str(best_dir/"vina_profile.csv"), "--cd-csv", str(cd_csv), "--dsd", args.dsd, "--out", str(out_root/"best_compare.svg")])
        print("Best:")
        print(best.to_string())


if __name__ == "__main__":
    main()
