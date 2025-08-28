#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import os
from pathlib import Path
import subprocess as sp
import pandas as pd
import numpy as np


def run(cmd: list[str]):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1] / "src") + os.pathsep + env.get("PYTHONPATH", "")
    cp = sp.run(cmd, capture_output=True, text=True, env=env)
    if cp.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{cp.stdout}\n{cp.stderr}")
    return cp.stdout


def com_rms(vina_csv: Path) -> float:
    df = pd.read_csv(vina_csv)
    if "com_dist" not in df.columns:
        return float("inf")
    d = df["com_dist"].dropna().to_numpy()
    if d.size == 0:
        return float("inf")
    return float(np.sqrt((d**2).mean()))


def main():
    ap = argparse.ArgumentParser(description="Tune tight-box + strict placement by minimizing COM RMS")
    ap.add_argument("--dsd", required=True)
    ap.add_argument("--receptor-pdbqt", required=True)
    ap.add_argument("--ligand-pdbqt", required=True)
    ap.add_argument("--vina-bin", default="vina")
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--cpu", type=int, default=4)
    # Sweep ranges
    ap.add_argument("--tight-margin", nargs="*", type=float, default=[0.25, 0.5, 0.75])
    ap.add_argument("--tight-min-box", nargs="*", type=float, default=[6.0, 8.0, 10.0])
    ap.add_argument("--tight-max-box", nargs="*", type=float, default=[10.0, 12.0])
    ap.add_argument("--strict-com-tol", nargs="*", type=float, default=[0.5, 0.75, 1.0])
    ap.add_argument("--strict-shell-tol", nargs="*", type=float, default=[0.25, 0.5, 0.75])
    args = ap.parse_args()

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    results = []
    for tm, tmin, tmax, ctol, stol in itertools.product(
        args.tight_margin, args.tight_min_box, args.tight_max_box, args.strict_com_tol, args.strict_shell_tol
    ):
        label = f"tm{tm}_min{tmin}_max{tmax}_ctol{ctol}_stol{stol}"
        out = out_root / label
        out.mkdir(parents=True, exist_ok=True)
        try:
            run([
                "python3", "scripts/run_lb_on_cd_discs.py",
                "--mode", "independent_tight",
                "--dsd", args.dsd,
                "--receptor-pdbqt", args.receptor_pdbqt,
                "--ligand-pdbqt", args.ligand_pdbqt,
                "--out-dir", str(out),
                "--vina-bin", args.vina_bin,
                "--cpu", str(args.cpu),
                "--tight-margin", str(tm),
                "--tight-min-box", str(tmin),
                "--tight-max-box", str(tmax),
                "--strict-com-tol", str(ctol),
                "--strict-shell-tol", str(stol),
            ])
            rms = com_rms(out / "vina_profile.csv")
            results.append({
                "label": label,
                "tight_margin": tm,
                "tight_min_box": tmin,
                "tight_max_box": tmax,
                "strict_com_tol": ctol,
                "strict_shell_tol": stol,
                "com_rms": rms,
            })
            print(f"{label}: COM_RMS={rms:.3f} Å")
        except Exception as e:
            print(f"{label}: failed: {e}")

    res = pd.DataFrame(results).sort_values("com_rms")
    res.to_csv(out_root / "com_tuning_results.csv", index=False)
    print(f"Wrote {out_root/'com_tuning_results.csv'}")
    if len(res):
        print("Best:")
        print(res.iloc[0].to_string())


if __name__ == "__main__":
    main()

