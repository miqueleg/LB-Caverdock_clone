#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def load_cd(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Best-effort column normalization
    dist_col = None
    for c in ["distance", "s", "step"]:
        if c in df.columns:
            dist_col = c
            break
    energy_col = None
    for c in ["upper_bound_energy", "upper", "upper_bound", "energy_upper", "energy"]:
        if c in df.columns:
            energy_col = c
            break
    if dist_col is None or energy_col is None:
        raise ValueError("Could not infer columns in CaverDock profile")
    return df.rename(columns={dist_col: "distance", energy_col: "energy"})[["distance", "energy"]]


def main():
    ap = argparse.ArgumentParser(description="Plot LinB WT and L177W comparison")
    ap.add_argument("--wt", required=True, help="WT vina_profile.csv")
    ap.add_argument("--wt-cd", required=True, help="WT CaverDock profile CSV")
    ap.add_argument("--mut", required=True, help="L177W vina_profile.csv")
    ap.add_argument("--mut-cd", required=True, help="L177W CaverDock profile CSV")
    ap.add_argument("--out", required=True, help="Output PNG path")
    args = ap.parse_args()

    wt_v = pd.read_csv(args.wt)
    wt_cd = load_cd(args.wt_cd)
    mu_v = pd.read_csv(args.mut)
    mu_cd = load_cd(args.mut_cd)

    fig, axs = plt.subplots(1, 2, figsize=(12, 4), sharey=True)

    # WT
    x = wt_v["distance"] if "distance" in wt_v.columns else wt_v["index"]
    axs[0].plot(x, wt_v["vina_score"], label="Lower-bound (Vina per-disc)", color="tab:blue")
    axs[0].plot(wt_cd["distance"], wt_cd["energy"], label="CaverDock upper-bound", color="tab:orange")
    axs[0].set_title("LinB WT + DBE")
    axs[0].set_xlabel("Distance (Å)")
    axs[0].set_ylabel("Energy (kcal/mol)")
    axs[0].grid(True, alpha=0.3)
    axs[0].legend(frameon=False)

    # L177W
    x = mu_v["distance"] if "distance" in mu_v.columns else mu_v["index"]
    axs[1].plot(x, mu_v["vina_score"], label="Lower-bound (Vina per-disc)", color="tab:blue")
    axs[1].plot(mu_cd["distance"], mu_cd["energy"], label="CaverDock upper-bound", color="tab:orange")
    axs[1].set_title("LinB L177W + DBE")
    axs[1].set_xlabel("Distance (Å)")
    axs[1].grid(True, alpha=0.3)

    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
