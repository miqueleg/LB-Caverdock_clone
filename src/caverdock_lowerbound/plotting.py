from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd


def plot_profiles(
    vina_profile_csv: str,
    out_png: str,
    caverdock_profile_csv: Optional[str] = None,
    caverdock_energy_column: str = "upper_bound_energy",
    caverdock_distance_column: str = "distance",
) -> None:
    vina_df = pd.read_csv(vina_profile_csv)
    fig, ax = plt.subplots(figsize=(8, 4))
    if "distance" in vina_df.columns and not vina_df["distance"].isna().all():
        x = vina_df["distance"]
        ax.set_xlabel("Distance along tunnel (Å)")
    else:
        x = vina_df["index"]
        ax.set_xlabel("Disc index")

    ax.plot(x, vina_df["vina_score"], label="Lower-bound (Vina per-disc)", color="tab:blue")

    if caverdock_profile_csv and Path(caverdock_profile_csv).exists():
        cdf = pd.read_csv(caverdock_profile_csv)
        if caverdock_distance_column not in cdf.columns:
            for cand in ["s", "distance", "step"]:
                if cand in cdf.columns:
                    caverdock_distance_column = cand
                    break
        if caverdock_energy_column not in cdf.columns:
            for cand in ["upper", "upper_bound", "energy_upper", "energy"]:
                if cand in cdf.columns:
                    caverdock_energy_column = cand
                    break
        ax.plot(
            cdf[caverdock_distance_column],
            cdf[caverdock_energy_column],
            label="CaverDock upper-bound",
            color="tab:orange",
        )

    ax.set_ylabel("Binding energy (kcal/mol)")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)

