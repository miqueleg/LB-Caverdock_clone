from __future__ import annotations

import json
from dataclasses import asdict, dataclass
import shutil
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import numpy as np

from . import caver as caver_mod
from . import discretizer as disc_mod
from .charge_providers import ChargeProvider, DiscContext, MGLToolsGasteigerChargeProvider
from . import vina as vina_mod
from .ligand_geom import read_pdbqt_coords, write_pdbqt_coords, recenter_ligand_to


@dataclass
class PipelineConfig:
    vina_bin: str = "vina"
    exhaustiveness: int = 8
    num_modes: int = 1
    cpu: int = 1
    margin: float = 2.0
    seed: Optional[int] = None
    lb_mode: str = "local"  # "local" = sequential local_only; "independent" = independent docks
    min_box_size: float = 8.0
    local_cycles: int = 3
    centering_fraction: float = 0.5
    local_box_margin: float = 0.5
    local_box_min: float = 6.0
    local_box_max: float = 12.0
    scoring_model: str = "ad4"


def discs_to_df(discs: List[object]) -> pd.DataFrame:
    return pd.DataFrame([asdict(d) for d in discs])


def run_pipeline(
    receptor_input: str,
    ligand_input: str,
    caver_output_dir: str,
    out_dir: str,
    tunnel_index: int = 1,
    dsd_path: Optional[str] = None,
    spacing: float = 0.5,
    charge_provider: Optional[ChargeProvider] = None,
    config: Optional[PipelineConfig] = None,
    caverdock_profile_csv: Optional[str] = None,
) -> Dict[str, object]:
    """
    Orchestrate loading discs, charge preparation and Vina docking per disc.

    Returns a dict with keys:
    - discs_df: DataFrame of discs
    - results_df: DataFrame with disc index, distance, radius, vina_score
    - out_dir: output directory used
    - caverdock_profile: optional DataFrame of provided profile
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    config = config or PipelineConfig()

    if dsd_path:
        discs_raw = disc_mod.parse_dsd(dsd_path)
        # Copy provided DSD into run directory
        try:
            shutil.copyfile(dsd_path, out / "tunnel.dsd")
        except Exception:
            pass
    else:
        # Load series from CAVER and discretize following CaverDock spacing
        s, (x, y, z), r = caver_mod.load_profile_series(caver_output_dir, tunnel_index)
        discs_raw = disc_mod.discretize_from_profile(
            distances=np.array(s),
            coords=(np.array(x), np.array(y), np.array(z)),
            radii=np.array(r),
            spacing=spacing,
        )
        # also write a dsd for reference
        disc_mod.write_dsd(discs_raw, out / "tunnel.dsd")

    # Convert to basic Disc-like dicts/DataFrame expected by rest of pipeline
    discs = [
        caver_mod.Disc(index=d.index, x=d.x, y=d.y, z=d.z, radius=d.radius, distance=d.distance)
        for d in discs_raw
    ]
    discs_df = discs_to_df(discs)
    discs_df.to_csv(out / "discs.csv", index=False)

    if charge_provider is None:
        charge_provider = MGLToolsGasteigerChargeProvider(
            receptor_input=receptor_input,
            ligand_input=ligand_input,
            out_dir=str(out / "prep"),
        )

    rows = []
    prev_out_lig: Optional[Path] = None
    # Prepare base ligand once
    base_lig = charge_provider.get_ligand_pdbqt(DiscContext(0, 0, 0, 0, 0, 0))
    base_lines, base_coords = read_pdbqt_coords(base_lig)

    for d in discs:
        ctx = DiscContext(
            index=d.index, x=d.x, y=d.y, z=d.z, radius=d.radius, distance=d.distance
        )
        rec_pdbqt = charge_provider.get_receptor_pdbqt(ctx)

        # Prepare input ligand for this disc
        if config.lb_mode == "local" and prev_out_lig is not None:
            # Start from previous optimized coords, re-center to current disc
            prev_lines, prev_coords = read_pdbqt_coords(prev_out_lig)
            lig_coords = recenter_ligand_to(prev_coords, (d.x, d.y, d.z))
            lig_input_path = out / f"lig_in_disc_{d.index}.pdbqt"
            write_pdbqt_coords(prev_lines, lig_coords, lig_input_path)
        else:
            # Start from base ligand re-centered to current disc
            lig_coords = recenter_ligand_to(base_coords, (d.x, d.y, d.z))
            lig_input_path = out / f"lig_in_disc_{d.index}.pdbqt"
            write_pdbqt_coords(base_lines, lig_coords, lig_input_path)

        # Box sizing with minimum size
        # Define box
        if config.lb_mode == "local":
            size = 2.0 * d.radius + config.local_box_margin
            size = max(config.local_box_min, min(size, config.local_box_max))
            box = vina_mod.VinaBox(d.x, d.y, d.z, size, size, size)
        else:
            box = vina_mod.box_from_sphere(d.x, d.y, d.z, d.radius, margin=config.margin)
            box.size_x = max(box.size_x, config.min_box_size)
            box.size_y = max(box.size_y, config.min_box_size)
            box.size_z = max(box.size_z, config.min_box_size)

        out_pdbqt = out / f"dock_disc_{d.index}.pdbqt"
        log_path = out / f"dock_disc_{d.index}.log"

        if config.lb_mode == "local":
            # Iterative COM restraint cycles
            current_input = Path(lig_input_path)
            for cyc in range(config.local_cycles):
                lines, coords = read_pdbqt_coords(current_input)
                if coords:
                    cx = sum(x for x, _, _ in coords) / len(coords)
                    cy = sum(y for _, y, _ in coords) / len(coords)
                    cz = sum(z for _, _, z in coords) / len(coords)
                else:
                    cx = cy = cz = 0.0
                frac = 1.0 if cyc == 0 else config.centering_fraction
                dx = (d.x - cx) * frac
                dy = (d.y - cy) * frac
                dz = (d.z - cz) * frac
                shifted = [(x + dx, y + dy, z + dz) for (x, y, z) in coords]
                tmp_in = out / f"lig_cycle{cyc}_disc_{d.index}.pdbqt"
                write_pdbqt_coords(lines, shifted, tmp_in)

                # local-only refinement toward local minimum
                _, _ = vina_mod.run_vina(
                    vina_bin=config.vina_bin,
                    receptor_pdbqt=rec_pdbqt,
                    ligand_pdbqt=str(tmp_in),
                    box=box,
                    out_pdbqt=str(out_pdbqt),
                    log_path=str(log_path),
                    cpu=config.cpu,
                    mode="local_only",
                )
                current_input = out_pdbqt

            # Score final pose; prefer AD4, but fall back to Vina if AD4 maps/flags unsupported
            try:
                score, _ = vina_mod.run_vina(
                    vina_bin=config.vina_bin,
                    receptor_pdbqt=rec_pdbqt,
                    ligand_pdbqt=str(out_pdbqt),
                    box=box,
                    out_pdbqt=str(out_pdbqt),
                    log_path=str(log_path),
                    cpu=config.cpu,
                    mode="score_only",
                    scoring=config.scoring_model,
                )
            except Exception:
                score, _ = vina_mod.run_vina(
                    vina_bin=config.vina_bin,
                    receptor_pdbqt=rec_pdbqt,
                    ligand_pdbqt=str(out_pdbqt),
                    box=box,
                    out_pdbqt=str(out_pdbqt),
                    log_path=str(log_path),
                    cpu=config.cpu,
                    mode="score_only",
                    scoring="vina",
                )
        else:
            score, _ = vina_mod.run_vina(
                vina_bin=config.vina_bin,
                receptor_pdbqt=rec_pdbqt,
                ligand_pdbqt=str(lig_input_path),
                box=box,
                out_pdbqt=str(out_pdbqt),
                log_path=str(log_path),
                exhaustiveness=config.exhaustiveness,
                num_modes=config.num_modes,
                seed=config.seed,
                cpu=config.cpu,
                mode="dock",
            )

        prev_out_lig = out_pdbqt

        rows.append({
            "index": d.index,
            "distance": d.distance,
            "radius": d.radius,
            "x": d.x,
            "y": d.y,
            "z": d.z,
            "vina_score": score,
        })

    # Always write a PDB tunnel for visualization
    try:
        caver_mod.write_centerline_pdb(discs, out / "tunnel.pdb")
    except Exception:
        pass

    results_df = pd.DataFrame(rows)
    results_df.sort_values("index", inplace=True)
    results_df.to_csv(out / "vina_profile.csv", index=False)

    caverdock_df = None
    if caverdock_profile_csv and Path(caverdock_profile_csv).exists():
        caverdock_df = pd.read_csv(caverdock_profile_csv)
        caverdock_df.to_csv(out / "caverdock_profile.csv", index=False)

    with open(out / "config.json", "w") as f:
        json.dump(
            {
                "pipeline": vars(config),
                "receptor_input": str(receptor_input),
                "ligand_input": str(ligand_input),
                "caver_output_dir": str(caver_output_dir),
                "tunnel_index": tunnel_index,
            },
            f,
            indent=2,
        )

    return {
        "discs_df": discs_df,
        "results_df": results_df,
        "out_dir": str(out),
        "caverdock_profile": caverdock_df,
    }
