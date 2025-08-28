from __future__ import annotations

import argparse
from pathlib import Path

from .charge_providers import MGLToolsGasteigerChargeProvider
from .caver import run_caver
from .pipeline import PipelineConfig, run_pipeline
from .plotting import plot_profiles


def main() -> None:
    p = argparse.ArgumentParser(
        prog="caverdock-lb",
        description="Recreate CaverDock-style lower-bound using Vina per-disc docking.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    cav_p = sub.add_parser("caver", help="Run CAVER 3 to compute tunnels")
    cav_p.add_argument("--pdb", required=True, help="Receptor PDB structure")
    cav_p.add_argument("--out-dir", required=True, help="CAVER output directory")
    cav_p.add_argument(
        "--caver-cmd",
        nargs=argparse.ONE_OR_MORE,
        default=["caver"],
        help="Command vector to run CAVER (e.g., 'caver' or 'java -jar Caver3.jar')",
    )
    cav_p.add_argument("--config", help="Optional CAVER config file")

    run_p = sub.add_parser("run", help="Run full pipeline: CAVER3 -> discretize -> prep -> Vina per-disc -> plot")
    run_p.add_argument("--receptor", required=True, help="Receptor PDB/PDBQT input for preparation")
    run_p.add_argument("--ligand", required=True, help="Ligand PDB/MOL2/PDBQT input for preparation")
    run_p.add_argument("--caver-output", required=True, help="CAVER output directory with tunnel data")
    run_p.add_argument("--tunnel-index", type=int, default=1)
    run_p.add_argument("--dsd", help="Use an existing tunnel.dsd (skip discretization)")
    run_p.add_argument("--spacing", type=float, default=0.5, help="Discretization spacing (Å)")
    run_p.add_argument("--out-dir", required=True)
    run_p.add_argument("--vina-bin", default="vina")
    run_p.add_argument("--exhaustiveness", type=int, default=8)
    run_p.add_argument("--num-modes", type=int, default=1)
    run_p.add_argument("--cpu", type=int, default=1)
    run_p.add_argument("--margin", type=float, default=2.0, help="Box margin beyond 2R")
    run_p.add_argument("--seed", type=int)
    run_p.add_argument("--lb-mode", choices=["local","independent"], default="local", help="Lower-bound mode: local (sequential local_only) or independent (separate docks)")
    run_p.add_argument("--min-box-size", type=float, default=8.0, help="Minimum Vina box size (Å)")
    # MGLTools
    run_p.add_argument("--prepare-receptor-script", default="prepare_receptor4.py")
    run_p.add_argument("--prepare-ligand-script", default="prepare_ligand4.py")
    run_p.add_argument("--mgltools-pythonsh", help="Path to MGLTools pythonsh interpreter")
    run_p.add_argument("--per-sphere-receptor", action="store_true")
    run_p.add_argument("--per-sphere-ligand", action="store_true")
    # Compare to CaverDock profile
    run_p.add_argument("--caverdock-profile", help="CSV with CaverDock upper-bound profile")

    args = p.parse_args()

    if args.cmd == "caver":
        run_caver(
            pdb_path=args.pdb,
            output_dir=args.out_dir,
            caver_cmd=args.caver_cmd,
            config_path=args.config,
        )
        print(f"CAVER outputs written to {args.out_dir}")

    elif args.cmd == "run":
        out = Path(args.out_dir)
        provider = MGLToolsGasteigerChargeProvider(
            receptor_input=args.receptor,
            ligand_input=args.ligand,
            out_dir=str(out / "prep"),
            prepare_receptor_script=args.prepare_receptor_script,
            prepare_ligand_script=args.prepare_ligand_script,
            mgltools_pythonsh=args.mgltools_pythonsh,
            per_sphere_recompute_receptor=args.per_sphere_receptor,
            per_sphere_recompute_ligand=args.per_sphere_ligand,
        )
        cfg = PipelineConfig(
            vina_bin=args.vina_bin,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            cpu=args.cpu,
            margin=args.margin,
            seed=args.seed,
            lb_mode=args.lb_mode,
            min_box_size=args.min_box_size,
        )
        result = run_pipeline(
            receptor_input=args.receptor,
            ligand_input=args.ligand,
            caver_output_dir=args.caver_output,
            out_dir=args.out_dir,
            tunnel_index=args.tunnel_index,
            dsd_path=args.dsd,
            spacing=args.spacing,
            charge_provider=provider,
            config=cfg,
            caverdock_profile_csv=args.caverdock_profile,
        )
        vina_csv = str(Path(result["out_dir"]) / "vina_profile.csv")
        plot_path = str(Path(result["out_dir"]) / "profile_compare.png")
        plot_profiles(
            vina_profile_csv=vina_csv,
            out_png=plot_path,
            caverdock_profile_csv=args.caverdock_profile,
        )
        print(f"Wrote: {vina_csv}\nPlot: {plot_path}")


if __name__ == "__main__":
    main()
