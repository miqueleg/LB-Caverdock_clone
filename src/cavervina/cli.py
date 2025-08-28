import argparse
import json
import os
import sys
from pathlib import Path

from .discs import load_discs, write_boxes_csv
from .prepare import prepare_receptor_pdbqt, prepare_ligand_pdbqt
from .vina import generate_disc_boxes, run_vina_for_discs, parse_energy_profile
from .profile import filter_poses_near_disc, export_profile_csv


def _path(p: str) -> Path:
    return Path(p).expanduser().resolve()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="cavervina",
        description=(
            "Recreate CaverDock upper-bound energies using CAVER3 tunnels and AutoDock Vina."
        ),
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # discs load/generate boxes
    p_discs = sub.add_parser("discs", help="Load discs and generate vina boxes")
    p_discs.add_argument("discs_file", help="Path to discs data from CAVER (txt/csv)")
    p_discs.add_argument(
        "--box-margin",
        type=float,
        default=3.0,
        help="Extra Å added to each disc radius for Vina box size.",
    )
    p_discs.add_argument(
        "--out",
        type=_path,
        required=True,
        help="Output directory for computed boxes (CSV).",
    )

    # prepare receptor/ligand
    p_prep = sub.add_parser("prep", help="Prepare receptor and ligand PDBQT")
    p_prep.add_argument("receptor", help="Receptor PDB/PDBQT input path")
    p_prep.add_argument("ligand", help="Ligand input (PDB/MOL2/SDF/PDBQT)")
    p_prep.add_argument("--outdir", type=_path, required=True)
    p_prep.add_argument("--adt-bin", type=_path, help="Path to AutoDockTools scripts dir")
    p_prep.add_argument("--obabel", type=_path, help="Path to obabel binary (optional)")

    # run vina per disc
    p_run = sub.add_parser("run", help="Run Vina for each disc (parallel)")
    p_run.add_argument("boxes_csv", help="CSV produced by 'discs' step")
    p_run.add_argument("receptor_pdbqt", help="Prepared receptor PDBQT")
    p_run.add_argument("ligand_pdbqt", help="Prepared ligand PDBQT")
    p_run.add_argument("--vina", type=_path, default=Path("vina"), help="vina binary")
    p_run.add_argument("--exhaustiveness", type=int, default=16)
    p_run.add_argument("--num-poses", type=int, default=20)
    p_run.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    p_run.add_argument("--outdir", type=_path, required=True)

    # filter and profile
    p_prof = sub.add_parser("profile", help="Filter poses near discs and export energies")
    p_prof.add_argument("boxes_csv", help="Boxes CSV (with disc centers/radii)")
    p_prof.add_argument("vina_outdir", help="Directory with Vina outputs per disc")
    p_prof.add_argument("--thickness", type=float, default=2.0, help="Disc half-thickness Å")
    p_prof.add_argument("--out", type=_path, required=True, help="Output CSV path")

    # compare to CaverDock
    p_cmp = sub.add_parser("compare", help="Compare profile to CaverDock energies")
    p_cmp.add_argument("our_csv", help="CSV from 'profile' step")
    p_cmp.add_argument("cdock_csv", help="CaverDock energy CSV")

    args = parser.parse_args(argv)

    if args.cmd == "discs":
        discs = load_discs(args.discs_file)
        boxes = generate_disc_boxes(discs, box_margin=args.box_margin)
        os.makedirs(args.out, exist_ok=True)
        out_csv = Path(args.out) / "boxes.csv"
        write_boxes_csv(boxes, out_csv)
        print(str(out_csv))
        return 0

    if args.cmd == "prep":
        os.makedirs(args.outdir, exist_ok=True)
        rec_out = Path(args.outdir) / "receptor.pdbqt"
        lig_out = Path(args.outdir) / "ligand.pdbqt"
        prepare_receptor_pdbqt(args.receptor, rec_out, adt_bin=args.adt_bin, obabel=args.obabel)
        prepare_ligand_pdbqt(args.ligand, lig_out, adt_bin=args.adt_bin, obabel=args.obabel)
        print(json.dumps({"receptor": str(rec_out), "ligand": str(lig_out)}))
        return 0

    if args.cmd == "run":
        os.makedirs(args.outdir, exist_ok=True)
        energies = run_vina_for_discs(
            boxes_csv=args.boxes_csv,
            receptor_pdbqt=args.receptor_pdbqt,
            ligand_pdbqt=args.ligand_pdbqt,
            vina_bin=str(args.vina),
            exhaustiveness=args.exhaustiveness,
            num_poses=args.num_poses,
            threads=args.threads,
            outdir=args.outdir,
        )
        print(json.dumps(energies, indent=2))
        return 0

    if args.cmd == "profile":
        export_profile_csv(
            boxes_csv=args.boxes_csv,
            vina_outdir=args.vina_outdir,
            disc_half_thickness=args.thickness,
            out_csv=args.out,
        )
        print(str(args.out))
        return 0

    if args.cmd == "compare":
        from .profile import compare_profiles

        diff = compare_profiles(args.our_csv, args.cdock_csv)
        print(json.dumps(diff, indent=2))
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())

