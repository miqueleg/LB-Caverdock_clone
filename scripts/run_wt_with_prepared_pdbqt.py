#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from caverdock_lowerbound.pipeline import run_pipeline, PipelineConfig
from caverdock_lowerbound.charge_providers import StaticPDBQTProvider


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor-pdbqt", required=True)
    ap.add_argument("--ligand-pdbqt", required=True)
    ap.add_argument("--caver-output", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--vina-bin", default="vina")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--num-modes", type=int, default=1)
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--margin", type=float, default=2.0)
    args = ap.parse_args()

    provider = StaticPDBQTProvider(args.receptor_pdbqt, args.ligand_pdbqt)
    cfg = PipelineConfig(
        vina_bin=args.vina_bin,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        cpu=args.cpu,
        margin=args.margin,
    )

    res = run_pipeline(
        receptor_input="",  # unused
        ligand_input="",  # unused
        caver_output_dir=args.caver_output,
        out_dir=args.out_dir,
        charge_provider=provider,
        config=cfg,
    )
    print(f"Wrote {Path(res['out_dir'])/'vina_profile.csv'}")


if __name__ == "__main__":
    main()
