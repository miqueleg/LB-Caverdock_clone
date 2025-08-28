#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import sqrt
from pathlib import Path

from caverdock_lowerbound import caver as caver_mod


@dataclass
class Disc:
    index: int
    x: float
    y: float
    z: float
    radius: float
    distance: float


def parse_dsd(path: str) -> list[Disc]:
    discs = []
    prev = None
    dist = 0.0
    with open(path, "r") as f:
        for i, line in enumerate(f, start=1):
            parts = line.split()
            if len(parts) < 7:
                continue
            x, y, z = map(float, parts[:3])
            r = float(parts[6])
            if prev is not None:
                dist += sqrt((x - prev[0]) ** 2 + (y - prev[1]) ** 2 + (z - prev[2]) ** 2)
            discs.append(Disc(i, x, y, z, r, dist))
            prev = (x, y, z)
    return discs


def discs_to_pdb(discs, out_path: str) -> None:
    # Reuse the library writer by mapping to its Disc type if available
    if discs and isinstance(discs[0], Disc):
        # Convert to library Disc type
        lib_discs = [caver_mod.Disc(d.index, d.x, d.y, d.z, d.radius, d.distance) for d in discs]
        caver_mod.write_centerline_pdb(lib_discs, out_path)
    else:
        caver_mod.write_centerline_pdb(discs, out_path)


def main():
    ap = argparse.ArgumentParser(description="Export tunnel centerline to PDB for comparison")
    sub = ap.add_subparsers(dest="cmd", required=True)

    c3 = sub.add_parser("from-caver3", help="Export from CAVER 3 output directory")
    c3.add_argument("--caver-output", required=True)
    c3.add_argument("--tunnel-index", type=int, default=1)
    c3.add_argument("--out", required=True)

    cd = sub.add_parser("from-dsd", help="Export from CaverDock tunnel.dsd file")
    cd.add_argument("--dsd", required=True)
    cd.add_argument("--out", required=True)

    args = ap.parse_args()

    if args.cmd == "from-caver3":
        discs = caver_mod.load_discs(args.caver_output, tunnel_index=args.tunnel_index)
        caver_mod.write_centerline_pdb(discs, args.out)
        print(f"Wrote {args.out} from CAVER 3 outputs")
    else:
        discs = parse_dsd(args.dsd)
        discs_to_pdb(discs, args.out)
        print(f"Wrote {args.out} from tunnel.dsd")


if __name__ == "__main__":
    main()

