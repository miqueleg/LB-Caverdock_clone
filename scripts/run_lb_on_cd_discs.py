#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import sqrt
from pathlib import Path

import pandas as pd
import numpy as np

from caverdock_lowerbound.charge_providers import StaticPDBQTProvider, DiscContext
from caverdock_lowerbound.vina import box_from_sphere, run_vina, VinaBox
from caverdock_lowerbound.ligand_geom import read_pdbqt_coords, write_pdbqt_coords, recenter_ligand_to


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


def main():
    ap = argparse.ArgumentParser(description="Run Vina per-disc using CaverDock tunnel.dsd discs")
    ap.add_argument("--dsd", required=True, help="Path to tunnel.dsd")
    ap.add_argument("--receptor-pdbqt", required=True)
    ap.add_argument("--ligand-pdbqt", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--vina-bin", default="vina")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--margin", type=float, default=2.0)
    ap.add_argument("--min-box-size", type=float, default=8.0)
    ap.add_argument("--scoring", choices=["vina","vinardo"], default="vina")
    ap.add_argument(
        "--mode",
        choices=["independent", "independent_tight", "centered_local"],
        default="independent",
        help=(
            "independent=default Vina per-disc docking; "
            "independent_tight=Vina per-disc with tight sphere-sized box; "
            "centered_local=recenter+local_only then score"
        ),
    )
    ap.add_argument("--tight-margin", type=float, default=0.5, help="Tight mode: margin added to 2R")
    ap.add_argument("--tight-min-box", type=float, default=6.0, help="Tight mode: min box size")
    ap.add_argument("--tight-max-box", type=float, default=12.0, help="Tight mode: max box size")
    ap.add_argument("--strict-com-tol", type=float, default=0.75, help="Max COM distance (Å) from sphere center")
    ap.add_argument("--strict-shell-tol", type=float, default=0.5, help="Allow atoms up to R+tol inside sphere")
    ap.add_argument("--local-cycles", type=int, default=2)
    ap.add_argument("--local-margin", type=float, default=0.5)
    ap.add_argument("--local-min-box", type=float, default=6.0)
    ap.add_argument("--local-max-box", type=float, default=12.0)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    discs = parse_dsd(args.dsd)
    pd.DataFrame([d.__dict__ for d in discs]).to_csv(out / "discs.csv", index=False)

    provider = StaticPDBQTProvider(args.receptor_pdbqt, args.ligand_pdbqt)
    base_lines, base_coords = read_pdbqt_coords(args.ligand_pdbqt)

    def parse_pdbqt_models(path: Path):
        models = []
        with open(path, "r") as f:
            lines = f.readlines()
        current = []
        energy = None
        in_model = False
        for ln in lines:
            if ln.startswith("MODEL"):
                current = [ln]
                energy = None
                in_model = True
            elif ln.startswith("ENDMDL"):
                current.append(ln)
                # extract coords
                coords = []
                for l in current:
                    if l.startswith("ATOM") or l.startswith("HETATM"):
                        try:
                            x = float(l[30:38]); y = float(l[38:46]); z = float(l[46:54])
                            coords.append((x, y, z))
                        except Exception:
                            pass
                    elif "REMARK VINA RESULT:" in l:
                        try:
                            parts = l.strip().split()
                            energy = float(parts[3])
                        except Exception:
                            energy = None
                if coords:
                    models.append((energy, coords, current))
                in_model = False
            else:
                if in_model:
                    current.append(ln)
        if not models:
            coords = []
            energy = None
            for l in lines:
                if l.startswith("ATOM") or l.startswith("HETATM"):
                    try:
                        x = float(l[30:38]); y = float(l[38:46]); z = float(l[46:54])
                        coords.append((x, y, z))
                    except Exception:
                        pass
                elif "REMARK VINA RESULT:" in l:
                    try:
                        parts = l.strip().split()
                        energy = float(parts[3])
                    except Exception:
                        energy = None
            if coords:
                models.append((energy, coords, lines))
        return models

    def com(coords):
        pts = np.array(coords)
        return pts.mean(axis=0)

    def coords_from_pdbqt_single(path: Path):
        coords = []
        try:
            with open(path, "r") as f:
                for l in f:
                    if l.startswith("ATOM") or l.startswith("HETATM"):
                        try:
                            x = float(l[30:38]); y = float(l[38:46]); z = float(l[46:54])
                            coords.append((x, y, z))
                        except Exception:
                            pass
        except Exception:
            pass
        return coords

    rows = []
    for d in discs:
        ctx = DiscContext(d.index, d.x, d.y, d.z, d.radius, d.distance)
        rec = provider.get_receptor_pdbqt(ctx)
        lig = provider.get_ligand_pdbqt(ctx)
        out_pdbqt = out / f"dock_disc_{d.index}.pdbqt"
        log_path = out / f"dock_disc_{d.index}.log"

        if args.mode == "independent":
            box = box_from_sphere(d.x, d.y, d.z, d.radius, margin=args.margin)
            box.size_x = max(box.size_x, args.min_box_size)
            box.size_y = max(box.size_y, args.min_box_size)
            box.size_z = max(box.size_z, args.min_box_size)
            score, _ = run_vina(
                vina_bin=args.vina_bin,
                receptor_pdbqt=rec,
                ligand_pdbqt=lig,
                box=box,
                out_pdbqt=str(out_pdbqt),
                log_path=str(log_path),
                exhaustiveness=args.exhaustiveness,
                num_modes=1,
                cpu=args.cpu,
                scoring=args.scoring,
            )
        elif args.mode == "independent_tight":
            # Default Vina search, but restrict to sphere vicinity so Vina chooses orientation
            size = max(args.tight_min_box, min(2.0 * d.radius + args.tight_margin, args.tight_max_box))
            box = VinaBox(d.x, d.y, d.z, size, size, size)
            score, _ = run_vina(
                vina_bin=args.vina_bin,
                receptor_pdbqt=rec,
                ligand_pdbqt=lig,
                box=box,
                out_pdbqt=str(out_pdbqt),
                log_path=str(log_path),
                exhaustiveness=args.exhaustiveness,
                num_modes=9,
                cpu=args.cpu,
                scoring=args.scoring,
            )
            # Post-filter by COM proximity and inside-sphere condition
            models = parse_pdbqt_models(out_pdbqt)
            center = np.array([d.x, d.y, d.z])
            rad = d.radius + args.strict_shell_tol
            best_idx = None
            best_e = None
            for i, (e, coords, mlines) in enumerate(models):
                c = com(coords)
                if np.linalg.norm(c - center) > args.strict_com_tol:
                    continue
                pts = np.array(coords)
                if np.any(np.linalg.norm(pts - center, axis=1) > rad):
                    continue
                if best_e is None or (e is not None and e < best_e):
                    best_e = e
                    best_idx = i
            if best_idx is not None:
                with open(out_pdbqt, "w") as f:
                    f.writelines(models[best_idx][2])
                score = models[best_idx][0]
            else:
                # Fallback: recenter and local refine then score
                lig_centered = out / f"lig_in_disc_{d.index}.pdbqt"
                centered_coords = recenter_ligand_to(base_coords, (d.x, d.y, d.z))
                write_pdbqt_coords(base_lines, centered_coords, lig_centered)
                _, _ = run_vina(
                    vina_bin=args.vina_bin,
                    receptor_pdbqt=rec,
                    ligand_pdbqt=str(lig_centered),
                    box=box,
                    out_pdbqt=str(out_pdbqt),
                    log_path=str(log_path),
                    cpu=args.cpu,
                    mode="local_only",
                    scoring=args.scoring,
                )
                score, _ = run_vina(
                    vina_bin=args.vina_bin,
                    receptor_pdbqt=rec,
                    ligand_pdbqt=str(out_pdbqt),
                    box=box,
                    out_pdbqt=str(out_pdbqt),
                    log_path=str(log_path),
                    cpu=args.cpu,
                    mode="score_only",
                    scoring=args.scoring,
                )
        else:
            # Recenter ligand to disc center and refine locally in a tight box, then score
            lig_centered = out / f"lig_in_disc_{d.index}.pdbqt"
            centered_coords = recenter_ligand_to(base_coords, (d.x, d.y, d.z))
            write_pdbqt_coords(base_lines, centered_coords, lig_centered)
            size = max(args.local_min_box, min(2.0 * d.radius + args.local_margin, args.local_max_box))
            box = VinaBox(d.x, d.y, d.z, size, size, size)
            # Compute tangent approx using neighbor discs
            if d.index == 1 and len(discs) > 1:
                tv = np.array([discs[1].x - d.x, discs[1].y - d.y, discs[1].z - d.z], float)
            elif d.index == len(discs) and len(discs) > 1:
                tv = np.array([d.x - discs[-2].x, d.y - discs[-2].y, d.z - discs[-2].z], float)
            elif len(discs) > 2:
                p0 = discs[d.index - 2]
                p2 = discs[d.index]
                tv = np.array([p2.x - p0.x, p2.y - p0.y, p2.z - p0.z], float)
            else:
                tv = np.array([1.0, 0.0, 0.0], float)
            nrm = np.linalg.norm(tv)
            if nrm < 1e-6:
                tv = np.array([1.0, 0.0, 0.0], float)
            tv = tv / np.linalg.norm(tv)

            # Align ligand principal axis to tangent
            def principal_axis(coords):
                pts = np.array(coords)
                cen = pts.mean(axis=0)
                pts0 = pts - cen
                cov = np.cov(pts0.T)
                w, v = np.linalg.eigh(cov)
                axis = v[:, np.argmax(w)]
                if axis.dot(tv) < 0:
                    axis = -axis
                return axis

            def rotation_matrix_from_vectors(a, b):
                a = a / np.linalg.norm(a)
                b = b / np.linalg.norm(b)
                v = np.cross(a, b)
                c = float(np.dot(a, b))
                if np.linalg.norm(v) < 1e-8:
                    return np.eye(3)
                vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
                R = np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))
                return R

            def apply_rotation(coords, R):
                pts = np.array(coords)
                cen = pts.mean(axis=0)
                pts0 = pts - cen
                ptsr = (R @ pts0.T).T + cen
                return [tuple(p) for p in ptsr]

            centered = recenter_ligand_to(base_coords, (d.x, d.y, d.z))
            axis0 = principal_axis(centered)
            R_align = rotation_matrix_from_vectors(axis0, tv)
            aligned = apply_rotation(centered, R_align)

            best_score = None
            best_pose_text = None
            # roll around tangent axis
            for roll in [0.0, 90.0, 180.0, 270.0]:
                theta = np.deg2rad(roll)
                k = tv
                K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
                R_roll = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
                rolled = apply_rotation(aligned, R_roll)
                tmp_in = out / f"lig_disc_{d.index}_roll_{int(roll)}.pdbqt"
                write_pdbqt_coords(base_lines, rolled, tmp_in)

                current = tmp_in
                for _ in range(max(1, args.local_cycles)):
                    _, _ = run_vina(
                        vina_bin=args.vina_bin,
                        receptor_pdbqt=rec,
                        ligand_pdbqt=str(current),
                        box=box,
                        out_pdbqt=str(out_pdbqt),
                        log_path=str(log_path),
                        cpu=args.cpu,
                        mode="local_only",
                        scoring=args.scoring,
                    )
                    current = out_pdbqt
                sc, _ = run_vina(
                    vina_bin=args.vina_bin,
                    receptor_pdbqt=rec,
                    ligand_pdbqt=str(out_pdbqt),
                    box=box,
                    out_pdbqt=str(out_pdbqt),
                    log_path=str(log_path),
                    cpu=args.cpu,
                    mode="score_only",
                    scoring=args.scoring,
                )
                if best_score is None or (sc is not None and sc < best_score):
                    best_score = sc
                    best_pose_text = out_pdbqt.read_text()
            if best_pose_text is not None:
                out_pdbqt.write_text(best_pose_text)
            score = best_score
        # Compute COM distance of final pose to disc center
        final_coords = coords_from_pdbqt_single(out_pdbqt)
        com_dist = None
        if final_coords:
            c = com(final_coords)
            com_dist = float(np.linalg.norm(c - np.array([d.x, d.y, d.z])))

        rows.append({
            "index": d.index,
            "distance": d.distance,
            "radius": d.radius,
            "x": d.x,
            "y": d.y,
            "z": d.z,
            "vina_score": score,
            "com_dist": com_dist,
        })

    pd.DataFrame(rows).to_csv(out / "vina_profile.csv", index=False)
    print(f"Wrote {out/'vina_profile.csv'}")


if __name__ == "__main__":
    main()
