#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


def read_csv_index_energy(path: str, idx_col: str, energy_col: str):
    data = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            try:
                idx = int(float(row[idx_col]))
                e = float(row[energy_col])
            except (KeyError, ValueError):
                continue
            data[idx] = e
    return data


def read_dsd_distances(path: str):
    xs, ys, zs = [], [], []
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 7:
                continue
            x, y, z = map(float, parts[:3])
            xs.append(x); ys.append(y); zs.append(z)
    dist = [0.0]
    import math
    for i in range(1, len(xs)):
        dx = xs[i] - xs[i-1]
        dy = ys[i] - ys[i-1]
        dz = zs[i] - zs[i-1]
        dist.append(dist[-1] + math.sqrt(dx*dx + dy*dy + dz*dz))
    # indices are 1-based
    return {i+1: d for i, d in enumerate(dist)}


def scale(value, a0, a1, b0, b1):
    if a1 == a0:
        return (b0 + b1) / 2
    t = (value - a0) / (a1 - a0)
    return b0 + t * (b1 - b0)


def polyline(points, stroke, width=2):
    d = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    return f'<polyline fill="none" stroke="{stroke}" stroke-width="{width}" points="{d}" />\n'


def main():
    import argparse

    ap = argparse.ArgumentParser(description="Plot Vina LB vs CaverDock LB as SVG (distance x-axis)")
    ap.add_argument("--vina-csv", required=True)
    ap.add_argument("--cd-csv", required=True)
    ap.add_argument("--dsd", required=True, help="CaverDock tunnel.dsd to define the x-axis distances")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # Read energies per index
    vina_e = read_csv_index_energy(args.vina_csv, "index", "vina_score")
    cd_e = read_csv_index_energy(args.cd_csv, "index", "cd_lb_energy")
    # Read distances from DSD and align
    idx_to_dist = read_dsd_distances(args.dsd)

    # Build aligned series by index
    common_idx = sorted(set(vina_e.keys()) & set(cd_e.keys()) & set(idx_to_dist.keys()))
    if not common_idx:
        raise SystemExit("No overlapping indices between inputs")
    vx = [idx_to_dist[i] for i in common_idx]
    vy = [vina_e[i] for i in common_idx]
    cx = [idx_to_dist[i] for i in common_idx]
    cy = [cd_e[i] for i in common_idx]

    xmin = min(vx)
    xmax = max(vx)
    ymin = min(min(vy), min(cy))
    ymax = max(max(vy), max(cy))

    width, height = 900, 420
    left, right, top, bottom = 70, 20, 20, 50
    plot_w = width - left - right
    plot_h = height - top - bottom

    def to_px(data_x, data_y):
        x = scale(data_x, xmin, xmax, left, left + plot_w)
        y = scale(data_y, ymin, ymax, top + plot_h, top)
        return x, y

    v_pts = [to_px(x, y) for x, y in zip(vx, vy)]
    c_pts = [to_px(x, y) for x, y in zip(cx, cy)]

    # Build SVG
    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">\n')
    svg.append('<rect width="100%" height="100%" fill="#ffffff"/>\n')
    # Axes
    svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#000"/>\n')
    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#000"/>\n')
    # Labels
    svg.append(f'<text x="{width/2}" y="{height-10}" text-anchor="middle" font-family="sans-serif" font-size="14">Distance along tunnel (Å)</text>\n')
    svg.append(f'<text x="15" y="{top+plot_h/2}" transform="rotate(-90 15,{top+plot_h/2})" text-anchor="middle" font-family="sans-serif" font-size="14">Energy (kcal/mol)</text>\n')
    # Title
    svg.append(f'<text x="{width/2}" y="{15}" text-anchor="middle" font-family="sans-serif" font-size="16">Lower-bound: Vina vs CaverDock</text>\n')
    # Curves
    svg.append(polyline(v_pts, "#1f77b4", 2))
    svg.append(polyline(c_pts, "#2ca02c", 2))
    # Legend
    lx, ly = left + 20, top + 20
    svg.append(f'<rect x="{lx-10}" y="{ly-15}" width="220" height="40" fill="#ffffff" stroke="#ccc"/>\n')
    svg.append(f'<line x1="{lx}" y1="{ly}" x2="{lx+30}" y2="{ly}" stroke="#1f77b4" stroke-width="3"/>\n')
    svg.append(f'<text x="{lx+40}" y="{ly+5}" font-family="sans-serif" font-size="12">Vina per-disc (LB)</text>\n')
    svg.append(f'<line x1="{lx}" y1="{ly+18}" x2="{lx+30}" y2="{ly+18}" stroke="#2ca02c" stroke-width="3"/>\n')
    svg.append(f'<text x="{lx+40}" y="{ly+23}" font-family="sans-serif" font-size="12">CaverDock LB</text>\n')
    svg.append('</svg>\n')

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text("".join(svg))
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
