from __future__ import annotations

"""
Lightweight helpers for working with CAVER 3 outputs.

This module documents expected inputs but does not execute CAVER itself.

Expected outputs after running CAVER 3 on a receptor:
- A set of tunnel centerlines with discs: a text/CSV file with columns
  x,y,z,radius[,nx,ny,nz] along the path from protein interior to surface.

Typical command (example; adjust paths/options):
  java -jar Caver.jar \
    -p receptor.pdb \
    -c config.txt \
    -o caver_output

Then export a specific tunnel’s centerline to a CSV for this tool.
"""

from pathlib import Path
from typing import List

from .discs import Disc


def save_discs_as_csv(discs: List[Disc], out_csv: str | Path) -> None:
    from .discs import write_boxes_csv, generate_disc_boxes

    # Reuse writer to keep a consistent CSV schema (includes normals)
    boxes = generate_disc_boxes(discs, box_margin=0.0)
    write_boxes_csv(boxes, out_csv)

