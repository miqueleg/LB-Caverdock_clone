### This project is just me experimenting and testing stuff. Do not use this, Use CaverDock for this type of calculations

Lower-Bound (Vina-based)
========================

Library and CLI to reproduce CaverDock-like lower-bound energy profiles using
AutoDock Vina + AutoDockTools (MGLTools) + CAVER 3. The pipeline docks a ligand
independently at each CAVER disc (centerline) with configurable per-sphere
charge parametrization hooks.

Key features
------------
- Compute tunnels with CAVER 3 or parse existing outputs
- Extract centerline discs (x, y, z, radius, distance)
- Prepare receptor/ligand PDBQT via MGLTools (Gasteiger), with pluggable charge providers
- Run Vina per disc to estimate a lower-bound trajectory
- Compile, plot, and compare the energy profile against CaverDock upper-bound outputs

Requirements
------------
- Python 3.9+
- External tools installed and on PATH (or configured):
  - CAVER 3 CLI (e.g., `caver` or `java -jar Caver.jar`)
  - AutoDock Vina (e.g., `vina`)
  - MGLTools prepare scripts (e.g., `prepare_receptor4.py`, `prepare_ligand4.py`)

Install
-------
```
pip install -e .[dev]
```

Quick start
-----------
```
caverdock-lb run \
  --receptor path/to/receptor.pdb \
  --ligand   path/to/ligand.pdb \
  --caver-output path/to/caver_out \
  --out-dir  runs/exp1 \
  --vina-bin vina \
  --prepare-receptor-script /path/to/prepare_receptor4.py \
  --prepare-ligand-script   /path/to/prepare_ligand4.py
```

You can pass `--caverdock-profile` to overlay/compare a known CaverDock upper-bound profile (CSV).

Extending charges per sphere
----------------------------
Implement `ChargeProvider` to compute new charges per disc (e.g., g-xtb, AMOEBA,
FDBβ electric field adjustments). Plug it via `--charge-provider` (Python import
path) or integrate in your scripts by constructing the provider directly.

Repository layout
-----------------
See `docs/README.md` and `configs/example.yaml` for configuration and more details.
