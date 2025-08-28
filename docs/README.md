# Documentation

Overview
--------
This library reproduces CaverDock-style lower-bound energy profiles by running
independent AutoDock Vina dockings at each CAVER disc along a tunnel.

Design
------
- `caver.py`: Run/parse CAVER outputs to extract discs (centerline + radii)
- `charge_providers.py`: Strategy interface for per-disc PDBQT preparation
  - Default: `MGLToolsGasteigerChargeProvider` (prepare_*4.py)
  - Extensible: plug XTB/AMOEBA/E-field providers in future
- `vina.py`: Build box from disc and run Vina
- `pipeline.py`: Orchestrate per-disc docking, collect results
- `plotting.py`: Plot Vina lower-bound vs. CaverDock upper-bound profiles
- `cli.py`: End-to-end CLI (`caverdock-ub run ...`)

Assumptions & Notes
-------------------
- CAVER outputs vary; `load_discs` tries common file names (`profile.csv`,
  `profile.txt`, `centerline.pdb`, `spheres.pdb`). Adjust if needed.
- The Vina search box is a cube of size `2*radius + margin`. Tune `margin`.
- Receptor/ligand preparation is factored for replaceable per-disc charge models.

Comparing to CaverDock
----------------------
Provide `--caverdock-profile` to overlay a known CaverDock upper-bound energy profile. You may
need to specify/rename energy and distance columns; the code guesses common names.

Extending Charges
-----------------
Implement `ChargeProvider` with methods returning PDBQT paths per disc context.
This allows inserting XTB-derived charges, AMOEBA polarization, or field-based
adjustments without changing the pipeline.

Testing
-------
Unit tests mock subprocess to avoid external dependencies. Run `make test`.
