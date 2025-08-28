# LinB WT and L177W with DBE

This guide reproduces CaverDock-style lower-bound (Vina per-disc) vs. CaverDock upper-bound comparisons for LinB WT and the
L177W variant using DBE as ligand.

Inputs expected
---------------
Place the following files in the repository:

- `data/dbe.pdb` (or `.mol2`): DBE ligand coordinates
- `data/linb_wt/receptor.pdb`: LinB WT structure (protonated as needed)
- `data/linb_wt/caver_out/`: CAVER 3 output for the selected tunnel
- `data/linb_wt/caverdock_profile.csv`: CaverDock upper-bound profile (distance, energy)
- `data/linb_l177w/receptor.pdb`: LinB L177W structure
- `data/linb_l177w/caver_out/`: CAVER 3 output for the selected tunnel
- `data/linb_l177w/caverdock_profile.csv`: CaverDock upper-bound profile

Running the benchmark
---------------------
Ensure external tools are installed and on PATH:

- AutoDock Vina (`vina`)
- MGLTools prepare scripts (`prepare_receptor4.py`, `prepare_ligand4.py`)

Then run:

```
bash scripts/benchmark_linb.sh
```

Outputs
-------
- `runs/linb_wt/vina_profile.csv` and `runs/linb_l177w/vina_profile.csv`:
  per-disc Vina energies with disc positions and radii
- `assets/plots/linb_wt_l177w_compare.png`:
  side-by-side plots comparing Vina lower-bound vs CaverDock upper-bound profiles

Notes and tuning
----------------
- Box size = `2*radius + margin` (default margin 2.0 Å). Adjust via config/CLI.
- Increase `exhaustiveness` for more thorough sampling; increases runtime.
- For per-sphere charge recalculation experiments, implement a new `ChargeProvider`
  and plug it into the pipeline.
