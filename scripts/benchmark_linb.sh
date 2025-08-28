#!/usr/bin/env bash
set -euo pipefail

# This script expects the following files/directories to exist:
# - data/dbe.pdb (or .mol2)
# - data/linb_wt/receptor.pdb
# - data/linb_wt/caver_out/  (CAVER outputs for the chosen tunnel)
# - data/linb_wt/caverdock_profile.csv (CaverDock upper bound profile)
# - data/linb_l177w/receptor.pdb
# - data/linb_l177w/caver_out/
# - data/linb_l177w/caverdock_profile.csv

MGL_PYTHONSH="/home/mestevez/Programs/Autodock4/mgltools_x86_64Linux2_1.5.7/bin/pythonsh"
PREP_RECEPTOR="/home/mestevez/Programs/Autodock4/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
PREP_LIGAND="/home/mestevez/Programs/Autodock4/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
VINA_BIN="/home/mestevez/Programs/Autodock4/vina_1.2.7_linux_x86_64"

echo "Running LinB WT with DBE..."
caverdock-lb run \
  --receptor data/linb_wt/receptor.pdb \
  --ligand   data/dbe.pdb \
  --caver-output data/linb_wt/caver_out \
  --out-dir runs/linb_wt \
  --vina-bin "$VINA_BIN" \
  --prepare-receptor-script "$PREP_RECEPTOR" \
  --prepare-ligand-script   "$PREP_LIGAND" \
  --mgltools-pythonsh "$MGL_PYTHONSH" \
  --caverdock-profile data/linb_wt/caverdock_profile.csv

echo "Running LinB L177W with DBE..."
caverdock-lb run \
  --receptor data/linb_l177w/receptor.pdb \
  --ligand   data/dbe.pdb \
  --caver-output data/linb_l177w/caver_out \
  --out-dir runs/linb_l177w \
  --vina-bin "$VINA_BIN" \
  --prepare-receptor-script "$PREP_RECEPTOR" \
  --prepare-ligand-script   "$PREP_LIGAND" \
  --mgltools-pythonsh "$MGL_PYTHONSH" \
  --caverdock-profile data/linb_l177w/caverdock_profile.csv

echo "Aggregating plots..."
python3 scripts/plot_linb_comparison.py \
  --wt runs/linb_wt/vina_profile.csv \
  --wt-cd data/linb_wt/caverdock_profile.csv \
  --mut runs/linb_l177w/vina_profile.csv \
  --mut-cd data/linb_l177w/caverdock_profile.csv \
  --out assets/plots/linb_wt_l177w_compare.png

echo "Done. See runs/ and assets/plots/."
