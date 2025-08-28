#!/usr/bin/env bash
set -euo pipefail

receptor=${1:-data/receptor.pdb}
ligand=${2:-data/ligand.pdb}
caver_out=${3:-data/caver_out}
out_dir=${4:-runs/example}

caverdock-lb run \
  --receptor "$receptor" \
  --ligand "$ligand" \
  --caver-output "$caver_out" \
  --out-dir "$out_dir" \
  --vina-bin vina \
  --prepare-receptor-script prepare_receptor4.py \
  --prepare-ligand-script prepare_ligand4.py

echo "Done. See $out_dir"
