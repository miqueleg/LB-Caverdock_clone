#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <receptor_dir> <out_dir> <config_file>" >&2
  exit 1
fi

RECEPTOR_DIR=$(realpath "$1")
OUT_DIR=$(realpath "$2")
CONFIG=$(realpath "$3")

JAVA_JAR="/home/mestevez/Programs/Caver3/caver/caver.jar"
CAVER_HOME="/home/mestevez/Programs/Caver3/caver"

mkdir -p "$OUT_DIR"
echo "Running CAVER: $RECEPTOR_DIR -> $OUT_DIR using $CONFIG"
java -jar "$JAVA_JAR" -home "$CAVER_HOME" -pdb "$RECEPTOR_DIR" -conf "$CONFIG" -out "$OUT_DIR"

echo "Done. Outputs in $OUT_DIR"
