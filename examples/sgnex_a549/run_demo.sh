#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXAMPLE_DIR="$ROOT/examples/sgnex_a549"
OUTPUT_DIR="$EXAMPLE_DIR/output"
PYTHON_BIN="${PYTHON_BIN:-python3}"

if [ -e "$OUTPUT_DIR" ]; then
  echo "Refusing to overwrite $OUTPUT_DIR; remove it before rerunning." >&2
  exit 1
fi
mkdir -p "$OUTPUT_DIR"
export MPLCONFIGDIR="$OUTPUT_DIR/.matplotlib"

"$PYTHON_BIN" "$ROOT/code/annotation_free_adaptive_depth.py" \
  --ssc "$EXAMPLE_DIR/SGNex_A549_directRNA_replicate6_run1.doradov300.sup.mt.sorted_flnc.ssc" \
  --outdir "$OUTPUT_DIR" \
  --sample-id SGNex_A549_directRNA_replicate6_run1 \
  --chrom MT \
  --sequence-length 16569 \
  --run-mode data \
  --adaptive-depth-mode fixed \
  --adaptive-raw-sqrt-divisor 10.0 \
  --adaptive-cluster-fraction 0.015 \
  --adaptive-display-fraction 0.035 \
  --adaptive-representative-fraction 0.006 \
  --adaptive-low-cluster-fraction 0.0 \
  --adaptive-low-display-fraction 0.0 \
  --adaptive-low-representative-fraction 0.01 \
  --adaptive-low-total-support-max 10000 \
  --adaptive-medium-total-support-max 100000 \
  --density-min-raw-pair-support 25 \
  --density-seed-min-raw-pair-support 5 \
  --density-backfill-max-raw-pair-support 4 \
  --density-connect-start-bp 20 \
  --density-connect-end-bp 20 \
  --density-min-cluster-support 200 \
  --display-min-cluster-support 1000 \
  --display-min-top-raw-fraction 0.1 \
  --display-min-representative-raw-support 392 \
  --display-max-start-span-bp 0 \
  --display-max-end-span-bp 0 \
  --compact-delta-start-bp 100 \
  --compact-delta-end-bp 100 \
  --shared-endpoint-delta-bp 40 \
  --shadow-other-endpoint-delta-bp 400 \
  --containment-min 0.78 \
  --plot-annotation-reference-format auto
