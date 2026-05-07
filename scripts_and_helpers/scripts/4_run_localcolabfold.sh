#!/bin/bash
set -e

# Grab the arguments, with safe fallbacks
NUM_RECYCLES=${1:-12}
TOLERANCE=${2:-0.5}
NUM_MODELS=${3:-5}

INPUT_DIR="../../outputs/colabfold_multimer_inputs"
OUTPUT_DIR="../../outputs/colabfold_multimer_results"

# Make the output directories if they don't exist
mkdir -p "$OUTPUT_DIR"

echo "--- Running ColabFold (AlphaFold2 Multimer) ---"
echo "Recycles: $NUM_RECYCLES | Early Stop Tolerance: $TOLERANCE"

# Run ColabFold
colabfold_batch "$INPUT_DIR" "$OUTPUT_DIR" \
    --model-type alphafold2_multimer_v3 \
    --num-recycle "$NUM_RECYCLES" \
    --num-models "$NUM_MODELS" \
    --recycle-early-stop-tolerance "$TOLERANCE"

echo "ColabFold batch complete! Results saved to $OUTPUT_DIR"