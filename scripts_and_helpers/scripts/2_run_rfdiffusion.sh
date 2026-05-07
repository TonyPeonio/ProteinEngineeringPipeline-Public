#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e 

# ==========================================
# ARGUMENT PARSING
# ==========================================
TARGET_PDB=$1
AA_RANGE=$2
BINDER_LENGTH=$3
TARGET_NUM_DESIGNS=$4
HOTSPOTS=$5

# ==========================================
# DYNAMIC PATH ANCHORING
# ==========================================
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/../../" && pwd)

RFDIFFUSION_DIR="$PROJECT_ROOT/RFdiffusion"
OUTPUT_DIR="$PROJECT_ROOT/outputs/rfdiffusion"
mkdir -p "$OUTPUT_DIR"

TARGET_BASENAME=$(basename "$TARGET_PDB" .pdb)
PREFIX="$OUTPUT_DIR/${TARGET_BASENAME}_binder"

# ==========================================
# AUTO-RESUME LOGIC
# ==========================================
EXISTING=$(ls ${PREFIX}_*.pdb 2>/dev/null | wc -l)

echo "Found $EXISTING existing designs."

REMAINING=$((TARGET_NUM_DESIGNS - EXISTING))

if [ "$REMAINING" -le 0 ]; then
    echo "Target of $TARGET_NUM_DESIGNS designs already reached. Nothing to do."
    exit 0
fi

echo "Generating the remaining $REMAINING designs targeting $TARGET_PDB..."

# Create a clean staging directory for the batch
STAGING_DIR="$OUTPUT_DIR/.staging"
mkdir -p "$STAGING_DIR"
STAGING_PREFIX="$STAGING_DIR/batch"

# ==========================================
# RUN INFERENCE
# ==========================================
python "$RFDIFFUSION_DIR/scripts/run_inference.py" \
    inference.output_prefix="$STAGING_PREFIX" \
    inference.model_directory_path="$RFDIFFUSION_DIR/models" \
    inference.input_pdb="$TARGET_PDB" \
    contigmap.contigs=["$AA_RANGE/0 $BINDER_LENGTH"] \
    inference.num_designs=$REMAINING \
    ppi.hotspot_res=["$HOTSPOTS"] \
    hydra.run.dir="$STAGING_DIR" \
    hydra.output_subdir=null

echo "Generation complete. Integrating files into main directory..."

# ==========================================
# FILE RENAMING AND CLEANUP
# ==========================================
for (( i=0; i<REMAINING; i++ )); do
    staged_pdb="${STAGING_PREFIX}_${i}.pdb"
    staged_trb="${STAGING_PREFIX}_${i}.trb"

    if [ -e "$staged_pdb" ]; then
        mv "$staged_pdb" "${PREFIX}_${EXISTING}.pdb"
        
        if [ -e "$staged_trb" ]; then
            mv "$staged_trb" "${PREFIX}_${EXISTING}.trb"
        fi
        
        EXISTING=$((EXISTING + 1))
    fi
done

# Clean up the staging directory
rm -rf "$STAGING_DIR"

if [ "$EXISTING" -ge "$TARGET_NUM_DESIGNS" ]; then
    echo "Success: All $TARGET_NUM_DESIGNS designs are finished and perfectly numbered."
else
    echo "Warning: RFdiffusion dropped $((TARGET_NUM_DESIGNS - EXISTING)) designs. Re-run to generate the remaining ones."
fi