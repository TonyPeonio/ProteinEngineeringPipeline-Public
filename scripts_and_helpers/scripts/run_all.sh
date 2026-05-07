#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e 

# ==========================================
# DYNAMIC PATH ANCHORING
# ==========================================
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/../../" && pwd)

# Define key paths relative to root
OUTPUT_DIR="$PROJECT_ROOT/outputs"
COLAB_DIR="$PROJECT_ROOT/localcolabfold/colabfold-conda" 

echo "INITIATING PIPELINE AT: $PROJECT_ROOT"

# ==========================================
# MASTER CONFIGURATION
# ==========================================
# PDB Fetching Config
PDB_ID="7ut8"
DESIGNABLE_CHAINS="B"

# RFdiffusion Config
BINDER_LENGTH="55-65"
NUM_DESIGNS=1
AA_RANGE="A43-123/0 A126-318"
HOTSPOTS="" 

# ProteinMPNN Config
SEQ_PER_BACKBONE=1
MPNN_TEMP=0.1
WILDCARDS=0

# ColabFold / Filtering Config
SERINE_PAINTING=1
CORE_PROTECTION=0
NUM_RECYCLES=12
RECYCLE_TOLERANCE=0.5
NUM_MODELS=5

# Ensure output directory exists using the dynamic path
mkdir -p "$OUTPUT_DIR"
echo "Configuration saved to $OUTPUT_DIR/config.txt"

cat <<EOF > "$OUTPUT_DIR/config.txt"
# Pipeline Configuration
PDB_ID="$PDB_ID"
PROJECT_ROOT="$PROJECT_ROOT"
BINDER_LENGTH="$BINDER_LENGTH"
NUM_DESIGNS="$NUM_DESIGNS"
AA_RANGE="$AA_RANGE"
HOTSPOTS="$HOTSPOTS"
SEQ_PER_BACKBONE="$SEQ_PER_BACKBONE"
MPNN_TEMP="$MPNN_TEMP"
WILDCARDS="$WILDCARDS"
NUM_RECYCLES="$NUM_RECYCLES"
SERINE_PAINTING="$SERINE_PAINTING"
EOF

# ==========================================
# ENVIRONMENT SETUP
# ==========================================
eval "$(conda shell.bash hook)"
export PYTHONWARNINGS="ignore"

# ==========================================
# PIPELINE EXECUTION
# ==========================================

echo "--- STEP 1: PDB Fetching (Pipeline Preparation) ---"
conda activate env_rfdiffusion
python 1_pdb_fetcher.py "$PDB_ID"

TARGET_PDB="$PROJECT_ROOT/scripts_and_helpers/pdb/${PDB_ID}_clean.pdb"

echo "--- STEP 2: RFdiffusion (Complex Shape Generation) ---"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
conda activate env_rfdiffusion
bash 2_run_rfdiffusion.sh "$TARGET_PDB" "$AA_RANGE" "$BINDER_LENGTH" "$NUM_DESIGNS" "$HOTSPOTS"

echo "--- STEP 3: ProteinMPNN (Complex Sequence Generation) ---"
conda activate env_mpnn
bash 3_run_mpnn.sh "$SEQ_PER_BACKBONE" "$MPNN_TEMP" "$DESIGNABLE_CHAINS"
conda activate env_rosetta
python 3-2_mpnn_sequence_filtering.py "$WILDCARDS" "$TARGET_PDB" "$SERINE_PAINTING" "$CORE_PROTECTION"

echo "--- STEP 4: ColabFold (Multimer Validation) ---"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform
export TF_FORCE_UNIFIED_MEMORY=1

conda activate "$COLAB_DIR"
bash 4_run_localcolabfold.sh "$NUM_RECYCLES" "$RECYCLE_TOLERANCE" "$NUM_MODELS"

echo "PIPELINE COMPLETE!"