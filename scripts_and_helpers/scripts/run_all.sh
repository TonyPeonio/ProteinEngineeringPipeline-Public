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
PDB_ID="1ycr"
DESIGNABLE_CHAINS="B"

# RFdiffusion Config
NUM_DESIGNS=1
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
# PIPELINE EXECUTION (EVOLUTIONARY LOOP)
# ==========================================
MAX_ROUNDS=30
EVO_LOG="$OUTPUT_DIR/evolution_log.csv"

echo "--- STEP 1: PDB Fetching (Pipeline Preparation) ---"
conda activate env_rfdiffusion
python 1_pdb_fetcher.py "$PDB_ID"

# Initialize the "Common Ancestor"
CURRENT_BEST_PDB="$PROJECT_ROOT/scripts_and_helpers/pdb/${PDB_ID}_clean.pdb"

# Create the log file for your graphs
echo "Round,Score_pLDDT,Score_PAE,Mutations_from_Ancestor" > "$EVO_LOG"
for ROUND in $(seq 1 $MAX_ROUNDS); do
    echo "=========================================="
    echo "       STARTING EVOLUTION ROUND $ROUND    "
    echo "=========================================="
    
    # --- THE MAGIC IF STATEMENT ---
    if [ "$ROUND" -eq 1 ]; then
        # Round 1: Original PDB numbering
        ACTIVE_TARGET="A25-109"
        ACTIVE_BINDER="B17-29"
    else
        # Round 2+: ColabFold renumbers chains starting from 1
        ACTIVE_TARGET="A1-85"
        ACTIVE_BINDER="B1-13"
    fi
    
    # Create round-specific output folders
    ROUND_DIR="$OUTPUT_DIR/round_$ROUND"
    mkdir -p "$ROUND_DIR"

    echo "--- STEP 2: RFdiffusion (Mutating the Backbone) ---"
    export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
    conda activate env_rfdiffusion
    
    # Pass ACTIVE_BINDER instead of a fixed length
    bash 2_run_rfdiffusion.sh "$CURRENT_BEST_PDB" "$ACTIVE_TARGET" "$ACTIVE_BINDER" "$NUM_DESIGNS" "$HOTSPOTS" "$ROUND_DIR"

    echo "--- STEP 3: ProteinMPNN (Translating to Sequence) ---"
    conda activate env_mpnn
    bash 3_run_mpnn.sh "$SEQ_PER_BACKBONE" "$MPNN_TEMP" "$DESIGNABLE_CHAINS" "$ROUND_DIR"
    conda activate env_rosetta
    python 3-2_mpnn_sequence_filtering.py "$WILDCARDS" "$CURRENT_BEST_PDB" "$SERINE_PAINTING" "$CORE_PROTECTION" "$ROUND_DIR"

    echo "--- STEP 4: ColabFold (The Oracle's Judgment) ---"
    export XLA_PYTHON_CLIENT_PREALLOCATE=false
    export XLA_PYTHON_CLIENT_ALLOCATOR=platform
    export TF_FORCE_UNIFIED_MEMORY=1
    conda activate "$COLAB_DIR"
    bash 4_run_localcolabfold.sh "$NUM_RECYCLES" "$RECYCLE_TOLERANCE" "$NUM_MODELS" "$ROUND_DIR"

    echo "--- STEP 5: Selection & Feedback (Closing the Loop) ---"
    conda activate base
    # This python script reads the ColabFold JSON, decides if the new design is better, 
    # updates the log for your graphs, and outputs the path to the PDB for the next round.
    CURRENT_BEST_PDB=$(python 5_evaluate_and_select.py --round $ROUND --round_dir "$ROUND_DIR" --log "$EVO_LOG" --current_best "$CURRENT_BEST_PDB")
    
    echo "Round $ROUND Complete. Current Best Seed: $CURRENT_BEST_PDB"
done

echo "EVOLUTIONARY PIPELINE COMPLETE!"