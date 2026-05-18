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

echo "INITIATING MULTI-AGENT ARMS RACE PIPELINE AT: $PROJECT_ROOT"

# ==========================================
# MASTER CONFIGURATION
# ==========================================
# PDB Fetching Config
PDB_ID="1ycr"
DESIGNABLE_CHAINS="B"

# RFdiffusion Config (Only used for Gen 0 baseline)
BINDER_LENGTH="15-25"
NUM_DESIGNS=1
AA_RANGE="A25-109"
HOTSPOTS="" 

# ProteinMPNN Config (Only used for Gen 0 baseline)
SEQ_PER_BACKBONE=1
MPNN_TEMP=0.1
WILDCARDS=0

# ColabFold Config
NUM_RECYCLES=12
RECYCLE_TOLERANCE=0.5
NUM_MODELS=1

# Agentic AI Config
NUM_GENERATIONS=20
NUM_MUTATIONS_CANCER=3
NUM_MUTATIONS_DRUG=3

# Ensure output directory exists using the dynamic path
mkdir -p "$OUTPUT_DIR"
mkdir -p "$PROJECT_ROOT/results"
echo "Configuration saved to $OUTPUT_DIR/config.txt"

# ==========================================
# ENVIRONMENT SETUP
# ==========================================
eval "$(conda shell.bash hook)"
export PYTHONWARNINGS="ignore"

# Setup ColabFold Environment Variables once
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform
export TF_FORCE_UNIFIED_MEMORY=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export XLA_PYTHON_CLIENT_FRACTION=0.7

# ==========================================
# PIPELINE INITIALIZATION & CHECKPOINTING
# ==========================================
cd "$PROJECT_ROOT/scripts_and_helpers/scripts"
CSV_PATH="$PROJECT_ROOT/results/arms_race_data.csv"

if [ -f "$CSV_PATH" ]; then
    COMPLETED_GENS=$(($(wc -l < "$CSV_PATH") - 1))
    START_GEN=$((COMPLETED_GENS + 1))
    echo "RESUMING PIPELINE: Starting at Generation $START_GEN"

    # --- THE CRITICAL RESUMABILITY FIX ---
    # If we are resuming from Gen 2 or later, we must restore the staging files
    # from the last fully completed generation's archive.
    if [ "$COMPLETED_GENS" -ge 1 ]; then
        echo "Restoring state from Generation $COMPLETED_GENS..."
        cp "$PROJECT_ROOT/results/generation_$COMPLETED_GENS/cancer_gen${COMPLETED_GENS}.fasta" "$OUTPUT_DIR/current_cancer.fasta"
        cp "$PROJECT_ROOT/results/generation_$COMPLETED_GENS/drug_gen${COMPLETED_GENS}.fasta" "$OUTPUT_DIR/current_drug.fasta"
    fi

else
    START_GEN=1
    echo "STARTING NEW PIPELINE: No history found. Establishing Generation 0 Baseline."
    
    echo "--- STEP 1: PDB Fetching ---"
    conda activate env_rfdiffusion
    python 1_pdb_fetcher.py "$PDB_ID"
    TARGET_PDB="$PROJECT_ROOT/scripts_and_helpers/pdb/${PDB_ID}_clean.pdb"

    echo "--- STEP 2: RFdiffusion (Baseline Complex Shape) ---"
    bash 2_run_rfdiffusion.sh "$TARGET_PDB" "$AA_RANGE" "$BINDER_LENGTH" "$NUM_DESIGNS" "$HOTSPOTS"

    echo "--- STEP 3: ProteinMPNN (Baseline Sequence) ---"
    conda activate env_mpnn
    bash 3_run_mpnn.sh "$SEQ_PER_BACKBONE" "$MPNN_TEMP" "$DESIGNABLE_CHAINS"
    conda activate env_rosetta
    python 3-2_mpnn_sequence_filtering.py "$WILDCARDS" "$TARGET_PDB" "0" "1"

    echo "--- STEP 4: ColabFold (Baseline Validation) ---"
    conda activate "$COLAB_DIR"
    bash 4_run_localcolabfold.sh "$NUM_RECYCLES" "$RECYCLE_TOLERANCE" "$NUM_MODELS"
    
    echo "--- STEP 5: Extract Baseline FASTAs ---"
    conda activate env_rosetta
    python 5_evaluate_drug.py --gen 0
fi

# ==========================================
# MULTI-AGENT ADVERSARIAL LOOP BEGINS
# ==========================================

for (( GEN=$START_GEN; GEN<=NUM_GENERATIONS; GEN++ )); do
    echo ""
    echo "=========================================================="
    echo "                 INITIATING GENERATION $GEN               "
    echo "=========================================================="
    echo ""
    cd "$PROJECT_ROOT/scripts_and_helpers/scripts"

    echo "==================================="
    echo "--- STEP 6: The Cancer AI's Turn ---"
    echo "==================================="
    conda activate env_rosetta
    python 6_cancer_mutator.py \
        --current_cancer "$OUTPUT_DIR/current_cancer.fasta" \
        --history_csv "$CSV_PATH" \
        --out_fasta "$OUTPUT_DIR/cancer_gen${GEN}.fasta" \
        --num_mutations "$NUM_MUTATIONS_CANCER"

    sleep 60
    
    echo "==================================="
    echo "--- STEP 7: The Pharmacologist AI's Turn ---"
    echo "==================================="
    # The new Agent!
    python 7_pharma_mutator.py \
        --current_drug "$OUTPUT_DIR/current_drug.fasta" \
        --history_csv "$CSV_PATH" \
        --out_fasta "$OUTPUT_DIR/drug_gen${GEN}.fasta" \
        --num_mutations "$NUM_MUTATIONS_DRUG"

    echo "==================================="
    echo "--- STEP 8: Stitch FASTAs for Arbiter ---"
    echo "==================================="
    # Create a staging folder for this generation's combat inputs
    mkdir -p "$OUTPUT_DIR/colabfold_inputs_gen${GEN}"
    
    python 8_fasta_stitcher.py \
        --cancer_fasta "$OUTPUT_DIR/cancer_gen${GEN}.fasta" \
        --drug_fasta "$OUTPUT_DIR/drug_gen${GEN}.fasta" \
        --out_dir "$OUTPUT_DIR/colabfold_inputs_gen${GEN}"

    echo "==================================="
    echo "--- STEP 9: Combat (ColabFold Validation) ---"
    echo "==================================="
    conda activate "$COLAB_DIR"
    # We pass the entire directory, so ColabFold folds both the Drug and p53 complexes!
    colabfold_batch "$OUTPUT_DIR/colabfold_inputs_gen${GEN}" "$OUTPUT_DIR/colabfold_combat_results" \
        --num-recycle "$NUM_RECYCLES" \
        --recycle-early-stop-tolerance "$RECYCLE_TOLERANCE" \
        --num-models "$NUM_MODELS"

    echo "==================================="
    echo "--- STEP 10: Scorekeeper ---"
    echo "==================================="
    conda activate env_rosetta
    
    python 10_scorekeeper.py \
        --gen "$GEN" \
        --results_dir "$OUTPUT_DIR/colabfold_combat_results" \
        --history_csv "$CSV_PATH" \
        --wt_cancer_fasta "$PROJECT_ROOT/results/generation_0/cancer_gen0.fasta"

    # --- CLOSING THE LOOP ---
    echo "--- Archiving Generation $GEN Results ---"
    GEN_DIR="$PROJECT_ROOT/results/generation_$GEN"
    mkdir -p "$GEN_DIR"
    
    # Save the sequences
    cp "$OUTPUT_DIR/cancer_gen${GEN}.fasta" "$GEN_DIR/"
    cp "$OUTPUT_DIR/drug_gen${GEN}.fasta" "$GEN_DIR/"
    
    # Set the current state for the next loop
    cp "$OUTPUT_DIR/cancer_gen${GEN}.fasta" "$OUTPUT_DIR/current_cancer.fasta"
    cp "$OUTPUT_DIR/drug_gen${GEN}.fasta" "$OUTPUT_DIR/current_drug.fasta"

    # --- CLEANUP ---
    echo "--- Cleaning Staging Directories for Next Generation ---"
    rm -rf "$OUTPUT_DIR/colabfold_combat_results/"*
    
    echo "==================================="
    echo "--- Creating Final Graph ---"
    echo "==================================="
    python 11_graph_results.py

    echo "GENERATION $GEN COMPLETE!"
done

echo ""
echo "=========================================================="
echo "         $NUM_GENERATIONS-GENERATION ARMS RACE COMPLETE!  "
echo "         Check the '$PROJECT_ROOT/results' folder.        "
echo "=========================================================="