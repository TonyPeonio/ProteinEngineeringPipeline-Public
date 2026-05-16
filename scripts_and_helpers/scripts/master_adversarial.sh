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

echo "INITIATING ARMS RACE PIPELINE AT: $PROJECT_ROOT"

# ==========================================
# MASTER CONFIGURATION
# ==========================================
# PDB Fetching Config
PDB_ID="1ycr"
DESIGNABLE_CHAINS="B"

# RFdiffusion Config
BINDER_LENGTH="15-25"
NUM_DESIGNS=1
AA_RANGE="A25-109"
BACKUP_AA_RANGE="A1-85"
HOTSPOTS="" 

# ProteinMPNN Config
SEQ_PER_BACKBONE=1
MPNN_TEMP=0.1
WILDCARDS=0

# ColabFold / Filtering Config
SERINE_PAINTING=0
CORE_PROTECTION=1
NUM_RECYCLES=12
RECYCLE_TOLERANCE=0.5
NUM_MODELS=1

# Generations to run
NUM_GENERATIONS=100

# Mutations Config
NUM_MUTANTS=1
NUM_MUTATIONS=3

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

# Setup ColabFold Environment Variables once
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform
export TF_FORCE_UNIFIED_MEMORY=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

# ==========================================
# PIPELINE EXECUTION - STEP 1 (RUNS ONCE)
# ==========================================

echo "--- STEP 1: PDB Fetching (Pipeline Preparation) ---"
cd "$PROJECT_ROOT/scripts_and_helpers/scripts"
conda activate env_rfdiffusion
python 1_pdb_fetcher.py "$PDB_ID"

# Set the initial target to the wildtype PDB for Generation 1
TARGET_PDB="$PROJECT_ROOT/scripts_and_helpers/pdb/${PDB_ID}_clean.pdb"
CSV_PATH="$PROJECT_ROOT/results/arms_race_data.csv"


if [ -f "$CSV_PATH" ]; then
    COMPLETED_GENS=$(($(wc -l < "$CSV_PATH") - 1))
    START_GEN=$((COMPLETED_GENS + 1))

    if ["$COMPLETED_GENS$" -eq 0 ]; then
        TARGET_PDB="$PROJECT_ROOT/scripts_and_helpers/pdb/${PDB_ID}_clean.pdb"
    else
        TARGET_PDB="$PROJECT_ROOT/results/generation_${COMPLETED_GENS}/mutant_target_gen${COMPLETED_GENS+1}.pdb"
        AA_RANGE="$BACKUP_AA_RANGE"
    fi
    echo "RESUMING PIPELINE: Starting at Generation $START_GEN"
    echo "Current Target: $TARGET_PDB"

else
    START_GEN=1
    TARGET_PDB="$PROJECT_ROOT/scripts_and_helpers/pdb/${PDB_ID}_clean.pdb"
    echo "STARTING NEW PIPELINE: Generation 1"
fi

# ==========================================
# ADVERSARIAL LOOP BEGINS
# ==========================================

for (( GEN=$START_GEN; GEN<=NUM_GENERATIONS; GEN++ )); do
    echo ""
    echo "=========================================================="
    echo "                 INITIATING GENERATION $GEN               "
    echo "                 Target: $(basename "$TARGET_PDB")        "
    echo "=========================================================="
    echo ""

    # Ensure we are in the scripts directory
    cd "$PROJECT_ROOT/scripts_and_helpers/scripts"

    # --- PHASE 1: THE PHARMACOLOGIST ---

    echo "==================================="
    echo "--- STEP 2: RFdiffusion (Complex Shape Generation) ---"
    echo "==================================="
    conda activate env_rfdiffusion
    bash 2_run_rfdiffusion.sh "$TARGET_PDB" "$AA_RANGE" "$BINDER_LENGTH" "$NUM_DESIGNS" "$HOTSPOTS"

    AA_RANGE="$BACKUP_AA_RANGE"

    echo "==================================="
    echo "--- STEP 3: ProteinMPNN (Complex Sequence Generation) ---"
    echo "==================================="
    conda activate env_mpnn
    bash 3_run_mpnn.sh "$SEQ_PER_BACKBONE" "$MPNN_TEMP" "$DESIGNABLE_CHAINS"
    conda activate env_rosetta
    python 3-2_mpnn_sequence_filtering.py "$WILDCARDS" "$TARGET_PDB" "$SERINE_PAINTING" "$CORE_PROTECTION"

    echo "==================================="
    echo "--- STEP 4: ColabFold (Drug Validation) ---"
    echo "==================================="
    conda activate "$COLAB_DIR"
    bash 4_run_localcolabfold.sh "$NUM_RECYCLES" "$RECYCLE_TOLERANCE" "$NUM_MODELS"
    
    echo "==================================="
    echo "--- STEP 5: Select Winning Drug ---"
    echo "==================================="
    conda activate env_rosetta
    python 5_evaluate_drug.py --gen "${GEN}"
    
    echo "==================================="
    echo "--- STEP 6: Mutate the Target ---"
    echo "==================================="
    conda activate env_rosetta
    
    LEAD_DRUG_PDB="$OUTPUT_DIR/evaluate_drug_results/lead_drug.pdb"
    
    # We now pass BOTH the true cancer target AND the new lead drug to the mutator!
    python 6_cancer_mutator.py \
        --target_mdm2 "$TARGET_PDB" \
        --lead_drug "$LEAD_DRUG_PDB" \
        --out_dir "$OUTPUT_DIR/phase2_fastas" \
        --num_mutants "$NUM_MUTANTS" \
        --num_mutations "$NUM_MUTATIONS"
    
    echo "==================================="
    echo "--- STEP 6b: ColabFold (Cancer Evaluation) ---"
    echo "==================================="
    conda activate "$COLAB_DIR"
    colabfold_batch "$OUTPUT_DIR/phase2_fastas" "$OUTPUT_DIR/colabfold_phase2_results" \
        --num-recycle "$NUM_RECYCLES" \
        --recycle-early-stop-tolerance "$RECYCLE_TOLERANCE" \
        --num-models "$NUM_MODELS"

    echo "==================================="
    echo "--- STEP 7: Select Winning Evasion Mutation ---"
    echo "==================================="
    conda activate env_rosetta
    python 7_evaluate_cancer.py --gen "$GEN"
    
    # --- PHASE 3: CLOSING THE LOOP & SAVING DATA ---
    echo "--- Archiving Generation $GEN Results ---"
    GEN_DIR="$PROJECT_ROOT/results/generation_$GEN"
    mkdir -p "$GEN_DIR"
    
    # Update these copy paths to pull from the new centralized outputs folder!
    cp "$OUTPUT_DIR/evaluate_drug_results/lead_drug.pdb" "$GEN_DIR/lead_drug_gen${GEN}.pdb"
    cp "$OUTPUT_DIR/evaluate_cancer_results/next_gen_target.pdb" "$GEN_DIR/mutant_target_gen${GEN}.pdb"

    # Swap the target for the next generation!
    TARGET_PDB="$OUTPUT_DIR/evaluate_cancer_results/next_gen_target.pdb"
    
    # --- PHASE 4: CLEANUP ---
    echo "--- Cleaning Staging Directories for Next Generation ---"
    rm -rf "$OUTPUT_DIR/colabfold_multimer_results/"*
    rm -rf "$OUTPUT_DIR/colabfold_phase2_results/"*
    rm -rf "$OUTPUT_DIR/phase2_fastas/"*
    rm -rf "$OUTPUT_DIR/rfdiffusion/"*
    rm -rf "$OUTPUT_DIR/mpnn/"* || true

    echo "GENERATION $GEN COMPLETE!"
done

echo "==================================="
echo "--- Creating Final Graph ---"
echo "==================================="
python 8_graph_results.py

echo ""
echo "=========================================================="
echo "         5-GENERATION ARMS RACE FULLY COMPLETE!           "
echo "         Check the '$PROJECT_ROOT/results' folder.        "
echo "=========================================================="