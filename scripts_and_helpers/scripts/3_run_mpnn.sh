#!/bin/bash
# ProteinMPNN Complex Sequence Design (SMART RESUME)
set -e 

# ==========================================
# ARGUMENT PARSING & PATHS
# ==========================================
if [ -z "$3" ]; then
    echo "Error: DESIGNABLE_CHAINS argument is missing."
    echo "Usage: $0 <NUM_SEQS> <TEMP> <DESIGNABLE_CHAINS>"
    echo "Example: $0 5 0.1 A"
    exit 1
fi

NUM_SEQS=${1:-5}
TEMP=${2:-0.1}
DESIGNABLE_CHAINS=${3}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/../../" && pwd)

MPNN_DIR="$PROJECT_ROOT/ProteinMPNN"
INPUT_DIR="$PROJECT_ROOT/outputs/rfdiffusion" 
OUTPUT_DIR="$PROJECT_ROOT/outputs/mpnn_results"
SEQS_DIR="$OUTPUT_DIR/seqs"
STAGING_DIR="$OUTPUT_DIR/staging_pdbs"

mkdir -p "$OUTPUT_DIR" "$SEQS_DIR" "$STAGING_DIR"

# Clean staging directory safely before starting
rm -f "$STAGING_DIR"/*.pdb || true
rm -f "$STAGING_DIR"/*.jsonl || true

# ==========================================
# SMART RESUME (PROCESS ONLY THE DELTA)
# ==========================================
echo "Checking for unprocessed RFdiffusion outputs..."
NEW_MODELS=0

# Loop through all PDBs in the rfdiffusion folder
for pdb_file in "$INPUT_DIR"/*.pdb; do
    [ -e "$pdb_file" ] || continue 
    base_name=$(basename "$pdb_file" .pdb)
    
    # If MPNN hasn't made a .fa file for this backbone yet, stage it
    if [ ! -f "$SEQS_DIR/$base_name.fa" ]; then
        # Using cp instead of ln -s to avoid any broken symlink headaches
        cp "$pdb_file" "$STAGING_DIR/"
        NEW_MODELS=$((NEW_MODELS + 1))
    fi
done

MASTER_PARSED="$OUTPUT_DIR/parsed_pdbs.jsonl"

# If there is no new work, just secure the handoff and exit cleanly
if [ "$NEW_MODELS" -eq 0 ]; then
    echo "All models have already been processed by ProteinMPNN. Skipping inference."
    
    # Ensure the master JSONL exists for the filtering script even if we skipped!
    if [ ! -f "$MASTER_PARSED" ]; then
        echo "Rebuilding master parsed_pdbs.jsonl for the next step..."
        python "$MPNN_DIR/helper_scripts/parse_multiple_chains.py" \
            --input_path="$INPUT_DIR" \
            --output_path="$MASTER_PARSED" > /dev/null 2>&1
    fi
    exit 0
fi

echo "Found $NEW_MODELS new model(s). Processing ONLY the new data..."

# ==========================================
# PREPARE INPUTS (ONLY FOR NEW MODELS)
# ==========================================
STAGING_PARSED="$STAGING_DIR/new_parsed.jsonl"
STAGING_ASSIGNED="$STAGING_DIR/new_assigned.jsonl"

echo "Parsing new complexes..."
python "$MPNN_DIR/helper_scripts/parse_multiple_chains.py" \
    --input_path="$STAGING_DIR" \
    --output_path="$STAGING_PARSED" > /dev/null

echo "Assigning designable chains..."
python "$MPNN_DIR/helper_scripts/assign_fixed_chains.py" \
    --input_path="$STAGING_PARSED" \
    --output_path="$STAGING_ASSIGNED" \
    --chain_list "$DESIGNABLE_CHAINS" > /dev/null

# ==========================================
# RUN PROTEINMPNN
# ==========================================
echo "Running ProteinMPNN..."
python "$MPNN_DIR/protein_mpnn_run.py" \
    --jsonl_path "$STAGING_PARSED" \
    --chain_id_jsonl "$STAGING_ASSIGNED" \
    --out_folder "$OUTPUT_DIR" \
    --num_seq_per_target "$NUM_SEQS" \
    --sampling_temp "$TEMP" \
    --seed 37 \
    --batch_size 1

# ==========================================
# REBUILD MASTER DATASET
# ==========================================
# Instead of risky text appending, we just parse the whole input folder in < 1 second.
# This guarantees a perfect JSONL file for the Python script every single time.
echo "Rebuilding complete master parsed_pdbs.jsonl..."
python "$MPNN_DIR/helper_scripts/parse_multiple_chains.py" \
    --input_path="$INPUT_DIR" \
    --output_path="$MASTER_PARSED" > /dev/null 2>&1

echo "ProteinMPNN successfully processed $NEW_MODELS new models."