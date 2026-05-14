#!/usr/bin/env python
import os
import glob
import sys
import json
import math

# ==========================================
# ARGUMENT PARSING & CONFIGURATION
# ==========================================
if len(sys.argv) < 6:
    print("Usage: python mpnn_sequence_filtering.py <WILDCARDS> <TARGET_PDB_PATH> <PAINT_FRACTION> <CORE_PROTECT_FRACTION> <ROUND_DIR>")
    sys.exit(1)

WILDCARDS_NEEDED = int(sys.argv[1])
WILD_TYPE_PDB_PATH = os.path.abspath(sys.argv[2])
SERINE_PAINT_FRACTION = float(sys.argv[3])
CORE_PROTECT_FRACTION = float(sys.argv[4])
ROUND_DIR = os.path.abspath(sys.argv[5])

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "../../"))

# Use the ROUND_DIR directly
INPUT_DIR = os.path.join(ROUND_DIR, "mpnn_results")
MPNN_OUTPUT_DIR = os.path.join(ROUND_DIR, "mpnn_results", "seqs")
PARSED_JSONL_PATH = os.path.join(ROUND_DIR, "mpnn_results", "parsed_pdbs.jsonl")
OUTPUT_FASTA_DIR = os.path.join(ROUND_DIR, "colabfold_multimer_inputs")

print("Extracting target sequences and filtering ProteinMPNN results...")
os.makedirs(OUTPUT_FASTA_DIR, exist_ok=True)

# ==========================================
# PDB PARSING & 3D MATH LOGIC
# ==========================================
def parse_pdb_ca_atoms(pdb_path, target_chain=None, ignore_chain=None):
    """Returns a list of dicts with Alpha-Carbon coordinates for specified chains."""
    atoms = []
    if not os.path.exists(pdb_path):
        print(f"Error: Cannot find PDB file at {pdb_path}.")
        sys.exit(1)
        
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                chain_id = line[21]
                
                if atom_name == "CA":
                    if target_chain and chain_id != target_chain: continue
                    if ignore_chain and chain_id == ignore_chain: continue
                    
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    atoms.append({
                        "chain": chain_id,
                        "coords": (x, y, z)
                    })
    return atoms

def calculate_structural_core(target_atoms, core_radius=10.0, protect_fraction=0.20):
    """
    Calculates internal contacts to find the 'load-bearing' structural backbone.
    Returns the indices of the top X% most highly connected residues.
    """
    if protect_fraction <= 0.0:
        return set()

    neighbor_counts = {}
    
    # Count how many neighbors each residue has within its own chain
    for i, atom_i in enumerate(target_atoms):
        count = 0
        xi, yi, zi = atom_i["coords"]
        for j, atom_j in enumerate(target_atoms):
            if i == j: continue
            xj, yj, zj = atom_j["coords"]
            dist = math.sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
            if dist <= core_radius:
                count += 1
        neighbor_counts[i] = count
        
    # Sort residues from highest internal connectivity to lowest
    sorted_by_core = sorted(neighbor_counts.items(), key=lambda item: item[1], reverse=True)
    
    # Grab the top X% to protect
    num_to_protect = int(math.ceil(len(sorted_by_core) * protect_fraction))
    protected_indices = set([item[0] for item in sorted_by_core[:num_to_protect]])
    
    print(f"Structural Protection: Identified the top {num_to_protect} 'load-bearing' core residues.")
    return protected_indices

def calculate_ranked_proximity_indices(pdb_path, target_chain="A", paint_fraction=1.0, protect_fraction=0.20, interface_radius=8.0):
    """
    Finds interface residues to paint, while mathematically protecting the structural core.
    """
    target_atoms = parse_pdb_ca_atoms(pdb_path, target_chain=target_chain)
    context_atoms = parse_pdb_ca_atoms(pdb_path, ignore_chain=target_chain)
    
    if not context_atoms:
        print(f"No contextual chains found in {pdb_path}. Serine painting skipped.")
        return set()

    # 1. Figure out which residues are structurally vital
    protected_core_indices = calculate_structural_core(target_atoms, protect_fraction=protect_fraction)

    # 2. Figure out distance to context chains
    residue_distances = {}
    for i, t_atom in enumerate(target_atoms):
        tx, ty, tz = t_atom["coords"]
        min_dist = float('inf')
        for c_atom in context_atoms:
            cx, cy, cz = c_atom["coords"]
            dist = math.sqrt((tx-cx)**2 + (ty-cy)**2 + (tz-cz)**2)
            if dist < min_dist:
                min_dist = dist
        residue_distances[i] = min_dist

    # Filter to interface pool and sort closest first
    interface_pool = {k: v for k, v in residue_distances.items() if v <= interface_radius}
    if not interface_pool:
        return set()

    sorted_interface = sorted(interface_pool.items(), key=lambda item: item[1])
    num_to_paint = int(math.ceil(len(sorted_interface) * paint_fraction))
    raw_indices_to_paint = set([item[0] for item in sorted_interface[:num_to_paint]])
    
    # 3. Filter out the structural core
    final_indices_to_paint = raw_indices_to_paint - protected_core_indices
    protected_count = len(raw_indices_to_paint) - len(final_indices_to_paint)
    
    print(f"Serine Painting: Identified {num_to_paint} target interface residues.")
    if protected_count > 0:
        print(f"-> Spared {protected_count} of them because they are vital to the structural core.")
    print(f"-> Final painted residues: {len(final_indices_to_paint)}")
    
    return final_indices_to_paint

# ==========================================
# EXTRACT DYNAMIC TARGET SEQUENCES
# ==========================================
target_sequences = {}
if os.path.exists(PARSED_JSONL_PATH):
    with open(PARSED_JSONL_PATH, 'r') as f:
        for line in f:
            data = json.loads(line.strip())
            name = data.get("name")
            seq_a = data.get("seq_chain_A", "")
            if name and seq_a:
                target_sequences[name] = seq_a
else:
    print(f"Error: Cannot find {PARSED_JSONL_PATH}. Run MPNN first.")
    sys.exit(1)

# Dynamically calculate the trap array
TRAP_INDICES = set()
if SERINE_PAINT_FRACTION > 0.0:
    print("Calculating distance-ranked Serine Paint array with Core Protection...")
    TRAP_INDICES = calculate_ranked_proximity_indices(
        WILD_TYPE_PDB_PATH, 
        target_chain="A", 
        paint_fraction=SERINE_PAINT_FRACTION,
        protect_fraction=CORE_PROTECT_FRACTION
    )

# ==========================================
# PARSING FUNCTION (MPNN to Multimer FASTA)
# ==========================================
def parse_mpnn_fasta(filepath):
    sequences = []
    backbone_name = os.path.basename(filepath).replace(".fa", "").replace(".fasta", "")
    
    target_seq = target_sequences.get(backbone_name, "")
    if not target_seq:
        return []

    # >>> DYNAMIC SERINE PAINTING <<<
    seq_list = list(target_seq)
    if SERINE_PAINT_FRACTION > 0.0:
        for idx in TRAP_INDICES:
            if idx < len(seq_list):
                seq_list[idx] = 'S'
                
    patched_target_seq = "".join(seq_list)

    with open(filepath, 'r') as f:
        content = [e for e in f.read().split(">") if e.strip()]
    
    for entry in content:
        lines = entry.strip().splitlines()
        header = lines[0]
        raw_seq = "".join(lines[1:]).strip()
        
        if "sample=" in header and "score=" in header:
            try:
                score_str = header.split("score=")[1].split(",")[0]
                score = float(score_str)
                binder_seq = raw_seq.split("/")[-1] if "/" in raw_seq else raw_seq
                
                # THE STITCH: Painted Target Chain : Binder
                multimer_seq = f"{patched_target_seq}:{binder_seq}"

                multimer_seq = multimer_seq.replace("-", "G")

                sequences.append({
                    "backbone": backbone_name,
                    "header": header,
                    "sequence": multimer_seq,
                    "score": score
                })
            except (IndexError, ValueError):
                continue
                
    return sequences

# ==========================================
# MAIN EXECUTION
# ==========================================
fasta_files = glob.glob(os.path.join(MPNN_OUTPUT_DIR, "*.fa"))

if not fasta_files:
    print(f"Error: No .fa files found in {MPNN_OUTPUT_DIR}/")
    sys.exit(1)

top_ones = []
remaining_pool = []
total_parsed = 0

for f in fasta_files:
    seqs = parse_mpnn_fasta(f)
    if not seqs: continue
    
    total_parsed += len(seqs)
    seqs.sort(key=lambda x: x["score"])
    
    best_seq = seqs.pop(0)
    best_seq["name"] = f"{best_seq['backbone']}_TOP1"
    top_ones.append(best_seq)
    remaining_pool.extend(seqs)

remaining_pool.sort(key=lambda x: x["score"])
wildcards_to_grab = min(WILDCARDS_NEEDED, len(remaining_pool))

wildcards = []
for i in range(wildcards_to_grab):
    w_seq = remaining_pool[i]
    w_seq["name"] = f"{w_seq['backbone']}_WILDCARD_{i+1}"
    wildcards.append(w_seq)

final_results = top_ones + wildcards

for item in final_results:
    out_path = os.path.join(OUTPUT_FASTA_DIR, f"{item['name']}.fasta")
    with open(out_path, 'w') as out_f:
        out_f.write(f">{item['name']} | score={item['score']:.4f}\n")
        out_f.write(item["sequence"] + "\n")

print(f"\nSuccess: Parsed {total_parsed} generated sequences.")
print(f"Saved {len(final_results)} multimer-formatted FASTAs to: {OUTPUT_FASTA_DIR}/")