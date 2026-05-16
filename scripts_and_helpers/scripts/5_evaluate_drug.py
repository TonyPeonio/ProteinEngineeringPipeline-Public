import os
import glob
import json
import shutil
import numpy as np
import re
import argparse

def evaluate_binders(colabfold_out_dir, output_dir, target_gen):
    """
    Evaluates all json files in the target directory to find the winning drug.
    Assumes the directory is wiped clean each generation by the Bash loop.
    """
    current_gen_files = glob.glob(os.path.join(colabfold_out_dir, "*_scores*.json"))

    if not current_gen_files:
        # Throw an actual Exception so Bash set -e catches it and halts the pipeline!
        raise FileNotFoundError(f"CRITICAL: No JSON score files found in {colabfold_out_dir}!")
    
    os.makedirs(output_dir, exist_ok=True)

    best_score = float('inf') 
    best_json = None
    best_metrics = {}

    print(f"Evaluating {len(current_gen_files)} designs for Generation {target_gen}...")

    for j_file in current_gen_files:
        with open(j_file, 'r') as f:
            data = json.load(f)
        
        mean_pae = np.mean(data.get('pae', 100)) 
        mean_plddt = np.mean(data.get('plddt', 0))
        iptm = data.get('iptm', 0)

        if mean_pae < best_score:
            best_score = mean_pae
            best_json = j_file
            best_metrics = {'pae': mean_pae, 'plddt': mean_plddt, 'iptm': iptm}

    # Save the winning drug PAE for the graphing script
    results_root = "/home/tonypeonio/ProteinDesignChallenge/results"
    os.makedirs(results_root, exist_ok=True)
    with open(os.path.join(results_root, "current_drug_pae.txt"), "w") as f:
        f.write(str(best_score))
    with open(os.path.join(results_root, "current_drug_plddt.txt"), "w") as f:
        f.write(str(best_metrics['plddt']))

    if best_json:
        # Match the PDB to the JSON using the rank/model identifier
        match = re.search(r'(rank_\d+_[a-zA-Z0-9_]+_model_\d+_seed_\d+)', best_json)
        
        if match:
            identifier = match.group(1)
            possible_pdbs = glob.glob(os.path.join(colabfold_out_dir, f"*{identifier}*.pdb"))
            
            if possible_pdbs:
                best_pdb = possible_pdbs[0]
                final_output_path = os.path.join(output_dir, "lead_drug.pdb")
                
                print("\n" + "="*40)
                print(f"WINNING DRUG SELECTED FOR GEN {target_gen}:")
                print(f"Source: {os.path.basename(best_pdb)}")
                print(f"Scores: PAE={best_metrics['pae']:.2f}, pLDDT={best_metrics['plddt']:.2f}")
                print(f"Copying to -> {final_output_path}")
                print("="*40)
                
                shutil.copy(best_pdb, final_output_path)
            else:
                raise FileNotFoundError(f"Error: Matching PDB not found for {identifier}")
    else:
        raise ValueError("Error: No winning JSON found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gen", type=int, required=True, help="Current generation number from bash")
    args = parser.parse_args()
    
    # Ensure these paths match your centralized outputs folder
    IN_DIR = "/home/tonypeonio/ProteinDesignChallenge/outputs/colabfold_multimer_results"
    OUT_DIR = "/home/tonypeonio/ProteinDesignChallenge/outputs/evaluate_drug_results"
    
    evaluate_binders(IN_DIR, OUT_DIR, args.gen)