import os
import glob
import json
import shutil
import numpy as np
import re
import argparse

def evaluate_binders(colabfold_out_dir, output_dir, target_gen):
    """
    Finds the model with the lowest PAE for the SPECIFIC generation.
    Handles the case where Gen 1 is untagged and Gen 2+ is 'genN'.
    """
    all_json_files = glob.glob(os.path.join(colabfold_out_dir, "*_scores*.json"))
    current_gen_files = []

    # Logic: If we are on Bash GEN 1, look for files WITHOUT 'gen' in the name.
    # If we are on Bash GEN 2, look for 'gen1'. 
    # Therefore, we look for 'gen' + (target_gen - 1)
    search_tag = f"gen{target_gen - 1}"

    for file in all_json_files:
        basename = os.path.basename(file)
        match = re.search(r'gen(\d+)', basename)
        
        if match:
            # File has a tag; does it match our search_tag?
            if f"gen{match.group(1)}" == search_tag:
                current_gen_files.append(file)
        else:
            # File has NO tag; it must be from the very first wildtype round
            if target_gen == 1:
                current_gen_files.append(file)

    if not current_gen_files:
        print(f"CRITICAL: No files found for Generation {target_gen} (Search tag: {search_tag})!")
        return
    
    os.makedirs(output_dir, exist_ok=True)

    best_score = float('inf') 
    best_json = None
    best_metrics = {}

    # !!! FIXED: Iterate over current_gen_files, NOT all_json_files !!!
    print(f"Evaluating {len(current_gen_files)} designs for Generation {target_gen}...")

    for j_file in current_gen_files:
        with open(j_file, 'r') as f:
            data = json.load(f)
        
        mean_pae = np.mean(data.get('pae', 100)) 
        mean_plddt = np.mean(data.get('plddt', 0))
        iptm = data.get('iptm', 0)

        # Optional: Print progress for each
        # print(f"  Checking {os.path.basename(j_file)} -> PAE: {mean_pae:.2f}")

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
                print(f"Error: Matching PDB not found for {identifier}")
    else:
        print("Error: No winning JSON found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gen", type=int, required=True, help="Current generation number from bash")
    args = parser.parse_args()
    
    IN_DIR = "/home/tonypeonio/ProteinDesignChallenge/outputs/colabfold_multimer_results"
    OUT_DIR = "/home/tonypeonio/ProteinDesignChallenge/evaluate_drug_results/"
    
    evaluate_binders(IN_DIR, OUT_DIR, args.gen)