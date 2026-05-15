import os
import glob
import json
import shutil
import numpy as np

def evaluate_binders(colabfold_out_dir, output_pdb_name):
    """
    Parses ColabFold JSON outputs, finds the model with the lowest PAE, 
    and isolates it as the current lead drug.
    """
    # Find all ColabFold score JSON files
    json_files = glob.glob(os.path.join(colabfold_out_dir, "*_scores*.json"))
    
    if not json_files:
        print(f"Error: No score JSON files found in {colabfold_out_dir}")
        return
    
    os.makedirs(output_pdb_name, exist_ok=True)

    best_score = float('inf') # We want the lowest PAE
    best_pdb = None
    best_metrics = {}

    print(f"Evaluating {len(json_files)} designs from Phase 1...")

    for j_file in json_files:
        with open(j_file, 'r') as f:
            data = json.load(f)
        
        # ColabFold stores PAE as a 2D matrix. We calculate the overall mean.
        # It also stores per-residue pLDDT, and usually a single iPTM float for multimers.
        mean_pae = np.mean(data.get('pae', 100)) 
        mean_plddt = np.mean(data.get('plddt', 0))
        iptm = data.get('iptm', 0)

        print(f"Design: {os.path.basename(j_file)}")
        print(f"  -> Mean PAE: {mean_pae:.2f} | Mean pLDDT: {mean_plddt:.2f} | iPTM: {iptm:.3f}")

        # MVP Logic: Pick the design with the lowest Mean PAE
        if mean_pae < best_score:
            best_score = mean_pae
            
            # The PDB file shares the exact same prefix as the JSON
            expected_pdb = j_file.replace("_scores.json", ".pdb")
            
            if os.path.exists(expected_pdb):
                best_pdb = expected_pdb
                best_metrics = {'pae': mean_pae, 'plddt': mean_plddt, 'iptm': iptm}
            else:
                print(f"  -> Warning: Expected PDB {expected_pdb} not found!")

    if best_pdb:
        print("\n" + "="*40)
        print(f"WINNING DRUG SELECTED:")
        print(f"File: {os.path.basename(best_pdb)}")
        print(f"Scores: PAE={best_metrics['pae']:.2f}, pLDDT={best_metrics['plddt']:.2f}, iPTM={best_metrics['iptm']:.3f}")
        print(f"Copying to -> {output_pdb_name}")
        print("="*40)
        
        shutil.copy(best_pdb, output_pdb_name)
    else:
        print("Failed to find a winning PDB. Check your ColabFold output directory.")

if __name__ == "__main__":
    evaluate_binders("/home/tonypeonio/ProteinDesignChallenge/outputs/colabfold_multimer_results", "/home/tonypeonio/ProteinDesignChallenge/evaluate_drug_results/")