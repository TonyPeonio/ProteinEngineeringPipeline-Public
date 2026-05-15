import os
import glob
import json
import shutil
import numpy as np
import re
from Bio import PDB

def extract_chain_a(input_pdb, output_pdb):
    """Extracts only Chain A (the mutated MDM2) to feed back into RFdiffusion."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("complex", input_pdb)
    io = PDB.PDBIO()
    
    class ChainSelect(PDB.Select):
        def accept_chain(self, chain):
            return chain.get_id() == 'A'
            
    io.set_structure(structure)
    io.save(output_pdb, ChainSelect())

def evaluate_cancer(colabfold_out_dir, output_dir, num_mutants=5):
    os.makedirs(output_dir, exist_ok=True)
    
    best_delta = -float('inf')
    winning_mutant = None
    winning_drug_pdb = None
    best_metrics = {}

    print("Evaluating Phase 2 Cancer Mutations...")

    for i in range(1, num_mutants + 1):
        drug_json_pattern = os.path.join(colabfold_out_dir, f"*mutant_{i}_vs_drug*scores*.json")
        p53_json_pattern = os.path.join(colabfold_out_dir, f"*mutant_{i}_vs_p53*scores*.json")
        
        drug_jsons = glob.glob(drug_json_pattern)
        p53_jsons = glob.glob(p53_json_pattern)
        
        if not drug_jsons or not p53_jsons:
            print(f"  -> Missing ColabFold outputs for mutant {i}. Skipping.")
            continue
            
        drug_json = sorted(drug_jsons)[0]
        p53_json = sorted(p53_jsons)[0]
        
        with open(drug_json, 'r') as f:
            drug_pae = np.mean(json.load(f).get('pae', 100))
        with open(p53_json, 'r') as f:
            p53_pae = np.mean(json.load(f).get('pae', 100))
            
        delta = drug_pae - p53_pae
        print(f"Mutant {i}: Drug PAE = {drug_pae:.2f} | p53 PAE = {p53_pae:.2f} | Fitness Delta = {delta:.2f}")
        
        if delta > best_delta:
            best_delta = delta
            winning_mutant = i
            best_metrics = {'drug_pae': drug_pae, 'p53_pae': p53_pae, 'delta': delta}
            
            match = re.search(r'(rank_\d+_[a-zA-Z0-9_]+_model_\d+_seed_\d+)', drug_json)
            if match:
                identifier = match.group(1)
                winning_drug_pdb = glob.glob(os.path.join(colabfold_out_dir, f"*{identifier}*.pdb"))[0]

    if winning_drug_pdb:
        final_output_path = os.path.join(output_dir, "next_gen_target.pdb")
        
        print("\n" + "="*40)
        print(f"WINNING EVASION MUTATION: Mutant {winning_mutant}")
        print(f"Fitness Delta: {best_metrics['delta']:.2f}")
        print(f"Extracting new MDM2 target to -> {final_output_path}")
        print("="*40)
        
        extract_chain_a(winning_drug_pdb, final_output_path)
        
        # ==========================================
        # AUTOMATIC CSV LOGGING
        # ==========================================
        results_dir = "/home/tonypeonio/ProteinDesignChallenge/results"
        os.makedirs(results_dir, exist_ok=True)
        
        csv_path = os.path.join(results_dir, "arms_race_data.csv")
        drug_pae_path = os.path.join(results_dir, "current_drug_pae.txt")
        
        pharmacologist_pae = 0.0
        if os.path.exists(drug_pae_path):
            with open(drug_pae_path, "r") as f:
                pharmacologist_pae = float(f.read().strip())
                
        # Create CSV and header if it doesn't exist
        if not os.path.exists(csv_path):
            with open(csv_path, "w") as f:
                f.write("Generation,Drug_PAE,Cancer_vs_Drug_PAE,Cancer_vs_p53_PAE,Fitness_Delta\n")
                
        # Calculate generation based on line count (Line 1 is header, Line 2 is Gen 1, etc.)
        with open(csv_path, "r") as f:
            gen = sum(1 for _ in f)
            
        # Append the new data
        with open(csv_path, "a") as f:
            f.write(f"{gen},{pharmacologist_pae:.2f},{best_metrics['drug_pae']:.2f},{best_metrics['p53_pae']:.2f},{best_metrics['delta']:.2f}\n")
            
        print(f"Successfully appended Generation {gen} data to arms_race_data.csv!")

    else:
        print("Failed to identify a winning mutant.")

if __name__ == "__main__":
    IN_DIR = "/home/tonypeonio/ProteinDesignChallenge/outputs/colabfold_phase2_results"
    OUT_DIR = "/home/tonypeonio/ProteinDesignChallenge/evaluate_cancer_results"
    
    evaluate_cancer(IN_DIR, OUT_DIR)