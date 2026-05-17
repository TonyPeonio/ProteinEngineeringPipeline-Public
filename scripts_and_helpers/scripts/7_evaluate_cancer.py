import os
import glob
import json
import shutil
import numpy as np
import re
import argparse
from Bio import PDB

# Add this dictionary so Biopython can convert the PDB residues to a string
THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_sequence_from_pdb(pdb_path, chain_id='A'):
    """Extracts a 1-letter amino acid sequence from a specific chain in a PDB."""
    if not os.path.exists(pdb_path):
        return None
        
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_path)
    model = structure[0]
    
    if chain_id not in model:
        return None
        
    chain = model[chain_id]
    seq = ""
    for residue in chain:
        resname = residue.get_resname().strip().upper()
        if resname in THREE_TO_ONE:
            seq += THREE_TO_ONE[resname]
    return seq

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

def get_sequence_from_fasta(fasta_path):
    if not os.path.exists(fasta_path):
        return None
    with open(fasta_path, 'r') as f:
        return "".join([line.strip() for line in f if not line.startswith(">")])

def evaluate_cancer(colabfold_out_dir, output_dir, current_gen, num_mutants=5):
    os.makedirs(output_dir, exist_ok=True)
    results_dir = "/home/tonypeonio/ProteinDesignChallenge_Agentic/results"
    os.makedirs(results_dir, exist_ok=True)
    
    best_delta = -float('inf')
    winning_mutant = None
    winning_drug_pdb = None
    best_metrics = {}

    print(f"--- Evaluating Phase 2 Cancer Mutations (Gen {current_gen}) ---")

    for i in range(1, num_mutants + 1):
        drug_json_pattern = os.path.join(colabfold_out_dir, f"*mutant_{i}_vs_drug*scores*.json")
        p53_json_pattern = os.path.join(colabfold_out_dir, f"*mutant_{i}_vs_p53*scores*.json")
        
        drug_jsons = glob.glob(drug_json_pattern)
        p53_jsons = glob.glob(p53_json_pattern)
        
        if not drug_jsons or not p53_jsons:
            continue
            
        # Select the best model for this specific mutant
        drug_json = sorted(drug_jsons)[0]
        p53_json = sorted(p53_jsons)[0]
        
        with open(drug_json, 'r') as f:
            drug_pae = np.mean(json.load(f).get('pae', 100))
        with open(p53_json, 'r') as f:
            p53_pae = np.mean(json.load(f).get('pae', 100))
            
        delta = drug_pae - p53_pae
        print(f"Mutant {i}: Drug PAE = {drug_pae:.2f} | p53 PAE = {p53_pae:.2f} | Delta = {delta:.2f}")
        
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
        extract_chain_a(winning_drug_pdb, final_output_path)
        
        # --- DIAGNOSTIC: Calculate Mutation Drift ---
        # Update this path to point exactly to your Gen 1 starting PDB
        wt_pdb = "/home/tonypeonio/ProteinDesignChallenge_Agentic/scripts_and_helpers/pdb/1ycr_clean.pdb"
        
        mutant_fasta_list = glob.glob(os.path.join("/home/tonypeonio/ProteinDesignChallenge_Agentic/outputs/phase2_fastas", f"*mutant_{winning_mutant}_vs_drug.fasta"))
        
        mutation_count = 0
        if os.path.exists(wt_pdb) and mutant_fasta_list:
            # Extract Chain A (MDM2) from the Wildtype PDB
            wt_seq = extract_sequence_from_pdb(wt_pdb, chain_id='A')
            mut_seq = get_sequence_from_fasta(mutant_fasta_list[0])
            
            if wt_seq and mut_seq:
                # Isolate just the MDM2 portion of the multimer FASTA
                mut_mdm2_only = mut_seq.split(":")[0] 
                
                # Compare the two strings, stopping at whichever is shorter just in case
                mutation_count = sum(1 for a, b in zip(wt_seq, mut_mdm2_only) if a != b)
                
                print(f"DIAGNOSTIC -> Wildtype Length: {len(wt_seq)} | Mutant Length: {len(mut_mdm2_only)}")
                print(f"DIAGNOSTIC -> Drift calculated: {mutation_count} residues")
        else:
            print("⚠️ WARNING: Could not calculate drift. Missing WT PDB or Mutant FASTA.")

        # --- DIAGNOSTIC: Load Pharmacologist Stats ---
        pharmacologist_pae = 0.0
        pharmacologist_plddt = 0.0
        
        pae_file = os.path.join(results_dir, "current_drug_pae.txt")
        plddt_file = os.path.join(results_dir, "current_drug_plddt.txt")
        
        if os.path.exists(pae_file):
            with open(pae_file, "r") as f: pharmacologist_pae = float(f.read().strip())
        if os.path.exists(plddt_file):
            with open(plddt_file, "r") as f: pharmacologist_plddt = float(f.read().strip())

        # --- LOG TO CSV ---
        csv_path = os.path.join(results_dir, "arms_race_data.csv")
        header = "Generation,Drug_PAE,Cancer_vs_Drug_PAE,Cancer_vs_p53_PAE,Fitness_Delta,Drug_pLDDT,Mutations_vs_WT\n"
        
        if not os.path.exists(csv_path):
            with open(csv_path, "w") as f: f.write(header)
                
        with open(csv_path, "a") as f:
            f.write(f"{current_gen},{pharmacologist_pae:.2f},{best_metrics['drug_pae']:.2f},"
                    f"{best_metrics['p53_pae']:.2f},{best_metrics['delta']:.2f},"
                    f"{pharmacologist_plddt:.2f},{mutation_count}\n")
            
        print("\n" + "="*40)
        print(f"GENERATION {current_gen} SUMMARY")
        print(f"Winning Mutant: {winning_mutant} | Total Drift: {mutation_count} residues")
        print(f"Fitness Delta: {best_metrics['delta']:.2f}")
        print(f"Exported target for Gen {current_gen + 1}")
        print("="*40)
    else:
        print("Error: Could not identify a winning mutant.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gen", type=int, required=True)
    args = parser.parse_args()
    
    IN_DIR = "../../outputs/colabfold_phase2_results"
    OUT_DIR = "../../outputs/evaluate_cancer_results"    
    evaluate_cancer(IN_DIR, OUT_DIR, args.gen)