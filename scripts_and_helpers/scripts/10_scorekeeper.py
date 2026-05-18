import os
import glob
import json
import argparse
import numpy as np
import pandas as pd

def get_interface_pae(pae_matrix, len_chain_a):
    """Calculates cross-chain PAE by averaging the off-diagonal quadrants."""
    quadrant_1 = np.mean(pae_matrix[0:len_chain_a, len_chain_a:])
    quadrant_2 = np.mean(pae_matrix[len_chain_a:, 0:len_chain_a])
    return (quadrant_1 + quadrant_2) / 2.0

def calculate_sequence_drift(current_fasta, wt_fasta):
    """Calculates the exact number of amino acid differences from the wildtype."""
    try:
        with open(current_fasta, "r") as f:
            current_seq = "".join([l.strip() for l in f if not l.startswith(">")])
        with open(wt_fasta, "r") as f:
            wt_seq = "".join([l.strip() for l in f if not l.startswith(">")])
            
        # Count mismatches (Hamming distance)
        # We use zip to pair them up; if lengths differ, it safely stops at the shorter length
        drift = sum(1 for a, b in zip(current_seq, wt_seq) if a != b)
        
        # Add any length difference as mutations (in case of insertions/deletions later)
        drift += abs(len(current_seq) - len(wt_seq))
        return drift
    except Exception as e:
        print(f"Warning: Could not calculate sequence drift ({e}). Defaulting to 0.")
        return 0

def extract_metrics(results_dir, output_dir):
    """Extracts metrics from both the Drug and p53 ColabFold JSONs."""
    drug_jsons = glob.glob(os.path.join(results_dir, "complex_vs_drug*rank_001*_scores.json"))
    p53_jsons = glob.glob(os.path.join(results_dir, "complex_vs_p53*rank_001*_scores.json"))
    
    if not drug_jsons or not p53_jsons:
        raise FileNotFoundError(f"Missing required scores.json in {results_dir}")
        
    cancer_fasta = os.path.join(output_dir, "current_cancer.fasta")
    with open(cancer_fasta, "r") as f:
        len_cancer = len("".join([l.strip() for l in f if not l.startswith(">")]))
    
    # Process Drug Matchup
    with open(drug_jsons[0], "r") as f:
        drug_data = json.load(f)
    drug_pae_matrix = np.array(drug_data["pae"])
    drug_plddt_array = np.array(drug_data["plddt"])
    
    drug_pae = np.mean(drug_pae_matrix[len_cancer:, len_cancer:])
    cancer_vs_drug_pae = get_interface_pae(drug_pae_matrix, len_cancer)
    drug_plddt = np.mean(drug_plddt_array[len_cancer:])
    
    # Process p53 Control
    with open(p53_jsons[0], "r") as f:
        p53_data = json.load(f)
    p53_pae_matrix = np.array(p53_data["pae"])
    
    cancer_vs_p53_pae = get_interface_pae(p53_pae_matrix, len_cancer)
    
    return drug_pae, cancer_vs_drug_pae, drug_plddt, cancer_vs_p53_pae

def update_history(gen, results_dir, csv_path, wt_cancer_fasta):
    output_dir = os.path.dirname(results_dir)
    current_cancer_fasta = os.path.join(output_dir, "current_cancer.fasta")
    
    try:
        drug_pae, cancer_vs_drug_pae, drug_plddt, cancer_vs_p53_pae = extract_metrics(results_dir, output_dir)
    except Exception as e:
        print(f"Error extracting metrics: {e}")
        drug_pae, cancer_vs_drug_pae, drug_plddt, cancer_vs_p53_pae = 30.0, 30.0, 40.0, 30.0 

    # --- THE FIX: Calculate exact mutations ---
    exact_mutations = calculate_sequence_drift(current_cancer_fasta, wt_cancer_fasta)
    
    fitness_delta = cancer_vs_drug_pae - cancer_vs_p53_pae
    
    new_data = {
        "Generation": gen,
        "Drug_PAE": round(drug_pae, 2),
        "Cancer_vs_Drug_PAE": round(cancer_vs_drug_pae, 2),
        "Cancer_vs_p53_PAE": round(cancer_vs_p53_pae, 2),
        "Fitness_Delta": round(fitness_delta, 2),
        "Drug_pLDDT": round(drug_plddt, 2),
        "Mutations_vs_WT": exact_mutations
    }
    
    df_new = pd.DataFrame([new_data])
    
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        df = df[df['Generation'] != gen] 
        df = pd.concat([df, df_new], ignore_index=True)
    else:
        df = df_new
        
    df.to_csv(csv_path, index=False)
    
    print(f"--- GENERATION {gen} SCORECARD ---")
    print(f"Cancer Evasion (Drug PAE): {new_data['Cancer_vs_Drug_PAE']}")
    print(f"Cancer Function (p53 PAE): {new_data['Cancer_vs_p53_PAE']}")
    print(f"Overall Fitness Delta:     {new_data['Fitness_Delta']}")
    print(f"Exact Mutations vs WT:     {new_data['Mutations_vs_WT']}")
    print(f"Metrics saved to {csv_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gen", type=int, required=True)
    parser.add_argument("--results_dir", type=str, required=True)
    parser.add_argument("--history_csv", type=str, required=True)
    parser.add_argument("--wt_cancer_fasta", type=str, required=True, help="Path to the Generation 0 baseline cancer fasta")
    args = parser.parse_args()
    
    update_history(args.gen, args.results_dir, args.history_csv, args.wt_cancer_fasta)