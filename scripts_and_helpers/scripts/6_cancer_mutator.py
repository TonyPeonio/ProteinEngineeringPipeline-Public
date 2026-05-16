import os
import random
import argparse
from Bio.PDB import PDBParser

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
NATIVE_P53_SEQ = "ETFSDLWKLLPENNV"

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_sequence(chain):
    seq = ""
    for residue in chain:
        resname = residue.get_resname().strip().upper()
        if resname in THREE_TO_ONE:
            seq += THREE_TO_ONE[resname]
    return seq

def mutate_sequence(seq, num_mutations=2, start_idx=40, end_idx=90):
    seq_list = list(seq)
    start = max(0, start_idx)
    end = min(len(seq) - 1, end_idx)
    
    mutated_positions = []
    for _ in range(num_mutations):
        idx = random.randint(start, end)
        current_aa = seq_list[idx]
        new_aa = random.choice([aa for aa in AMINO_ACIDS if aa != current_aa])
        seq_list[idx] = new_aa
        mutated_positions.append(f"{current_aa}{idx+1}{new_aa}")
        
    return "".join(seq_list), mutated_positions

def generate_mutants(target_mdm2_pdb, lead_drug_pdb, output_dir, num_mutants, num_mutations):
    os.makedirs(output_dir, exist_ok=True)
    parser = PDBParser(QUIET=True)
    
    # 1. EXTRACT THE TRUE PATHOGEN SEQUENCE from the safe archive!
    print(f"Loading true Pathogen sequence from: {target_mdm2_pdb}")
    target_structure = parser.get_structure("target", target_mdm2_pdb)
    mdm2_seq = extract_sequence(target_structure[0]['A'])
    
    # 2. EXTRACT THE NEW DRUG SEQUENCE from the pharmacologist output!
    print(f"Loading new Drug sequence from: {lead_drug_pdb}")
    drug_structure = parser.get_structure("drug", lead_drug_pdb)
    chains = list(drug_structure[0].get_chains())
    drug_chain = min(chains, key=lambda c: len(list(c.get_residues())))
    drug_seq = extract_sequence(drug_chain)
    
    print(f"Pathogen Length: {len(mdm2_seq)} | Lead Drug Length: {len(drug_seq)}")
    
    for i in range(1, num_mutants + 1):
        mutant_mdm2_seq, mutations = mutate_sequence(mdm2_seq, num_mutations)
        mut_string = "_".join(mutations)
        
        with open(os.path.join(output_dir, f"mutant_{i}_vs_drug.fasta"), "w") as f:
            f.write(f">mutant_{i}_vs_drug | {mut_string}\n{mutant_mdm2_seq}:{drug_seq}\n")
            
        with open(os.path.join(output_dir, f"mutant_{i}_vs_p53.fasta"), "w") as f:
            f.write(f">mutant_{i}_vs_p53 | {mut_string}\n{mutant_mdm2_seq}:{NATIVE_P53_SEQ}\n")

    print(f"\nGenerated {num_mutants * 2} FASTA files in '{output_dir}'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--target_mdm2", type=str, required=True, help="The currently evolved pathogen target")
    parser.add_argument("--lead_drug", type=str, required=True, help="The newly designed drug")
    parser.add_argument("--out_dir", type=str, required=True)
    parser.add_argument("--num_mutants", type=int, required=True)
    parser.add_argument("--num_mutations", type=int, required=True)
    args = parser.parse_args()
    
    generate_mutants(args.target_mdm2, args.lead_drug, args.out_dir, args.num_mutants, args.num_mutations)