import os
import glob
import random
from Bio.PDB import PDBParser

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
NATIVE_P53_SEQ = "ETFSDLWKLLPENNV"

# Bulletproof dictionary to avoid Biopython version errors
THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_sequence(chain):
    """Extracts a 1-letter amino acid sequence safely using a standard dictionary."""
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

def generate_mutants(input_dir, output_dir, num_mutants=5):
    os.makedirs(output_dir, exist_ok=True)
    
    # AUTOMATION: Automatically find the lead drug PDB
    pdbs = glob.glob(os.path.join(input_dir, "*.pdb"))
    if not pdbs:
        print(f"Error: No PDB files found in {input_dir}")
        return
        
    latest_pdb = max(pdbs, key=os.path.getctime)
    print(f"Auto-detected lead drug: {latest_pdb}")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", latest_pdb)
    
    model = structure[0]
    mdm2_seq = extract_sequence(model['A'])
    drug_seq = extract_sequence(model['B'])
    
    print(f"Wildtype MDM2 Length: {len(mdm2_seq)} | Lead Drug Length: {len(drug_seq)}")
    
    for i in range(1, num_mutants + 1):
        mutant_mdm2_seq, mutations = mutate_sequence(mdm2_seq)
        mut_string = "_".join(mutations)
        
        # Test A: Evasion (Mutant vs Drug)
        with open(os.path.join(output_dir, f"mutant_{i}_vs_drug.fasta"), "w") as f:
            f.write(f">mutant_{i}_vs_drug | {mut_string}\n{mutant_mdm2_seq}:{drug_seq}\n")
            
        # Test B: Constraint (Mutant vs p53)
        with open(os.path.join(output_dir, f"mutant_{i}_vs_p53.fasta"), "w") as f:
            f.write(f">mutant_{i}_vs_p53 | {mut_string}\n{mutant_mdm2_seq}:{NATIVE_P53_SEQ}\n")

    print(f"\nGenerated {num_mutants * 2} FASTA files in '{output_dir}'")

if __name__ == "__main__":
    # Hardcoded paths matching your specific environment structure
    IN_DIR = "/home/tonypeonio/ProteinDesignChallenge/evaluate_drug_results/"
    OUT_DIR = "/home/tonypeonio/ProteinDesignChallenge/phase2_fastas/"
    
    generate_mutants(IN_DIR, OUT_DIR)