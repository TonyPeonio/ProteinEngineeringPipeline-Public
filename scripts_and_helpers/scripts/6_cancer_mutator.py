import os
import random
import argparse
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_one, is_aa

# Standard amino acids for mutation
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
# The native p53 peptide sequence from 1YCR that MDM2 MUST still bind
NATIVE_P53_SEQ = "ETFSDLWKLLPENNV"

def extract_sequence(chain):
    """Extracts a 1-letter amino acid sequence from a PDB chain."""
    seq = ""
    for residue in chain:
        if is_aa(residue, standard=True):
            try:
                seq += three_to_one(residue.get_resname())
            except KeyError:
                pass
    return seq

def mutate_sequence(seq, num_mutations=2, start_idx=40, end_idx=90):
    """Randomly mutates amino acids in the target binding region."""
    seq_list = list(seq)
    # Ensure indices are within bounds
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

def generate_mutants(input_pdb, output_dir, num_mutants=5):
    os.makedirs(output_dir, exist_ok=True)
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", input_pdb)
    
    # Assuming Chain A is MDM2 and Chain B is the Drug
    model = structure[0]
    mdm2_seq = extract_sequence(model['A'])
    drug_seq = extract_sequence(model['B'])
    
    print(f"Wildtype MDM2 Length: {len(mdm2_seq)}")
    print(f"Lead Drug Length: {len(drug_seq)}")
    
    for i in range(1, num_mutants + 1):
        # 1. Mutate MDM2
        mutant_mdm2_seq, mutations = mutate_sequence(mdm2_seq)
        mut_string = "_".join(mutations)
        print(f"Generating Mutant {i}: {mut_string}")
        
        # 2. Write Mutant vs Drug FASTA (Test A: Evasion)
        drug_fasta_path = os.path.join(output_dir, f"mutant_{i}_vs_drug.fasta")
        with open(drug_fasta_path, "w") as f:
            f.write(f">mutant_{i}_vs_drug | {mut_string}\n")
            f.write(f"{mutant_mdm2_seq}:{drug_seq}\n")
            
        # 3. Write Mutant vs Native p53 FASTA (Test B: Constraint)
        p53_fasta_path = os.path.join(output_dir, f"mutant_{i}_vs_p53.fasta")
        with open(p53_fasta_path, "w") as f:
            f.write(f">mutant_{i}_vs_p53 | {mut_string}\n")
            f.write(f"{mutant_mdm2_seq}:{NATIVE_P53_SEQ}\n")

    print(f"\nGenerated {num_mutants * 2} FASTA files in '{output_dir}'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_pdb", help="Path to the winning drug PDB")
    parser.add_argument("--outdir", default="./phase2_fastas", help="Output directory for fastas")
    
    args = parser.parse_args()
    generate_mutants(args.input_pdb, args.outdir)