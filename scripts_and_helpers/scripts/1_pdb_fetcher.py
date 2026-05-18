import sys
import os
import urllib.request

# Standard amino acid 3-to-1 letter translation
AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def main():
    # Ensure the user provided a PDB ID
    if len(sys.argv) != 2:
        print("Usage: python fetch_and_clean.py <PDB_ID>")
        print("Example: python fetch_and_clean.py 7ut8")
        sys.exit(1)

    pdb_id = sys.argv[1].lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    # Define paths
    base_dir = "../pdb"
    
    raw_pdb_path = os.path.join(base_dir, f"{pdb_id}.pdb")
    clean_pdb_path = os.path.join(base_dir, f"{pdb_id}_clean.pdb")
    fasta_path = os.path.join(base_dir, f"../../outputs/current_cancer.fasta")

    # Create directories if they do not exist
    os.makedirs(base_dir, exist_ok=True)

    # Step 1: Fetch the PDB from RCSB
    print(f"Fetching PDB {pdb_id.upper()} from RCSB...")
    try:
        urllib.request.urlretrieve(url, raw_pdb_path)
    except Exception as e:
        print(f"Error fetching PDB {pdb_id.upper()}: {e}")
        print("Please check your internet connection or verify the PDB ID.")
        sys.exit(1)

    # Step 2: Clean the PDB & Extract FASTA
    print(f"Cleaning {pdb_id.upper()} and extracting sequence...")
    
    seq = []
    seen_residues = set()
    
    with open(raw_pdb_path, 'r') as infile, open(clean_pdb_path, 'w') as outfile:
        for line in infile:
            if line.startswith('ATOM  '):
                outfile.write(line) # Save the clean ATOM line
                
                # Extract FASTA logic
                res_name = line[17:20].strip()
                res_num = line[22:26].strip()
                chain_id = line[21]
                
                residue_id = f"{chain_id}_{res_num}"
                
                # Only grab the amino acid once per residue (not every atom)
                if residue_id not in seen_residues:
                    seen_residues.add(residue_id)
                    seq.append(AA_MAP.get(res_name, 'X'))

    # Save the FASTA
    with open(fasta_path, 'w') as f:
        f.write(f">cancer_gen0_baseline\n{''.join(seq)}\n")

    print(f"Success! Cleaned PDB saved to: {clean_pdb_path}")
    print(f"Success! Cancer FASTA saved to: {fasta_path}")

if __name__ == "__main__":
    main()