import sys
import os
import urllib.request

def main():
    # Ensure the user provided a PDB ID
    if len(sys.argv) != 2:
        print("Usage: python fetch_and_clean.py <PDB_ID>")
        print("Example: python fetch_and_clean.py 7ut8")
        sys.exit(1)

    pdb_id = sys.argv[1].lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    # Define paths (adjust base_dir if your folder structure requires it)
    base_dir = "../pdb"
    
    raw_pdb_path = os.path.join(base_dir, f"{pdb_id}.pdb")
    clean_pdb_path = os.path.join(base_dir, f"{pdb_id}_clean.pdb")

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

    # Step 2: Clean the PDB (Keep only ATOM records)
    print(f"Cleaning {pdb_id.upper()}...")
    with open(raw_pdb_path, 'r') as infile, open(clean_pdb_path, 'w') as outfile:
        for line in infile:
            if line.startswith('ATOM  '):
                outfile.write(line)

    print(f"Success! Cleaned PDB saved to: {clean_pdb_path}")

if __name__ == "__main__":
    main()