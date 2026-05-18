import argparse
import os

NATIVE_P53_SEQ = "ETFSDLWKLLPENNV"

def stitch_fastas(cancer_fasta, drug_fasta, out_dir):
    # Read cancer sequence
    with open(cancer_fasta, "r") as f:
        cancer_seq = "".join([line.strip() for line in f.readlines() if not line.startswith(">")])
        
    # Read drug sequence
    with open(drug_fasta, "r") as f:
        drug_seq = "".join([line.strip() for line in f.readlines() if not line.startswith(">")])
        
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. Cancer vs Drug (The Arms Race)
    with open(os.path.join(out_dir, "complex_vs_drug.fasta"), "w") as f:
        f.write(f">complex_vs_drug\n{cancer_seq}:{drug_seq}\n")
        
    # 2. Cancer vs p53 (The Functional Control)
    with open(os.path.join(out_dir, "complex_vs_p53.fasta"), "w") as f:
        f.write(f">complex_vs_p53\n{cancer_seq}:{NATIVE_P53_SEQ}\n")
        
    print(f"Stitched combat files saved to {out_dir}")
    print(f"Cancer Length: {len(cancer_seq)} | Drug Length: {len(drug_seq)} | p53 Length: {len(NATIVE_P53_SEQ)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer_fasta", type=str, required=True)
    parser.add_argument("--drug_fasta", type=str, required=True)
    parser.add_argument("--out_dir", type=str, required=True)
    args = parser.parse_args()
    
    stitch_fastas(args.cancer_fasta, args.drug_fasta, args.out_dir)