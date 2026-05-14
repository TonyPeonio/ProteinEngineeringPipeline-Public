import argparse
import glob
import json
import os
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("--round", type=int, required=True)
parser.add_argument("--round_dir", type=str, required=True)
parser.add_argument("--log", type=str, required=True)
parser.add_argument("--current_best", type=str, required=True)
args = parser.parse_args()

# 1. Find the ColabFold JSON and PDB output in the correct sub-directory
colabfold_dir = os.path.join(args.round_dir, "colabfold_multimer_results")

# Look specifically for the Rank 1 (best) outputs
json_files = glob.glob(f"{colabfold_dir}/*_scores_rank_001*.json")
pdb_files = glob.glob(f"{colabfold_dir}/*_unrelaxed_rank_001*.pdb")

if not json_files or not pdb_files:
    # If AF2 failed to produce a Rank 1 output, return the old best PDB and don't log
    print(args.current_best) 
    exit()

best_json = json_files[0]
best_pdb = pdb_files[0]

# 2. Extract Metrics
with open(best_json, 'r') as f:
    data = json.load(f)
    plddt = sum(data["plddt"]) / len(data["plddt"])
    # If doing binder design, PAE or pTM is the real metric of success
    pae = data.get("ptm", 0) # Adjust based on your specific ColabFold output format

# 3. Simple Selection Logic (Greedy Algorithm)
# In a real setup, you'd compare this to the previous round's score. 
# For now, we assume it survives and becomes the new parent.
new_seed_path = f"{args.round_dir}/selected_seed_round_{args.round}.pdb"
shutil.copy(best_pdb, new_seed_path)

# 4. Log for your Advisor's Graphs
with open(args.log, 'a') as f:
    # Round, pLDDT, PAE, Mutation tracker (placeholder)
    f.write(f"{args.round},{plddt:.2f},{pae:.4f},N/A\n")

# 5. Print the new PDB path so bash can capture it for the next loop
print(new_seed_path)