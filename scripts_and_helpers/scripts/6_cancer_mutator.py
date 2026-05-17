import os
import random
import argparse
import json
import pandas as pd
from Bio.PDB import PDBParser
import google.generativeai as genai

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

def get_arms_race_history():
    """Reads the CSV to give the LLM memory of past generations."""
    csv_path = "../../results/arms_race_data.csv"
    if not os.path.exists(csv_path):
        return "No history available. This is Generation 1."
    
    try:
        df = pd.read_csv(csv_path)
        # Give the LLM the last 5 generations to keep the prompt concise
        recent_history = df.tail(5)[['Generation', 'Fitness_Delta', 'Drug_PAE', 'Cancer_vs_Drug_PAE', 'Mutations_vs_WT']].to_dict(orient='records')
        return json.dumps(recent_history, indent=2)
    except Exception as e:
        return f"Could not read history: {e}"

def mutate_sequence_random(seq, num_mutations=2, start_idx=40, end_idx=90):
    """The original stochastic mutator, used as a fallback if the AI fails."""
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

def get_agentic_mutations(seq, num_mutants, num_mutations, start_idx=40, end_idx=90):
    """Calls Gemini to strategically select mutations based on evolutionary history."""
    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        print("WARNING: GEMINI_API_KEY not found in environment. Falling back to random mutations.")
        return [mutate_sequence_random(seq, num_mutations, start_idx, end_idx) for _ in range(num_mutants)]

    genai.configure(api_key=api_key)
    
    history = get_arms_race_history()
    
    prompt = f"""
    You are an evolving cancer protein (MDM2) locked in an evolutionary arms race with a pharmacologist's drug.
    Your goal is to mutate to EVADE the drug (maximize Drug PAE) while maintaining your function (minimize p53 PAE).
    CRITICAL REMINDER: PAE is an error metric. To EVADE the drug, you must MAXIMIZE the Drug PAE. To MAINTAIN function, you must MINIMIZE the p53 PAE.
    
    Here is the recent history of the arms race:
    {history}
    
    Your current amino acid sequence is:
    {seq}
    
    Task:
    Provide exactly {num_mutants} different mutant designs.
    For each design, select exactly {num_mutations} mutations.
    - You must only select 0-indexed positions between {start_idx} and {end_idx}.
    - The new amino acid must be a 1-letter uppercase code from: {AMINO_ACIDS}.
    - Do not pick the amino acid that is already at that position.
    
    Respond STRICTLY with a JSON object matching this schema:
    {{
        "thought_process": "Explain your evolutionary strategy here. What did you learn from the history? Why are you targeting these specific positions?",
        "designs": [
            [
                {{"position": 45, "new_aa": "A"}},
                {{"position": 52, "new_aa": "Y"}}
            ]
        ]
    }}
    """
    
    try:
        print("Consulting Gemini for strategic mutations...")
        # Force strict JSON output so the code doesn't crash on conversational text
        model = genai.GenerativeModel('gemini-2.5-flash', generation_config={"response_mime_type": "application/json"})
        response = model.generate_content(prompt)
        
        data = json.loads(response.text)

        log_path = "../../results/ai_strategy_log.txt"
        thought_process = data.get("thought_process", "No thought process provided.")

        try:
            with open(log_path, "a") as log_file:
                log_file.write(f"=== STRATEGY UPDATE ===\n")
                log_file.write(f"{thought_process}\n\n")
        except Exception as e:
            print(f"Warning: Could not write to log ({e})")
        
        results = []
        for design in data.get("designs", [])[:num_mutants]:
            seq_list = list(seq)
            mutated_positions = []
            
            for mut in design[:num_mutations]:
                idx = mut["position"]
                new_aa = mut["new_aa"]
                
                # Safety checks
                if start_idx <= idx <= end_idx and new_aa in AMINO_ACIDS:
                    current_aa = seq_list[idx]
                    seq_list[idx] = new_aa
                    mutated_positions.append(f"{current_aa}{idx+1}{new_aa}")
            
            results.append(("".join(seq_list), mutated_positions))
            
        # Ensure we got enough designs back; if not, fill with random
        while len(results) < num_mutants:
            results.append(mutate_sequence_random(seq, num_mutations, start_idx, end_idx))
            
        print("Successfully generated AI-driven mutations!")
        return results

    except Exception as e:
        print(f"Gemini API failed ({e}). Falling back to random mutations.")
        return [mutate_sequence_random(seq, num_mutations, start_idx, end_idx) for _ in range(num_mutants)]

def generate_mutants(target_mdm2_pdb, lead_drug_pdb, output_dir, num_mutants, num_mutations):
    os.makedirs(output_dir, exist_ok=True)
    parser = PDBParser(QUIET=True)
    
    print(f"Loading true Pathogen sequence from: {target_mdm2_pdb}")
    target_structure = parser.get_structure("target", target_mdm2_pdb)
    mdm2_seq = extract_sequence(target_structure[0]['A'])
    
    print(f"Loading new Drug sequence from: {lead_drug_pdb}")
    drug_structure = parser.get_structure("drug", lead_drug_pdb)
    chains = list(drug_structure[0].get_chains())
    drug_chain = min(chains, key=lambda c: len(list(c.get_residues())))
    drug_seq = extract_sequence(drug_chain)
    
    print(f"Pathogen Length: {len(mdm2_seq)} | Lead Drug Length: {len(drug_seq)}")
    
    # Get all mutations at once (Agentic or Fallback)
    mutant_data = get_agentic_mutations(mdm2_seq, num_mutants, num_mutations)
    
    for i, (mutant_mdm2_seq, mutations) in enumerate(mutant_data, start=1):
        mut_string = "_".join(mutations) if mutations else "NO_MUT"
        
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