import os
import random
import argparse
import json
import pandas as pd
import google.generativeai as genai
from dotenv import load_dotenv

load_dotenv(dotenv_path="../../.env")

gemini_api_key = os.getenv("GEMINI_API_KEY")

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

def get_arms_race_history(csv_path):
    """Reads the CSV to give the LLM memory of past generations."""
    if not os.path.exists(csv_path):
        return "No history available. This is Generation 1."
    
    try:
        df = pd.read_csv(csv_path)
        # Give the LLM the last 5 generations to keep the prompt concise
        recent_history = df.tail(5)[['Generation', 'Cancer_vs_Drug_PAE']].to_dict(orient='records')
        return json.dumps(recent_history, indent=2)
    except Exception as e:
        return f"Could not read history: {e}"

def mutate_sequence_random(seq, num_mutations):
    """The original stochastic mutator, used as a fallback if the AI fails."""
    seq_list = list(seq)
    end = len(seq) - 1
    
    mutated_positions = []
    for _ in range(num_mutations):
        idx = random.randint(0, end)
        current_aa = seq_list[idx]
        new_aa = random.choice([aa for aa in AMINO_ACIDS if aa != current_aa])
        seq_list[idx] = new_aa
        mutated_positions.append(f"{current_aa}{idx+1}{new_aa}")
        
    return "".join(seq_list), mutated_positions

def get_agentic_mutations(seq, num_mutations, csv_path):
    """Calls Gemini to strategically select drug mutations."""
    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        print("WARNING: GEMINI_API_KEY not found. Falling back to random mutations.")
        return mutate_sequence_random(seq, num_mutations)

    genai.configure(api_key=api_key)
    history = get_arms_race_history(csv_path)
    
    start_idx = 0
    end_idx = len(seq) - 1
    
    prompt = f"""
    You are a pharmacologist evolving a peptide drug to fight in an evolutionary arms race with a rapidly mutating cancer protein (MDM2).
    Your goal is to mutate your drug sequence to tightly bind to the cancer protein (you want to MINIMIZE the Cancer_vs_Drug_PAE score).
    
    Here is the recent history of the arms race PAE scores:
    {history}
    
    Your current drug amino acid sequence is:
    {seq}
    
    Task:
    Select exactly {num_mutations} strategic point mutations to make to your drug.
    - You must only select 0-indexed positions between {start_idx} and {end_idx}.
    - The new amino acid must be a 1-letter uppercase code from: {AMINO_ACIDS}.
    - Do not pick the amino acid that is already at that position.
    
    Respond STRICTLY with a JSON object matching this schema:
    {{
        "thought_process": "Explain your evolutionary strategy here. What did you learn from the history? Why are you targeting these specific positions?",
        "design": [
            {{"position": 2, "new_aa": "A"}},
            {{"position": 12, "new_aa": "Y"}},
            {{"position": 15, "new_aa": "W"}}
        ]
    }}
    """
    
    try:
        print("Consulting Gemini (Pharmacologist AI) for strategic mutations...")
        model = genai.GenerativeModel('gemini-2.5-flash', generation_config={"response_mime_type": "application/json"})
        response = model.generate_content(prompt)
        
        data = json.loads(response.text)

        log_path = "../../results/pharma_strategy_log.txt"
        thought_process = data.get("thought_process", "No thought process provided.")

        try:
            with open(log_path, "a") as log_file:
                log_file.write(f"=== STRATEGY UPDATE ===\n")
                log_file.write(f"{thought_process}\n\n")
        except Exception as e:
            print(f"Warning: Could not write to log ({e})")
        
        seq_list = list(seq)
        mutated_positions = []
        design = data.get("design", [])
        
        for mut in design[:num_mutations]:
            idx = mut["position"]
            new_aa = mut["new_aa"]
            
            # Safety checks to prevent crashes
            if start_idx <= idx <= end_idx and new_aa in AMINO_ACIDS:
                current_aa = seq_list[idx]
                seq_list[idx] = new_aa
                mutated_positions.append(f"{current_aa}{idx+1}{new_aa}")
                
        print("Successfully generated AI-driven drug mutations!")
        return "".join(seq_list), mutated_positions

    except Exception as e:
        print(f"Gemini API failed ({e}). Falling back to random mutations.")
        return mutate_sequence_random(seq, num_mutations)

def mutate_drug(input_fasta, history_csv, output_fasta, num_mutations):
    # Read the current drug sequence
    with open(input_fasta, "r") as f:
        lines = f.readlines()
        # Join lines ignoring the > header
        drug_seq = "".join([line.strip() for line in lines if not line.startswith(">")])
        
    print(f"Current Drug Sequence: {drug_seq} (Length: {len(drug_seq)})")
    
    # Mutate
    new_seq, mutations = get_agentic_mutations(drug_seq, num_mutations, history_csv)
    mut_string = "_".join(mutations) if mutations else "NO_MUT"
    
    # Write the new FASTA
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    with open(output_fasta, "w") as f:
        f.write(f">pharmacologist_drug | {mut_string}\n{new_seq}\n")
        
    print(f"Saved mutated drug to {output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--current_drug", type=str, required=True, help="Path to current drug FASTA")
    parser.add_argument("--history_csv", type=str, required=True, help="Path to arms race history")
    parser.add_argument("--out_fasta", type=str, required=True, help="Where to save the new drug FASTA")
    parser.add_argument("--num_mutations", type=int, required=True)
    args = parser.parse_args()
    
    mutate_drug(args.current_drug, args.history_csv, args.out_fasta, args.num_mutations)