import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_arms_race(csv_path, output_image_path):
    if not os.path.exists(csv_path):
        print(f"Error: Could not find {csv_path}. Make sure you created the CSV!")
        return

    # Load the data
    df = pd.read_csv(csv_path)

    # Set up the figure with two subplots (stacked vertically)
    plt.style.use('ggplot') # Clean, professional styling
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    fig.suptitle('AI Protein Arms Race: Drug vs. Cancer Target', fontsize=16, fontweight='bold')

    # ==========================================
    # TOP PANEL: The PAE Battles
    # ==========================================
    # Drug catching the target (Lower is better)
    ax1.plot(df['Generation'], df['Drug_PAE'], marker='o', color='blue', linewidth=2, label='Pharmacologist (Drug PAE)')
    
    # Cancer evading the drug (Higher is better for cancer)
    ax1.plot(df['Generation'], df['Cancer_vs_Drug_PAE'], marker='s', color='red', linewidth=2, linestyle='--', label='Cancer Evasion (Mutant vs Drug PAE)')
    
    # Cancer maintaining function (Lower is better for cancer)
    ax1.plot(df['Generation'], df['Cancer_vs_p53_PAE'], marker='^', color='green', linewidth=2, label='Cancer Function (Mutant vs Native p53 PAE)')

    ax1.set_ylabel('Predicted Aligned Error (PAE)', fontsize=12)
    ax1.set_title('Binding Affinities Over Generations', fontsize=14)
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.invert_yaxis() # Invert PAE so "better binding" is visually higher on the graph for Drug/p53
    ax1.set_ylim(35, 0) # Standard ColabFold PAE max is roughly 31. Inverting 35 to 0.

    # ==========================================
    # BOTTOM PANEL: The Fitness Delta
    # ==========================================
    # This is the cancer's ultimate objective function
    ax2.plot(df['Generation'], df['Fitness_Delta'], marker='D', color='purple', linewidth=2, markersize=8)
    
    # Fill the area under the curve for visual impact
    ax2.fill_between(df['Generation'], df['Fitness_Delta'], color='purple', alpha=0.2)

    ax2.set_xlabel('Generation', fontsize=12)
    ax2.set_ylabel('Fitness Delta\n(Evasion - Function)', fontsize=12)
    ax2.set_title('Overall Cancer Evolutionary Fitness', fontsize=14)
    ax2.set_xticks(df['Generation'])

    # Clean up layout and save
    plt.tight_layout()
    plt.savefig(output_image_path, dpi=300, bbox_inches='tight')
    print(f"Success! Graph saved to {output_image_path}")

if __name__ == "__main__":
    # Define paths
    ROOT_DIR = "/home/tonypeonio/ProteinDesignChallenge"
    CSV_FILE = os.path.join(ROOT_DIR, "results", "arms_race_data.csv")
    OUTPUT_IMAGE = os.path.join(ROOT_DIR, "results", "arms_race_plot.png")
    
    plot_arms_race(CSV_FILE, OUTPUT_IMAGE)