import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_arms_race(csv_path, output_image_path):
    if not os.path.exists(csv_path):
        print(f"Error: Could not find {csv_path}.")
        return

    # Load the data
    df = pd.read_csv(csv_path)

    # Check if we have the advanced diagnostic columns
    has_advanced = 'Drug_pLDDT' in df.columns and 'Mutations_vs_WT' in df.columns
    
    # Dynamically set up 2 or 4 subplots
    num_panels = 4 if has_advanced else 2
    fig, axes = plt.subplots(num_panels, 1, figsize=(12, 5 * num_panels), sharex=True)
    fig.suptitle('AI Protein Arms Race: Diagnostic Dashboard', fontsize=18, fontweight='bold', y=0.98)
    
    # Handle array of axes consistently
    if num_panels == 2:
        ax1, ax2 = axes
    else:
        ax1, ax2, ax3, ax4 = axes

    # ==========================================
    # PANEL 1: The PAE Battles
    # ==========================================
    ax1.plot(df['Generation'], df['Drug_PAE'], marker='o', color='#0047AB', linewidth=2.5, label='Pharmacologist (Drug PAE)')
    ax1.plot(df['Generation'], df['Cancer_vs_Drug_PAE'], marker='s', color='#D22B2B', linewidth=2.5, linestyle='--', label='Cancer Evasion (Mutant vs Drug)')
    ax1.plot(df['Generation'], df['Cancer_vs_p53_PAE'], marker='^', color='#228B22', linewidth=2.5, label='Cancer Function (Mutant vs Native p53)')

    ax1.set_ylabel('Predicted Aligned Error (PAE)', fontsize=12)
    ax1.set_title('Binding Affinities Over Generations', fontsize=14)
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.invert_yaxis()
    ax1.set_ylim(35, 0)
    ax1.grid(True, linestyle=':', alpha=0.7)

    # ==========================================
    # PANEL 2: Overall Evolutionary Fitness
    # ==========================================
    ax2.plot(df['Generation'], df['Fitness_Delta'], marker='D', color='#800080', linewidth=2.5, markersize=8)
    ax2.fill_between(df['Generation'], df['Fitness_Delta'], color='#800080', alpha=0.15)
    ax2.axhline(0, color='black', linewidth=1, linestyle='--') # Baseline

    if not has_advanced:
        ax2.set_xlabel('Generation', fontsize=12)
        
    ax2.set_ylabel('Fitness Delta\n(Evasion - Function)', fontsize=12)
    ax2.set_title('Cancer Evolutionary Fitness Trajectory', fontsize=14)
    ax2.set_xticks(df['Generation'])
    ax2.grid(True, linestyle=':', alpha=0.7)

    # ==========================================
    # ADVANCED PANELS (If data is available)
    # ==========================================
    if has_advanced:
        # PANEL 3: Drug Structural Confidence (pLDDT)
        ax3.plot(df['Generation'], df['Drug_pLDDT'], marker='o', color='#FF8C00', linewidth=2.5)
        ax3.set_ylabel('pLDDT (0-100)', fontsize=12)
        ax3.set_title('Pharmacologist Structural Confidence', fontsize=14)
        ax3.set_ylim(40, 100)
        ax3.grid(True, linestyle=':', alpha=0.7)
        
        # Add a confidence threshold zone
        ax3.axhspan(0, 70, color='red', alpha=0.1)
        ax3.text(df['Generation'].iloc[0], 65, ' Low Confidence Zone', color='red', fontsize=10)

        # PANEL 4: Cumulative Sequence Drift
        ax4.plot(df['Generation'], df['Mutations_vs_WT'], marker='X', color='#008080', linewidth=2.5, markersize=10)
        ax4.fill_between(df['Generation'], df['Mutations_vs_WT'], color='#008080', alpha=0.15)
        ax4.set_xlabel('Generation', fontsize=12)
        ax4.set_ylabel('Cumulative Amino Acid Changes', fontsize=12)
        ax4.set_title('Sequence Drift (Distance from Wildtype MDM2)', fontsize=14)
        ax4.grid(True, linestyle=':', alpha=0.7)
        # Force integer y-ticks for mutations
        ax4.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    plt.tight_layout()
    plt.savefig(output_image_path, dpi=300, bbox_inches='tight')
    print(f"Success! Advanced Graph saved to {output_image_path}")

if __name__ == "__main__":
    ROOT_DIR = "/home/tonypeonio/ProteinDesignChallenge"
    CSV_FILE = os.path.join(ROOT_DIR, "results", "arms_race_data.csv")
    OUTPUT_IMAGE = os.path.join(ROOT_DIR, "results", "arms_race_dashboard.png")
    
    plot_arms_race(CSV_FILE, OUTPUT_IMAGE)