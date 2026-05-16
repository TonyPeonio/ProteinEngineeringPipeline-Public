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
    
    # ==========================================
    # DYNAMIC GRID SETUP (Optimized for Slides)
    # ==========================================
    if has_advanced:
        # Create a 2x2 grid, standard 16:9-ish aspect ratio
        fig, axes = plt.subplots(2, 2, figsize=(16, 10))
        fig.suptitle('AI Protein Arms Race: Diagnostic Dashboard', fontsize=20, fontweight='bold', y=0.98)
        
        # Flatten the 2x2 axes matrix into individual variables
        ax1, ax2 = axes[0, 0], axes[0, 1]
        ax3, ax4 = axes[1, 0], axes[1, 1]
    else:
        # Create a 1x2 wide grid if only doing standard tracking
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        fig.suptitle('AI Protein Arms Race: Core Tracking', fontsize=20, fontweight='bold', y=0.98)
        ax1, ax2 = axes[0], axes[1]

    # ==========================================
    # PANEL 1: The PAE Battles (Top Left)
    # ==========================================
    ax1.plot(df['Generation'], df['Drug_PAE'], marker='o', color='#0047AB', linewidth=2.5, label='Pharmacologist (Drug PAE)')
    ax1.plot(df['Generation'], df['Cancer_vs_Drug_PAE'], marker='s', color='#D22B2B', linewidth=2.5, linestyle='--', label='Cancer Evasion')
    ax1.plot(df['Generation'], df['Cancer_vs_p53_PAE'], marker='^', color='#228B22', linewidth=2.5, label='Cancer Function')

    ax1.set_ylabel('Predicted Aligned Error (PAE)', fontsize=12)
    ax1.set_title('Binding Affinities Over Generations', fontsize=14)
    # Legend tucked inside the box to prevent overlapping Panel 2
    ax1.legend(loc='best', fontsize=10) 
    ax1.invert_yaxis()
    ax1.set_ylim(35, 0)
    ax1.grid(True, linestyle=':', alpha=0.7)

    # ==========================================
    # PANEL 2: Overall Evolutionary Fitness (Top Right)
    # ==========================================
    ax2.plot(df['Generation'], df['Fitness_Delta'], marker='D', color='#800080', linewidth=2.5, markersize=8)
    ax2.fill_between(df['Generation'], df['Fitness_Delta'], color='#800080', alpha=0.15)
    ax2.axhline(0, color='black', linewidth=1, linestyle='--') # Baseline

    ax2.set_ylabel('Fitness Delta\n(Evasion - Function)', fontsize=12)
    ax2.set_title('Cancer Evolutionary Fitness Trajectory', fontsize=14)
    ax2.set_xticks(df['Generation'])
    ax2.grid(True, linestyle=':', alpha=0.7)

    # Label X axis on top row only if there is no bottom row
    if not has_advanced:
        ax1.set_xlabel('Generation', fontsize=12)
        ax2.set_xlabel('Generation', fontsize=12)

    # ==========================================
    # ADVANCED PANELS (Bottom Row)
    # ==========================================
    if has_advanced:
        # PANEL 3: Drug Structural Confidence (pLDDT) - Bottom Left
        ax3.plot(df['Generation'], df['Drug_pLDDT'], marker='o', color='#FF8C00', linewidth=2.5)
        ax3.set_ylabel('pLDDT (0-100)', fontsize=12)
        ax3.set_xlabel('Generation', fontsize=12)
        ax3.set_title('Pharmacologist Structural Confidence', fontsize=14)
        ax3.set_ylim(40, 100)
        ax3.grid(True, linestyle=':', alpha=0.7)
        
        # Add a confidence threshold zone
        ax3.axhspan(0, 70, color='red', alpha=0.1)
        ax3.text(df['Generation'].iloc[0], 65, ' Low Confidence Zone', color='red', fontsize=10)

        # PANEL 4: Cumulative Sequence Drift - Bottom Right
        ax4.plot(df['Generation'], df['Mutations_vs_WT'], marker='X', color='#008080', linewidth=2.5, markersize=10)
        ax4.fill_between(df['Generation'], df['Mutations_vs_WT'], color='#008080', alpha=0.15)
        ax4.set_xlabel('Generation', fontsize=12)
        ax4.set_ylabel('Cumulative Amino Acid Changes', fontsize=12)
        ax4.set_title('Sequence Drift (Distance from Wildtype MDM2)', fontsize=14)
        ax4.grid(True, linestyle=':', alpha=0.7)
        # Force integer y-ticks for mutations
        ax4.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Apply tight layout, leaving room at the top for the Suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_image_path, dpi=300, bbox_inches='tight')
    print(f"Success! Advanced Graph saved to {output_image_path}")

if __name__ == "__main__":
    ROOT_DIR = "/home/tonypeonio/ProteinDesignChallenge"
    CSV_FILE = os.path.join(ROOT_DIR, "results", "arms_race_data.csv")
    OUTPUT_IMAGE = os.path.join(ROOT_DIR, "results", "arms_race_dashboard.png")
    
    plot_arms_race(CSV_FILE, OUTPUT_IMAGE)