# 🧬 Protein Engineering Pipeline Quick-Reference

- [🧬 Protein Engineering Pipeline Quick-Reference](#-protein-engineering-pipeline-quick-reference)
    - [1. RFdiffusion (Backbone Generation)](#1-rfdiffusion-backbone-generation)
    - [2. ProteinMPNN (Sequence Design)](#2-proteinmpnn-sequence-design)
    - [3. Protein Folding](#3-protein-folding)
    - [4. PyMOL Structural Analysis](#4-pymol-structural-analysis)

### 1. RFdiffusion (Backbone Generation)

**Environment Setup:**

```bash
conda activate rfdiffusion
```

**Iteration 1 (Shorter Binder, Site 1):**

```bash
python RFdiffusion/scripts/run_inference.py \
    inference.output_prefix=outputs/rfdiffusion/nitrogenase_binder \
    'contigmap.contigs=[A4-480/0 50-80]' \
    inference.input_pdb=pdb/cleaned_pdb/nitrogenase_clean.pdb \
    inference.num_designs=24 \
    ppi.hotspot_res=[A356,A357,A360,A361,A364,A367,A368]
```

**Iteration 2 (Longer Binder, Site 2):**

```bash
python RFdiffusion/scripts/run_inference.py \
    inference.output_prefix=outputs/rfdiffusion/nitrogenase_binder \
    'contigmap.contigs=[A4-480/0 60-100]' \
    inference.input_pdb=pdb/cleaned_pdb/nitrogenase_clean.pdb \
    inference.num_designs=100 \
    ppi.hotspot_res=[A310,A313,A314,A318,A319,A322,A325,A326] \
    hydra.run.dir=outputs/rfdiffusion/.logs
```

---

### 2. ProteinMPNN (Sequence Design)

Choose **one** of the following execution methods:

* **Option A: Sequential (Run after RFdiffusion finishes)**

```bash
bash scripts/run_mpnn.sh
```

* **Option B: Parallel (Run simultaneously with RFdiffusion)**

```bash
python scripts/watchdog_mpnn.py
```

> **Sequence Processing & Filtering: Run all cells in mpnn_sequence_filtering.ipynb**

**Clean and Move the Output:**

```bash
# Clean the Output and Put Into outputs/esmfold_inputs/cleaned:
python scripts/folding_prep.py
```

---

### 3. Protein Folding

**ESMFold (Fast Folding for Funneling):**

```bash
conda deactivate
conda activate esmfold_hf
```

```bash
#1. Run batch_esmfold.py
python scripts/batch_esmfold.py
```

> **Sequence Processing & Filtering: Run all cells in esmfold_sequence_filtering.ipynb**

**ColabFold / AlphaFold2 Multimer (High-Accuracy):**

```bash
# 1. Activate environment
source localcolabfold/conda/bin/activate localcolabfold/colabfold-conda

# 2. Prevent GPU memory allocation errors
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

# 3. Run ColabFold in the background
colabfold_batch outputs/colabfold_inputs/ outputs/colabfold_results/ \
    --model-type alphafold2_multimer_v3 \
    --num-recycle 12 \
    --num-models 5 \
    --recycle-early-stop-tolerance 0.5 &
```

### 4. PyMOL Structural Analysis

**Exporting Files from WSL:**

1. Open the file explorer in your current WSL directory:

```bash
explorer.exe .
```

2. Move the generated `.pdb` files to a Windows-accessible folder (e.g., Downloads).

**PyMOL Alignment Workflow:**

1. Import the native nitrogenase target (`3u7q` or `7ut8`).
2. Import your designed binder `.pdb`.
3. Run the following command in the PyMOL console to align the target chains *(replace the `#` with your specific binder numbers)*:

```pymol
align nitrogenase_binder_##_sample_#_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000 and chain A, 7ut8 and chain A
```
