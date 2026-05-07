# Quick Reference: Basic Commands

This cheat sheet covers the standard commands needed to navigate the environment, version control, and background processing for this pipeline. 

## 1. Conda / Mamba (Environment Management)

*Note: I highly recommend using `mamba` for faster dependency solving, but `conda` can be used interchangeably for all commands below.*

**Create a new environment from a YAML file (provided in scripts_and_helpers/envs):**
```bash
mamba env create -f environment.yml
```

**Activate the environment:**
```bash
mamba activate protein-design
```

**Update an existing environment (if the YAML file changes):**
```bash
mamba env update -n protein-design -f environment.yml --prune
```

**Deactivate the current environment:**
```bash
mamba deactivate
```

---

## 2. Git (Code Version Control)

**Pushing updates to GitHub:**
```bash
git add .                             # Stage all modified files
git commit -m "Briefly explain edits" # Commit the changes
git push origin main                  # Push to remote (may be 'master' on older repos)
```

**Pulling updates from GitHub:**
```bash
git pull origin main
```

**Fixing a blocked pull:**
*(If Git blocks your pull because you have unsaved local changes, you can "stash" them aside temporarily).*
```bash
git stash
git pull origin main
```

---

## 3. DVC (Data Version Control)

*DVC handles the heavy structural files (PDBs, model weights) so Git doesn't crash. Never commit raw outputs directly to Git.*

**Tracking new heavy data:**
```bash
# 1. Tell DVC to track the directory (moves data to cache, creates a .dvc file)
dvc add outputs/ old_results/ pyrosetta_outputs/

# 2. Tell Git to track the lightweight pointer file and ignore the heavy data
git add outputs.dvc old_results.dvc pyrosetta_outputs.dvc .gitignore

# 3. Commit the pointer to Git
git commit -m "Tracked new outputs batch via DVC"
```

**Syncing data:**
```bash
dvc push   # Uploads heavy files to remote storage (if configured)
dvc pull   # Downloads heavy files associated with your current Git commit
```

---

## 4. Tmux (Background Processing)

*Tmux allows your scripts to keep running on the server even if you close your SSH connection or terminal window.*

**Start a new named session (Highly Recommended):**
```bash
tmux new -s my_pipeline
```

**Detach from the session (Leave it running in the background):**
Press `Ctrl + B`, release both keys, then press `D`.

**View all running sessions:**
```bash
tmux ls
```

**Reattach to a session:**
```bash
tmux attach -t my_pipeline
```