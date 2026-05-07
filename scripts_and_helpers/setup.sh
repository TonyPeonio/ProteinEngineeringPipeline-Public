#!/bin/bash
set -e 

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
ENVS_DIR="$SCRIPT_DIR/envs"

echo "Starting automated pipeline setup in $PROJECT_ROOT..."
eval "$(conda shell.bash hook)"

# 1. Version Control
cd "$PROJECT_ROOT"
[ ! -d ".git" ] && git init
if ! command -v dvc &> /dev/null; then pip install dvc; fi
[ ! -d ".dvc" ] && dvc init && git commit -m "Initialize DVC" || true

# 2. Clone Repositories
[ ! -d "RFdiffusion" ] && git clone https://github.com/RosettaCommons/RFdiffusion.git
[ ! -d "ProteinMPNN" ] && git clone https://github.com/dauparas/ProteinMPNN.git
if [ ! -d "localcolabfold" ]; then
    mkdir -p localcolabfold
    cd localcolabfold
    wget -qN https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
    bash install_colabbatch_linux.sh
    cd ..
fi

# 3. Setup Basic Environments
echo "Setting up basic environments..."
[ -f "$ENVS_DIR/env_mpnn.yml" ] && conda env create -f "$ENVS_DIR/env_mpnn.yml" || true
[ -f "$ENVS_DIR/env_rosetta.yml" ] && conda env create -f "$ENVS_DIR/env_rosetta.yml" || true
[ -f "$ENVS_DIR/env_colabfold.yml" ] && conda env create -f "$ENVS_DIR/env_colabfold.yml" || true

# 4. Download RFdiffusion Weights
echo "Downloading RFdiffusion model weights..."
cd "$PROJECT_ROOT/RFdiffusion"
bash scripts/download_models.sh models

# 5. Create Base RFdiffusion Environment
echo "Creating base Python environment for RFdiffusion..."
conda create -n env_rfdiffusion python=3.9 -y

echo ""
echo "Automated setup complete! Please check the README.md to finish configuring the RFdiffusion environment for your specific hardware."