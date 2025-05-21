#!/bin/bash
# fix_conda_environment.sh - Script to fix conda environment issues for CryoProtect

set -e

echo "=== CryoProtect Environment Fix Tool ==="
echo "This script will fix issues with the conda environment for CryoProtect."

# Function to log messages
log() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] $1"
}

# Function to log success messages
log_success() {
    echo -e "[$(date +"%Y-%m-%d %H:%M:%S")] \033[32m$1\033[0m"
}

# Function to log warning messages
log_warning() {
    echo -e "[$(date +"%Y-%m-%d %H:%M:%S")] \033[33m$1\033[0m"
}

# Function to log error messages
log_error() {
    echo -e "[$(date +"%Y-%m-%d %H:%M:%S")] \033[31m$1\033[0m"
}

# Check if conda is properly set up
log "Checking conda installation..."
if [ ! -f ~/miniconda3/bin/conda ]; then
    log_error "Conda not found at expected location ~/miniconda3/bin/conda"
    log "Please install Miniconda first. Visit: https://docs.conda.io/projects/miniconda/en/latest/"
    exit 1
fi

# Initialize conda for shell
log "Initializing conda for current shell..."
eval "$(~/miniconda3/bin/conda shell.bash hook)"
if [ $? -ne 0 ]; then
    log_error "Failed to initialize conda for shell"
    exit 1
fi
log_success "Conda initialized successfully"

# Check existing environments
log "Checking existing conda environments..."
if conda env list | grep -q "cryoprotect"; then
    log_warning "Existing 'cryoprotect' environment found"
    read -p "Do you want to remove and recreate the environment? (y/n): " recreate_env
    if [[ "$recreate_env" =~ ^[Yy]$ ]]; then
        log "Removing existing environment..."
        conda env remove -n cryoprotect
        log_success "Environment removed"
    else
        log "Keeping existing environment"
    fi
fi

# Create/update environment from file
log "Creating/updating conda environment from environment.yml..."
conda env update -f environment.yml --prune
if [ $? -ne 0 ]; then
    log_error "Failed to create/update environment from environment.yml"
    log "Trying alternative approach with mamba..."
    
    # Try to install mamba if it fails
    conda install -y -c conda-forge mamba
    if [ $? -ne 0 ]; then
        log_error "Failed to install mamba. Please try manual installation."
        exit 1
    fi
    
    # Try with mamba
    mamba env update -f environment.yml --prune
    if [ $? -ne 0 ]; then
        log_error "Failed to create environment with mamba."
        exit 1
    fi
fi
log_success "Environment created/updated successfully"

# Activate the environment
log "Activating 'cryoprotect' environment..."
conda activate cryoprotect
if [ $? -ne 0 ]; then
    log_error "Failed to activate 'cryoprotect' environment"
    exit 1
fi
log_success "Environment activated"

# Verify RDKit installation
log "Verifying RDKit installation..."
python -c "from rdkit import Chem; print('RDKit verification successful')"
if [ $? -ne 0 ]; then
    log_error "RDKit verification failed"
    log "Attempting to reinstall RDKit..."
    conda install -y -c conda-forge rdkit=2024.03.4
    if [ $? -ne 0 ]; then
        log_error "Failed to reinstall RDKit. Please check the conda-forge channel."
        exit 1
    fi
    
    # Verify again
    python -c "from rdkit import Chem; print('RDKit verification successful')"
    if [ $? -ne 0 ]; then
        log_error "RDKit still not working after reinstallation."
        exit 1
    fi
fi
log_success "RDKit verified successfully"

# Install pip packages
log "Installing pip packages from requirements.txt..."
pip install -r requirements.txt
if [ $? -ne 0 ]; then
    log_error "Failed to install pip packages"
    exit 1
fi
log_success "Pip packages installed successfully"

# Create shell activation script
log "Creating activation helper script..."
cat > activate_cryoprotect.sh << 'EOL'
#!/bin/bash
# Helper script to activate the cryoprotect environment

# Store current directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Initialize conda
eval "$(~/miniconda3/bin/conda shell.bash hook)"

# Activate the environment
conda activate cryoprotect

# Print success message
echo -e "\033[32mCryoProtect environment activated successfully!\033[0m"
echo "You are now using Python: $(which python)"
echo "Run the application with: ./run_app.sh"
echo "Run tests with: ./run_tests.sh"
echo "Deactivate environment with: conda deactivate"
EOL

chmod +x activate_cryoprotect.sh
log_success "Created activate_cryoprotect.sh script"

# Fix run_app.sh if it exists
if [ -f run_app.sh ]; then
    log "Fixing run_app.sh script..."
    # Make backup
    cp run_app.sh run_app.sh.bak
    
    # Modify the script
    sed -i 's/source $(conda info --base)\/etc\/profile.d\/conda.sh/eval "$(~/miniconda3\/bin\/conda shell.bash hook)"/g' run_app.sh
    
    chmod +x run_app.sh
    log_success "Fixed run_app.sh script (backup at run_app.sh.bak)"
fi

# Summary
log_success "==========================================="
log_success "Environment setup completed successfully!"
log_success "==========================================="
echo ""
echo "To use the environment:"
echo "1. Run: source ./activate_cryoprotect.sh"
echo "2. Run application: ./run_app.sh"
echo ""
echo "Current Python path: $(which python)"
echo "Current environment: $(conda info --envs | grep '*' || echo 'None')"
echo ""
echo "If you continue to have issues, please consult the project documentation"
echo "or run 'conda init bash' and restart your shell before trying again."