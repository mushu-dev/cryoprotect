#!/bin/bash
# RDKit Property Calculator wrapper script
# Provides a convenient way to run the property calculator

set -e  # Exit on error

# Check for Python environment
if command -v conda &> /dev/null; then
    # If using conda, activate the environment if it exists
    if conda env list | grep -q "cryoprotect"; then
        echo "Activating cryoprotect conda environment..."
        conda activate cryoprotect
    fi
fi

# Display help if requested
if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "RDKit Property Calculator for CryoProtect"
    echo ""
    echo "Usage: ./run_rdkit_calculator.sh [options]"
    echo ""
    echo "Options:"
    echo "  --known              Process only known cryoprotectants"
    echo "  --sample N           Process a random sample of N molecules"
    echo "  --limit N            Limit processing to N molecules"
    echo "  --url URL            Supabase URL (overrides environment variable)"
    echo "  --key KEY            Supabase service role key (overrides environment variable)"
    echo "  --help, -h           Show this help message"
    echo ""
    echo "Example: ./run_rdkit_calculator.sh --known --limit 10"
    exit 0
fi

# Check if .env file exists
if [[ -f .env ]]; then
    echo "Using credentials from .env file"
else
    echo "Note: No .env file found. You'll be prompted for credentials if needed."
    echo "You can create a .env file based on .env-example to avoid being prompted."
fi

# Check if RDKit is installed
if python -c "import rdkit" 2>/dev/null; then
    echo "RDKit is installed and available."
else
    echo "Warning: RDKit is not installed. Will use mock implementation."
    echo "For more accurate results, install RDKit with: conda install -c conda-forge rdkit"
fi

# Run the calculator
echo "Starting RDKit Property Calculator..."
python rdkit_property_calculator.py "$@"