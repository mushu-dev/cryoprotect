#!/bin/bash
# CryoProtect v2 - Robust Server Startup Script for Linux (Conda + RDKit)

# Activate the conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cryoprotect

# Diagnostic: Print Python executable and installed packages
python --version
which python
python -m pip list

# Set RDKit base path if needed (usually not required if installed via conda)
# export RDBASE=$CONDA_PREFIX/share/rdkit

# Set Flask app entry point
export FLASK_APP=app.py

# Optionally set Flask environment (development/production)
export FLASK_ENV=development

# Start the Flask server (development)
python -m flask run --host=0.0.0.0 --port=5000

# For production, use Waitress (uncomment below after installing waitress)
# waitress-serve --port=5000 app:app