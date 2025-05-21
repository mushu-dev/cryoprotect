#!/bin/bash
# Setup mock RDKit for testing when RDKit is not available

set -e

echo "Creating mock RDKit module for testing..."
python /home/mushu/Projects/CryoProtect/mock_rdkit.py

# Get the mock module path
MOCK_PATH=$(python -c "import sys; print(f'/tmp/mock_modules:\{\":\".\".join(sys.path)}')")

echo "Mock RDKit created successfully."
echo
echo "To use the mock RDKit in the container, run:"
echo "./run_in_cryoprotect.sh \"export PYTHONPATH= && python your_script.py\""
