# Environment Setup and SciPy Import Error Fix

## Issue

We're experiencing SciPy import errors when running tests. The issue stems from using the virtual environment (`.venv`) for running tests instead of the Conda environment where SciPy is properly installed.

## Analysis

1. Our project has two Python environments:
   - A Conda environment (`cryoprotect`) with scientific packages including scipy
   - A `.venv` virtual environment created by some tools automatically

2. The test runner is using the `.venv` Python, which doesn't have scipy installed properly.

## Solution: Setup Test Environment Properly

### Option 1: Fix the `.venv` Environment (Preferred)

1. Install scipy directly into the `.venv` environment:

```bash
# Activate the virtual environment
source .venv/bin/activate  # On Linux/Mac
# OR
.venv\Scripts\activate  # On Windows

# Install scipy
pip install scipy

# Verify installation
python -c "import scipy; print(scipy.__version__)"

# Run the tests to confirm it works
python tests/run_tests.py
```

### Option 2: Ensure Tests Use Conda Environment

1. Create a test runner script that explicitly uses the Conda environment:

```bash
# Create a run_tests_conda.bat file for Windows
echo @echo off > run_tests_conda.bat
echo call conda activate cryoprotect >> run_tests_conda.bat
echo python tests/run_tests.py %* >> run_tests_conda.bat

# Create a run_tests_conda.sh file for Linux/Mac
echo '#!/bin/bash' > run_tests_conda.sh
echo 'eval "$(conda shell.bash hook)"' >> run_tests_conda.sh
echo 'conda activate cryoprotect' >> run_tests_conda.sh
echo 'python tests/run_tests.py "$@"' >> run_tests_conda.sh
chmod +x run_tests_conda.sh

# Then run tests using:
# ./run_tests_conda.sh  # On Linux/Mac
# run_tests_conda.bat  # On Windows
```

### Option 3: Consolidate to One Environment Approach

For longer-term solution, standardize on either:

1. **Conda-only**: Remove `.venv` and ensure all team members activate the conda environment before working
2. **Pip-only**: Export conda packages to requirements.txt and use only pip with virtual environment

## Implementation Plan

1. Try Option 1 first since it's simple and non-disruptive
2. If that doesn't work, implement Option 2 as a failsafe
3. Consider Option 3 for long-term standardization once current issues are resolved

## Update CLAUDE.md

After fixing the environment issues, update CLAUDE.md to reflect the correct environment setup instructions:

```markdown
## Environment Setup

- **Conda environment**: `conda env create -f environment.yml`  
- **Activate environment**: `conda activate cryoprotect`
- **Alternative pip environment**: If using virtual environment instead of conda, run:
  ```
  python -m venv .venv
  source .venv/bin/activate  # Linux/Mac
  .venv\Scripts\activate  # Windows
  pip install -r requirements.txt
  ```

# Running Tests
- **With conda**: `python tests/run_tests.py` (ensure conda environment is activated)
- **With venv**: `python tests/run_tests.py` (ensure venv is activated and scipy is installed)
```