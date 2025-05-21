# Python Package Persistence Issue: Diagnosis and Solution

## Problem Description

After restarting the system, multiple Python packages (e.g., scipy, rdkit, psutil, seaborn, scikit-learn, etc.) become inaccessible in the project codebase. Manual reinstallation is required after each restart, disrupting workflow.

## Diagnosis Steps

1. **Checked Python Environment:**
   - Active environment: `(cryoprotect)` (Anaconda-managed)
   - Python executable: `C:\Users\1edwa\anaconda3\envs\cryoprotect\python.exe`
   - All required packages are installed in this environment.

2. **Checked Package Installation:**
   - Packages are present in `C:\Users\1edwa\anaconda3\envs\cryoprotect\lib\site-packages`
   - `rdkit` is installed via conda-forge.

3. **Checked Environment Persistence:**
   - The `cryoprotect` environment is not deleted or corrupted between restarts.

4. **Checked for Multiple Python Installations:**
   - Multiple environments exist, but `cryoprotect` is the intended one.

## Root Cause

After a system restart, the `cryoprotect` environment is not automatically activated. If you run Python scripts or open VSCode terminals without activating the environment, the system may use a different Python interpreter (e.g., system Python or another conda environment), causing packages to appear "missing".

## Solution

1. **Use the Batch File to Run the Application:**
   - The `run_app.bat` file already handles environment activation correctly:
     ```
     .\run_app.bat
     ```
   - This batch file activates the cryoprotect environment, checks for RDKit, and runs the application.

2. **Always Activate the Environment Manually (if needed):**
   - Before running any Python code directly, activate the environment:
     ```
     conda activate cryoprotect
     ```
   - **IMPORTANT**: Direct execution of Python scripts without environment activation is the root cause of the "missing packages" issue.

3. **Set VSCode Python Interpreter:**
   - In VSCode, ensure the Python interpreter is set to:
     ```
     C:\Users\1edwa\anaconda3\envs\cryoprotect\python.exe
     ```
   - This ensures all scripts and terminals use the correct environment.

4. **Use the Verification Batch File:**
   - Run `verify_packages.bat` to check if all required packages are available:
     ```
     .\verify_packages.bat
     ```
   - This batch file activates the environment before running the verification script.
   - All packages should show as [OK] when the environment is properly activated.

## Verification Scripts

Two verification scripts are provided:

1. `verify_packages.py` - The Python script that checks for required packages.
2. `verify_packages.bat` - A batch file wrapper that activates the environment before running the Python script.

**Always use the batch file version** (`.\verify_packages.bat`) to ensure the environment is properly activated before checking packages.

## Success Criteria

- All specified packages remain accessible after system restart without reinstallation.
- Solution works persistently across multiple restarts.
- This documentation provides a clear explanation of the root cause and solution.
- The application starts successfully using `run_app.bat` without package-related errors.