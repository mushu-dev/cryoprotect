# CryoProtect Package Management

## Overview

This document explains how to manage Python package dependencies for the CryoProtect project, including how to resolve package persistence issues that may occur after system restarts.

## Package Persistence Issue

Some users may experience an issue where Python packages become inaccessible after system restart. This typically happens due to one of the following reasons:

1. **Virtual Environment Path Issues**: The Python interpreter may lose track of the virtual environment path after restart.
2. **Package Installation Location**: Packages might be installed in a temporary location that doesn't persist across restarts.
3. **Environment Variable Changes**: System environment variables (like PYTHONPATH) may change during restart.
4. **Multiple Python Installations**: Conflicts between different Python installations can cause packages to become inaccessible.

## Solution

We've implemented an automated solution to ensure all required packages are properly installed and accessible:

1. **Package Installation Script**: The `install_packages.py` script automatically installs all required packages.
2. **Startup Scripts**: The `startup_packages.bat` (Windows) and `startup_packages.sh` (Unix/Linux/Mac) scripts can be run after system startup to ensure all packages are properly installed.
3. **Package Verification**: The `verify_packages.py` script checks if all required packages are accessible.

## Required Packages

The following packages are required for CryoProtect:

- scipy
- rdkit
- psutil
- xlsxwriter
- seaborn
- scikit-learn (imported as sklearn)
- python-json-logger (imported as pythonjsonlogger.jsonlogger)
- ecs-logger (imported as ecs_logger)
- prometheus_client
- PyYAML (imported as yaml)

## Usage Instructions

### Automatic Package Management

The application's startup scripts (`run_app.bat` and `run_app.sh`) now include automatic package verification and installation:

1. When you run the application using `run_app.bat` (Windows) or `run_app.sh` (Unix/Linux/Mac), the script will:
   - Verify that all required packages are installed and accessible
   - Automatically run the package installation script if any packages are missing
   - Continue with starting the application once all packages are available

This ensures that all necessary packages are always available when running the application, even after system restarts.

### Manual Package Installation

If you prefer to manually install packages:

1. Run the appropriate startup script for your operating system:
   - Windows: Double-click `startup_packages.bat` or run it from the command line
   - Unix/Linux/Mac: Run `./startup_packages.sh` from the terminal (you may need to make it executable with `chmod +x startup_packages.sh`)

2. The script will:
   - Activate the virtual environment if it exists
   - Use the system Python installation if no virtual environment is found
   - Install all required packages
   - Confirm successful installation

### Verifying Package Installation

To verify that all packages are properly installed and accessible:

```bash
python verify_packages.py
```

This will generate a report showing the status of each required package.

### Manual Package Installation

If you prefer to install packages manually, you can:

```bash
# Activate the virtual environment
# Windows
.venv\Scripts\activate
# Unix/Linux/Mac
source .venv/bin/activate

# Install packages
python -m pip install scipy rdkit psutil xlsxwriter seaborn scikit-learn python-json-logger ecs-logger prometheus_client PyYAML
```

## Troubleshooting

If you continue to experience package persistence issues:

1. **Check Python Path**: Verify the Python path using:
   ```python
   import sys
   print(sys.path)
   ```

2. **Verify Virtual Environment**: Ensure you're using the correct virtual environment:
   ```bash
   # Windows
   where python
   # Unix/Linux/Mac
   which python
   ```

3. **Check Package Installation Location**: Verify where packages are installed:
   ```python
   import package_name
   print(package_name.__file__)
   ```

4. **Reinstall Virtual Environment**: If problems persist, consider recreating the virtual environment:
   ```bash
   # Windows
   rmdir /s /q .venv
   python -m venv .venv
   # Unix/Linux/Mac
   rm -rf .venv
   python -m venv .venv
   ```

## Additional Information

For more detailed information about specific packages:

- RDKit: See `README_RDKit_Troubleshooting.md`
- Authentication: See `README_Authentication.md`
- Scoring: See `README_Scoring.md`