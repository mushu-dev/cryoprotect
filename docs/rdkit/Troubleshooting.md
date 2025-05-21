# CryoProtect Analyzer - RDKit Installation Troubleshooting Guide

This guide provides solutions for common RDKit installation issues and alternative installation methods to ensure the CryoProtect Analyzer application runs properly.

## Common Installation Issues

### Conda Solver Getting Stuck

The most common issue when installing RDKit via conda is the solver getting stuck at the "Solving environment" step:

```
conda install -c conda-forge rdkit
Retrieving notices: ...working... done
Collecting package metadata (current_repodata.json): done
Solving environment: unsuccessful initial attempt using frozen solve. Retrying with flexible solve.
Solving environment: unsuccessful attempt using repodata from current_repodata.json, retrying with next repodata source.
Collecting package metadata (repodata.json): done
Solving environment: \|
```

This happens because:
1. RDKit has many dependencies that conda needs to resolve
2. Conflicts with existing packages in your environment
3. Network issues when retrieving package metadata

## Alternative Installation Methods

### Method 1: Clean Conda Environment (Recommended)

Creating a fresh conda environment is the most reliable approach:

```bash
# Create a new conda environment with Python 3.9
conda create -n cryoprotect python=3.9
conda activate cryoprotect

# Install RDKit with specific version
conda install -c conda-forge rdkit=2023.9.1

# Install other dependencies
pip install -r requirements.txt
```

### Method 2: Conda with Mamba Solver

Mamba is a faster drop-in replacement for conda that often resolves dependency issues more efficiently:

```bash
# Install mamba in base environment
conda install -c conda-forge mamba

# Create environment and install RDKit using mamba
mamba create -n cryoprotect python=3.9
mamba activate cryoprotect
mamba install -c conda-forge rdkit

# Install other dependencies
pip install -r requirements.txt
```

### Method 3: Pip Installation

For some systems, pip installation may work better:

```bash
# Create a virtual environment
python -m venv cryoprotect-env
source cryoprotect-env/bin/activate  # On Windows: cryoprotect-env\Scripts\activate

# Install RDKit via pip
pip install rdkit

# Install other dependencies
pip install -r requirements.txt
```

### Method 4: Docker Container

Using a pre-built Docker container with RDKit already installed:

```bash
# Pull the RDKit container
docker pull informaticsmatters/rdkit-python3-debian

# Run the container with the CryoProtect directory mounted
docker run -it --rm -v /path/to/cryoprotect:/cryoprotect -p 5000:5000 informaticsmatters/rdkit-python3-debian

# Inside the container, install requirements and run the app
cd /cryoprotect
pip install -r requirements.txt
python app.py
```

## Comprehensive Environment Setup Script

Save the following script as `setup_environment.bat` (Windows) or `setup_environment.sh` (Linux/Mac) in your project directory:

### Windows (setup_environment.bat)

```batch
@echo off
echo Setting up CryoProtect Analyzer environment...

:: Try to create a new conda environment
call conda create -n cryoprotect python=3.9 -y
if %ERRORLEVEL% neq 0 (
    echo Failed to create conda environment. Trying with mamba...
    call conda install -c conda-forge mamba -y
    call mamba create -n cryoprotect python=3.9 -y
    if %ERRORLEVEL% neq 0 (
        echo Failed to create environment with mamba. Exiting.
        exit /b 1
    )
)

:: Activate the environment
call conda activate cryoprotect

:: Try to install RDKit with conda
echo Installing RDKit...
call conda install -c conda-forge rdkit=2023.9.1 -y
if %ERRORLEVEL% neq 0 (
    echo Failed to install RDKit with conda. Trying with mamba...
    call mamba install -c conda-forge rdkit=2023.9.1 -y
    if %ERRORLEVEL% neq 0 (
        echo Failed to install RDKit with mamba. Trying with pip...
        call pip install rdkit
        if %ERRORLEVEL% neq 0 (
            echo Failed to install RDKit with pip. Please try manual installation.
            exit /b 1
        )
    )
)

:: Install other dependencies
echo Installing other dependencies...
call pip install -r requirements.txt

:: Verify RDKit installation
echo Verifying RDKit installation...
python -c "from rdkit import Chem; print('RDKit installation successful!')"
if %ERRORLEVEL% neq 0 (
    echo RDKit verification failed. Please check the installation.
    exit /b 1
)

echo Environment setup complete! Activate with: conda activate cryoprotect
```

### Linux/Mac (setup_environment.sh)

```bash
#!/bin/bash
echo "Setting up CryoProtect Analyzer environment..."

# Try to create a new conda environment
conda create -n cryoprotect python=3.9 -y
if [ $? -ne 0 ]; then
    echo "Failed to create conda environment. Trying with mamba..."
    conda install -c conda-forge mamba -y
    mamba create -n cryoprotect python=3.9 -y
    if [ $? -ne 0 ]; then
        echo "Failed to create environment with mamba. Exiting."
        exit 1
    fi
fi

# Activate the environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cryoprotect

# Try to install RDKit with conda
echo "Installing RDKit..."
conda install -c conda-forge rdkit=2023.9.1 -y
if [ $? -ne 0 ]; then
    echo "Failed to install RDKit with conda. Trying with mamba..."
    mamba install -c conda-forge rdkit=2023.9.1 -y
    if [ $? -ne 0 ]; then
        echo "Failed to install RDKit with mamba. Trying with pip..."
        pip install rdkit
        if [ $? -ne 0 ]; then
            echo "Failed to install RDKit with pip. Please try manual installation."
            exit 1
        fi
    fi
fi

# Install other dependencies
echo "Installing other dependencies..."
pip install -r requirements.txt

# Verify RDKit installation
echo "Verifying RDKit installation..."
python -c "from rdkit import Chem; print('RDKit installation successful!')"
if [ $? -ne 0 ]; then
    echo "RDKit verification failed. Please check the installation."
    exit 1
fi

echo "Environment setup complete! Activate with: conda activate cryoprotect"
```

## Application Launch Troubleshooting

If the application fails to launch even after successful RDKit installation, try the following:

### 1. Check for Import Errors

Run the following command to check for import errors:

```bash
python -c "from rdkit import Chem; print('RDKit import successful')"
```

### 2. Debug Flask Application

Create a debug script (`debug_app.py`) to identify Flask application issues:

```python
#!/usr/bin/env python3
import os
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG)

# Import Flask app with error handling
try:
    from app import create_app
    app = create_app()
    print("App created successfully!")
except Exception as e:
    logging.exception("Error creating Flask app:")
    raise

# Try to run the app
if __name__ == '__main__':
    try:
        port = int(os.environ.get('PORT', 5000))
        app.run(host='0.0.0.0', port=port, debug=True)
    except Exception as e:
        logging.exception("Error running Flask app:")
        raise
```

Run with: `python debug_app.py`

### 3. Check Environment Variables

Ensure all required environment variables are set in your `.env` file:

```
FLASK_APP=app.py
FLASK_ENV=development
SUPABASE_URL=your_supabase_url
SUPABASE_KEY=your_supabase_key
```

## Dependency Management

### Updated requirements.txt with Pinned Versions

Create an updated `requirements.txt` with pinned versions that are known to work together:

```
# CryoProtect Analyzer API Requirements

# Flask and extensions
Flask==2.3.3
Flask-RESTful==0.3.10
Flask-Cors==4.0.0
flask-apispec==0.11.4

# API documentation
apispec==6.3.0
marshmallow==3.20.1

# Supabase
supabase==2.0.3

# Environment variables
python-dotenv==1.0.0

# Date handling
python-dateutil==2.8.2

# Utilities
requests==2.31.0
numpy==1.26.0
Pillow==10.0.1

# RDKit - Note: Install via conda with: conda install -c conda-forge rdkit=2023.9.1
# rdkit==2023.9.1
```

### Conda Environment File

Create an `environment.yml` file for more reliable conda environment creation:

```yaml
name: cryoprotect
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.9
  - rdkit=2023.9.1
  - pip=23.1.2
  - pip:
    - flask==2.3.3
    - flask-restful==0.3.10
    - flask-cors==4.0.0
    - flask-apispec==0.11.4
    - apispec==6.3.0
    - marshmallow==3.20.1
    - supabase==2.0.3
    - python-dotenv==1.0.0
    - python-dateutil==2.8.2
    - requests==2.31.0
    - numpy==1.26.0
    - pillow==10.0.1
```

Use with: `conda env create -f environment.yml`

## Application Launch Script

Create a reliable application launch script (`run_app.bat` for Windows or `run_app.sh` for Linux/Mac):

### Windows (run_app.bat)

```batch
@echo off
echo Starting CryoProtect Analyzer...

:: Activate the conda environment
call conda activate cryoprotect

:: Check if RDKit is installed
python -c "from rdkit import Chem; print('RDKit check passed')" > nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo ERROR: RDKit is not installed or not working properly.
    echo Please run setup_environment.bat first.
    exit /b 1
)

:: Run the Flask application with error handling
echo Starting Flask application...
python app.py
if %ERRORLEVEL% neq 0 (
    echo Application exited with error code %ERRORLEVEL%
    echo Check the logs above for details.
)
```

### Linux/Mac (run_app.sh)

```bash
#!/bin/bash
echo "Starting CryoProtect Analyzer..."

# Activate the conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cryoprotect

# Check if RDKit is installed
python -c "from rdkit import Chem; print('RDKit check passed')" > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: RDKit is not installed or not working properly."
    echo "Please run setup_environment.sh first."
    exit 1
fi

# Run the Flask application with error handling
echo "Starting Flask application..."
python app.py
if [ $? -ne 0 ]; then
    echo "Application exited with error code $?"
    echo "Check the logs above for details."
fi
```

## Conclusion

By following this troubleshooting guide, you should be able to resolve RDKit installation issues and get the CryoProtect Analyzer application running properly. If you continue to experience problems, please check the RDKit documentation or seek help from the RDKit community.