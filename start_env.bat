@echo off
REM CryoProtect Analyzer - Automated Environment Setup and Startup Script

REM 1. Check if conda is available
where conda >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo ERROR: Conda is not installed or not in your PATH. Please install Anaconda or Miniconda and try again.
    exit /b 1
)

REM 2. Check if the "cryoprotect" environment exists
conda info --envs | findstr /C:"cryoprotect" >nul
if %ERRORLEVEL% neq 0 (
    echo Creating conda environment "cryoprotect" with Python 3.9...
    conda create -n cryoprotect python=3.9 -y
    if %ERRORLEVEL% neq 0 (
        echo ERROR: Failed to create conda environment.
        exit /b 1
    )
) else (
    echo Conda environment "cryoprotect" already exists.
)

REM 3. Install RDKit in the environment
echo Installing RDKit in "cryoprotect" environment...
conda install -n cryoprotect -c conda-forge rdkit=2023.9.1 -y
if %ERRORLEVEL% neq 0 (
    echo ERROR: Failed to install RDKit.
    exit /b 1
)

REM 4. Install all other dependencies
echo Installing Python dependencies from requirements.txt...
conda run -n cryoprotect pip install -r requirements.txt
if %ERRORLEVEL% neq 0 (
    echo ERROR: Failed to install Python dependencies.
    exit /b 1
)

REM 5. Verify RDKit installation
echo Verifying RDKit installation...
conda run -n cryoprotect python -c "from rdkit import Chem; print('RDKit installation successful!')"
if %ERRORLEVEL% neq 0 (
    echo ERROR: RDKit verification failed. Please check the installation.
    exit /b 1
)

REM 6. (Optional) Start the API server
echo.
echo To activate the environment and start the API server, run:
echo   conda activate cryoprotect
echo   python app.py
echo.
echo Environment setup complete!