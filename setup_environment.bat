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