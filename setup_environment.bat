@echo off
echo Setting up CryoProtect v2 environment...

:: Check for .env file and create from template if it doesn't exist
if not exist .env (
    echo Creating .env file from template...
    copy .env.template .env
    echo ✅ Created .env file from template
    echo ⚠️ IMPORTANT: You need to update the .env file with your credentials
    echo    Please edit the .env file and fill in the required values marked with [YOUR-*]
) else (
    echo ✅ .env file already exists
)

:: Try to create a new conda environment
echo Creating conda environment...
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

:: Print next steps
echo.
echo 🎉 Environment setup complete! 🎉
echo.
echo Next steps:
echo 1. Activate the environment with: conda activate cryoprotect
echo 2. Edit the .env file to fill in your Supabase credentials and other required values
echo 3. Run the application with: run_app.bat
echo 4. For debugging, use: python debug_app.py
echo.
echo For RDKit troubleshooting, see: README_RDKit_Troubleshooting.md
echo For API documentation, see: README_API.md
echo For authentication setup, see: README_Authentication.md