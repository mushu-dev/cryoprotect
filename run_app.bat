@echo off
echo Starting CryoProtect Analyzer...

:: Activate the conda environment
call conda activate cryoprotect

:: Verify and install required packages
echo Verifying required packages...
python verify_packages.py > nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo Some packages are missing or not accessible.
    echo Running package installation script...
    call startup_packages.bat
)

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