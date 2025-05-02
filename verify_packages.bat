@echo off
echo Verifying Python packages in cryoprotect environment...

:: Activate the conda environment
call conda activate cryoprotect

:: Run the verification script
python verify_packages.py

:: Pause to see the results
pause