@echo off
REM CryoProtect v2 - Robust Server Startup Script for Windows (Conda + RDKit)

REM Activate the conda environment
call activate cryoprotect

REM Diagnostic: Print Python executable and installed packages
python --version
where python
python -m pip list

REM Set RDKit base path if needed (usually not required if installed via conda)
REM set RDBASE=%CONDA_PREFIX%\Library\share\rdkit

REM Set Flask app entry point
set FLASK_APP=app.py

REM Optionally set Flask environment (development/production)
set FLASK_ENV=development

REM Start the Flask server (development)
REM Use the conda environment's python.exe to run Flask
C:\Users\1edwa\anaconda3\envs\cryoprotect\python.exe -m flask run --host=0.0.0.0 --port=5000

REM For production, use Waitress (uncomment below after installing waitress)
REM waitress-serve --port=5000 app:app

pause