@echo off
echo Starting CryoProtect Analyzer with Service Role Authentication Fix...

:: Check if auth_config.py exists
if not exist auth_config.py (
    echo Service role authentication fix not found.
    echo Running fix_auth_simple.py to apply the fix...
    python fix_auth_simple.py
    if %ERRORLEVEL% neq 0 (
        echo Failed to apply service role authentication fix.
        echo Please check the error messages above.
        pause
        exit /b 1
    )
) else (
    echo Service role authentication fix already applied.
)

:: Activate the conda environment
call conda activate cryoprotect

:: Check if RDKit is installed
python -c "from rdkit import Chem; print('RDKit check passed')" > nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo ERROR: RDKit is not installed or not working properly.
    echo Please run setup_environment.bat first.
    pause
    exit /b 1
)

:: Test the authentication
echo Testing Supabase authentication with service role approach...
python test_service_role_auth.py
if %ERRORLEVEL% neq 0 (
    echo Authentication test failed with error code %ERRORLEVEL%
    echo Please check the error messages above.
    pause
    exit /b 1
)

:: Run the Flask application with error handling
echo Starting Flask application...
python app.py
if %ERRORLEVEL% neq 0 (
    echo Application exited with error code %ERRORLEVEL%
    echo Check the logs above for details.
    pause
)
