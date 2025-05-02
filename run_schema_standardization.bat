@echo off
echo ===================================================
echo CryoProtect Database Schema Standardization
echo ===================================================
echo.

REM Check if Python is installed
where python >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo ERROR: Python is not installed or not in PATH.
    echo Please install Python 3.6+ and try again.
    goto :end
)

REM Check Python version
for /f "tokens=2" %%I in ('python --version 2^>^&1') do set PYTHON_VERSION=%%I
echo Detected Python version: %PYTHON_VERSION%
echo.

REM Create a virtual environment if it doesn't exist
if not exist .venv (
    echo Creating virtual environment...
    python -m venv .venv
    if %ERRORLEVEL% neq 0 (
        echo ERROR: Failed to create virtual environment.
        goto :end
    )
)

REM Activate the virtual environment
echo Activating virtual environment...
call .venv\Scripts\activate.bat
if %ERRORLEVEL% neq 0 (
    echo ERROR: Failed to activate virtual environment.
    goto :end
)

REM Install required packages
echo Installing required packages...
pip install -q requests
if %ERRORLEVEL% neq 0 (
    echo ERROR: Failed to install required packages.
    goto :end
)

REM Run the test script first
echo.
echo ===================================================
echo Testing Supabase connection...
echo ===================================================
python test_supabase_connection.py
if %ERRORLEVEL% neq 0 (
    echo.
    echo Connection test failed. Please fix the issues before continuing.
    echo.
    set /p CONTINUE=Do you want to continue anyway? (y/n): 
    if /i "%CONTINUE%" neq "y" goto :end
)

echo.
echo ===================================================
echo IMPORTANT: This script will standardize the database schema
echo by converting singular table names to plural and adding
echo proper constraints. This operation is potentially destructive.
echo ===================================================
echo.
echo Please make sure you have:
echo  1. Backed up your database
echo  2. Reviewed the standardize_schema.py script
echo  3. Tested your connection with test_supabase_connection.py
echo.
set /p CONFIRM=Are you sure you want to continue? (y/n): 

if /i "%CONFIRM%" neq "y" (
    echo Operation cancelled by user.
    goto :end
)

echo.
echo ===================================================
echo Running schema standardization...
echo ===================================================
python standardize_schema.py

if %ERRORLEVEL% equ 0 (
    echo.
    echo ===================================================
    echo Schema standardization completed successfully!
    echo ===================================================
) else (
    echo.
    echo ===================================================
    echo Schema standardization completed with errors.
    echo Please check the log file for details.
    echo ===================================================
)

:end
echo.
echo Press any key to exit...
pause >nul