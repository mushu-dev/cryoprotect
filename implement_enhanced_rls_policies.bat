@echo off
REM Script to implement enhanced RLS policies for CryoProtect v2

REM Ensure we're in the project root directory
cd /d "%~dp0"

REM Check if Python is installed
where python >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo Python is not installed or not in PATH. Please install Python and try again.
    exit /b 1
)

REM Check if virtual environment exists
if exist .venv\Scripts\activate.bat (
    echo Activating virtual environment...
    call .venv\Scripts\activate.bat
)

REM Check if .env file exists
if not exist .env (
    echo Warning: .env file not found. Using default configuration.
    if exist .env.template (
        echo Consider copying .env.template to .env and updating the values.
    )
)

REM Create logs directory if it doesn't exist
if not exist logs mkdir logs
if not exist reports\security mkdir reports\security

REM Run the Python script
echo Implementing enhanced RLS policies...
python implement_enhanced_rls_policies.py %*

REM Check if the script executed successfully
if %ERRORLEVEL% equ 0 (
    echo ✅ Enhanced RLS policies implemented successfully!
    echo Check the reports/security directory for verification reports.
) else (
    echo ❌ Error implementing enhanced RLS policies. Check logs for details.
    exit /b 1
)