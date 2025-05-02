@echo off
echo ===================================================
echo CryoProtect v2 Database Remediation Verification
echo ===================================================
echo.

REM Activate virtual environment if it exists
if exist .venv\Scripts\activate (
    call .venv\Scripts\activate
) else (
    echo WARNING: Virtual environment not found. Using system Python.
)

REM Check if Python is available
python --version >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Python not found. Please install Python 3.8 or higher.
    exit /b 1
)

REM Check if required packages are installed
echo Checking required packages...
python -c "import supabase" >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Installing supabase package...
    pip install supabase
)

python -c "import requests" >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Installing requests package...
    pip install requests
)

python -c "import dotenv" >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Installing python-dotenv package...
    pip install python-dotenv
)

REM Create reports directory if it doesn't exist
if not exist reports mkdir reports

echo.
echo Running database verification tests...
echo.

REM Run the verification script
python verify_database_remediation.py --verbose

echo.
echo Verification complete. Check the reports directory for detailed results.
echo.

REM Pause to keep the window open
pause