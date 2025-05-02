@echo off
REM CryoProtect v2 - Test Runner Script
REM This script runs all the tests and generates a comprehensive test report.

echo CryoProtect v2 - Running Tests
echo =============================
echo.

REM Check if Python is installed
python --version >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python is not installed or not in PATH
    exit /b 1
)

REM Check if the virtual environment exists
if not exist ".venv" (
    echo Creating virtual environment...
    python -m venv .venv
    
    if %ERRORLEVEL% NEQ 0 (
        echo Error: Failed to create virtual environment
        exit /b 1
    )
)

REM Activate the virtual environment
echo Activating virtual environment...
if exist ".venv\Scripts\activate.bat" (
    call .venv\Scripts\activate.bat
) else (
    echo Error: Could not find activation script for virtual environment
    exit /b 1
)

REM Install required packages
echo Installing required packages...
pip install -r requirements.txt

if %ERRORLEVEL% NEQ 0 (
    echo Error: Failed to install required packages
    exit /b 1
)

REM Run the test runner
echo Running tests...
python run_all_tests.py

REM Check the exit code
if %ERRORLEVEL% EQU 0 (
    echo All tests passed!
    echo See TEST_RESULTS_REPORT_*.md for detailed results
) else (
    echo Some tests failed. See TEST_RESULTS_REPORT_*.md for detailed results
)

REM Deactivate the virtual environment
call deactivate

echo.
echo Testing complete