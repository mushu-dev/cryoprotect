@echo off
echo Running CryoProtect v2 API Integration Tests...
echo.

REM Check if Python is installed
where python >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo Python is not installed or not in the PATH.
    echo Please install Python and try again.
    exit /b 1
)

REM Check if the virtual environment exists
if not exist .venv (
    echo Virtual environment not found.
    echo Creating virtual environment...
    python -m venv .venv
    if %ERRORLEVEL% neq 0 (
        echo Failed to create virtual environment.
        exit /b 1
    )
)

REM Activate the virtual environment
call .venv\Scripts\activate.bat

REM Install required packages
echo Installing required packages...
pip install requests
if %ERRORLEVEL% neq 0 (
    echo Failed to install required packages.
    exit /b 1
)

REM Run the test script
echo.
echo Starting API integration tests...
echo.
python test_api_integration.py
set TEST_RESULT=%ERRORLEVEL%

REM Deactivate the virtual environment
call .venv\Scripts\deactivate.bat

REM Return the test result
exit /b %TEST_RESULT%