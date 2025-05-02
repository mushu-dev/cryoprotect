@echo off
echo Running API integration fixes for CryoProtect v2...
echo.

REM Check if Python is installed
python --version >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python is not installed or not in PATH.
    echo Please install Python and try again.
    exit /b 1
)

REM Run the fix script
python fix_api_integration.py

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo Error: Failed to apply API fixes.
    echo Please check the logs for details.
    exit /b 1
) else (
    echo.
    echo API integration fixes applied successfully!
    echo The API now works with the standardized database schema.
    echo.
    echo Changes made:
    echo 1. Updated API resource classes to use plural table names
    echo 2. Fixed endpoint duplication issues
    echo 3. Implemented consistent error handling
    echo 4. Added retry logic for API calls
    echo.
    echo Backups of the original files were created with .bak extension.
)

pause