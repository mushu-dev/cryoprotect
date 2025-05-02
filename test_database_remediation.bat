@echo off
echo ===============================================================
echo CryoProtect v2 - Test Database Remediation
echo ===============================================================
echo.

REM Check if Python is installed
where python >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python is not installed or not in the PATH.
    echo Please install Python and try again.
    exit /b 1
)

echo This script will test the database remediation process in a safe environment.
echo It will:
echo 1. Create a test schema with sample data that mimics the issues in production
echo 2. Run the remediation process on this test schema
echo 3. Verify that the remediation was successful
echo 4. Clean up the test schema
echo.
echo This is a safe way to test the remediation process without affecting your production data.
echo.

set /p CONFIRM=Do you want to continue? (y/n): 
if /i "%CONFIRM%" NEQ "y" (
    echo Operation cancelled.
    exit /b 0
)

echo.
echo Select an option:
echo 1. Run test and clean up afterward
echo 2. Run test and keep the test schema for inspection
echo.

set /p OPTION=Enter option (1-2): 

if "%OPTION%"=="1" (
    echo.
    echo Running test and cleaning up afterward...
    python test_database_remediation.py
) else if "%OPTION%"=="2" (
    echo.
    echo Running test and keeping the test schema...
    python test_database_remediation.py --keep-schema
) else (
    echo Invalid option.
    exit /b 1
)

echo.
echo Test process completed.
echo Check the log file for details.
echo.

pause