@echo off
echo ===============================================================
echo CryoProtect v2 - Verify Database Remediation
echo ===============================================================
echo.

REM Check if Python is installed
where python >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python is not installed or not in the PATH.
    echo Please install Python and try again.
    exit /b 1
)

echo This script will verify that the database remediation was successful.
echo It will check:
echo.
echo 1. RLS is enabled on all tables
echo 2. Anonymous access is properly restricted
echo 3. All tables have appropriate RLS policies
echo 4. All foreign keys have corresponding indexes
echo 5. Application roles are created correctly
echo 6. Table names are standardized to plural form
echo 7. Junction tables are created to fix fan traps
echo.

set /p CONFIRM=Do you want to continue? (y/n): 
if /i "%CONFIRM%" NEQ "y" (
    echo Operation cancelled.
    exit /b 0
)

echo.
echo Running verification...
python verify_database_remediation.py

echo.
echo Verification process completed.
echo Check the log file for details.
echo.

pause