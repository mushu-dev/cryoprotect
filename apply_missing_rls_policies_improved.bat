@echo off
echo.
echo ============================================================
echo CryoProtect Analyzer - Apply Enhanced RLS Policies
echo ============================================================
echo This script will apply improved RLS policies to missing tables/views:
echo 1. experiment_with_results
echo 2. migrations
echo 3. mixture_with_components
echo 4. molecule_with_properties
echo.
echo Enhancements include:
echo - Transaction support for atomic operations
echo - Enhanced verification and effectiveness testing
echo - Performance benchmarking for RLS impact
echo - Comprehensive audit trail for scientific data
echo - Optimized RLS policy conditions for better performance
echo.
echo Prerequisites:
echo - Python installed with psycopg2 package
echo - Supabase connection configured in .env file
echo.

REM Check for Python
where python >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Python is not installed or not in PATH.
    echo Please install Python and try again.
    exit /b 1
)

REM Ensure backup before proceeding
echo Creating database backup before applying RLS policies...
python create_database_backup.py --description "Pre-RLS-enhancement backup"

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo WARNING: Could not create database backup. This is recommended before proceeding.
    echo          Do you want to continue anyway? (Y/N)
    set /p CONTINUE=
    if /i "%CONTINUE%" NEQ "Y" (
        echo Operation cancelled by user.
        exit /b 0
    )
)

REM Check for required packages
python -c "import psycopg2" 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Installing required packages...
    pip install psycopg2-binary python-dotenv
)

echo Running enhanced RLS policy application...
python apply_missing_rls_policies_improved.py

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo ERROR: RLS policy application failed. Check the logs for details.
    exit /b 1
) else (
    echo.
    echo SUCCESS: Enhanced RLS policies have been applied successfully.
    echo          A detailed report has been generated in the reports/security directory.
    echo.
)

REM Optional: Open the report
echo Do you want to open the latest report? (Y/N)
set /p OPEN_REPORT=
if /i "%OPEN_REPORT%" EQU "Y" (
    for /f "delims=" %%i in ('dir /b /od /a-d "reports\security\rls_implementation_summary_*.md"') do set LATEST_REPORT=%%i
    start "" "reports\security\%LATEST_REPORT%"
)

exit /b 0