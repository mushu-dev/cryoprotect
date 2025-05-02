@echo off
echo.
echo ============================================================
echo CryoProtect Analyzer - Apply Missing RLS Policies
echo ============================================================
echo This script will apply RLS policies to missing tables/views:
echo 1. experiment_with_results
echo 2. migrations
echo 3. mixture_with_components
echo 4. molecule_with_properties
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

REM Check for required packages
python -c "import psycopg2" 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Installing required packages...
    pip install psycopg2-binary python-dotenv
)

echo Running RLS policy application...
python apply_missing_rls_policies.py

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo ERROR: RLS policy application failed. Check the logs for details.
    exit /b 1
) else (
    echo.
    echo SUCCESS: RLS policies have been applied successfully.
    echo.
)

exit /b 0