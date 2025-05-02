@echo off
echo ===============================================================
echo CryoProtect v2 - Test exec_sql Function
echo ===============================================================
echo.

REM Check if Python is installed
where python >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python is not installed or not in the PATH.
    echo Please install Python and try again.
    exit /b 1
)

REM Check if .env file exists
if not exist .env (
    echo Error: .env file not found.
    echo Please create a .env file with your Supabase credentials.
    echo You can use .env.template as a reference.
    exit /b 1
)

echo This script will test the exec_sql function in your Supabase database.
echo It will create a temporary test table, insert data, query it, and then drop it.
echo.

set /p CONFIRM=Do you want to continue? (y/n): 
if /i "%CONFIRM%" NEQ "y" (
    echo Operation cancelled.
    exit /b 0
)

echo.
echo Running test script...
python test_exec_sql_function.py

echo.
echo Test process completed.
echo.

pause