@echo off
echo ===============================================================
echo CryoProtect v2 - Setup Database Remediation Environment
echo ===============================================================
echo.

REM Check if Python is installed
where python >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python is not installed or not in the PATH.
    echo Please install Python and try again.
    exit /b 1
)

echo This script will set up the environment for the database remediation process.
echo It will:
echo 1. Install required Python packages
echo 2. Create a .env file for your Supabase credentials
echo.

set /p CONFIRM=Do you want to continue? (y/n): 
if /i "%CONFIRM%" NEQ "y" (
    echo Operation cancelled.
    exit /b 0
)

echo.
echo Installing required Python packages...
pip install supabase python-dotenv

if %ERRORLEVEL% NEQ 0 (
    echo Error: Failed to install required packages.
    echo Please check your internet connection and try again.
    exit /b 1
)

echo.
echo Creating .env file...

if exist .env (
    echo .env file already exists.
    set /p OVERWRITE=Do you want to overwrite it? (y/n): 
    if /i "%OVERWRITE%" NEQ "y" (
        echo Keeping existing .env file.
        goto :setup_complete
    )
)

echo # CryoProtect v2 - Database Remediation Environment Variables > .env
echo. >> .env

echo Please enter your Supabase credentials:
echo.

set /p SUPABASE_URL=Supabase URL (e.g., https://abcdefghijklm.supabase.co): 
set /p SUPABASE_KEY=Supabase Service Role Key: 
set /p SUPABASE_PROJECT_ID=Supabase Project ID (optional): 

echo. >> .env
echo # Supabase Project URL >> .env
echo SUPABASE_URL=%SUPABASE_URL% >> .env
echo. >> .env
echo # Supabase Service Role Key >> .env
echo SUPABASE_KEY=%SUPABASE_KEY% >> .env

if not "%SUPABASE_PROJECT_ID%"=="" (
    echo. >> .env
    echo # Supabase Project ID >> .env
    echo SUPABASE_PROJECT_ID=%SUPABASE_PROJECT_ID% >> .env
)

echo .env file created successfully.

:setup_complete
echo.
echo Environment setup completed successfully.
echo You can now run the database remediation process using:
echo   run_database_remediation.bat
echo.

pause