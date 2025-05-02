@echo off
echo CryoProtect v2 - Database Population Script
echo ==========================================
echo.

REM Check if Python is installed
where python >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo Error: Python is not installed or not in PATH.
    echo Please install Python 3.8 or higher and try again.
    exit /b 1
)

REM Check if .env file exists
if not exist .env (
    echo Warning: .env file not found.
    echo Creating a template .env file. Please edit it with your Supabase credentials.
    echo SUPABASE_URL=your_supabase_url> .env
    echo SUPABASE_KEY=your_supabase_key>> .env
    echo SUPABASE_USER=your_supabase_user_email>> .env
    echo SUPABASE_PASSWORD=your_supabase_user_password>> .env
    echo.
    echo Template .env file created. Please edit it and run this script again.
    exit /b 1
)

REM Check if required packages are installed
echo Checking required packages...
pip show python-dotenv supabase >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo Installing required packages...
    pip install python-dotenv supabase
)

REM Run the population script
echo.
echo Running database population script...
echo.
python populate_database_supabase.py %*

REM Check if the script ran successfully
if %ERRORLEVEL% neq 0 (
    echo.
    echo Error: Database population failed.
    exit /b 1
)

echo.
echo Database population completed successfully.
echo.
echo Running verification script...
echo.
python verify_database_population.py

echo.
echo Process completed.
echo See database_population.log and database_verification.log for details.
echo.