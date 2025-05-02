@echo off
REM Run CryoProtect Analyzer with Collaboration Features

echo CryoProtect Analyzer - Starting with Collaboration Features

REM Check if .env file exists
if not exist .env (
    echo Error: .env file not found.
    echo Please create a .env file with your Supabase credentials.
    echo Example:
    echo SUPABASE_URL=https://your-project.supabase.co
    echo SUPABASE_KEY=your-api-key
    echo SUPABASE_USER=your-email@example.com
    echo SUPABASE_PASSWORD=your-password
    echo FLASK_ENV=development
    exit /b 1
)

REM Apply the teams migration if needed
echo Checking if teams migration needs to be applied...
node migrations/apply_teams_migration.js

REM Check if the migration was successful
if %ERRORLEVEL% neq 0 (
    echo Error applying teams migration.
    exit /b 1
)

REM Run the application
echo Starting CryoProtect Analyzer...
python app.py

echo CryoProtect Analyzer has been stopped.