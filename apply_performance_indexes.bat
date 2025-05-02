@echo off
echo ===============================================================
echo CryoProtect v2 - Apply Performance Indexes Migration
echo ===============================================================
echo.

REM Check if Node.js is installed
where node >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: Node.js is not installed or not in the PATH.
    echo Please install Node.js and try again.
    exit /b 1
)

echo Applying performance indexes migration...
echo.

REM Run the migration script
node migrations\apply_performance_indexes_migration.js

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo Error: Migration failed.
    echo Please check the error messages above.
    exit /b 1
)

echo.
echo Migration completed successfully!
echo.