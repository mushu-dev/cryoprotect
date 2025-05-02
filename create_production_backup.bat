@echo off
echo ===============================================================
echo CryoProtect v2 - Create Production Database Backup
echo ===============================================================
echo.

REM Check if Python is installed
where python >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python is not installed or not in the PATH.
    echo Please install Python and try again.
    exit /b 1
)

REM Check if the script exists
if not exist create_production_backup.py (
    echo Error: create_production_backup.py not found.
    echo Please make sure the script is in the current directory.
    exit /b 1
)

echo Creating production database backup...
echo This may take several minutes depending on the database size.
echo.

REM Run the backup script
python create_production_backup.py --format both

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Production database backup created successfully!
    echo Please check the production_backups directory for the backup files.
    exit /b 0
) else if %ERRORLEVEL% EQU 1 (
    echo.
    echo Production database backup created with warnings.
    echo Please check the log file for details.
    exit /b 1
) else (
    echo.
    echo Error: Failed to create production database backup.
    echo Please check the log file for details.
    exit /b 2
)