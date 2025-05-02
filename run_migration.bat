@echo off
REM CryoProtect v2 - Database Migration Runner Script
REM This script provides a convenient way to run the migration script with different options

echo CryoProtect v2 - Database Migration to Plural Tables
echo ====================================================
echo.

REM Check if Python is installed
where python >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo Error: Python is not installed or not in PATH
    exit /b 1
)

REM Check if the migration script exists
if not exist migrate_to_plural_tables.py (
    echo Error: migrate_to_plural_tables.py not found
    exit /b 1
)

REM Parse command line arguments
set DRY_RUN=
set ROLLBACK=
set REPORT=
set TEST=

:parse_args
if "%~1"=="" goto :end_parse_args
if /i "%~1"=="--dry-run" set DRY_RUN=--dry-run
if /i "%~1"=="--rollback" set ROLLBACK=--rollback
if /i "%~1"=="--report" set REPORT=--report %~2 & shift
if /i "%~1"=="--test" set TEST=1
shift
goto :parse_args
:end_parse_args

REM Run the test script if requested
if defined TEST (
    echo Running migration test script...
    echo.
    python test_migration_script.py
    if %ERRORLEVEL% neq 0 (
        echo.
        echo Migration test failed.
        exit /b 1
    )
    echo.
    echo Migration test completed successfully.
    exit /b 0
)

REM Run the migration script with the specified options
if defined ROLLBACK (
    if not defined REPORT (
        echo Error: --report is required with --rollback
        echo Usage: run_migration.bat --rollback --report migration_report_YYYYMMDD_HHMMSS.json
        exit /b 1
    )
    echo Running migration rollback...
    echo.
    python migrate_to_plural_tables.py --rollback %REPORT%
) else if defined DRY_RUN (
    echo Running migration in dry run mode...
    echo.
    python migrate_to_plural_tables.py --dry-run
) else (
    echo.
    echo WARNING: You are about to run the actual migration.
    echo This will modify the database structure and should be done during a maintenance window.
    echo.
    set /p CONFIRM=Are you sure you want to continue? (y/n): 
    if /i not "%CONFIRM%"=="y" (
        echo Migration cancelled.
        exit /b 0
    )
    
    echo.
    echo Running migration...
    echo.
    python migrate_to_plural_tables.py
)

if %ERRORLEVEL% neq 0 (
    echo.
    echo Migration failed. Check the log file for details.
    exit /b 1
)

echo.
echo Migration completed successfully.
exit /b 0