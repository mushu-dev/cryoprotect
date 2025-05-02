@echo off
echo ===============================================================
echo CryoProtect v2 - Database Remediation
echo ===============================================================
echo.

REM Check if Python is installed
where python >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python is not installed or not in the PATH.
    echo Please install Python and try again.
    exit /b 1
)

echo This script will run the complete database remediation process.
echo It will address the following critical issues:
echo.
echo 1. SECURITY: Enable Row Level Security (RLS)
echo 2. STRUCTURE: Standardize Schema ^& Fix Relationships
echo 3. PERFORMANCE: Add Missing Indexes
echo 4. ROLES: Create Application-Specific Roles
echo 5. DATA: Consolidate Duplicate Tables
echo.
echo WARNING: This will make significant changes to your database.
echo Make sure you have a backup before proceeding.
echo.

set /p CONFIRM=Do you want to continue? (y/n): 
if /i "%CONFIRM%" NEQ "y" (
    echo Operation cancelled.
    exit /b 0
)

echo.
echo Select an option:
echo 1. Run all phases
echo 2. Run in dry-run mode (no changes will be made)
echo 3. Run a specific phase
echo.

set /p OPTION=Enter option (1-3): 

if "%OPTION%"=="1" (
    echo.
    echo Running all phases...
    python complete_database_remediation.py
) else if "%OPTION%"=="2" (
    echo.
    echo Running in dry-run mode...
    python complete_database_remediation.py --dry-run
) else if "%OPTION%"=="3" (
    echo.
    echo Select a phase to run:
    echo 1. SECURITY: Enable Row Level Security (RLS)
    echo 2. STRUCTURE: Standardize Schema ^& Fix Relationships
    echo 3. PERFORMANCE: Add Missing Indexes
    echo 4. ROLES: Create Application-Specific Roles
    echo 5. DATA: Consolidate Duplicate Tables
    echo.
    
    set /p PHASE=Enter phase (1-5): 
    
    if "%PHASE%" GEQ "1" if "%PHASE%" LEQ "5" (
        echo.
        echo Running phase %PHASE%...
        python complete_database_remediation.py --phase %PHASE%
    ) else (
        echo Invalid phase number.
        exit /b 1
    )
) else (
    echo Invalid option.
    exit /b 1
)

echo.
echo Database remediation process completed.
echo Check the log file for details.
echo.

pause