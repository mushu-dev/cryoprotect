@echo off
echo CryoProtect v2 - Foreign Key Relationship Fixer
echo ===============================================
echo.

REM Activate the virtual environment if it exists
if exist .\.venv\Scripts\activate (
    call .\.venv\Scripts\activate
) else (
    echo Warning: Virtual environment not found. Make sure dependencies are installed.
)

echo.
echo Available options:
echo 1. Run with dry-run (show what would be done without making changes)
echo 2. Verify only (check for issues without making changes)
echo 3. Apply all fixes
echo 4. Rollback to previous state
echo.

set /p option="Enter option (1-4): "

if "%option%"=="1" (
    echo.
    echo Running in dry-run mode...
    python fix_foreign_key_relationships.py --dry-run
) else if "%option%"=="2" (
    echo.
    echo Verifying foreign key constraints...
    python fix_foreign_key_relationships.py --verify-only
) else if "%option%"=="3" (
    echo.
    echo Applying all foreign key fixes...
    python fix_foreign_key_relationships.py
) else if "%option%"=="4" (
    echo.
    echo Rolling back to previous state...
    python fix_foreign_key_relationships.py --rollback
) else (
    echo.
    echo Invalid option. Please run the script again and select a valid option.
    exit /b 1
)

echo.
echo Script execution completed.
pause