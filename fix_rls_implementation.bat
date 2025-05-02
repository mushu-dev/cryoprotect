@echo off
echo Running RLS Implementation Fix...

REM Pass all command-line arguments to the Python script
python fix_rls_implementation.py %*

if %ERRORLEVEL% NEQ 0 (
    echo Error: RLS Implementation Fix failed.
    exit /b 1
)
echo RLS Implementation Fix completed successfully.

echo.
echo Usage examples:
echo   fix_rls_implementation.bat                   - Run the fix with default settings
echo   fix_rls_implementation.bat --dry-run         - Show what would be done without making changes
echo   fix_rls_implementation.bat --verify-only     - Only verify the current RLS implementation
echo   fix_rls_implementation.bat --rollback        - Rollback changes to the original state
echo   fix_rls_implementation.bat --schema public   - Run on a different schema (default: remediation_test)