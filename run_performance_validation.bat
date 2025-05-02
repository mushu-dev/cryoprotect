@echo off
echo ===================================================
echo CryoProtect v2 - Database Performance Validation
echo ===================================================
echo.
echo This script will run performance tests to validate:
echo 1. Query performance on tables with RLS enabled
echo 2. Query performance on tables with foreign key relationships
echo 3. Query performance on junction tables
echo.
echo Press Ctrl+C to cancel or any key to continue...
pause > nul

echo.
echo Setting up environment...
call setup_environment.bat > nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo Failed to set up environment. Please run setup_environment.bat manually first.
    exit /b 1
)

echo.
echo Running performance validation tests...
python test_database_performance_remediation.py
if %ERRORLEVEL% neq 0 (
    echo Performance validation failed with error code %ERRORLEVEL%.
    exit /b %ERRORLEVEL%
)

echo.
echo Performance validation completed successfully.
echo Results are available in:
echo - database_performance_validation_report.txt
echo - database_performance_validation_report.json
echo.
echo Press any key to exit...
pause > nul