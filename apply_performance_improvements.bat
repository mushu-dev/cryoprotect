@echo off
echo ===============================================================
echo CryoProtect v2 - Apply Critical Performance Improvements
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
if not exist apply_performance_improvements.py (
    echo Error: apply_performance_improvements.py not found.
    echo Please make sure the script is in the current directory.
    exit /b 1
)

echo Applying critical performance improvements...
echo.

REM Run the script with backup option
python apply_performance_improvements.py --backup

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Performance improvements applied successfully!
    exit /b 0
) else if %ERRORLEVEL% EQU 1 (
    echo.
    echo Performance improvements applied with warnings.
    echo Please check the log file for details.
    exit /b 1
) else (
    echo.
    echo Error: Failed to apply performance improvements.
    echo Please check the log file for details.
    exit /b 2
)