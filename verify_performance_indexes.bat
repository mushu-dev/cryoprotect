@echo off
echo ===============================================================
echo CryoProtect v2 - Verify Performance Indexes
echo ===============================================================
echo.

python verify_performance_indexes.py

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo Error: Performance indexes verification failed.
    echo Please check the log file for details.
    exit /b 1
)

echo.
echo Performance indexes verification completed successfully!
echo.