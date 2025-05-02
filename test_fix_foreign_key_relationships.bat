@echo off
echo CryoProtect v2 - Test Foreign Key Relationship Fixes
echo ===================================================
echo.

REM Activate the Python environment if it exists
if exist .venv\Scripts\activate.bat (
    echo Activating Python environment...
    call .venv\Scripts\activate.bat
) else (
    echo Warning: Python virtual environment not found.
    echo Using system Python installation.
)

REM Run the test script
echo.
echo Running test for foreign key relationship fixes...
python test_fix_foreign_key_relationships.py

REM Check if the script ran successfully
if %ERRORLEVEL% EQU 0 (
    echo.
    echo Test completed successfully.
) else (
    echo.
    echo Test failed.
    echo Please check the log file for details.
)

echo.
echo Press any key to exit...
pause > nul