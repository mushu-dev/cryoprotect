@echo off
echo Running Supabase API Integration Tests...
echo.

REM Activate virtual environment if it exists
if exist .venv\Scripts\activate.bat (
    call .venv\Scripts\activate.bat
)

REM Run the tests
python tests\run_supabase_api_tests.py

REM Store the exit code
set TEST_RESULT=%ERRORLEVEL%

REM Deactivate virtual environment if it was activated
if exist .venv\Scripts\deactivate.bat (
    call .venv\Scripts\deactivate.bat
)

REM Exit with the test result
exit /b %TEST_RESULT%