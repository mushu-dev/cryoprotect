@echo off
echo ===============================================================
echo CryoProtect v2 Database Performance Testing
echo ===============================================================
echo.

REM Check if Python is installed
python --version >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python is not installed or not in the PATH.
    echo Please install Python 3.7+ and try again.
    exit /b 1
)

REM Check if required packages are installed
echo Checking required packages...
python -c "import supabase" >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Installing required packages...
    pip install supabase psutil python-dotenv
)

REM Check if .env file exists
if not exist .env (
    echo Warning: .env file not found.
    echo Creating a template .env file. Please edit it with your Supabase credentials.
    echo SUPABASE_URL=your_supabase_url> .env
    echo SUPABASE_KEY=your_supabase_service_role_key>> .env
    echo.
    echo Please edit the .env file with your Supabase credentials and run this script again.
    exit /b 1
)

echo Starting performance tests...
echo This may take several minutes to complete.
echo.

REM Run the main test script
python run_all_performance_tests.py

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo Error: Performance tests failed.
    echo Please check the logs for more information.
    exit /b 1
)

echo.
echo ===============================================================
echo Performance tests completed successfully!
echo ===============================================================
echo.
echo Reports are available in:
echo - comprehensive_performance_report.txt (Summary report)
echo - comprehensive_performance_report.json (Detailed JSON report)
echo - database_performance_report.txt (Database performance)
echo - query_plan_analysis_report.txt (Query plan analysis)
echo - concurrent_load_test_report.txt (Concurrent load test)
echo.
echo See README_Performance_Testing.md for information on interpreting the results.
echo.

REM Ask if user wants to open the report
set /p open_report=Would you like to open the summary report now? (y/n): 

if /i "%open_report%"=="y" (
    start notepad comprehensive_performance_report.txt
)

echo.
echo Thank you for using CryoProtect v2 Database Performance Testing.
echo.