@echo off
echo Running API Performance Benchmark...

REM Check if Python is installed
where python >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo Python is not installed or not in PATH. Please install Python and try again.
    exit /b 1
)

REM Check if psutil is installed
python -c "import psutil" >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo Installing psutil...
    pip install psutil
)

REM Run the benchmark script
python benchmark_api_endpoints.py %*

if %ERRORLEVEL% neq 0 (
    echo Benchmark failed with error code %ERRORLEVEL%
    exit /b %ERRORLEVEL%
)

echo Benchmark completed successfully.
echo Report saved to API_Performance_Benchmark_Report.md