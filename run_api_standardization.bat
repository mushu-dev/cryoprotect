@echo off
REM CryoProtect Analyzer API - Standardization Runner (Windows)
REM
REM This script runs the API audit and applies standardization to all API endpoints.
REM It provides a convenient way to standardize the API in one step.

echo Running API standardization process...

REM Check if Python is installed
where python >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo Python is not installed or not in PATH. Please install Python and try again.
    exit /b 1
)

REM Parse command line arguments
set DRY_RUN=
set VERBOSE=
set SKIP_AUDIT=
set SKIP_APPLY=

:parse_args
if "%~1"=="" goto :run
if /i "%~1"=="--dry-run" set DRY_RUN=--dry-run
if /i "%~1"=="--verbose" set VERBOSE=--verbose
if /i "%~1"=="--skip-audit" set SKIP_AUDIT=--skip-audit
if /i "%~1"=="--skip-apply" set SKIP_APPLY=--skip-apply
shift
goto :parse_args

:run
REM Run the standardization process
python run_api_standardization.py %DRY_RUN% %VERBOSE% %SKIP_AUDIT% %SKIP_APPLY%

if %ERRORLEVEL% neq 0 (
    echo API standardization process failed.
    exit /b 1
)

echo API standardization process completed successfully.
exit /b 0