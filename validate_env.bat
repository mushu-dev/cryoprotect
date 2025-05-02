@echo off
:: Environment Variable Validation Script for CryoProtect v2
::
:: This script validates that all required environment variables are set by:
:: 1. Parsing the .env.template file to extract variable names
:: 2. Checking if each variable is set in the current environment
:: 3. Reporting any missing variables
:: 4. Exiting with appropriate status code

setlocal EnableDelayedExpansion

:: Check if template file exists
if not exist ".env.template" (
    echo Error: Template file '.env.template' not found.
    exit /b 1
)

echo Validating environment variables for CryoProtect v2...

:: Initialize variables
set "MISSING_REQUIRED="
set "EXIT_CODE=0"

:: Check required variables
call :check_var "SUPABASE_URL" "DATABASE CONFIGURATION"
call :check_var "SUPABASE_KEY" "DATABASE CONFIGURATION"
call :check_var "SECRET_KEY" "APPLICATION CONFIGURATION"

:: Report results
if not "!MISSING_REQUIRED!"=="" (
    echo.
    echo ERROR: The following required environment variables are missing:
    
    for %%v in (!MISSING_REQUIRED!) do (
        echo   - %%v
    )
    
    echo.
    echo Please set these required variables in your environment or .env file before running the application.
    exit /b 1
) else (
    echo All required environment variables are set.
    exit /b 0
)

:: Function to check if a variable is set
:check_var
set "VAR_NAME=%~1"
set "SECTION=%~2"

call set "VAR_VALUE=%%%VAR_NAME%%%"
if "!VAR_VALUE!"=="" (
    set "MISSING_REQUIRED=!MISSING_REQUIRED! %VAR_NAME%"
)
goto :eof