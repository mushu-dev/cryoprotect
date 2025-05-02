@echo off
:: CryoProtect v2 - Database Operations
:: This script provides shortcuts for common database operations

echo CryoProtect v2 - Database Operations

if "%1" == "" (
    echo Usage: run_database [operation] [options]
    echo.
    echo Operations:
    echo   populate    - Populate database with data
    echo   verify      - Verify database integrity
    echo   fix         - Fix database issues
    echo.
    echo Example: run_database populate --molecules-file molecules.json
    exit /b 1
)

python database_cli.py %*