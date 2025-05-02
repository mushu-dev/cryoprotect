@echo off
setlocal enabledelayedexpansion

echo ===============================================================
echo CryoProtect v2 - Database Remediation Quick Start Guide
echo ===============================================================
echo.

:menu
echo Please select an option:
echo.
echo 1. Run complete database remediation
echo 2. Run in dry-run mode (no changes will be made)
echo 3. Run verification only
echo 4. Run a specific phase
echo 5. Test the remediation manager
echo 6. View documentation
echo 7. Exit
echo.

set /p choice=Enter your choice (1-7): 

if "%choice%"=="1" (
    echo.
    echo Running complete database remediation...
    echo.
    call run_database_remediation.bat
    goto end
)

if "%choice%"=="2" (
    echo.
    echo Running in dry-run mode (no changes will be made)...
    echo.
    call run_database_remediation.bat --dry-run
    goto end
)

if "%choice%"=="3" (
    echo.
    echo Running verification only...
    echo.
    call run_database_remediation.bat --verify-only
    goto end
)

if "%choice%"=="4" (
    echo.
    echo Available phases:
    echo 1. Security Lockdown
    echo 2. Schema Standardization
    echo 3. Relationship Remediation
    echo 4. Data Migration
    echo 5. Performance Optimization
    echo.
    set /p phase=Enter phase number (1-5): 
    
    if "!phase!"=="1" (
        echo.
        echo Running Phase 1: Security Lockdown...
        echo.
        call run_database_remediation.bat --phase 1
        goto end
    )
    
    if "!phase!"=="2" (
        echo.
        echo Running Phase 2: Schema Standardization...
        echo.
        call run_database_remediation.bat --phase 2
        goto end
    )
    
    if "!phase!"=="3" (
        echo.
        echo Running Phase 3: Relationship Remediation...
        echo.
        call run_database_remediation.bat --phase 3
        goto end
    )
    
    if "!phase!"=="4" (
        echo.
        echo Running Phase 4: Data Migration...
        echo.
        call run_database_remediation.bat --phase 4
        goto end
    )
    
    if "!phase!"=="5" (
        echo.
        echo Running Phase 5: Performance Optimization...
        echo.
        call run_database_remediation.bat --phase 5
        goto end
    )
    
    echo Invalid phase number.
    goto menu
)

if "%choice%"=="5" (
    echo.
    echo Testing the remediation manager...
    echo.
    call test_database_remediation.bat
    goto end
)

if "%choice%"=="6" (
    echo.
    echo Opening documentation...
    echo.
    start "" notepad.exe README_Database_Remediation_Manager.md
    goto menu
)

if "%choice%"=="7" (
    echo.
    echo Exiting...
    goto exit
)

echo Invalid choice. Please try again.
goto menu

:end
echo.
echo Press any key to return to the menu...
pause > nul
goto menu

:exit
endlocal