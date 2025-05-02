@echo off
:: scan_docker_image.bat - Security scanning for Docker images on Windows
:: 
:: This script scans Docker images for vulnerabilities using Trivy
:: and can fail the build if critical vulnerabilities are found.
::
:: Usage: scan_docker_image.bat [options] IMAGE_NAME
::
:: Options:
::   --format FORMAT       Output format (table, json, sarif, cyclonedx) [default: table]
::   --output FILE         Output file name [default: trivy-results.{format}]
::   --severity SEVERITY   Severities to scan for (comma-separated) [default: CRITICAL,HIGH,MEDIUM]
::   --exit-code CODE      Exit code when vulnerabilities are found [default: 0]
::   --fail-on SEVERITY    Fail on specific severity (CRITICAL, HIGH, MEDIUM, LOW, UNKNOWN) [default: CRITICAL]
::   --help                Show this help message and exit

setlocal enabledelayedexpansion

:: Default values
set "FORMAT=table"
set "OUTPUT="
set "SEVERITY=CRITICAL,HIGH,MEDIUM"
set "EXIT_CODE=0"
set "FAIL_ON=CRITICAL"
set "IMAGE_NAME="

:: Parse arguments
:parse_args
if "%~1"=="" goto :check_args
if "%~1"=="--format" (
    set "FORMAT=%~2"
    shift /1
    shift /1
    goto :parse_args
)
if "%~1"=="--output" (
    set "OUTPUT=%~2"
    shift /1
    shift /1
    goto :parse_args
)
if "%~1"=="--severity" (
    set "SEVERITY=%~2"
    shift /1
    shift /1
    goto :parse_args
)
if "%~1"=="--exit-code" (
    set "EXIT_CODE=%~2"
    shift /1
    shift /1
    goto :parse_args
)
if "%~1"=="--fail-on" (
    set "FAIL_ON=%~2"
    shift /1
    shift /1
    goto :parse_args
)
if "%~1"=="--help" (
    echo Usage: scan_docker_image.bat [options] IMAGE_NAME
    echo.
    echo Options:
    echo   --format FORMAT       Output format (table, json, sarif, cyclonedx) [default: table]
    echo   --output FILE         Output file name [default: trivy-results.{format}]
    echo   --severity SEVERITY   Severities to scan for (comma-separated) [default: CRITICAL,HIGH,MEDIUM]
    echo   --exit-code CODE      Exit code when vulnerabilities are found [default: 0]
    echo   --fail-on SEVERITY    Fail on specific severity (CRITICAL, HIGH, MEDIUM, LOW, UNKNOWN) [default: CRITICAL]
    echo   --help                Show this help message and exit
    exit /b 0
)
set "IMAGE_NAME=%~1"
shift /1
goto :parse_args

:check_args
:: Check if image name is provided
if "%IMAGE_NAME%"=="" (
    echo Error: IMAGE_NAME is required
    echo Run 'scan_docker_image.bat --help' for usage information
    exit /b 1
)

:: Set default output file if not provided
if "%OUTPUT%"=="" (
    set "OUTPUT=trivy-results.%FORMAT%"
)

echo Scanning Docker image: %IMAGE_NAME%
echo Scanning for vulnerabilities with severity: %SEVERITY%
echo Output format: %FORMAT%
echo Output file: %OUTPUT%

:: Check if Trivy is installed
trivy --version >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo Error: Trivy is not installed
    echo Please install Trivy: https://aquasecurity.github.io/trivy/latest/getting-started/installation/
    exit /b 1
)

:: Run Trivy scan
trivy image --format "%FORMAT%" --output "%OUTPUT%" --severity "%SEVERITY%" "%IMAGE_NAME%"

:: Check for vulnerabilities based on FAIL_ON severity
if not "%FAIL_ON%"=="NONE" (
    if "%FORMAT%"=="json" (
        :: For JSON format, use findstr
        findstr /C:"\"Severity\": \"%FAIL_ON%\"" "%OUTPUT%" >nul
        if not errorlevel 1 (
            echo Found %FAIL_ON% severity vulnerabilities!
            if not "%EXIT_CODE%"=="0" (
                echo Failing build as requested (exit code %EXIT_CODE%)
                exit /b %EXIT_CODE%
            ) else (
                echo Warning: Vulnerabilities found but continuing build (exit code 0)
            )
        ) else (
            echo No %FAIL_ON% severity vulnerabilities found.
        )
    ) else (
        :: For other formats, use findstr
        findstr /C:"%FAIL_ON%" "%OUTPUT%" >nul
        if not errorlevel 1 (
            echo Found %FAIL_ON% severity vulnerabilities!
            if not "%EXIT_CODE%"=="0" (
                echo Failing build as requested (exit code %EXIT_CODE%)
                exit /b %EXIT_CODE%
            ) else (
                echo Warning: Vulnerabilities found but continuing build (exit code 0)
            )
        ) else (
            echo No %FAIL_ON% severity vulnerabilities found.
        )
    )
)

echo Scan completed successfully.
exit /b 0