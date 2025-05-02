@echo off
REM scan_dependency_check.bat - OWASP Dependency-Check scanner for Windows
REM 
REM This script scans project dependencies for known vulnerabilities using OWASP Dependency-Check.
REM It can be run as a standalone script or integrated into CI/CD pipelines.
REM
REM Usage: scan_dependency_check.bat [options]
REM
REM Options:
REM   --path PATH           Path to scan (default: current directory)
REM   --format FORMAT       Output format (HTML, XML, CSV, JSON, JUNIT, SARIF) (default: HTML)
REM   --output DIR          Output directory (default: .\dependency-check-reports)
REM   --name NAME           Project name (default: CryoProtect)
REM   --exit-on-critical    Exit with code 1 if critical vulnerabilities found (CVSS >= 7.0)
REM   --nvd-api-key KEY     NVD API key for better performance
REM   --help                Show this help message and exit

setlocal enabledelayedexpansion

REM Default values
set SCAN_PATH=.
set FORMAT=HTML
set OUTPUT_DIR=.\dependency-check-reports
set PROJECT_NAME=CryoProtect
set EXIT_ON_CRITICAL=0
set NVD_API_KEY=

REM Parse arguments
:parse_args
if "%~1"=="" goto :end_parse_args
if "%~1"=="--path" (
    set SCAN_PATH=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--format" (
    set FORMAT=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--output" (
    set OUTPUT_DIR=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--name" (
    set PROJECT_NAME=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--exit-on-critical" (
    set EXIT_ON_CRITICAL=1
    shift
    goto :parse_args
)
if "%~1"=="--nvd-api-key" (
    set NVD_API_KEY=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--help" (
    echo Usage: scan_dependency_check.bat [options]
    echo.
    echo Options:
    echo   --path PATH           Path to scan (default: current directory)
    echo   --format FORMAT       Output format (HTML, XML, CSV, JSON, JUNIT, SARIF) (default: HTML)
    echo   --output DIR          Output directory (default: .\dependency-check-reports)
    echo   --name NAME           Project name (default: CryoProtect)
    echo   --exit-on-critical    Exit with code 1 if critical vulnerabilities found (CVSS ^>= 7.0)
    echo   --nvd-api-key KEY     NVD API key for better performance
    echo   --help                Show this help message and exit
    exit /b 0
)
echo Unknown option: %~1
echo Run 'scan_dependency_check.bat --help' for usage information
exit /b 1

:end_parse_args

REM Create output directory
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

REM Generate timestamp for reports
for /f "tokens=2 delims==" %%a in ('wmic OS Get localdatetime /value') do set "dt=%%a"
set "TIMESTAMP=%dt:~0,8%%dt:~8,6%"
set "REPORT_NAME=%PROJECT_NAME%-%TIMESTAMP%"

REM Check if Docker is available
where docker >nul 2>nul
if %ERRORLEVEL% equ 0 (
    echo Using Docker to run OWASP Dependency-Check...
    
    REM Prepare Docker command
    set "DOCKER_CMD=docker run --rm"
    set "DOCKER_CMD=!DOCKER_CMD! -v "%cd%\%SCAN_PATH%:/src""
    set "DOCKER_CMD=!DOCKER_CMD! -v "%cd%\%OUTPUT_DIR%:/report""
    set "DOCKER_CMD=!DOCKER_CMD! -v "%cd%\dependency-check-data:/usr/share/dependency-check/data""
    
    REM Add NVD API key if provided
    if not "%NVD_API_KEY%"=="" (
        set "DOCKER_CMD=!DOCKER_CMD! -e "NVD_API_KEY=%NVD_API_KEY%""
    )
    
    REM Complete the command
    set "DOCKER_CMD=!DOCKER_CMD! owasp/dependency-check:latest"
    set "DOCKER_CMD=!DOCKER_CMD! --scan /src"
    set "DOCKER_CMD=!DOCKER_CMD! --format %FORMAT%"
    set "DOCKER_CMD=!DOCKER_CMD! --project "%PROJECT_NAME%""
    set "DOCKER_CMD=!DOCKER_CMD! --out /report"
    set "DOCKER_CMD=!DOCKER_CMD! --enableExperimental"
    
    REM Run the scan
    echo Running OWASP Dependency-Check scan on %SCAN_PATH%...
    %DOCKER_CMD%
    
    REM Rename the report files with timestamp
    if exist "%OUTPUT_DIR%\dependency-check-report.html" (
        move "%OUTPUT_DIR%\dependency-check-report.html" "%OUTPUT_DIR%\%REPORT_NAME%.html"
    )
    
    if exist "%OUTPUT_DIR%\dependency-check-report.json" (
        move "%OUTPUT_DIR%\dependency-check-report.json" "%OUTPUT_DIR%\%REPORT_NAME%.json"
    )
    
    if exist "%OUTPUT_DIR%\dependency-check-report.xml" (
        move "%OUTPUT_DIR%\dependency-check-report.xml" "%OUTPUT_DIR%\%REPORT_NAME%.xml"
    )
    
    if exist "%OUTPUT_DIR%\dependency-check-report.csv" (
        move "%OUTPUT_DIR%\dependency-check-report.csv" "%OUTPUT_DIR%\%REPORT_NAME%.csv"
    )
    
    if exist "%OUTPUT_DIR%\dependency-check-junit.xml" (
        move "%OUTPUT_DIR%\dependency-check-junit.xml" "%OUTPUT_DIR%\%REPORT_NAME%-junit.xml"
    )
    
    if exist "%OUTPUT_DIR%\dependency-check-sarif.json" (
        move "%OUTPUT_DIR%\dependency-check-sarif.json" "%OUTPUT_DIR%\%REPORT_NAME%-sarif.json"
    )
    
    echo Scan completed. Reports saved to %OUTPUT_DIR%
    
    REM Check for critical vulnerabilities if requested
    if %EXIT_ON_CRITICAL% equ 1 (
        if exist "%OUTPUT_DIR%\%REPORT_NAME%.json" (
            echo Checking for critical vulnerabilities (CVSS ^>= 7.0)...
            
            REM Use findstr to check for critical vulnerabilities
            findstr /C:"\"baseScore\": [7-9]" "%OUTPUT_DIR%\%REPORT_NAME%.json" > nul
            if !ERRORLEVEL! equ 0 (
                echo Found critical vulnerabilities (CVSS ^>= 7.0)!
                echo See the report for details: %OUTPUT_DIR%\%REPORT_NAME%.html
                exit /b 1
            ) else (
                echo No critical vulnerabilities found.
            )
        )
    )
    
) else (
    REM Check if dependency-check script is installed
    where dependency-check.bat >nul 2>nul
    if %ERRORLEVEL% equ 0 (
        set DEPENDENCY_CHECK_CMD=dependency-check.bat
    ) else if exist ".\dependency-check\bin\dependency-check.bat" (
        set DEPENDENCY_CHECK_CMD=.\dependency-check\bin\dependency-check.bat
    ) else (
        echo Error: OWASP Dependency-Check not found and Docker is not available.
        echo Please install OWASP Dependency-Check or Docker.
        echo Installation instructions: https://jeremylong.github.io/DependencyCheck/dependency-check-cli/index.html
        exit /b 1
    )
    
    REM Prepare command
    set "CMD=%DEPENDENCY_CHECK_CMD%"
    set "CMD=!CMD! --scan "%SCAN_PATH%""
    set "CMD=!CMD! --format %FORMAT%"
    set "CMD=!CMD! --project "%PROJECT_NAME%""
    set "CMD=!CMD! --out "%OUTPUT_DIR%""
    set "CMD=!CMD! --enableExperimental"
    
    REM Add NVD API key if provided
    if not "%NVD_API_KEY%"=="" (
        set "CMD=!CMD! --nvdApiKey "%NVD_API_KEY%""
    )
    
    REM Run the scan
    echo Running OWASP Dependency-Check scan on %SCAN_PATH%...
    %CMD%
    
    REM Rename the report files with timestamp
    if exist "%OUTPUT_DIR%\dependency-check-report.html" (
        move "%OUTPUT_DIR%\dependency-check-report.html" "%OUTPUT_DIR%\%REPORT_NAME%.html"
    )
    
    if exist "%OUTPUT_DIR%\dependency-check-report.json" (
        move "%OUTPUT_DIR%\dependency-check-report.json" "%OUTPUT_DIR%\%REPORT_NAME%.json"
    )
    
    if exist "%OUTPUT_DIR%\dependency-check-report.xml" (
        move "%OUTPUT_DIR%\dependency-check-report.xml" "%OUTPUT_DIR%\%REPORT_NAME%.xml"
    )
    
    if exist "%OUTPUT_DIR%\dependency-check-report.csv" (
        move "%OUTPUT_DIR%\dependency-check-report.csv" "%OUTPUT_DIR%\%REPORT_NAME%.csv"
    )
    
    if exist "%OUTPUT_DIR%\dependency-check-junit.xml" (
        move "%OUTPUT_DIR%\dependency-check-junit.xml" "%OUTPUT_DIR%\%REPORT_NAME%-junit.xml"
    )
    
    if exist "%OUTPUT_DIR%\dependency-check-sarif.json" (
        move "%OUTPUT_DIR%\dependency-check-sarif.json" "%OUTPUT_DIR%\%REPORT_NAME%-sarif.json"
    )
    
    echo Scan completed. Reports saved to %OUTPUT_DIR%
    
    REM Check for critical vulnerabilities if requested
    if %EXIT_ON_CRITICAL% equ 1 (
        if exist "%OUTPUT_DIR%\%REPORT_NAME%.json" (
            echo Checking for critical vulnerabilities (CVSS ^>= 7.0)...
            
            REM Use findstr to check for critical vulnerabilities
            findstr /C:"\"baseScore\": [7-9]" "%OUTPUT_DIR%\%REPORT_NAME%.json" > nul
            if !ERRORLEVEL! equ 0 (
                echo Found critical vulnerabilities (CVSS ^>= 7.0)!
                echo See the report for details: %OUTPUT_DIR%\%REPORT_NAME%.html
                exit /b 1
            ) else (
                echo No critical vulnerabilities found.
            )
        )
    )
)

REM Generate summary file
set "SUMMARY_FILE=%OUTPUT_DIR%\%REPORT_NAME%-summary.json"

echo { > "%SUMMARY_FILE%"
echo   "scanner": "owasp-dependency-check", >> "%SUMMARY_FILE%"
echo   "timestamp": "%date% %time%", >> "%SUMMARY_FILE%"
echo   "project": "%PROJECT_NAME%", >> "%SUMMARY_FILE%"
echo   "scan_path": "%SCAN_PATH%", >> "%SUMMARY_FILE%"
echo   "report_files": { >> "%SUMMARY_FILE%"

REM Add report files to summary
set "COMMA="
if exist "%OUTPUT_DIR%\%REPORT_NAME%.html" (
    echo     "html": "%REPORT_NAME%.html"%COMMA% >> "%SUMMARY_FILE%"
    set "COMMA=,"
)

if exist "%OUTPUT_DIR%\%REPORT_NAME%.json" (
    echo     "json": "%REPORT_NAME%.json"%COMMA% >> "%SUMMARY_FILE%"
    set "COMMA=,"
)

if exist "%OUTPUT_DIR%\%REPORT_NAME%.xml" (
    echo     "xml": "%REPORT_NAME%.xml"%COMMA% >> "%SUMMARY_FILE%"
    set "COMMA=,"
)

if exist "%OUTPUT_DIR%\%REPORT_NAME%.csv" (
    echo     "csv": "%REPORT_NAME%.csv"%COMMA% >> "%SUMMARY_FILE%"
    set "COMMA=,"
)

if exist "%OUTPUT_DIR%\%REPORT_NAME%-junit.xml" (
    echo     "junit": "%REPORT_NAME%-junit.xml"%COMMA% >> "%SUMMARY_FILE%"
    set "COMMA=,"
)

if exist "%OUTPUT_DIR%\%REPORT_NAME%-sarif.json" (
    echo     "sarif": "%REPORT_NAME%-sarif.json" >> "%SUMMARY_FILE%"
)

echo   } >> "%SUMMARY_FILE%"
echo } >> "%SUMMARY_FILE%"

echo Summary saved to %SUMMARY_FILE%

exit /b 0