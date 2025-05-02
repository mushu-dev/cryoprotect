@echo off
:: generate_sbom.bat - Generate Software Bill of Materials (SBOM) for Docker images on Windows
:: 
:: This script generates an SBOM for Docker images using Trivy
:: and can store it in various formats and locations.
::
:: Usage: generate_sbom.bat [options] IMAGE_NAME
::
:: Options:
::   --format FORMAT       Output format (cyclonedx, spdx, json) [default: cyclonedx]
::   --output FILE         Output file name [default: sbom.{format}]
::   --store-dir DIR       Directory to store SBOM files [default: .\sbom]
::   --version VERSION     Version tag for the SBOM file [default: latest]
::   --help                Show this help message and exit

setlocal enabledelayedexpansion

:: Default values
set "FORMAT=cyclonedx"
set "OUTPUT="
set "STORE_DIR=.\sbom"
set "VERSION=latest"
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
if "%~1"=="--store-dir" (
    set "STORE_DIR=%~2"
    shift /1
    shift /1
    goto :parse_args
)
if "%~1"=="--version" (
    set "VERSION=%~2"
    shift /1
    shift /1
    goto :parse_args
)
if "%~1"=="--help" (
    echo Usage: generate_sbom.bat [options] IMAGE_NAME
    echo.
    echo Options:
    echo   --format FORMAT       Output format (cyclonedx, spdx, json) [default: cyclonedx]
    echo   --output FILE         Output file name [default: sbom.{format}]
    echo   --store-dir DIR       Directory to store SBOM files [default: .\sbom]
    echo   --version VERSION     Version tag for the SBOM file [default: latest]
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
    echo Run 'generate_sbom.bat --help' for usage information
    exit /b 1
)

:: Set default output file if not provided
if "%OUTPUT%"=="" (
    if "%FORMAT%"=="cyclonedx" (
        set "OUTPUT=sbom.json"
    ) else if "%FORMAT%"=="spdx" (
        set "OUTPUT=sbom.spdx.json"
    ) else if "%FORMAT%"=="json" (
        set "OUTPUT=sbom.json"
    ) else (
        set "OUTPUT=sbom.json"
    )
)

echo Generating SBOM for Docker image: %IMAGE_NAME%
echo Output format: %FORMAT%
echo Output file: %OUTPUT%
echo Version: %VERSION%

:: Check if Trivy is installed
trivy --version >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo Error: Trivy is not installed
    echo Please install Trivy: https://aquasecurity.github.io/trivy/latest/getting-started/installation/
    exit /b 1
)

:: Create storage directory if it doesn't exist
if not exist "%STORE_DIR%" mkdir "%STORE_DIR%"

:: Generate SBOM
trivy image --format "%FORMAT%" --output "%OUTPUT%" "%IMAGE_NAME%"

:: Create a versioned copy in the storage directory
for /f "tokens=2 delims==" %%a in ('wmic os get localdatetime /format:list') do set "DATETIME=%%a"
set "TIMESTAMP=%DATETIME:~0,14%"

:: Get filename and extension
for %%i in ("%OUTPUT%") do (
    set "FILENAME=%%~ni"
    set "EXTENSION=%%~xi"
)

:: Create versioned filename (replace / and : with _)
set "VERSIONED_FILENAME=%IMAGE_NAME:/=_%_%VERSION%_%TIMESTAMP%%EXTENSION%"
set "VERSIONED_FILENAME=%VERSIONED_FILENAME::=_%"

:: Copy to versioned file
copy /Y "%OUTPUT%" "%STORE_DIR%\%VERSIONED_FILENAME%"

:: Create a latest version
set "LATEST_FILENAME=%IMAGE_NAME:/=_%_latest%EXTENSION%"
set "LATEST_FILENAME=%LATEST_FILENAME::=_%"

if exist "%STORE_DIR%\%LATEST_FILENAME%" del "%STORE_DIR%\%LATEST_FILENAME%"
copy /Y "%OUTPUT%" "%STORE_DIR%\%LATEST_FILENAME%"

echo SBOM generated successfully.
echo Stored in: %STORE_DIR%\%VERSIONED_FILENAME%
echo Latest version: %STORE_DIR%\%LATEST_FILENAME%

:: Generate metadata file
set "METADATA_FILENAME=%VERSIONED_FILENAME:.json=.metadata.json%"
set "METADATA_FILE=%STORE_DIR%\%METADATA_FILENAME%"

echo { > "%METADATA_FILE%"
echo   "image": "%IMAGE_NAME%", >> "%METADATA_FILE%"
echo   "version": "%VERSION%", >> "%METADATA_FILE%"
echo   "timestamp": "%TIMESTAMP:~0,4%-%TIMESTAMP:~4,2%-%TIMESTAMP:~6,2%T%TIMESTAMP:~8,2%:%TIMESTAMP:~10,2%:%TIMESTAMP:~12,2%Z", >> "%METADATA_FILE%"
echo   "format": "%FORMAT%", >> "%METADATA_FILE%"
echo   "sbom_file": "%VERSIONED_FILENAME%" >> "%METADATA_FILE%"
echo } >> "%METADATA_FILE%"

echo Metadata stored in: %METADATA_FILE%
exit /b 0