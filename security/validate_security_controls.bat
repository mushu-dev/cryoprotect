@echo off
:: CryoProtect v2 Security Controls Validation Script (Windows)
::
:: This script runs the comprehensive security validation process:
:: 1. Executes automated security tests
:: 2. Runs penetration tests
:: 3. Generates validation reports
::
:: Usage: validate_security_controls.bat [--url URL] [--output-dir DIR] [--html] [--pdf]

setlocal enabledelayedexpansion

:: Default values
set URL=http://localhost:5000
set OUTPUT_DIR=reports\security
set GENERATE_HTML=false
set GENERATE_PDF=false
for /f "tokens=2-4 delims=/ " %%a in ('date /t') do (set DATESTAMP=%%c%%a%%b)
for /f "tokens=1-2 delims=: " %%a in ('time /t') do (set TIMESTAMP=%%a%%b)
set REPORT_PREFIX=security_validation_%DATESTAMP%_%TIMESTAMP%

:: Parse command line arguments
:parse_args
if "%~1"=="" goto :end_parse_args
if "%~1"=="--url" (
    set URL=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--output-dir" (
    set OUTPUT_DIR=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--html" (
    set GENERATE_HTML=true
    shift
    goto :parse_args
)
if "%~1"=="--pdf" (
    set GENERATE_PDF=true
    shift
    goto :parse_args
)
if "%~1"=="--help" (
    echo Usage: %0 [--url URL] [--output-dir DIR] [--html] [--pdf]
    echo.
    echo Options:
    echo   --url URL         Base URL of the application to test (default: http://localhost:5000)
    echo   --output-dir DIR  Directory to save reports (default: reports\security)
    echo   --html            Generate HTML report in addition to JSON
    echo   --pdf             Generate PDF report in addition to JSON
    echo   --help            Show this help message
    exit /b 0
)
echo Unknown option: %~1
echo Use --help for usage information
exit /b 1

:end_parse_args

:: Create output directory if it doesn't exist
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"
echo Output directory: %OUTPUT_DIR%

:: Function to print section headers
:print_section
echo.
echo ================================================================================
echo   %~1
echo ================================================================================
echo.
goto :eof

:: Start validation process
call :print_section "Starting CryoProtect v2 Security Controls Validation"
echo Timestamp: %date% %time%
echo Target URL: %URL%
echo Output directory: %OUTPUT_DIR%

:: Step 1: Run automated security tests
call :print_section "Running Automated Security Tests"
set TEST_OUTPUT_FILE=%OUTPUT_DIR%\%REPORT_PREFIX%_test_results.txt
python -m pytest tests/test_security_controls.py -v > "%TEST_OUTPUT_FILE%" 2>&1
set TEST_EXIT_STATUS=%ERRORLEVEL%

echo Test results saved to: %TEST_OUTPUT_FILE%
if %TEST_EXIT_STATUS% EQU 0 (
    echo Test status: PASSED
) else (
    echo Test status: FAILED (Exit code: %TEST_EXIT_STATUS%)
)

:: Step 2: Run penetration tests
call :print_section "Running Penetration Tests"
set PENTEST_ARGS=--url %URL% --output-dir %OUTPUT_DIR%
if "%GENERATE_HTML%"=="true" set PENTEST_ARGS=%PENTEST_ARGS% --html
if "%GENERATE_PDF%"=="true" set PENTEST_ARGS=%PENTEST_ARGS% --pdf

echo Running: python security/run_pentest.py %PENTEST_ARGS%
python security/run_pentest.py %PENTEST_ARGS%
set PENTEST_EXIT_STATUS=%ERRORLEVEL%

if %PENTEST_EXIT_STATUS% EQU 0 (
    echo Penetration test status: PASSED
) else (
    echo Penetration test status: FAILED (Exit code: %PENTEST_EXIT_STATUS%)
)

:: Step 3: Generate validation summary
call :print_section "Generating Validation Summary"
set SUMMARY_FILE=%OUTPUT_DIR%\%REPORT_PREFIX%_summary.md

:: Create summary file
echo # CryoProtect v2 Security Controls Validation Summary > "%SUMMARY_FILE%"
echo. >> "%SUMMARY_FILE%"
echo **Date:** %date% %time% >> "%SUMMARY_FILE%"
echo **Target:** %URL% >> "%SUMMARY_FILE%"
echo. >> "%SUMMARY_FILE%"
echo ## Validation Results >> "%SUMMARY_FILE%"
echo. >> "%SUMMARY_FILE%"
echo ^| Component ^| Status ^| Details ^| >> "%SUMMARY_FILE%"
echo ^|-----------|--------|---------|  >> "%SUMMARY_FILE%"

:: Add test results
if %TEST_EXIT_STATUS% EQU 0 (
    echo ^| Automated Tests ^| ✅ PASSED ^| See [Test Results](%REPORT_PREFIX%_test_results.txt) ^| >> "%SUMMARY_FILE%"
) else (
    echo ^| Automated Tests ^| ❌ FAILED ^| See [Test Results](%REPORT_PREFIX%_test_results.txt) ^| >> "%SUMMARY_FILE%"
)

:: Add penetration test results
if %PENTEST_EXIT_STATUS% EQU 0 (
    echo ^| Penetration Tests ^| ✅ PASSED ^| See [Pentest Report](pentest_summary.md) ^| >> "%SUMMARY_FILE%"
) else (
    echo ^| Penetration Tests ^| ❌ FAILED ^| See [Pentest Report](pentest_summary.md) ^| >> "%SUMMARY_FILE%"
)

:: Add overall validation status
if %TEST_EXIT_STATUS% EQU 0 if %PENTEST_EXIT_STATUS% EQU 0 (
    echo ^| Overall Validation ^| ✅ PASSED ^| ^| >> "%SUMMARY_FILE%"
) else (
    echo ^| Overall Validation ^| ❌ FAILED ^| ^| >> "%SUMMARY_FILE%"
)

echo. >> "%SUMMARY_FILE%"
echo ## Security Controls Status >> "%SUMMARY_FILE%"
echo. >> "%SUMMARY_FILE%"
echo ^| Security Control ^| Status ^| Validation Method ^| >> "%SUMMARY_FILE%"
echo ^|------------------|--------|------------------|  >> "%SUMMARY_FILE%"

:: Add security controls status
if %TEST_EXIT_STATUS% EQU 0 if %PENTEST_EXIT_STATUS% EQU 0 (
    echo ^| CSRF Protection ^| ✅ Validated ^| Automated Tests, Penetration Testing ^| >> "%SUMMARY_FILE%"
    echo ^| Security Headers ^| ✅ Validated ^| Automated Tests, Penetration Testing ^| >> "%SUMMARY_FILE%"
    echo ^| Cookie Security ^| ✅ Validated ^| Automated Tests, Penetration Testing ^| >> "%SUMMARY_FILE%"
    echo ^| Encryption at Rest ^| ✅ Validated ^| Automated Tests, Penetration Testing ^| >> "%SUMMARY_FILE%"
    echo ^| Vulnerability Scanning ^| ✅ Validated ^| Automated Tests, Penetration Testing ^| >> "%SUMMARY_FILE%"
) else (
    echo ^| CSRF Protection ^| ⚠️ Needs Review ^| Automated Tests, Penetration Testing ^| >> "%SUMMARY_FILE%"
    echo ^| Security Headers ^| ⚠️ Needs Review ^| Automated Tests, Penetration Testing ^| >> "%SUMMARY_FILE%"
    echo ^| Cookie Security ^| ⚠️ Needs Review ^| Automated Tests, Penetration Testing ^| >> "%SUMMARY_FILE%"
    echo ^| Encryption at Rest ^| ⚠️ Needs Review ^| Automated Tests, Penetration Testing ^| >> "%SUMMARY_FILE%"
    echo ^| Vulnerability Scanning ^| ⚠️ Needs Review ^| Automated Tests, Penetration Testing ^| >> "%SUMMARY_FILE%"
)

echo. >> "%SUMMARY_FILE%"
echo ## Next Steps >> "%SUMMARY_FILE%"
echo. >> "%SUMMARY_FILE%"

if %TEST_EXIT_STATUS% EQU 0 if %PENTEST_EXIT_STATUS% EQU 0 (
    echo All security controls have been successfully validated. The application is ready for production use. >> "%SUMMARY_FILE%"
) else (
    echo Some security controls need review. Please check the detailed reports and address any issues before proceeding to production. >> "%SUMMARY_FILE%"
)

echo. >> "%SUMMARY_FILE%"
echo ## References >> "%SUMMARY_FILE%"
echo. >> "%SUMMARY_FILE%"
echo - [Security Validation Report](security_validation_report.md) >> "%SUMMARY_FILE%"
echo - [Test Results](%REPORT_PREFIX%_test_results.txt) >> "%SUMMARY_FILE%"
echo - [Penetration Test Report](pentest_summary.md) >> "%SUMMARY_FILE%"

echo Validation summary saved to: %SUMMARY_FILE%

:: Step 4: Copy the security validation report to the output directory
call :print_section "Finalizing Reports"
if exist "reports\security\security_validation_report.md" (
    copy "reports\security\security_validation_report.md" "%OUTPUT_DIR%\%REPORT_PREFIX%_report.md" > nul
    echo Security validation report copied to: %OUTPUT_DIR%\%REPORT_PREFIX%_report.md
) else (
    echo WARNING: Security validation report not found at reports\security\security_validation_report.md
)

:: Step 5: Determine overall validation status
call :print_section "Validation Complete"
if %TEST_EXIT_STATUS% EQU 0 if %PENTEST_EXIT_STATUS% EQU 0 (
    echo ✅ SUCCESS: All security controls have been successfully validated.
    set VALIDATION_STATUS=0
) else (
    echo ❌ FAILURE: Some security controls need review. Please check the detailed reports.
    set VALIDATION_STATUS=1
)

echo.
echo Reports generated:
echo - Validation Summary: %SUMMARY_FILE%
echo - Test Results: %TEST_OUTPUT_FILE%
echo - Penetration Test Reports: See %OUTPUT_DIR%\pentest_*.json
if "%GENERATE_HTML%"=="true" echo - HTML Reports: See %OUTPUT_DIR%\pentest_*.html
if "%GENERATE_PDF%"=="true" echo - PDF Reports: See %OUTPUT_DIR%\pentest_*.pdf

echo.
echo Validation process completed at %date% %time%
exit /b %VALIDATION_STATUS%