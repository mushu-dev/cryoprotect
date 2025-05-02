@echo off
REM setup_codex_cli.bat - Configures the environment for Codex CLI usage

echo Setting up Codex CLI environment...

REM Check if Node.js is installed
where node >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Node.js is not installed or not in PATH.
    echo Please download and install Node.js from https://nodejs.org/
    exit /b 1
)

REM Check if npm is installed
where npm >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: npm is not installed or not in PATH.
    echo Please ensure npm is installed with Node.js.
    exit /b 1
)

REM Check if Codex CLI is installed
where codex >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Codex CLI not found. Installing globally...
    npm install -g openai-codex
    
    REM Check installation success
    where codex >nul 2>nul
    if %ERRORLEVEL% NEQ 0 (
        echo ERROR: Failed to install Codex CLI.
        exit /b 1
    )
)

REM Check if OPENAI_API_KEY is set
if "%OPENAI_API_KEY%"=="" (
    echo WARNING: OPENAI_API_KEY environment variable is not set.
    
    REM Check if .env file exists and has the key
    if exist .env (
        findstr /C:"OPENAI_API_KEY" .env >nul
        if %ERRORLEVEL% EQU 0 (
            echo OPENAI_API_KEY found in .env file. Loading from .env...
            for /F "tokens=*" %%a in ('findstr /C:"OPENAI_API_KEY" .env') do set %%a
        ) else (
            echo No OPENAI_API_KEY found in .env file.
            echo Please set your OpenAI API key in the .env file or as an environment variable.
            exit /b 1
        )
    ) else (
        echo No .env file found.
        echo Please set your OpenAI API key as an environment variable or create an .env file.
        exit /b 1
    )
)

REM Test Codex CLI with a simple prompt
echo Testing Codex CLI...
codex "Write 'Hello, World!' in Python." > nul
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Codex CLI test failed.
    echo Please verify your API key and internet connection.
    exit /b 1
)

echo Codex CLI setup successful!
echo.
echo Usage:
echo   1. Run codex_analysis_helper.py to analyze your codebase:
echo      python codex_analysis_helper.py analyze-codebase --output analysis_report.md
echo.
echo   2. Create an improvement plan based on the analysis:
echo      python codex_analysis_helper.py create-plan --input analysis_report.md --output improvement_plan.json
echo.
echo   3. Delegate tasks to Roocode agents:
echo      python codex_analysis_helper.py delegate-tasks --input improvement_plan.json
echo.
echo   4. Or use the Codex Analysis mode directly in Roocode:
echo      roo "Analyze the codebase and suggest improvements" --mode codex-analysis
echo.
echo You can also use the Codex CLI directly:
echo   codex "Your prompt here"
echo.

exit /b 0
