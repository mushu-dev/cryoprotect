@echo off
REM CryoProtect v2 Secret Management Test Script for Windows
REM This script tests the secret management system to ensure it's working correctly

setlocal enabledelayedexpansion

REM ANSI color codes for Windows 10+
set RED=[91m
set GREEN=[92m
set YELLOW=[93m
set BLUE=[94m
set NC=[0m

REM Function to log messages with timestamp
call :log_message "%BLUE%CryoProtect v2 Secret Management Test Script for Windows%NC%"
call :log_message "This script will test the secret management system to ensure it's working correctly."
echo.

REM Check if Docker is installed
where docker >nul 2>&1
if %ERRORLEVEL% neq 0 (
    call :log_message "%RED%Error: Docker is not installed or not in PATH%NC%"
    exit /b 1
)

REM Check if Docker Compose is installed
where docker-compose >nul 2>&1
if %ERRORLEVEL% neq 0 (
    call :log_message "%RED%Error: Docker Compose is not installed or not in PATH%NC%"
    exit /b 1
)

REM Create a temporary directory for test secrets
set TEST_DIR=%TEMP%\cryoprotect-secret-test-%RANDOM%
mkdir "%TEST_DIR%"
call :log_message "Created temporary directory for test secrets: %TEST_DIR%"

REM Register cleanup on exit
REM Note: This doesn't work the same as in bash, but we'll clean up at the end

REM Create test secrets
call :log_message "Creating test secrets..."
echo test-value-1> "%TEST_DIR%\test_secret_1"
echo test-value-2> "%TEST_DIR%\test_secret_2"
for /f "tokens=*" %%a in ('powershell -Command "Get-Date -UFormat %%s"') do set timestamp=%%a
echo %timestamp%> "%TEST_DIR%\test_secret_rotation"

REM Create Docker secrets
call :log_message "Creating Docker secrets..."
type "%TEST_DIR%\test_secret_1" | docker secret create test_secret_1 - 2>nul
if %ERRORLEVEL% neq 0 (
    call :log_message "%YELLOW%Warning: Could not create Docker secret. Are you running in swarm mode?%NC%"
)
type "%TEST_DIR%\test_secret_2" | docker secret create test_secret_2 - 2>nul
if %ERRORLEVEL% neq 0 (
    call :log_message "%YELLOW%Warning: Could not create Docker secret. Are you running in swarm mode?%NC%"
)
type "%TEST_DIR%\test_secret_rotation" | docker secret create test_secret_rotation - 2>nul
if %ERRORLEVEL% neq 0 (
    call :log_message "%YELLOW%Warning: Could not create Docker secret. Are you running in swarm mode?%NC%"
)

REM Create a test docker-compose.yml file
echo version: '3.8'> "%TEST_DIR%\docker-compose.test.yml"
echo.>> "%TEST_DIR%\docker-compose.test.yml"
echo services:>> "%TEST_DIR%\docker-compose.test.yml"
echo   secret-test:>> "%TEST_DIR%\docker-compose.test.yml"
echo     image: alpine:latest>> "%TEST_DIR%\docker-compose.test.yml"
echo     command: sh -c 'cat /run/secrets/TEST_SECRET_1; echo; cat /run/secrets/TEST_SECRET_2; echo; cat /run/secrets/SECRET_ROTATION_TIMESTAMP; echo; sleep 5'>> "%TEST_DIR%\docker-compose.test.yml"
echo     secrets:>> "%TEST_DIR%\docker-compose.test.yml"
echo       - source: test_secret_1>> "%TEST_DIR%\docker-compose.test.yml"
echo         target: TEST_SECRET_1>> "%TEST_DIR%\docker-compose.test.yml"
echo       - source: test_secret_2>> "%TEST_DIR%\docker-compose.test.yml"
echo         target: TEST_SECRET_2>> "%TEST_DIR%\docker-compose.test.yml"
echo       - source: test_secret_rotation>> "%TEST_DIR%\docker-compose.test.yml"
echo         target: SECRET_ROTATION_TIMESTAMP>> "%TEST_DIR%\docker-compose.test.yml"
echo.>> "%TEST_DIR%\docker-compose.test.yml"
echo secrets:>> "%TEST_DIR%\docker-compose.test.yml"
echo   test_secret_1:>> "%TEST_DIR%\docker-compose.test.yml"
echo     external: true>> "%TEST_DIR%\docker-compose.test.yml"
echo   test_secret_2:>> "%TEST_DIR%\docker-compose.test.yml"
echo     external: true>> "%TEST_DIR%\docker-compose.test.yml"
echo   test_secret_rotation:>> "%TEST_DIR%\docker-compose.test.yml"
echo     external: true>> "%TEST_DIR%\docker-compose.test.yml"

REM Create a test docker-compose-dev.yml file for file-based secrets
mkdir "%TEST_DIR%\.secrets"
copy "%TEST_DIR%\test_secret_1" "%TEST_DIR%\.secrets\TEST_SECRET_1" >nul
copy "%TEST_DIR%\test_secret_2" "%TEST_DIR%\.secrets\TEST_SECRET_2" >nul
copy "%TEST_DIR%\test_secret_rotation" "%TEST_DIR%\.secrets\SECRET_ROTATION_TIMESTAMP" >nul

echo version: '3.8'> "%TEST_DIR%\docker-compose-dev.yml"
echo.>> "%TEST_DIR%\docker-compose-dev.yml"
echo services:>> "%TEST_DIR%\docker-compose-dev.yml"
echo   secret-test-dev:>> "%TEST_DIR%\docker-compose-dev.yml"
echo     image: alpine:latest>> "%TEST_DIR%\docker-compose-dev.yml"
echo     command: sh -c 'cat /run/secrets/TEST_SECRET_1; echo; cat /run/secrets/TEST_SECRET_2; echo; cat /run/secrets/SECRET_ROTATION_TIMESTAMP; echo; sleep 5'>> "%TEST_DIR%\docker-compose-dev.yml"
echo     secrets:>> "%TEST_DIR%\docker-compose-dev.yml"
echo       - source: test_secret_1>> "%TEST_DIR%\docker-compose-dev.yml"
echo         target: TEST_SECRET_1>> "%TEST_DIR%\docker-compose-dev.yml"
echo       - source: test_secret_2>> "%TEST_DIR%\docker-compose-dev.yml"
echo         target: TEST_SECRET_2>> "%TEST_DIR%\docker-compose-dev.yml"
echo       - source: test_secret_rotation>> "%TEST_DIR%\docker-compose-dev.yml"
echo         target: SECRET_ROTATION_TIMESTAMP>> "%TEST_DIR%\docker-compose-dev.yml"
echo.>> "%TEST_DIR%\docker-compose-dev.yml"
echo secrets:>> "%TEST_DIR%\docker-compose-dev.yml"
echo   test_secret_1:>> "%TEST_DIR%\docker-compose-dev.yml"
echo     file: ./.secrets/TEST_SECRET_1>> "%TEST_DIR%\docker-compose-dev.yml"
echo   test_secret_2:>> "%TEST_DIR%\docker-compose-dev.yml"
echo     file: ./.secrets/TEST_SECRET_2>> "%TEST_DIR%\docker-compose-dev.yml"
echo   test_secret_rotation:>> "%TEST_DIR%\docker-compose-dev.yml"
echo     file: ./.secrets/SECRET_ROTATION_TIMESTAMP>> "%TEST_DIR%\docker-compose-dev.yml"

REM Create a test entrypoint script
echo #!/bin/sh> "%TEST_DIR%\test-entrypoint.sh"
echo set -e>> "%TEST_DIR%\test-entrypoint.sh"
echo.>> "%TEST_DIR%\test-entrypoint.sh"
echo # Function to log messages with timestamp>> "%TEST_DIR%\test-entrypoint.sh"
echo log_message() {>> "%TEST_DIR%\test-entrypoint.sh"
echo   echo "[$(date -Iseconds)] $1">> "%TEST_DIR%\test-entrypoint.sh"
echo }>> "%TEST_DIR%\test-entrypoint.sh"
echo.>> "%TEST_DIR%\test-entrypoint.sh"
echo # Function to load secrets into environment variables>> "%TEST_DIR%\test-entrypoint.sh"
echo load_secrets() {>> "%TEST_DIR%\test-entrypoint.sh"
echo   log_message "Loading secrets...">> "%TEST_DIR%\test-entrypoint.sh"
echo.>> "%TEST_DIR%\test-entrypoint.sh"
echo   # Define required secrets>> "%TEST_DIR%\test-entrypoint.sh"
echo   REQUIRED_SECRETS=(>> "%TEST_DIR%\test-entrypoint.sh"
echo     "TEST_SECRET_1">> "%TEST_DIR%\test-entrypoint.sh"
echo     "TEST_SECRET_2">> "%TEST_DIR%\test-entrypoint.sh"
echo   )>> "%TEST_DIR%\test-entrypoint.sh"
echo.>> "%TEST_DIR%\test-entrypoint.sh"
echo   # Load each secret if it exists>> "%TEST_DIR%\test-entrypoint.sh"
echo   for SECRET in "${REQUIRED_SECRETS[@]}"; do>> "%TEST_DIR%\test-entrypoint.sh"
echo     SECRET_FILE="/run/secrets/${SECRET}">> "%TEST_DIR%\test-entrypoint.sh"
echo     if [ -f "$SECRET_FILE" ]; then>> "%TEST_DIR%\test-entrypoint.sh"
echo       # Read the secret value and set it as an environment variable>> "%TEST_DIR%\test-entrypoint.sh"
echo       export "$SECRET"="$(cat "$SECRET_FILE")">> "%TEST_DIR%\test-entrypoint.sh"
echo       log_message "Loaded secret: $SECRET = ${!SECRET}">> "%TEST_DIR%\test-entrypoint.sh"
echo     else>> "%TEST_DIR%\test-entrypoint.sh"
echo       log_message "Error: Required secret $SECRET not found">> "%TEST_DIR%\test-entrypoint.sh"
echo       exit 1>> "%TEST_DIR%\test-entrypoint.sh"
echo     fi>> "%TEST_DIR%\test-entrypoint.sh"
echo   done>> "%TEST_DIR%\test-entrypoint.sh"
echo.>> "%TEST_DIR%\test-entrypoint.sh"
echo   # Check for rotation timestamp>> "%TEST_DIR%\test-entrypoint.sh"
echo   ROTATION_FILE="/run/secrets/SECRET_ROTATION_TIMESTAMP">> "%TEST_DIR%\test-entrypoint.sh"
echo   if [ -f "$ROTATION_FILE" ]; then>> "%TEST_DIR%\test-entrypoint.sh"
echo     ROTATION_TIMESTAMP=$(cat "$ROTATION_FILE")>> "%TEST_DIR%\test-entrypoint.sh"
echo     CURRENT_TIMESTAMP=$(date +%%s)>> "%TEST_DIR%\test-entrypoint.sh"
echo.>> "%TEST_DIR%\test-entrypoint.sh"
echo     log_message "Secret rotation timestamp: $ROTATION_TIMESTAMP">> "%TEST_DIR%\test-entrypoint.sh"
echo.>> "%TEST_DIR%\test-entrypoint.sh"
echo     # Calculate age of secrets in days>> "%TEST_DIR%\test-entrypoint.sh"
echo     SECRET_AGE_SECONDS=$((CURRENT_TIMESTAMP - ROTATION_TIMESTAMP))>> "%TEST_DIR%\test-entrypoint.sh"
echo     SECRET_AGE_DAYS=$((SECRET_AGE_SECONDS / 86400))>> "%TEST_DIR%\test-entrypoint.sh"
echo.>> "%TEST_DIR%\test-entrypoint.sh"
echo     log_message "Secret age: $SECRET_AGE_DAYS days">> "%TEST_DIR%\test-entrypoint.sh"
echo   else>> "%TEST_DIR%\test-entrypoint.sh"
echo     log_message "Warning: No rotation timestamp found">> "%TEST_DIR%\test-entrypoint.sh"
echo   fi>> "%TEST_DIR%\test-entrypoint.sh"
echo.>> "%TEST_DIR%\test-entrypoint.sh"
echo   log_message "Secret loading complete">> "%TEST_DIR%\test-entrypoint.sh"
echo }>> "%TEST_DIR%\test-entrypoint.sh"
echo.>> "%TEST_DIR%\test-entrypoint.sh"
echo # Load secrets>> "%TEST_DIR%\test-entrypoint.sh"
echo load_secrets>> "%TEST_DIR%\test-entrypoint.sh"
echo.>> "%TEST_DIR%\test-entrypoint.sh"
echo # Execute the command passed to the container>> "%TEST_DIR%\test-entrypoint.sh"
echo log_message "Executing: $*">> "%TEST_DIR%\test-entrypoint.sh"
echo exec "$@">> "%TEST_DIR%\test-entrypoint.sh"

REM Create a test Dockerfile
echo FROM alpine:latest> "%TEST_DIR%\Dockerfile.test"
echo.>> "%TEST_DIR%\Dockerfile.test"
echo # Create non-root user>> "%TEST_DIR%\Dockerfile.test"
echo RUN adduser -D -u 1000 appuser ^&^& \>> "%TEST_DIR%\Dockerfile.test"
echo     mkdir -p /app /run/secrets ^&^& \>> "%TEST_DIR%\Dockerfile.test"
echo     chown -R appuser:appuser /app ^&^& \>> "%TEST_DIR%\Dockerfile.test"
echo     chown -R appuser:appuser /run/secrets ^&^& \>> "%TEST_DIR%\Dockerfile.test"
echo     chmod 750 /app ^&^& \>> "%TEST_DIR%\Dockerfile.test"
echo     chmod 700 /run/secrets>> "%TEST_DIR%\Dockerfile.test"
echo.>> "%TEST_DIR%\Dockerfile.test"
echo # Copy entrypoint script>> "%TEST_DIR%\Dockerfile.test"
echo COPY test-entrypoint.sh /app/>> "%TEST_DIR%\Dockerfile.test"
echo RUN chmod +x /app/test-entrypoint.sh>> "%TEST_DIR%\Dockerfile.test"
echo.>> "%TEST_DIR%\Dockerfile.test"
echo # Switch to non-root user>> "%TEST_DIR%\Dockerfile.test"
echo USER appuser>> "%TEST_DIR%\Dockerfile.test"
echo.>> "%TEST_DIR%\Dockerfile.test"
echo WORKDIR /app>> "%TEST_DIR%\Dockerfile.test"
echo.>> "%TEST_DIR%\Dockerfile.test"
echo ENTRYPOINT ["/app/test-entrypoint.sh"]>> "%TEST_DIR%\Dockerfile.test"
echo CMD ["sh", "-c", "echo 'Test completed successfully'"]>> "%TEST_DIR%\Dockerfile.test"

REM Test 1: Run with Docker Compose (swarm mode)
call :log_message "%BLUE%Test 1: Running with Docker Compose (swarm mode)...%NC%"
docker info 2>nul | findstr "Swarm: active" >nul
if %ERRORLEVEL% equ 0 (
    cd /d "%TEST_DIR%"
    docker-compose -f docker-compose.test.yml up --abort-on-container-exit
    if %ERRORLEVEL% equ 0 (
        call :log_message "%GREEN%Test 1 passed: Docker Compose with swarm secrets works correctly%NC%"
    ) else (
        call :log_message "%RED%Test 1 failed: Docker Compose with swarm secrets failed%NC%"
    )
) else (
    call :log_message "%YELLOW%Skipping Test 1: Docker swarm is not active%NC%"
    call :log_message "To enable swarm mode, run: docker swarm init"
)

REM Test 2: Run with Docker Compose (file-based secrets for development)
call :log_message "%BLUE%Test 2: Running with Docker Compose (file-based secrets)...%NC%"
cd /d "%TEST_DIR%"
docker-compose -f docker-compose-dev.yml up --abort-on-container-exit
if %ERRORLEVEL% equ 0 (
    call :log_message "%GREEN%Test 2 passed: Docker Compose with file-based secrets works correctly%NC%"
) else (
    call :log_message "%RED%Test 2 failed: Docker Compose with file-based secrets failed%NC%"
)

REM Test 3: Build and run with custom entrypoint script
call :log_message "%BLUE%Test 3: Building and running with custom entrypoint script...%NC%"
cd /d "%TEST_DIR%"
docker build -t cryoprotect-secret-test -f Dockerfile.test .
if %ERRORLEVEL% equ 0 (
    call :log_message "Docker image built successfully"
    
    REM Create a temporary directory for mounting secrets
    set SECRETS_MOUNT=%TEMP%\cryoprotect-secrets-%RANDOM%
    mkdir "%SECRETS_MOUNT%"
    copy "%TEST_DIR%\test_secret_1" "%SECRETS_MOUNT%\TEST_SECRET_1" >nul
    copy "%TEST_DIR%\test_secret_2" "%SECRETS_MOUNT%\TEST_SECRET_2" >nul
    copy "%TEST_DIR%\test_secret_rotation" "%SECRETS_MOUNT%\SECRET_ROTATION_TIMESTAMP" >nul
    
    REM Run the container with mounted secrets
    docker run --rm -v "%SECRETS_MOUNT%:/run/secrets:ro" cryoprotect-secret-test
    if %ERRORLEVEL% equ 0 (
        call :log_message "%GREEN%Test 3 passed: Custom entrypoint script loads secrets correctly%NC%"
    ) else (
        call :log_message "%RED%Test 3 failed: Custom entrypoint script failed to load secrets%NC%"
    )
    
    REM Clean up
    rmdir /s /q "%SECRETS_MOUNT%"
) else (
    call :log_message "%RED%Test 3 failed: Could not build Docker image%NC%"
)

REM Summary
echo.
call :log_message "%BLUE%Secret Management Test Summary:%NC%"
call :log_message "1. Docker Compose with swarm secrets: %YELLOW%See above results%NC%"
call :log_message "2. Docker Compose with file-based secrets: %YELLOW%See above results%NC%"
call :log_message "3. Custom entrypoint script: %YELLOW%See above results%NC%"
echo.
call :log_message "%GREEN%Tests completed. Check the results above to ensure all tests passed.%NC%"
call :log_message "%BLUE%If any tests failed, check the error messages and fix the issues.%NC%"
echo.
call :log_message "For more information on secret management, see docs/secret_management.md"

REM Clean up
call :log_message "Cleaning up temporary files..."
docker secret rm test_secret_1 test_secret_2 test_secret_rotation 2>nul
rmdir /s /q "%TEST_DIR%"

goto :eof

:log_message
echo [%date% %time%] %~1
goto :eof