@echo off
REM CryoProtect v2 Secret Management Script for Windows
REM This script helps create, update, and manage Docker secrets for CryoProtect v2

setlocal enabledelayedexpansion

REM Default values
set SECRET_PREFIX=cryoprotect
set SECRET_DIR=.\.secrets
set MODE=help
set SECRET_NAME=
set SECRET_VALUE=
set ENVIRONMENT=
set DEPLOYMENT_COLOR=
set LIST_SECRETS=false
set ROTATE_SECRETS=false
set EXTERNAL_PROVIDER=
set DOCKER_COMPOSE_FILE=docker-compose.yml

REM ANSI color codes for Windows 10+
set RED=[91m
set GREEN=[92m
set YELLOW=[93m
set BLUE=[94m
set NC=[0m

REM Function to display help
:show_help
echo %BLUE%CryoProtect v2 Secret Management Script for Windows%NC%
echo.
echo Usage: %0 [options]
echo.
echo Options:
echo   -m, --mode MODE           Operation mode: create, update, delete, dev, swarm, k8s
echo   -n, --name NAME           Secret name (e.g., SUPABASE_URL)
echo   -v, --value VALUE         Secret value (or - to read from stdin)
echo   -e, --environment ENV     Environment prefix: staging, production
echo   -c, --color COLOR         Deployment color: blue, green
echo   -p, --prefix PREFIX       Secret name prefix (default: cryoprotect)
echo   -d, --dir DIRECTORY       Directory for dev secrets (default: .\.secrets)
echo   -f, --file FILE           Docker compose file (default: docker-compose.yml)
echo   -l, --list                List all secrets
echo   -r, --rotate              Rotate secrets (update timestamp)
echo   -x, --external PROVIDER   Configure external provider: aws, vault, azure, gcp
echo   -h, --help                Show this help message
echo.
echo Examples:
echo   # Create a Docker Swarm secret
echo   %0 --mode swarm --name SUPABASE_URL --value "https://example.supabase.co"
echo.
echo   # Create a secret for staging environment
echo   %0 --mode swarm --name SUPABASE_URL --value "https://staging.supabase.co" --environment staging
echo.
echo   # Create a secret for blue deployment
echo   %0 --mode swarm --name SUPABASE_URL --value "https://blue.supabase.co" --color blue
echo.
echo   # Set up development secrets
echo   %0 --mode dev --name SUPABASE_URL --value "https://dev.supabase.co"
echo.
echo   # List all secrets
echo   %0 --list
echo.
echo   # Rotate secrets (update timestamp)
echo   %0 --rotate
echo.
echo   # Configure AWS Secrets Manager
echo   %0 --external aws --name AWS_ACCESS_KEY_ID --value "your-access-key"
echo.
goto :eof

REM Function to validate input
:validate_input
if "%MODE%" neq "help" if "%MODE%" neq "list" if "%MODE%" neq "rotate" if "%SECRET_NAME%"=="" (
    echo %RED%Error: Secret name is required%NC%
    exit /b 1
)

if "%MODE%"=="create" goto :check_value
if "%MODE%"=="update" goto :check_value
if "%MODE%"=="swarm" goto :check_value
if "%MODE%"=="dev" goto :check_value
if "%MODE%"=="k8s" goto :check_value
goto :eof

:check_value
if "%SECRET_VALUE%"=="" if "%SECRET_VALUE%" neq "-" (
    echo %RED%Error: Secret value is required%NC%
    exit /b 1
)
goto :eof

REM Function to get full secret name
:get_full_secret_name
set base_name=%~1
set full_name=%SECRET_PREFIX%_%base_name%

if not "%ENVIRONMENT%"=="" (
    set full_name=%SECRET_PREFIX%_%ENVIRONMENT%_%base_name%
)

if not "%DEPLOYMENT_COLOR%"=="" (
    set full_name=%SECRET_PREFIX%_%DEPLOYMENT_COLOR%_%base_name%
)

set result=%full_name%
goto :eof

REM Function to create a Docker Swarm secret
:create_swarm_secret
set name=%~1
set value=%~2
call :get_full_secret_name "%name%"
set full_name=!result!

REM Check if secret already exists
docker secret inspect "%full_name%" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    echo %YELLOW%Secret %full_name% already exists. Use update mode to change it.%NC%
    exit /b 1
)

REM Create the secret
if "%value%"=="-" (
    echo %BLUE%Enter secret value for %name% (input will be visible):
    set /p value=
    echo !value!| docker secret create "%full_name%" -
) else (
    echo %value%| docker secret create "%full_name%" -
)

echo %GREEN%Created Docker Swarm secret: %full_name%%NC%
goto :eof

REM Function to update a Docker Swarm secret
:update_swarm_secret
set name=%~1
set value=%~2
call :get_full_secret_name "%name%"
set full_name=!result!

REM Check if secret exists
docker secret inspect "%full_name%" >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo %YELLOW%Secret %full_name% does not exist. Creating it instead.%NC%
    call :create_swarm_secret "%name%" "%value%"
    goto :eof
)

REM Docker Swarm secrets cannot be updated directly, so we need to delete and recreate
docker secret rm "%full_name%"

REM Create the secret
if "%value%"=="-" (
    echo %BLUE%Enter new secret value for %name% (input will be visible):
    set /p value=
    echo !value!| docker secret create "%full_name%" -
) else (
    echo %value%| docker secret create "%full_name%" -
)

echo %GREEN%Updated Docker Swarm secret: %full_name%%NC%
goto :eof

REM Function to delete a Docker Swarm secret
:delete_swarm_secret
set name=%~1
call :get_full_secret_name "%name%"
set full_name=!result!

REM Check if secret exists
docker secret inspect "%full_name%" >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo %YELLOW%Secret %full_name% does not exist.%NC%
    exit /b 1
)

REM Delete the secret
docker secret rm "%full_name%"

echo %GREEN%Deleted Docker Swarm secret: %full_name%%NC%
goto :eof

REM Function to create a development secret
:create_dev_secret
set name=%~1
set value=%~2

REM Create the secrets directory if it doesn't exist
if not exist "%SECRET_DIR%" mkdir "%SECRET_DIR%"

REM Add to .gitignore if it's not already there
if exist ".gitignore" (
    findstr /C:"%SECRET_DIR%" .gitignore >nul
    if %ERRORLEVEL% neq 0 (
        echo %SECRET_DIR%>> .gitignore
        echo %BLUE%Added %SECRET_DIR% to .gitignore%NC%
    )
) else (
    echo %SECRET_DIR%> .gitignore
    echo %BLUE%Created .gitignore with %SECRET_DIR%%NC%
)

REM Create the secret file
set secret_file=%SECRET_DIR%\%name%

if exist "%secret_file%" (
    echo %YELLOW%Secret file %secret_file% already exists. Overwriting.%NC%
)

if "%value%"=="-" (
    echo %BLUE%Enter secret value for %name% (input will be visible):
    set /p value=
    echo !value!> "%secret_file%"
) else (
    echo %value%> "%secret_file%"
)

REM Set proper permissions (Windows equivalent - restrict to current user)
icacls "%secret_file%" /inheritance:r /grant:r "%USERNAME%:(F)" >nul

echo %GREEN%Created development secret: %secret_file%%NC%

REM Update docker-compose.yml if it exists
if exist "%DOCKER_COMPOSE_FILE%" (
    echo %BLUE%To use this secret in development, add the following to your %DOCKER_COMPOSE_FILE%:%NC%
    echo.
    echo secrets:
    echo   %name%:
    echo     file: %secret_file%
    echo.
    echo services:
    echo   cryoprotect-dev:
    echo     secrets:
    echo       - source: %name%
    echo         target: %name%
    echo         mode: 0400
)
goto :eof

REM Function to create a Kubernetes secret
:create_k8s_secret
set name=%~1
set value=%~2
if defined K8S_NAMESPACE (
    set namespace=%K8S_NAMESPACE%
) else (
    set namespace=default
)

REM Check if kubectl is available
where kubectl >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo %RED%Error: kubectl is not installed or not in PATH%NC%
    exit /b 1
)

REM Check if secret already exists
kubectl get secret "%SECRET_PREFIX%" -n "%namespace%" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    REM Secret exists, update it
    if "%value%"=="-" (
        echo %BLUE%Enter secret value for %name% (input will be visible):
        set /p value=
    )
    
    REM Base64 encode the value (Windows PowerShell method)
    for /f "usebackq delims=" %%a in (`powershell -Command "[Convert]::ToBase64String([System.Text.Encoding]::UTF8.GetBytes('%value%'))"`) do set encoded_value=%%a
    
    kubectl patch secret "%SECRET_PREFIX%" -n "%namespace%" -p "{\"data\":{\"%name%\":\"%encoded_value%\"}}"
    echo %GREEN%Updated Kubernetes secret: %name% in %SECRET_PREFIX%%NC%
) else (
    REM Secret doesn't exist, create it
    if "%value%"=="-" (
        echo %BLUE%Enter secret value for %name% (input will be visible):
        set /p value=
    )
    
    kubectl create secret generic "%SECRET_PREFIX%" -n "%namespace%" --from-literal="%name%=%value%"
    echo %GREEN%Created Kubernetes secret: %name% in %SECRET_PREFIX%%NC%
)

echo %BLUE%To use this secret in Kubernetes, set the following environment variables:%NC%
echo K8S_SECRET_NAMESPACE=%namespace%
echo K8S_SECRET_NAME=%SECRET_PREFIX%
goto :eof

REM Function to list secrets
:list_secrets
echo %BLUE%Listing secrets:%NC%

REM Check for Docker Swarm secrets
where docker >nul 2>&1
if %ERRORLEVEL% equ 0 (
    echo %BLUE%Docker Swarm secrets:%NC%
    docker secret ls | findstr "%SECRET_PREFIX%" || echo No Docker Swarm secrets found
    echo.
)

REM Check for development secrets
if exist "%SECRET_DIR%" (
    echo %BLUE%Development secrets:%NC%
    dir /b "%SECRET_DIR%" || echo No development secrets found
    echo.
)

REM Check for Kubernetes secrets
where kubectl >nul 2>&1
if %ERRORLEVEL% equ 0 (
    echo %BLUE%Kubernetes secrets:%NC%
    kubectl get secret "%SECRET_PREFIX%" --all-namespaces -o wide 2>nul || echo No Kubernetes secrets found
    echo.
)
goto :eof

REM Function to rotate secrets
:rotate_secrets
for /f "tokens=*" %%a in ('powershell -Command "Get-Date -UFormat %%s"') do set timestamp=%%a

echo %BLUE%Rotating secrets (updating timestamp to %timestamp%)%NC%

REM Update Docker Swarm secret
where docker >nul 2>&1
if %ERRORLEVEL% equ 0 (
    set rotation_secret=%SECRET_PREFIX%_SECRET_ROTATION_TIMESTAMP
    
    REM Check if secret already exists
    docker secret inspect "%rotation_secret%" >nul 2>&1
    if %ERRORLEVEL% equ 0 (
        docker secret rm "%rotation_secret%"
    )
    
    REM Create the secret
    echo %timestamp%| docker secret create "%rotation_secret%" -
    echo %GREEN%Updated Docker Swarm secret rotation timestamp: %rotation_secret%%NC%
)

REM Update development secret
if exist "%SECRET_DIR%" (
    set rotation_file=%SECRET_DIR%\SECRET_ROTATION_TIMESTAMP
    echo %timestamp%> "%rotation_file%"
    icacls "%rotation_file%" /inheritance:r /grant:r "%USERNAME%:(F)" >nul
    echo %GREEN%Updated development secret rotation timestamp: %rotation_file%%NC%
)

REM Update Kubernetes secret
where kubectl >nul 2>&1
if %ERRORLEVEL% equ 0 (
    if defined K8S_NAMESPACE (
        set namespace=%K8S_NAMESPACE%
    ) else (
        set namespace=default
    )
    
    kubectl get secret "%SECRET_PREFIX%" -n "%namespace%" >nul 2>&1
    if %ERRORLEVEL% equ 0 (
        REM Base64 encode the value (Windows PowerShell method)
        for /f "usebackq delims=" %%a in (`powershell -Command "[Convert]::ToBase64String([System.Text.Encoding]::UTF8.GetBytes('%timestamp%'))"`) do set encoded_timestamp=%%a
        
        kubectl patch secret "%SECRET_PREFIX%" -n "%namespace%" -p "{\"data\":{\"SECRET_ROTATION_TIMESTAMP\":\"%encoded_timestamp%\"}}"
        echo %GREEN%Updated Kubernetes secret rotation timestamp in %SECRET_PREFIX%%NC%
    ) else (
        kubectl create secret generic "%SECRET_PREFIX%" -n "%namespace%" --from-literal="SECRET_ROTATION_TIMESTAMP=%timestamp%"
        echo %GREEN%Created Kubernetes secret rotation timestamp in %SECRET_PREFIX%%NC%
    )
)
goto :eof

REM Function to configure external provider
:configure_external_provider
set provider=%~1
set name=%~2
set value=%~3

if "%provider%"=="aws" (
    echo %BLUE%Configuring AWS Secrets Manager integration%NC%
    if "%name%"=="" if "%value%"=="" (
        echo %YELLOW%For AWS integration, you need to create the following secrets:%NC%
        echo   - AWS_ACCESS_KEY_ID
        echo   - AWS_SECRET_ACCESS_KEY
        echo   - AWS_REGION (optional)
        echo.
        echo %YELLOW%And set the following environment variable:%NC%
        echo   - AWS_SECRET_PREFIX (e.g., /cryoprotect/)
        goto :eof
    )
) else if "%provider%"=="vault" (
    echo %BLUE%Configuring HashiCorp Vault integration%NC%
    if "%name%"=="" if "%value%"=="" (
        echo %YELLOW%For Vault integration, you need to create the following secrets:%NC%
        echo   - VAULT_TOKEN
        echo   - VAULT_ADDR (optional)
        echo.
        echo %YELLOW%And set the following environment variable:%NC%
        echo   - VAULT_SECRET_PATH (e.g., secret/cryoprotect)
        goto :eof
    )
) else if "%provider%"=="azure" (
    echo %BLUE%Configuring Azure Key Vault integration%NC%
    if "%name%"=="" if "%value%"=="" (
        echo %YELLOW%For Azure integration, you need to create the following secrets:%NC%
        echo   - AZURE_CLIENT_ID
        echo   - AZURE_CLIENT_SECRET
        echo   - AZURE_TENANT_ID
        echo.
        echo %YELLOW%And set the following environment variable:%NC%
        echo   - AZURE_KEYVAULT_NAME
        goto :eof
    )
) else if "%provider%"=="gcp" (
    echo %BLUE%Configuring Google Cloud Secret Manager integration%NC%
    if "%name%"=="" if "%value%"=="" (
        echo %YELLOW%For GCP integration, you need to create the following secret:%NC%
        echo   - GOOGLE_APPLICATION_CREDENTIALS (service account key file content)
        echo.
        echo %YELLOW%And set the following environment variables:%NC%
        echo   - GCP_PROJECT_ID
        echo   - GCP_SECRET_PREFIX
        goto :eof
    )
) else (
    echo %RED%Error: Unknown external provider: %provider%%NC%
    echo %YELLOW%Supported providers: aws, vault, azure, gcp%NC%
    exit /b 1
)

REM Create the secret based on the current mode
if "%MODE%"=="swarm" (
    call :create_swarm_secret "%name%" "%value%"
) else if "%MODE%"=="dev" (
    call :create_dev_secret "%name%" "%value%"
) else if "%MODE%"=="k8s" (
    call :create_k8s_secret "%name%" "%value%"
) else (
    echo %RED%Error: Please specify a mode (swarm, dev, k8s) when configuring external provider%NC%
    exit /b 1
)

echo %BLUE%To use %provider%, set EXTERNAL_SECRET_PROVIDER=%provider% in your environment%NC%
goto :eof

REM Parse command line arguments
:parse_args
if "%~1"=="" goto :execute

if "%~1"=="-m" (
    set MODE=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--mode" (
    set MODE=%~2
    shift
    shift
    goto :parse_args
)

if "%~1"=="-n" (
    set SECRET_NAME=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--name" (
    set SECRET_NAME=%~2
    shift
    shift
    goto :parse_args
)

if "%~1"=="-v" (
    set SECRET_VALUE=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--value" (
    set SECRET_VALUE=%~2
    shift
    shift
    goto :parse_args
)

if "%~1"=="-e" (
    set ENVIRONMENT=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--environment" (
    set ENVIRONMENT=%~2
    shift
    shift
    goto :parse_args
)

if "%~1"=="-c" (
    set DEPLOYMENT_COLOR=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--color" (
    set DEPLOYMENT_COLOR=%~2
    shift
    shift
    goto :parse_args
)

if "%~1"=="-p" (
    set SECRET_PREFIX=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--prefix" (
    set SECRET_PREFIX=%~2
    shift
    shift
    goto :parse_args
)

if "%~1"=="-d" (
    set SECRET_DIR=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--dir" (
    set SECRET_DIR=%~2
    shift
    shift
    goto :parse_args
)

if "%~1"=="-f" (
    set DOCKER_COMPOSE_FILE=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--file" (
    set DOCKER_COMPOSE_FILE=%~2
    shift
    shift
    goto :parse_args
)

if "%~1"=="-l" (
    set LIST_SECRETS=true
    shift
    goto :parse_args
)
if "%~1"=="--list" (
    set LIST_SECRETS=true
    shift
    goto :parse_args
)

if "%~1"=="-r" (
    set ROTATE_SECRETS=true
    shift
    goto :parse_args
)
if "%~1"=="--rotate" (
    set ROTATE_SECRETS=true
    shift
    goto :parse_args
)

if "%~1"=="-x" (
    set EXTERNAL_PROVIDER=%~2
    shift
    shift
    goto :parse_args
)
if "%~1"=="--external" (
    set EXTERNAL_PROVIDER=%~2
    shift
    shift
    goto :parse_args
)

if "%~1"=="-h" (
    call :show_help
    exit /b 0
)
if "%~1"=="--help" (
    call :show_help
    exit /b 0
)

echo %RED%Error: Unknown option: %~1%NC%
call :show_help
exit /b 1

REM Execute the requested operation
:execute

REM Handle list mode
if "%LIST_SECRETS%"=="true" (
    call :list_secrets
    exit /b 0
)

REM Handle rotate mode
if "%ROTATE_SECRETS%"=="true" (
    call :rotate_secrets
    exit /b 0
)

REM Handle external provider configuration
if not "%EXTERNAL_PROVIDER%"=="" (
    call :configure_external_provider "%EXTERNAL_PROVIDER%" "%SECRET_NAME%" "%SECRET_VALUE%"
    exit /b 0
)

REM Handle help mode
if "%MODE%"=="help" (
    call :show_help
    exit /b 0
)

REM Validate input
call :validate_input
if %ERRORLEVEL% neq 0 exit /b %ERRORLEVEL%

REM Execute the requested operation
if "%MODE%"=="create" (
    call :create_swarm_secret "%SECRET_NAME%" "%SECRET_VALUE%"
) else if "%MODE%"=="swarm" (
    call :create_swarm_secret "%SECRET_NAME%" "%SECRET_VALUE%"
) else if "%MODE%"=="update" (
    call :update_swarm_secret "%SECRET_NAME%" "%SECRET_VALUE%"
) else if "%MODE%"=="delete" (
    call :delete_swarm_secret "%SECRET_NAME%"
) else if "%MODE%"=="dev" (
    call :create_dev_secret "%SECRET_NAME%" "%SECRET_VALUE%"
) else if "%MODE%"=="k8s" (
    call :create_k8s_secret "%SECRET_NAME%" "%SECRET_VALUE%"
) else (
    echo %RED%Error: Unknown mode: %MODE%%NC%
    call :show_help
    exit /b 1
)

exit /b 0

REM Start execution
call :parse_args %*