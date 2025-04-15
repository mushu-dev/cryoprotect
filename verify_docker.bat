@echo off
REM CryoProtect Analyzer - Docker Verification Script for Windows
REM This script builds and runs the Docker container to verify that the fixes work.

echo Starting Docker verification...

REM Build the Docker image
echo Building Docker image...
docker-compose build

REM Check if the build was successful
if %ERRORLEVEL% neq 0 (
    echo Error: Docker build failed.
    exit /b 1
)

echo Docker image built successfully.

REM Run the Docker container in detached mode
echo Starting Docker container...
docker-compose up -d

REM Check if the container is running
if %ERRORLEVEL% neq 0 (
    echo Error: Failed to start Docker container.
    exit /b 1
)

echo Docker container started successfully.

REM Wait for the container to initialize
echo Waiting for container to initialize (10 seconds)...
timeout /t 10 /nobreak > nul

REM Get the container ID
for /f "tokens=*" %%i in ('docker-compose ps -q cryoprotect') do set CONTAINER_ID=%%i

if "%CONTAINER_ID%"=="" (
    echo Error: Could not find container ID.
    docker-compose down
    exit /b 1
)

echo Container ID: %CONTAINER_ID%

REM Run the verification script inside the container
echo Running RDKit verification script inside the container...
docker exec -it %CONTAINER_ID% conda run -n cryoprotect python verify_rdkit.py

REM Check if the verification was successful
if %ERRORLEVEL% neq 0 (
    echo Error: RDKit verification failed.
    docker-compose down
    exit /b 1
)

echo RDKit verification completed successfully.

REM Check the Flask application
echo Checking Flask application...
docker exec -it %CONTAINER_ID% conda run -n cryoprotect curl -s http://localhost:5000/health

REM Check if the health check was successful
if %ERRORLEVEL% neq 0 (
    echo Error: Flask application health check failed.
    docker-compose down
    exit /b 1
)

echo Flask application is running correctly.

REM Stop the Docker container
echo Stopping Docker container...
docker-compose down

echo Docker verification completed successfully.
echo The Docker container is now working correctly with RDKit integration.