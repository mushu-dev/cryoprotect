#!/bin/bash
# CryoProtect Analyzer - Docker Verification Script
# This script builds and runs the Docker container to verify that the fixes work.

echo "Starting Docker verification..."

# Build the Docker image
echo "Building Docker image..."
docker-compose build

# Check if the build was successful
if [ $? -ne 0 ]; then
    echo "Error: Docker build failed."
    exit 1
fi

echo "Docker image built successfully."

# Run the Docker container in detached mode
echo "Starting Docker container..."
docker-compose up -d

# Check if the container is running
if [ $? -ne 0 ]; then
    echo "Error: Failed to start Docker container."
    exit 1
fi

echo "Docker container started successfully."

# Wait for the container to initialize
echo "Waiting for container to initialize (10 seconds)..."
sleep 10

# Get the container ID
CONTAINER_ID=$(docker-compose ps -q cryoprotect)

if [ -z "$CONTAINER_ID" ]; then
    echo "Error: Could not find container ID."
    docker-compose down
    exit 1
fi

echo "Container ID: $CONTAINER_ID"

# Run the verification script inside the container
echo "Running RDKit verification script inside the container..."
docker exec -it $CONTAINER_ID conda run -n cryoprotect python verify_rdkit.py

# Check if the verification was successful
if [ $? -ne 0 ]; then
    echo "Error: RDKit verification failed."
    docker-compose down
    exit 1
fi

echo "RDKit verification completed successfully."

# Check the Flask application
echo "Checking Flask application..."
docker exec -it $CONTAINER_ID conda run -n cryoprotect curl -s http://localhost:5000/health

# Check if the health check was successful
if [ $? -ne 0 ]; then
    echo "Error: Flask application health check failed."
    docker-compose down
    exit 1
fi

echo "Flask application is running correctly."

# Stop the Docker container
echo "Stopping Docker container..."
docker-compose down

echo "Docker verification completed successfully."
echo "The Docker container is now working correctly with RDKit integration."