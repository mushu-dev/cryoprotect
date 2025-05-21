#!/bin/bash
# Simple Podman runner script for CryoProtect without Miniconda
# This script builds and runs the application with Podman directly

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Simple Podman Runner${NC}"
echo "This script builds and runs CryoProtect directly with Podman"
echo "==========================================================================="

# Load environment variables
if [ -f .env ]; then
  source .env
  echo -e "${GREEN}Loaded environment variables from .env${NC}"
else
  echo -e "${YELLOW}Warning: No .env file found. Using default environment variables.${NC}"
  export FLASK_APP=app.py
  export FLASK_ENV=development
  export SECRET_KEY=dev-secret-key-please-change-in-production
  export LOG_LEVEL=DEBUG
fi

# Check if cryoprotect-net network exists, create if not
if ! podman network ls | grep -q "cryoprotect-net"; then
  echo -e "${YELLOW}Network 'cryoprotect-net' not found. Creating...${NC}"
  podman network create --ipv6=false cryoprotect-net
  echo -e "${GREEN}Created network 'cryoprotect-net'${NC}"
fi

# Create necessary directories
mkdir -p logs
mkdir -p data
mkdir -p backup/data
mkdir -p cache

# Build a simplified Docker image
echo -e "${BLUE}Building simplified CryoProtect container...${NC}"

# Create a temporary Dockerfile without Miniconda
cat > Dockerfile.simple << EOF
FROM python:3.10-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && \\
    apt-get install -y --no-install-recommends \\
        curl \\
        build-essential \\
        libpq-dev \\
        tini && \\
    apt-get clean && \\
    rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements_updated.txt /app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --no-cache-dir psycopg2-binary flask flask-restful flask-cors psycopg2-binary rdkit

# Copy application code
COPY . /app/

# Create non-root user
RUN useradd -m -s /bin/bash -u 1000 appuser && \\
    mkdir -p /app /app/logs /app/cache && \\
    chown -R appuser:appuser /app && \\
    chmod 750 /app

# Switch to non-root user
USER appuser

# Expose port
EXPOSE 5000

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \\
    PYTHONUNBUFFERED=1 \\
    FLASK_APP=app.py \\
    FLASK_DEBUG=1

# Use tini as init
ENTRYPOINT ["/usr/bin/tini", "--", "python", "-m", "flask", "run", "--host=0.0.0.0", "--port=5000"]
EOF

# Build the image
podman build -t cryoprotect-simple:latest -f Dockerfile.simple .
build_status=$?

if [ $build_status -ne 0 ]; then
  echo -e "${RED}Container build failed. Trying with minimal requirements...${NC}"
  
  # Try with minimal requirements
  cat > minimal_requirements.txt << EOF
flask==3.0.2
flask-restful==0.3.10
flask-cors==4.0.0
psycopg2-binary==2.9.9
python-dotenv==1.0.1
requests==2.31.0
EOF

  # Update the Dockerfile to use minimal requirements
  sed -i 's|COPY requirements_updated.txt /app/requirements.txt|COPY minimal_requirements.txt /app/requirements.txt|' Dockerfile.simple
  
  # Try building again
  podman build -t cryoprotect-simple:latest -f Dockerfile.simple .
  build_status=$?
  
  if [ $build_status -ne 0 ]; then
    echo -e "${RED}Container build failed again. Please check the build errors above.${NC}"
    exit 1
  fi
fi

echo -e "${GREEN}Container built successfully!${NC}"

# Run the container
echo -e "${BLUE}Starting CryoProtect container...${NC}"
podman run --rm -it \
  --name cryoprotect-simple \
  --network=cryoprotect-net \
  --dns=1.1.1.1 --dns=8.8.8.8 \
  -p 5000:5000 \
  -v .:/app:z \
  -v ./logs:/app/logs:z \
  -v ./data:/app/data:z \
  -v ./cache:/app/cache:z \
  -e FLASK_APP=app.py \
  -e FLASK_ENV=development \
  -e FLASK_DEBUG=1 \
  -e SUPABASE_URL="${SUPABASE_URL}" \
  -e SUPABASE_KEY="${SUPABASE_KEY}" \
  -e SECRET_KEY="${SECRET_KEY}" \
  -e LOG_LEVEL="${LOG_LEVEL:-DEBUG}" \
  -e LOG_TO_FILE="${LOG_TO_FILE:-1}" \
  -e LOG_TO_CONSOLE="${LOG_TO_CONSOLE:-1}" \
  --security-opt label=disable \
  cryoprotect-simple:latest

# Clean up the temporary Dockerfile
rm -f Dockerfile.simple
rm -f minimal_requirements.txt