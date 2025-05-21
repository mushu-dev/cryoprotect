#!/bin/bash
# Ultra-minimal CryoProtect runner with no volume mounts

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Minimal Runner${NC}"
echo "Running a simplified version with no volume mounts"
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

# Create a temporary Dockerfile with all files copied in
echo -e "${BLUE}Creating minimal Dockerfile...${NC}"

cat > Dockerfile.minimal << EOF
FROM python:3.10-slim

WORKDIR /app

# Install minimal dependencies
RUN pip install --no-cache-dir flask==3.0.2 flask-restful==0.3.10 flask-cors==4.0.0 \
    requests==2.31.0 python-dotenv==1.0.1 apispec==6.3.0 flask-apispec==0.11.4 \
    marshmallow==3.21.0

# Copy application code
COPY . /app/

# Expose port
EXPOSE 5000

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \\
    PYTHONUNBUFFERED=1 \\
    FLASK_APP=app.py \\
    FLASK_DEBUG=1

CMD ["python", "-m", "flask", "run", "--host=0.0.0.0", "--port=5000"]
EOF

# Build the image
echo -e "${BLUE}Building minimal image...${NC}"
podman build -t cryoprotect-minimal:latest -f Dockerfile.minimal .
build_status=$?

if [ $build_status -ne 0 ]; then
  echo -e "${RED}Container build failed.${NC}"
  exit 1
fi

echo -e "${GREEN}Container built successfully!${NC}"

# Run with absolutely no volume mounts
echo -e "${BLUE}Starting minimal container...${NC}"
podman run --rm -it \
  --name cryoprotect-minimal \
  --network=cryoprotect-net \
  --dns=1.1.1.1 --dns=8.8.8.8 \
  -p 5000:5000 \
  -e FLASK_APP=app.py \
  -e FLASK_ENV=development \
  -e FLASK_DEBUG=1 \
  -e SUPABASE_URL="${SUPABASE_URL}" \
  -e SUPABASE_KEY="${SUPABASE_KEY}" \
  -e SECRET_KEY="${SECRET_KEY}" \
  -e LOG_LEVEL="${LOG_LEVEL:-DEBUG}" \
  -e LOG_TO_FILE="${LOG_TO_FILE:-1}" \
  -e LOG_TO_CONSOLE="${LOG_TO_CONSOLE:-1}" \
  cryoprotect-minimal:latest

# Clean up the temporary Dockerfile
rm -f Dockerfile.minimal