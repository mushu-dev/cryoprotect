#!/bin/bash
# CryoProtect v2 Secret Management Test Script
# This script tests the secret management system to ensure it's working correctly

set -e

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to log messages with timestamp
log_message() {
  echo -e "[$(date -Iseconds)] $1"
}

log_message "${BLUE}CryoProtect v2 Secret Management Test Script${NC}"
log_message "This script will test the secret management system to ensure it's working correctly."
echo ""

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
  log_message "${RED}Error: Docker is not installed or not in PATH${NC}"
  exit 1
fi

# Check if Docker Compose is installed
if ! command -v docker-compose &> /dev/null; then
  log_message "${RED}Error: Docker Compose is not installed or not in PATH${NC}"
  exit 1
fi

# Create a temporary directory for test secrets
TEST_DIR=$(mktemp -d)
log_message "Created temporary directory for test secrets: ${TEST_DIR}"

# Clean up on exit
trap 'log_message "Cleaning up temporary files..."; rm -rf ${TEST_DIR}; docker secret rm test_secret_1 test_secret_2 test_secret_rotation 2>/dev/null || true' EXIT

# Create test secrets
log_message "Creating test secrets..."
echo "test-value-1" > "${TEST_DIR}/test_secret_1"
echo "test-value-2" > "${TEST_DIR}/test_secret_2"
echo "$(date +%s)" > "${TEST_DIR}/test_secret_rotation"

# Create Docker secrets
log_message "Creating Docker secrets..."
cat "${TEST_DIR}/test_secret_1" | docker secret create test_secret_1 - || log_message "${YELLOW}Warning: Could not create Docker secret. Are you running in swarm mode?${NC}"
cat "${TEST_DIR}/test_secret_2" | docker secret create test_secret_2 - || log_message "${YELLOW}Warning: Could not create Docker secret. Are you running in swarm mode?${NC}"
cat "${TEST_DIR}/test_secret_rotation" | docker secret create test_secret_rotation - || log_message "${YELLOW}Warning: Could not create Docker secret. Are you running in swarm mode?${NC}"

# Create a test docker-compose.yml file
cat > "${TEST_DIR}/docker-compose.test.yml" << EOF
version: '3.8'

services:
  secret-test:
    image: alpine:latest
    command: sh -c 'cat /run/secrets/TEST_SECRET_1; echo; cat /run/secrets/TEST_SECRET_2; echo; cat /run/secrets/SECRET_ROTATION_TIMESTAMP; echo; sleep 5'
    secrets:
      - source: test_secret_1
        target: TEST_SECRET_1
      - source: test_secret_2
        target: TEST_SECRET_2
      - source: test_secret_rotation
        target: SECRET_ROTATION_TIMESTAMP

secrets:
  test_secret_1:
    external: true
  test_secret_2:
    external: true
  test_secret_rotation:
    external: true
EOF

# Create a test docker-compose-dev.yml file for file-based secrets
mkdir -p "${TEST_DIR}/.secrets"
cp "${TEST_DIR}/test_secret_1" "${TEST_DIR}/.secrets/TEST_SECRET_1"
cp "${TEST_DIR}/test_secret_2" "${TEST_DIR}/.secrets/TEST_SECRET_2"
cp "${TEST_DIR}/test_secret_rotation" "${TEST_DIR}/.secrets/SECRET_ROTATION_TIMESTAMP"

cat > "${TEST_DIR}/docker-compose-dev.yml" << EOF
version: '3.8'

services:
  secret-test-dev:
    image: alpine:latest
    command: sh -c 'cat /run/secrets/TEST_SECRET_1; echo; cat /run/secrets/TEST_SECRET_2; echo; cat /run/secrets/SECRET_ROTATION_TIMESTAMP; echo; sleep 5'
    secrets:
      - source: test_secret_1
        target: TEST_SECRET_1
      - source: test_secret_2
        target: TEST_SECRET_2
      - source: test_secret_rotation
        target: SECRET_ROTATION_TIMESTAMP

secrets:
  test_secret_1:
    file: ./.secrets/TEST_SECRET_1
  test_secret_2:
    file: ./.secrets/TEST_SECRET_2
  test_secret_rotation:
    file: ./.secrets/SECRET_ROTATION_TIMESTAMP
EOF

# Create a test entrypoint script
cat > "${TEST_DIR}/test-entrypoint.sh" << 'EOF'
#!/bin/sh
set -e

# Function to log messages with timestamp
log_message() {
  echo "[$(date -Iseconds)] $1"
}

# Function to load secrets into environment variables
load_secrets() {
  log_message "Loading secrets..."
  
  # Define required secrets
  REQUIRED_SECRETS=(
    "TEST_SECRET_1"
    "TEST_SECRET_2"
  )
  
  # Load each secret if it exists
  for SECRET in "${REQUIRED_SECRETS[@]}"; do
    SECRET_FILE="/run/secrets/${SECRET}"
    if [ -f "$SECRET_FILE" ]; then
      # Read the secret value and set it as an environment variable
      export "$SECRET"="$(cat "$SECRET_FILE")"
      log_message "Loaded secret: $SECRET = ${!SECRET}"
    else
      log_message "Error: Required secret $SECRET not found"
      exit 1
    fi
  done
  
  # Check for rotation timestamp
  ROTATION_FILE="/run/secrets/SECRET_ROTATION_TIMESTAMP"
  if [ -f "$ROTATION_FILE" ]; then
    ROTATION_TIMESTAMP=$(cat "$ROTATION_FILE")
    CURRENT_TIMESTAMP=$(date +%s)
    
    log_message "Secret rotation timestamp: $ROTATION_TIMESTAMP"
    
    # Calculate age of secrets in days
    SECRET_AGE_SECONDS=$((CURRENT_TIMESTAMP - ROTATION_TIMESTAMP))
    SECRET_AGE_DAYS=$((SECRET_AGE_SECONDS / 86400))
    
    log_message "Secret age: $SECRET_AGE_DAYS days"
  else
    log_message "Warning: No rotation timestamp found"
  fi
  
  log_message "Secret loading complete"
}

# Load secrets
load_secrets

# Execute the command passed to the container
log_message "Executing: $*"
exec "$@"
EOF

chmod +x "${TEST_DIR}/test-entrypoint.sh"

# Create a test Dockerfile
cat > "${TEST_DIR}/Dockerfile.test" << EOF
FROM alpine:latest

# Create non-root user
RUN adduser -D -u 1000 appuser && \\
    mkdir -p /app /run/secrets && \\
    chown -R appuser:appuser /app && \\
    chown -R appuser:appuser /run/secrets && \\
    chmod 750 /app && \\
    chmod 700 /run/secrets

# Copy entrypoint script
COPY test-entrypoint.sh /app/
RUN chmod +x /app/test-entrypoint.sh

# Switch to non-root user
USER appuser

WORKDIR /app

ENTRYPOINT ["/app/test-entrypoint.sh"]
CMD ["sh", "-c", "echo 'Test completed successfully'"]
EOF

# Test 1: Run with Docker Compose (swarm mode)
log_message "${BLUE}Test 1: Running with Docker Compose (swarm mode)...${NC}"
if docker info 2>/dev/null | grep -q "Swarm: active"; then
  cd "${TEST_DIR}"
  if docker-compose -f docker-compose.test.yml up --abort-on-container-exit; then
    log_message "${GREEN}Test 1 passed: Docker Compose with swarm secrets works correctly${NC}"
  else
    log_message "${RED}Test 1 failed: Docker Compose with swarm secrets failed${NC}"
  fi
else
  log_message "${YELLOW}Skipping Test 1: Docker swarm is not active${NC}"
  log_message "To enable swarm mode, run: docker swarm init"
fi

# Test 2: Run with Docker Compose (file-based secrets for development)
log_message "${BLUE}Test 2: Running with Docker Compose (file-based secrets)...${NC}"
cd "${TEST_DIR}"
if docker-compose -f docker-compose-dev.yml up --abort-on-container-exit; then
  log_message "${GREEN}Test 2 passed: Docker Compose with file-based secrets works correctly${NC}"
else
  log_message "${RED}Test 2 failed: Docker Compose with file-based secrets failed${NC}"
fi

# Test 3: Build and run with custom entrypoint script
log_message "${BLUE}Test 3: Building and running with custom entrypoint script...${NC}"
cd "${TEST_DIR}"
if docker build -t cryoprotect-secret-test -f Dockerfile.test .; then
  log_message "Docker image built successfully"
  
  # Create a temporary directory for mounting secrets
  SECRETS_MOUNT=$(mktemp -d)
  cp "${TEST_DIR}/test_secret_1" "${SECRETS_MOUNT}/TEST_SECRET_1"
  cp "${TEST_DIR}/test_secret_2" "${SECRETS_MOUNT}/TEST_SECRET_2"
  cp "${TEST_DIR}/test_secret_rotation" "${SECRETS_MOUNT}/SECRET_ROTATION_TIMESTAMP"
  
  # Run the container with mounted secrets
  if docker run --rm -v "${SECRETS_MOUNT}:/run/secrets:ro" cryoprotect-secret-test; then
    log_message "${GREEN}Test 3 passed: Custom entrypoint script loads secrets correctly${NC}"
  else
    log_message "${RED}Test 3 failed: Custom entrypoint script failed to load secrets${NC}"
  fi
  
  # Clean up
  rm -rf "${SECRETS_MOUNT}"
else
  log_message "${RED}Test 3 failed: Could not build Docker image${NC}"
fi

# Summary
echo ""
log_message "${BLUE}Secret Management Test Summary:${NC}"
log_message "1. Docker Compose with swarm secrets: ${YELLOW}See above results${NC}"
log_message "2. Docker Compose with file-based secrets: ${YELLOW}See above results${NC}"
log_message "3. Custom entrypoint script: ${YELLOW}See above results${NC}"
echo ""
log_message "${GREEN}Tests completed. Check the results above to ensure all tests passed.${NC}"
log_message "${BLUE}If any tests failed, check the error messages and fix the issues.${NC}"
echo ""
log_message "For more information on secret management, see docs/secret_management.md"