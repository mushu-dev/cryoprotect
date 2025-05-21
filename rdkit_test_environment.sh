#!/bin/bash
# RDKit Test Environment Manager for CryoProtect
# This script sets up and manages a dedicated environment for RDKit testing

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

CONTAINER_NAME="CryoProtect-RDKit-Conda"
TEST_DIR="/app/rdkit_tests"

show_help() {
  echo -e "${BLUE}RDKit Test Environment Manager${NC}"
  echo "Usage: $0 [command]"
  echo
  echo "Commands:"
  echo "  setup           Set up the RDKit test environment"
  echo "  run SCRIPT      Run a specific test script in the RDKit environment"
  echo "  shell           Open an interactive shell in the RDKit container"
  echo "  status          Show the status of the RDKit test environment"
  echo "  stop            Stop the RDKit container"
  echo "  clean           Clean up the RDKit test environment"
  echo "  help            Show this help message"
  echo
  echo "Examples:"
  echo "  $0 setup"
  echo "  $0 run verify_rdkit.py"
  echo "  $0 shell"
}

setup_test_environment() {
  echo -e "${BLUE}Setting up RDKit test environment...${NC}"
  
  # Check if the RDKit container exists and is running
  if podman container exists "$CONTAINER_NAME"; then
    if ! podman container inspect "$CONTAINER_NAME" --format '{{.State.Running}}' | grep -q "true"; then
      echo "Starting existing RDKit container..."
      podman start "$CONTAINER_NAME"
    else
      echo "RDKit container is already running"
    fi
  else
    echo "Creating new RDKit container..."
    ./quick_conda_container.sh
  fi
  
  # Set up the test directory in the container
  echo "Creating test directory in container..."
  podman exec -it "$CONTAINER_NAME" bash -c "mkdir -p $TEST_DIR"
  
  # Copy all test scripts to the test directory
  echo "Copying RDKit test scripts to container..."
  for script in *rdkit*.py; do
    if [ -f "$script" ]; then
      podman cp "$script" "$CONTAINER_NAME:$TEST_DIR/"
      echo "  Copied $script"
    fi
  done
  
  echo -e "${GREEN}RDKit test environment is ready${NC}"
  echo "You can run tests with: $0 run [script_name]"
}

run_test_script() {
  if [ -z "$1" ]; then
    echo -e "${RED}Error: No test script specified${NC}"
    echo "Usage: $0 run [script_name]"
    return 1
  fi
  
  SCRIPT="$1"
  
  echo -e "${BLUE}Running $SCRIPT in RDKit environment...${NC}"
  
  # Check if the RDKit container is running
  if ! podman container exists "$CONTAINER_NAME" || \
     ! podman container inspect "$CONTAINER_NAME" --format '{{.State.Running}}' | grep -q "true"; then
    echo "RDKit container is not running. Setting up test environment..."
    setup_test_environment
  fi
  
  # Check if the script exists in the container
  if ! podman exec -it "$CONTAINER_NAME" bash -c "test -f $TEST_DIR/$SCRIPT"; then
    # If the script doesn't exist in the container's test directory, copy it
    if [ -f "$SCRIPT" ]; then
      echo "Copying $SCRIPT to container..."
      podman cp "$SCRIPT" "$CONTAINER_NAME:$TEST_DIR/"
    else
      echo -e "${RED}Error: Script $SCRIPT not found${NC}"
      return 1
    fi
  fi
  
  # Make the script executable in the container
  podman exec -it "$CONTAINER_NAME" bash -c "chmod +x $TEST_DIR/$SCRIPT"
  
  # Run the script in the container with the RDKit environment
  podman exec -it "$CONTAINER_NAME" bash -c "
    cd $TEST_DIR && \
    python3 -c 'import sys; print(\"Python path:\", sys.path)' && \
    python3 -c 'import rdkit; print(\"RDKit version:\", rdkit.__version__)' && \
    python3 $SCRIPT
  "
  
  echo -e "${GREEN}Test execution completed${NC}"
}

open_shell() {
  echo -e "${BLUE}Opening interactive shell in RDKit environment...${NC}"
  
  # Check if the RDKit container is running
  if ! podman container exists "$CONTAINER_NAME" || \
     ! podman container inspect "$CONTAINER_NAME" --format '{{.State.Running}}' | grep -q "true"; then
    echo "RDKit container is not running. Setting up test environment..."
    setup_test_environment
  fi
  
  # Open an interactive shell in the container
  podman exec -it "$CONTAINER_NAME" bash -c "
    cd $TEST_DIR && \
    echo -e \"\\n${GREEN}RDKit Test Environment Shell${NC}\" && \
    echo -e \"Current directory: $TEST_DIR\\n\" && \
    python3 -c 'import rdkit; print(\"RDKit version:\", rdkit.__version__)' && \
    bash
  "
}

show_status() {
  echo -e "${BLUE}RDKit Test Environment Status:${NC}"
  
  # Check if the RDKit container exists
  if podman container exists "$CONTAINER_NAME"; then
    # Check if the container is running
    if podman container inspect "$CONTAINER_NAME" --format '{{.State.Running}}' | grep -q "true"; then
      echo -e "Container status: ${GREEN}Running${NC}"
      
      # Get available test scripts
      echo -e "${BLUE}Available test scripts:${NC}"
      podman exec -it "$CONTAINER_NAME" bash -c "
        ls $TEST_DIR/*rdkit*.py 2>/dev/null || echo 'No test scripts found'
      "
      
      # Show RDKit version
      echo -e "${BLUE}RDKit information:${NC}"
      podman exec -it "$CONTAINER_NAME" bash -c "
        python3 -c 'import rdkit; print(\"RDKit version:\", rdkit.__version__)'
      "
    else
      echo -e "Container status: ${YELLOW}Stopped${NC}"
    fi
  else
    echo -e "Container status: ${RED}Not Created${NC}"
  fi
}

stop_container() {
  echo -e "${BLUE}Stopping RDKit container...${NC}"
  
  # Check if the RDKit container exists and is running
  if podman container exists "$CONTAINER_NAME" && \
     podman container inspect "$CONTAINER_NAME" --format '{{.State.Running}}' | grep -q "true"; then
    podman stop "$CONTAINER_NAME"
    echo -e "${GREEN}RDKit container stopped${NC}"
  else
    echo "RDKit container is not running"
  fi
}

clean_environment() {
  echo -e "${BLUE}Cleaning up RDKit test environment...${NC}"
  
  # Ask for confirmation
  read -p "This will remove the RDKit container and all test files. Continue? (y/n): " -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Operation cancelled"
    return
  fi
  
  # Stop and remove the RDKit container
  if podman container exists "$CONTAINER_NAME"; then
    if podman container inspect "$CONTAINER_NAME" --format '{{.State.Running}}' | grep -q "true"; then
      podman stop "$CONTAINER_NAME"
    fi
    podman rm "$CONTAINER_NAME"
    echo "Removed RDKit container"
  fi
  
  # Clean up temporary directories
  echo "Cleaning up temporary directories..."
  rm -rf /tmp/mock_modules /tmp/cryoprotect-test
  
  echo -e "${GREEN}Cleanup complete${NC}"
}

# Main command dispatcher
case "$1" in
  "setup")
    setup_test_environment
    ;;
  "run")
    run_test_script "$2"
    ;;
  "shell")
    open_shell
    ;;
  "status")
    show_status
    ;;
  "stop")
    stop_container
    ;;
  "clean")
    clean_environment
    ;;
  "help"|"--help"|"-h")
    show_help
    ;;
  *)
    echo -e "${RED}Error: Unknown command '$1'${NC}"
    show_help
    exit 1
    ;;
esac

exit 0