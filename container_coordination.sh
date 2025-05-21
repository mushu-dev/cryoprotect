#!/bin/bash
# Container Coordination Script for CryoProtect
# This script helps manage and coordinate multiple containers for the CryoProtect application

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

show_help() {
  echo -e "${BLUE}CryoProtect Container Coordination Script${NC}"
  echo "Usage: $0 [command]"
  echo
  echo "Commands:"
  echo "  status              Show status of all containers"
  echo "  start-all           Start all containers"
  echo "  stop-all            Stop all containers"
  echo "  restart-all         Restart all containers"
  echo "  setup               Set up container network and configuration"
  echo "  start-rdkit         Start only the RDKit container"
  echo "  start-app           Start only the application container"
  echo "  clean               Remove all containers and temporary files"
  echo "  run-test [test]     Run a specific test in the appropriate container"
  echo "  help                Show this help message"
  echo
  echo "Examples:"
  echo "  $0 status"
  echo "  $0 start-all"
  echo "  $0 run-test verify_rdkit.py"
}

# Function to check if a container exists
container_exists() {
  podman container exists "$1"
}

# Function to check if a container is running
container_running() {
  podman container inspect "$1" --format '{{.State.Running}}' 2>/dev/null | grep -q "true"
}

# Function to show the status of all containers
show_status() {
  echo -e "${BLUE}Container Status:${NC}"
  
  # RDKit Container
  echo -n "RDKit Container (CryoProtect-RDKit-Conda): "
  if container_exists "CryoProtect-RDKit-Conda"; then
    if container_running "CryoProtect-RDKit-Conda"; then
      echo -e "${GREEN}Running${NC}"
    else
      echo -e "${YELLOW}Stopped${NC}"
    fi
  else
    echo -e "${RED}Not Created${NC}"
  fi
  
  # Main CryoProtect Container
  echo -n "Main Container (CryoProtect): "
  if container_exists "CryoProtect"; then
    if container_running "CryoProtect"; then
      echo -e "${GREEN}Running${NC}"
    else
      echo -e "${YELLOW}Stopped${NC}"
    fi
  else
    echo -e "${RED}Not Created${NC}"
  fi
  
  # Development Container
  echo -n "Development Container (cryoprotect-dev): "
  if container_exists "cryoprotect-dev"; then
    if container_running "cryoprotect-dev"; then
      echo -e "${GREEN}Running${NC}"
    else
      echo -e "${YELLOW}Stopped${NC}"
    fi
  else
    echo -e "${RED}Not Created${NC}"
  fi
  
  # Network information
  echo -e "${BLUE}Network Information:${NC}"
  podman network ls
}

# Function to set up container network and configuration
setup_environment() {
  echo -e "${BLUE}Setting up container environment...${NC}"
  
  # Create network if it doesn't exist
  if ! podman network ls | grep -q "cryoprotect-net"; then
    echo -e "${YELLOW}Network 'cryoprotect-net' not found. Creating...${NC}"
    podman network create --ipv6=false cryoprotect-net
    echo -e "${GREEN}Created network 'cryoprotect-net'${NC}"
  else
    echo -e "${GREEN}Network 'cryoprotect-net' already exists${NC}"
  fi
  
  # Create directories if they don't exist
  mkdir -p logs
  mkdir -p data
  mkdir -p backup/data
  mkdir -p cache
  
  echo -e "${GREEN}Environment setup complete${NC}"
}

# Function to start all containers
start_all_containers() {
  echo -e "${BLUE}Starting all containers...${NC}"
  
  # Start RDKit container if it exists
  if container_exists "CryoProtect-RDKit-Conda"; then
    if ! container_running "CryoProtect-RDKit-Conda"; then
      echo "Starting RDKit container..."
      podman start CryoProtect-RDKit-Conda
    else
      echo "RDKit container is already running"
    fi
  else
    echo -e "${YELLOW}RDKit container does not exist. Creating...${NC}"
    ./quick_conda_container.sh
  fi
  
  # Start main container if it exists
  if container_exists "CryoProtect"; then
    if ! container_running "CryoProtect"; then
      echo "Starting main container..."
      podman start CryoProtect
    else
      echo "Main container is already running"
    fi
  else
    echo -e "${YELLOW}Main container does not exist. Please run create_cryoprotect_container.sh${NC}"
  fi
  
  # Start dev container if it exists (using podman-compose)
  if container_exists "cryoprotect-dev"; then
    if ! container_running "cryoprotect-dev"; then
      echo "Starting development container..."
      podman-compose -f podman-compose.minimal.yml up -d
    else
      echo "Development container is already running"
    fi
  else
    echo -e "${YELLOW}Development container does not exist. You can create it with podman-compose${NC}"
  fi
  
  echo -e "${GREEN}All containers started${NC}"
}

# Function to stop all containers
stop_all_containers() {
  echo -e "${BLUE}Stopping all containers...${NC}"
  
  # Stop RDKit container if it exists and is running
  if container_exists "CryoProtect-RDKit-Conda" && container_running "CryoProtect-RDKit-Conda"; then
    echo "Stopping RDKit container..."
    podman stop CryoProtect-RDKit-Conda
  fi
  
  # Stop main container if it exists and is running
  if container_exists "CryoProtect" && container_running "CryoProtect"; then
    echo "Stopping main container..."
    podman stop CryoProtect
  fi
  
  # Stop dev container if it exists and is running
  if container_exists "cryoprotect-dev" && container_running "cryoprotect-dev"; then
    echo "Stopping development container..."
    podman-compose -f podman-compose.minimal.yml down
  fi
  
  echo -e "${GREEN}All containers stopped${NC}"
}

# Function to restart all containers
restart_all_containers() {
  echo -e "${BLUE}Restarting all containers...${NC}"
  stop_all_containers
  start_all_containers
}

# Function to start only the RDKit container
start_rdkit_container() {
  echo -e "${BLUE}Starting RDKit container...${NC}"
  
  if container_exists "CryoProtect-RDKit-Conda"; then
    if ! container_running "CryoProtect-RDKit-Conda"; then
      echo "Starting RDKit container..."
      podman start CryoProtect-RDKit-Conda
    else
      echo "RDKit container is already running"
    fi
  else
    echo -e "${YELLOW}RDKit container does not exist. Creating...${NC}"
    ./quick_conda_container.sh
  fi
  
  echo -e "${GREEN}RDKit container ready${NC}"
}

# Function to start only the application container
start_app_container() {
  echo -e "${BLUE}Starting application container...${NC}"
  
  # Try to start the development container first, if it exists
  if container_exists "cryoprotect-dev"; then
    if ! container_running "cryoprotect-dev"; then
      echo "Starting development container..."
      podman-compose -f podman-compose.minimal.yml up -d
    else
      echo "Development container is already running"
    fi
  # Otherwise, try the main container
  elif container_exists "CryoProtect"; then
    if ! container_running "CryoProtect"; then
      echo "Starting main container..."
      podman start CryoProtect
    else
      echo "Main container is already running"
    fi
  else
    echo -e "${YELLOW}No application container found. Please create one first.${NC}"
  fi
}

# Function to clean up all containers and temporary files
clean_environment() {
  echo -e "${BLUE}Cleaning up containers and temporary files...${NC}"
  
  # Ask for confirmation
  read -p "This will remove all containers and temporary files. Continue? (y/n): " -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Operation cancelled"
    return
  fi
  
  # Stop and remove all containers
  stop_all_containers
  
  # Remove containers
  if container_exists "CryoProtect-RDKit-Conda"; then
    echo "Removing RDKit container..."
    podman rm CryoProtect-RDKit-Conda
  fi
  
  if container_exists "CryoProtect"; then
    echo "Removing main container..."
    podman rm CryoProtect
  fi
  
  if container_exists "cryoprotect-dev"; then
    echo "Removing development container..."
    podman rm cryoprotect-dev
  fi
  
  # Clean up temporary directories
  echo "Cleaning up temporary directories..."
  rm -rf /tmp/mock_modules /tmp/cryoprotect-data /tmp/cryoprotect-app /tmp/cryoprotect-test
  
  echo -e "${GREEN}Cleanup complete${NC}"
}

# Function to run a test in the appropriate container
run_test() {
  if [ -z "$1" ]; then
    echo -e "${RED}Error: No test specified${NC}"
    echo "Usage: $0 run-test [test_script.py]"
    exit 1
  fi
  
  TEST_SCRIPT="$1"
  
  if [[ "$TEST_SCRIPT" == *"rdkit"* ]]; then
    echo -e "${BLUE}Running RDKit test in RDKit container...${NC}"
    if ! container_exists "CryoProtect-RDKit-Conda" || ! container_running "CryoProtect-RDKit-Conda"; then
      start_rdkit_container
    fi
    ./run_with_rdkit.sh "$TEST_SCRIPT"
  else
    echo -e "${BLUE}Running general test in application container...${NC}"
    # Try to run in the development container first
    if container_exists "cryoprotect-dev" && container_running "cryoprotect-dev"; then
      podman exec -it cryoprotect-dev python "$TEST_SCRIPT"
    # Otherwise, try the main container
    elif container_exists "CryoProtect" && container_running "CryoProtect"; then
      ./run_in_cryoprotect.sh "python $TEST_SCRIPT"
    else
      echo -e "${RED}Error: No suitable container is running${NC}"
      exit 1
    fi
  fi
}

# Main command dispatcher
case "$1" in
  "status")
    show_status
    ;;
  "start-all")
    setup_environment
    start_all_containers
    ;;
  "stop-all")
    stop_all_containers
    ;;
  "restart-all")
    restart_all_containers
    ;;
  "setup")
    setup_environment
    ;;
  "start-rdkit")
    start_rdkit_container
    ;;
  "start-app")
    start_app_container
    ;;
  "clean")
    clean_environment
    ;;
  "run-test")
    run_test "$2"
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