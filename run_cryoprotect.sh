#!/bin/bash
# CryoProtect runner script with multiple security options

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Display header
echo -e "${BLUE}"
echo "  _____ _____          _____  _____   ____  _______ ______ _____ _______ "
echo " / ____|  __ \        |  __ \|  __ \ / __ \|__   __|  ____|  __ \__   __|"
echo "| |    | |__) |_   _  | |__) | |__) | |  | |  | |  | |__  | |  | | | |   "
echo "| |    |  _  /| | | | |  ___/|  _  /| |  | |  | |  |  __| | |  | | | |   "
echo "| |____| | \ \| |_| | | |    | | \ \| |__| |  | |  | |____| |__| | | |   "
echo " \_____|_|  \_\\__, | |_|    |_|  \_\\____/   |_|  |______|_____/  |_|   "
echo "                __/ |                                                     "
echo "               |___/           Fedora Runner Script                      "
echo -e "${NC}"

# Load environment variables
if [ -f .env ]; then
  source .env
  echo -e "${GREEN}✓ Loaded environment variables from .env file${NC}"
else
  echo -e "${YELLOW}⚠ No .env file found. Using default environment variables.${NC}"
  export FLASK_APP=app.py
  export FLASK_ENV=development
  export SECRET_KEY=dev-secret-key-please-change-in-production
  export LOG_LEVEL=DEBUG
fi

# Check Supabase URL and Key
if [ -z "$SUPABASE_URL" ] || [ -z "$SUPABASE_KEY" ]; then
  echo -e "${YELLOW}⚠ Supabase URL or key not set in .env file.${NC}"
  echo "Your application may not be able to connect to the database."
  echo
fi

# Function to confirm actions
confirm() {
  read -p "$(echo -e $CYAN"$1 [y/N]: "$NC)" choice
  case "$choice" in
    y|Y|yes|Yes|YES ) return 0;;
    * ) return 1;;
  esac
}

# Create necessary directories
mkdir -p logs
mkdir -p data
mkdir -p backup/data
mkdir -p cache

# Check if cryoprotect-net network exists, create if not
if ! podman network ls | grep -q "cryoprotect-net"; then
  echo -e "${YELLOW}Network 'cryoprotect-net' not found. Creating...${NC}"
  podman network create --ipv6=false cryoprotect-net
  echo -e "${GREEN}✓ Created network 'cryoprotect-net'${NC}"
fi

# Main menu
echo -e "${BLUE}=================================${NC}"
echo -e "${BLUE}     SELECT SECURITY MODE        ${NC}"
echo -e "${BLUE}=================================${NC}"
echo
echo "1. Run with SELinux enabled (requires setting contexts with sudo)"
echo "2. Run with SELinux disabled (easier but less secure)"
echo "3. Exit"
echo
read -p "$(echo -e $CYAN"Enter your choice [1-3]: "$NC)" mode_choice

case $mode_choice in
  1)
    echo -e "${BLUE}Running with SELinux enabled${NC}"
    
    # Check if contexts are already set
    if [ "$(ls -lZ ./logs | grep -c container_file_t)" -eq 0 ]; then
      echo -e "${YELLOW}⚠ SELinux contexts not set for volume directories.${NC}"
      echo "You need to run the setup_selinux_sudo.sh script with sudo first."
      
      if confirm "Would you like to run sudo ./setup_selinux_sudo.sh now?"; then
        sudo ./setup_selinux_sudo.sh
        if [ $? -ne 0 ]; then
          echo -e "${RED}✗ Failed to set SELinux contexts.${NC}"
          exit 1
        fi
      else
        echo -e "${YELLOW}⚠ Continuing without setting SELinux contexts.${NC}"
        echo "This may cause permission errors when mounting volumes."
      fi
    else
      echo -e "${GREEN}✓ SELinux contexts already set.${NC}"
    fi
    
    # Run with SELinux enabled
    echo -e "${BLUE}Starting CryoProtect container with SELinux enabled...${NC}"
    podman run --rm -it \
      --name cryoprotect-app \
      --network=cryoprotect-net \
      --dns=1.1.1.1 --dns=8.8.8.8 \
      -p 5000:5000 \
      -v .:/app:Z \
      -v ./logs:/app/logs:Z \
      -v ./data:/app/data:Z \
      -v ./cache:/app/cache:Z \
      -e FLASK_APP=app.py \
      -e FLASK_ENV=development \
      -e FLASK_DEBUG=1 \
      -e SUPABASE_URL="${SUPABASE_URL}" \
      -e SUPABASE_KEY="${SUPABASE_KEY}" \
      -e SECRET_KEY="${SECRET_KEY}" \
      -e LOG_LEVEL="${LOG_LEVEL:-DEBUG}" \
      -e LOG_TO_FILE="${LOG_TO_FILE:-1}" \
      -e LOG_TO_CONSOLE="${LOG_TO_CONSOLE:-1}" \
      --security-opt label=type:container_file_t \
      python:3.10-slim python -m flask run --host=0.0.0.0 --port=5000
    ;;
    
  2)
    echo -e "${BLUE}Running with SELinux disabled${NC}"
    
    # Run with SELinux disabled
    echo -e "${BLUE}Starting CryoProtect container with SELinux disabled...${NC}"
    podman run --rm -it \
      --name cryoprotect-app \
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
      python:3.10-slim python -m flask run --host=0.0.0.0 --port=5000
    ;;
    
  3)
    echo "Exiting..."
    exit 0
    ;;
    
  *)
    echo -e "${RED}✗ Invalid choice.${NC}"
    exit 1
    ;;
esac