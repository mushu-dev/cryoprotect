#!/bin/bash
# Setup SELinux contexts for CryoProtect directories (requires sudo)

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect SELinux Context Setup${NC}"
echo "This script will set the proper SELinux contexts for volume directories"
echo "==========================================================================="

# Check if running as root
if [ "$EUID" -ne 0 ]; then
  echo -e "${RED}Error: This script must be run with sudo${NC}"
  echo "Please run: sudo ./setup_selinux_sudo.sh"
  exit 1
fi

echo -e "${BLUE}Creating required directories if they don't exist...${NC}"
mkdir -p logs
mkdir -p backup/data
mkdir -p data
mkdir -p cache

echo -e "${BLUE}Setting SELinux contexts...${NC}"

# Set context for logs directory
echo "Setting context for logs directory..."
chcon -Rt container_file_t ./logs
echo -e "${GREEN}✓ Logs directory context set${NC}"

# Set context for backup directory and subdirectories
echo "Setting context for backup directory..."
chcon -R -t container_file_t ./backup
echo -e "${GREEN}✓ Backup directory context set${NC}"

# Set context for data directory
echo "Setting context for data directory..."
chcon -Rt container_file_t ./data
echo -e "${GREEN}✓ Data directory context set${NC}"

# Set context for cache directory
echo "Setting context for cache directory..."
chcon -Rt container_file_t ./cache
echo -e "${GREEN}✓ Cache directory context set${NC}"

# Set context for database directory if it exists
if [ -d "./database" ]; then
  echo "Setting context for database directory..."
  chcon -Rt container_file_t ./database
  echo -e "${GREEN}✓ Database directory context set${NC}"
fi

echo -e "\n${GREEN}SELinux contexts set successfully!${NC}"
echo -e "${YELLOW}Now you can run: ./simple_podman_run.sh${NC}"