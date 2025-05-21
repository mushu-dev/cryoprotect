#!/bin/bash
# Helper script to quickly install missing dependencies
# Author: Claude

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if package name was provided
if [ -z "$1" ]; then
    echo -e "${RED}Error: No package name provided${NC}"
    echo "Usage: $0 <package_name>"
    echo "Example: $0 numpy"
    exit 1
fi

PACKAGE=$1

echo -e "${BLUE}Installing missing dependency: ${PACKAGE}${NC}"
echo "==========================================================================="

# Activate virtual environment
echo -e "${BLUE}Activating virtual environment...${NC}"
source quick_env/bin/activate

# Install the package
echo -e "${BLUE}Installing ${PACKAGE}...${NC}"
pip install ${PACKAGE}

# Check if installation was successful
if [ $? -eq 0 ]; then
    echo -e "${GREEN}Successfully installed ${PACKAGE}${NC}"
else
    echo -e "${RED}Failed to install ${PACKAGE}${NC}"
    echo "Try specifying a different version or check for compatibility issues."
fi

# Deactivate virtual environment
deactivate

echo -e "${BLUE}You can now try running the application again.${NC}"