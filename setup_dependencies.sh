#!/bin/bash
# Install all dependencies for CryoProtect

# ANSI colors
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}Installing CryoProtect dependencies...${NC}"

# System dependencies
echo -e "${YELLOW}Installing system dependencies...${NC}"
sudo dnf install -y postgresql-devel python3-devel gcc

# Create virtual environment
echo -e "${BLUE}Creating virtual environment...${NC}"
python -m venv venv
source venv/bin/activate

# Install Python dependencies
echo -e "${BLUE}Installing Python dependencies...${NC}"
pip install --upgrade pip
pip install -r requirements_full.txt

echo -e "${GREEN}Dependencies installed successfully!${NC}"
echo "You can now run the application with: ./run_local.sh"
