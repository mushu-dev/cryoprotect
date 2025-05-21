#!/bin/bash
# Runs the minimal Flask application without containers

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Minimal App Runner${NC}"
echo "Running minimal Flask application (no dependencies on RDKit, etc.)"
echo "==========================================================================="

# Create and activate virtual environment if it doesn't exist
if [ ! -d "minimal_venv" ]; then
    echo -e "${BLUE}Creating virtual environment...${NC}"
    python3 -m venv minimal_venv
    echo -e "${GREEN}Virtual environment created.${NC}"
fi

echo -e "${BLUE}Activating virtual environment...${NC}"
source minimal_venv/bin/activate

# Install only the absolute minimal dependencies
echo -e "${BLUE}Installing minimal dependencies...${NC}"
pip install flask python-dotenv requests

# Load environment variables
if [ -f .env ]; then
    source .env
    echo -e "${GREEN}Loaded environment variables from .env file${NC}"
else
    echo -e "${YELLOW}Warning: No .env file found. Using default environment variables.${NC}"
    export FLASK_APP=minimal_app.py
    export FLASK_ENV=development
fi

# Run the minimal application
echo -e "${BLUE}Starting minimal Flask application...${NC}"
echo -e "${GREEN}The application will be available at http://localhost:5000${NC}"
python minimal_app.py

# Deactivate virtual environment when done
deactivate