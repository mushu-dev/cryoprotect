#!/bin/bash
# Run CryoProtect application in the enriched environment
# Author: Claude

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Application Runner${NC}"
echo "This script runs the CryoProtect application in the enriched virtual environment"
echo "==========================================================================="

# Activate virtual environment
echo -e "${BLUE}Activating virtual environment...${NC}"
source quick_env/bin/activate

# Load environment variables
if [ -f .env ]; then
    source .env
    echo -e "${GREEN}Loaded environment variables from .env file${NC}"
else
    echo -e "${RED}No .env file found. This is required.${NC}"
    deactivate
    exit 1
fi

# Choose which app to run
echo "Select which app to run:"
echo "1. Minimal test app (minimal_app.py)"
echo "2. Full application (app.py)"
read -p "Enter choice [1-2]: " app_choice

if [ "$app_choice" == "2" ]; then
    # Run full application
    echo -e "${BLUE}Starting full Flask application...${NC}"
    echo -e "${YELLOW}Note: This might fail if some dependencies are still missing.${NC}"
    echo -e "${GREEN}The application will be available at http://localhost:5000${NC}"
    export FLASK_APP=app.py
    export FLASK_ENV=development
    echo -e "${BLUE}Running command: python -m flask run --host=0.0.0.0 --port=5000${NC}"
    python -m flask run --host=0.0.0.0 --port=5000
else
    # Run minimal application
    echo -e "${BLUE}Starting minimal Flask application...${NC}"
    echo -e "${GREEN}The application will be available at http://localhost:5000${NC}"
    echo -e "${BLUE}Running command: python minimal_app.py${NC}"
    python minimal_app.py
fi

# Deactivate virtual environment when done
deactivate