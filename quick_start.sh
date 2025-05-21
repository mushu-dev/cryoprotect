#!/bin/bash
# Quick start script for CryoProtect that focuses on the minimal app
# Author: Claude

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Quick Start${NC}"
echo "This script runs the minimal app with only essential dependencies"
echo "==========================================================================="

# Create quick virtual environment
echo -e "${BLUE}Creating quick virtual environment...${NC}"
python -m venv quick_env
source quick_env/bin/activate

# Install only essential dependencies for minimal app
echo -e "${BLUE}Installing minimal dependencies...${NC}"
pip install --upgrade pip
pip install flask python-dotenv requests

# Check if .env file exists and create if needed
if [ ! -f .env ]; then
    echo -e "${YELLOW}No .env file found. Creating a template...${NC}"
    cat > .env << EOF
# Supabase credentials (required for application to run)
SUPABASE_URL=https://your-project.supabase.co
SUPABASE_KEY=your-supabase-anon-key
SECRET_KEY=dev-secret-key-please-change-in-production

# Development settings
FLASK_APP=minimal_app.py
FLASK_ENV=development
LOG_LEVEL=DEBUG
LOG_TO_FILE=1
LOG_TO_CONSOLE=1
EOF
    echo -e "${YELLOW}Please edit the .env file with your Supabase credentials.${NC}"
    read -p "Edit .env now? (y/n): " edit_env
    if [[ "$edit_env" == "y" || "$edit_env" == "Y" ]]; then
        editor .env 2>/dev/null || nano .env 2>/dev/null || vi .env
    fi
fi

# Load environment variables
source .env

# Check for valid Supabase URL and key
if [[ "$SUPABASE_URL" == "https://your-project.supabase.co" || "$SUPABASE_KEY" == "your-supabase-anon-key" ]]; then
    echo -e "${RED}Warning: You need to update your Supabase credentials in .env${NC}"
    read -p "Continue anyway? (y/n): " continue_anyway
    if [[ "$continue_anyway" != "y" && "$continue_anyway" != "Y" ]]; then
        echo "Aborted."
        deactivate
        exit 1
    fi
fi

# Run the minimal app
echo -e "${BLUE}Starting minimal Flask application...${NC}"
echo -e "${GREEN}The application will be available at http://localhost:5000${NC}"
python minimal_app.py

# Deactivate virtual environment when done
deactivate