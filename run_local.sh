#!/bin/bash
# CryoProtect Local Runner (no containers)
# Runs the application in a Python virtual environment without containers

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Local Runner${NC}"
echo "Running directly with Python virtual environment (no containers)"
echo "==========================================================================="

# Check for Python
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Python 3 is not installed. Please install Python 3 first.${NC}"
    exit 1
fi

# Create and activate virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo -e "${BLUE}Creating virtual environment...${NC}"
    python3 -m venv venv
    echo -e "${GREEN}Virtual environment created.${NC}"
fi

echo -e "${BLUE}Activating virtual environment...${NC}"
source venv/bin/activate

# Install dependencies
echo -e "${BLUE}Installing minimal dependencies...${NC}"
pip install flask==3.0.2 flask-restful==0.3.10 flask-cors==4.0.0 python-dotenv==1.0.1 requests==2.31.0 \
    apispec==6.3.0 flask-apispec==0.11.4 marshmallow==3.21.0 flask-mail==0.9.1 \
    postgrest supabase pyjwt cryptography prometheus-client \
    numpy scipy pandas matplotlib scikit-learn

# Load environment variables
if [ -f .env ]; then
    source .env
    echo -e "${GREEN}Loaded environment variables from .env file${NC}"
else
    echo -e "${YELLOW}Warning: No .env file found. Using default environment variables.${NC}"
    export FLASK_APP=app.py
    export FLASK_ENV=development
    export SECRET_KEY=dev-secret-key-please-change-in-production
    export LOG_LEVEL=DEBUG
fi

# Create needed directories
mkdir -p logs
mkdir -p data
mkdir -p cache

# Run the Flask application
echo -e "${BLUE}Starting Flask application...${NC}"
echo -e "${GREEN}The application will be available at http://localhost:5000${NC}"
python -m flask run --host=0.0.0.0 --port=5000

# Deactivate virtual environment when done
deactivate