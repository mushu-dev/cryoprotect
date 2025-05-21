#!/bin/bash
# Setup essential dependencies for CryoProtect

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Essential Setup${NC}"
echo "This script installs only the essential dependencies needed to run the minimal app"
echo "==========================================================================="

# Create and activate virtual environment
if [ ! -d "essential_venv" ]; then
    echo -e "${BLUE}Creating virtual environment...${NC}"
    python3 -m venv essential_venv
    echo -e "${GREEN}Virtual environment created.${NC}"
fi

echo -e "${BLUE}Activating virtual environment...${NC}"
source essential_venv/bin/activate

# Install essential dependencies
echo -e "${BLUE}Installing essential dependencies...${NC}"
pip install --upgrade pip
pip install -r requirements_essential.txt

# Create a new run script for the minimal app
cat > run_essential.sh << 'EOF'
#!/bin/bash
# Run CryoProtect with essential dependencies

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Essential Runner${NC}"
echo "Running minimal Flask application with essential dependencies"
echo "==========================================================================="

# Activate virtual environment
source essential_venv/bin/activate

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
EOF

chmod +x run_essential.sh

echo -e "${GREEN}Setup complete!${NC}"
echo "You can now run the essential app with: ./run_essential.sh"
echo "This will start the minimal app with all essential dependencies installed"

# Deactivate virtual environment
deactivate