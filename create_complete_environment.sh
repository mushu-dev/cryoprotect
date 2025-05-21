#!/bin/bash
# Creates a complete Python environment for CryoProtect with all dependencies

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
echo "               |___/           Complete Environment Setup                 "
echo -e "${NC}"

# Confirm before proceeding
read -p "$(echo -e $CYAN"This will create a new Python virtual environment with all required dependencies. Continue? [y/N]: "$NC)" choice
case "$choice" in
  y|Y|yes|Yes|YES ) echo "Proceeding with setup...";;
  * ) echo "Setup cancelled."; exit 1;;
esac

# Check for system dependencies
echo -e "${BLUE}Checking system dependencies...${NC}"
if ! command -v dnf &> /dev/null; then
    echo -e "${YELLOW}Warning: DNF package manager not found. This may not be a Fedora system.${NC}"
else
    echo -e "${YELLOW}Installing required system packages...${NC}"
    echo "You may be asked for your sudo password."
    sudo dnf install -y postgresql-devel python3-devel gcc
fi

# Create virtual environment
echo -e "${BLUE}Creating virtual environment...${NC}"
python3 -m venv complete_env
echo -e "${GREEN}Virtual environment created: complete_env${NC}"

# Activate virtual environment
echo -e "${BLUE}Activating virtual environment...${NC}"
source complete_env/bin/activate

# Upgrade pip
echo -e "${BLUE}Upgrading pip...${NC}"
pip install --upgrade pip

# Install dependencies
echo -e "${BLUE}Installing dependencies...${NC}"

# Core Flask and extensions
echo -e "${YELLOW}Installing Flask and extensions...${NC}"
pip install \
    flask==3.0.2 \
    flask-restful==0.3.10 \
    flask-cors==4.0.0 \
    flask-apispec==0.11.4 \
    flask-mail==0.9.1 \
    python-dotenv==1.0.1

# API documentation and serialization
echo -e "${YELLOW}Installing API documentation and serialization...${NC}"
pip install \
    apispec==6.3.0 \
    marshmallow==3.21.0

# Database related
echo -e "${YELLOW}Installing database packages...${NC}"
pip install \
    supabase==2.1.0 \
    postgrest==0.13.2 \
    psycopg2-binary==2.9.9 \
    alembic==1.13.1 \
    sqlalchemy==2.0.27

# Security 
echo -e "${YELLOW}Installing security packages...${NC}"
pip install \
    PyJWT==2.8.0 \
    cryptography==42.0.5

# HTTP and communications
echo -e "${YELLOW}Installing HTTP and communication libraries...${NC}"
pip install \
    requests==2.31.0 \
    httpx==0.24.1

# Monitoring and observability
echo -e "${YELLOW}Installing monitoring and observability...${NC}"
pip install \
    prometheus-client==0.19.0

# Scientific computing
echo -e "${YELLOW}Installing scientific computing packages...${NC}"
pip install \
    numpy \
    scipy \
    pandas \
    matplotlib \
    scikit-learn

# Check if RDKit is needed
echo -e "${BLUE}Do you want to install RDKit for molecular modeling?${NC}"
echo "Note: This is a large package that can take time to install."
read -p "$(echo -e $CYAN"Install RDKit? [y/N]: "$NC)" install_rdkit
if [[ "$install_rdkit" =~ ^[Yy] ]]; then
    echo -e "${YELLOW}Installing RDKit...${NC}"
    pip install rdkit-pypi
fi

# Create run script
echo -e "${BLUE}Creating run script...${NC}"
cat > run_complete.sh << 'EOF'
#!/bin/bash
# Run CryoProtect using the complete environment

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Runner (Complete Environment)${NC}"
echo "Running application with complete Python environment"
echo "==========================================================================="

# Activate virtual environment
source complete_env/bin/activate

# Load environment variables
if [ -f .env ]; then
    source .env
    echo -e "${GREEN}Loaded environment variables from .env file${NC}"
else
    echo -e "${YELLOW}Warning: No .env file found.${NC}"
fi

# Choose which app to run
echo "Select which app to run:"
echo "1. Minimal test app (minimal_app.py)"
echo "2. Full application (app.py)"
read -p "Enter choice [1-2]: " app_choice

if [ "$app_choice" == "2" ]; then
    # Run full application
    echo -e "${BLUE}Starting full Flask application...${NC}"
    echo -e "${GREEN}The application will be available at http://localhost:5000${NC}"
    export FLASK_APP=app.py
    export FLASK_ENV=development
    python -m flask run --host=0.0.0.0 --port=5000
else
    # Run minimal application
    echo -e "${BLUE}Starting minimal Flask application...${NC}"
    echo -e "${GREEN}The application will be available at http://localhost:5000${NC}"
    python minimal_app.py
fi

# Deactivate virtual environment when done
deactivate
EOF

chmod +x run_complete.sh

echo -e "${GREEN}Setup complete!${NC}"
echo "You can now run the application with: ./run_complete.sh"

# Deactivate virtual environment
deactivate