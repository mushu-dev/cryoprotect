#!/bin/bash
# Fix Python environment for CryoProtect

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Python Environment Fixer${NC}"
echo "This script creates a clean Python virtual environment with all required dependencies"
echo "==========================================================================="

# Ask before proceeding
read -p "This will create a new virtual environment. Continue? (y/n): " confirm
if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
    echo "Aborted."
    exit 1
fi

# Create a new virtual environment
echo -e "${BLUE}Creating new virtual environment...${NC}"
python -m venv fresh_env
source fresh_env/bin/activate

# Upgrade pip
echo -e "${BLUE}Upgrading pip...${NC}"
pip install --upgrade pip

# Install dependencies one group at a time
echo -e "${BLUE}Installing Flask and extensions...${NC}"
pip install flask==3.0.2 flask-restful==0.3.10 flask-cors==4.0.0
if [ $? -ne 0 ]; then echo -e "${RED}Error installing Flask packages${NC}"; exit 1; fi

echo -e "${BLUE}Installing API documentation packages...${NC}"
pip install apispec==6.3.0 marshmallow==3.21.0 flask-apispec==0.11.4
if [ $? -ne 0 ]; then echo -e "${RED}Error installing API documentation packages${NC}"; exit 1; fi

echo -e "${BLUE}Installing utility packages...${NC}"
pip install python-dotenv==1.0.1 requests==2.31.0 flask-mail==0.9.1
if [ $? -ne 0 ]; then echo -e "${RED}Error installing utility packages${NC}"; exit 1; fi

echo -e "${BLUE}Installing security packages...${NC}"
pip install PyJWT==2.8.0 cryptography==42.0.5
if [ $? -ne 0 ]; then echo -e "${RED}Error installing security packages${NC}"; exit 1; fi

echo -e "${BLUE}Installing Supabase packages...${NC}"
pip install supabase==2.1.0 postgrest==0.13.2
if [ $? -ne 0 ]; then echo -e "${RED}Error installing Supabase packages${NC}"; exit 1; fi

echo -e "${BLUE}Installing numerical/scientific packages...${NC}"
pip install numpy scipy pandas matplotlib scikit-learn
if [ $? -ne 0 ]; then echo -e "${RED}Error installing scientific packages${NC}"; exit 1; fi

echo -e "${BLUE}Installing database packages...${NC}"
pip install psycopg2-binary==2.9.9
if [ $? -ne 0 ]; then echo -e "${RED}Error installing database packages${NC}"; exit 1; fi

echo -e "${BLUE}Installing monitoring packages...${NC}"
pip install prometheus-client==0.19.0
if [ $? -ne 0 ]; then echo -e "${RED}Error installing monitoring packages${NC}"; exit 1; fi

# Ask about RDKit installation
read -p "Do you want to install RDKit? This is a large package. (y/n): " install_rdkit
if [[ "$install_rdkit" == "y" || "$install_rdkit" == "Y" ]]; then
    echo -e "${BLUE}Installing RDKit...${NC}"
    pip install rdkit-pypi
    if [ $? -ne 0 ]; then echo -e "${RED}Error installing RDKit${NC}"; fi
fi

# Create a runner script
echo -e "${BLUE}Creating runner script...${NC}"
cat > run_fresh.sh << 'EOF'
#!/bin/bash
# Run CryoProtect with fresh Python environment

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Fresh Environment Runner${NC}"
echo "Running with clean Python environment"
echo "==========================================================================="

# Activate virtual environment
source fresh_env/bin/activate

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

chmod +x run_fresh.sh

echo -e "${GREEN}Environment setup complete!${NC}"
echo "You can now run the application with: ./run_fresh.sh"

# Deactivate virtual environment
deactivate