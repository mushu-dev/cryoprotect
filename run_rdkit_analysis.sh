#!/bin/bash
# Run RDKit analysis for CryoProtect

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_header() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

# Check if we have the RDKit wrapper and calculator
if [ ! -f "rdkit_wrapper.py" ]; then
    echo -e "${RED}Error: rdkit_wrapper.py not found${NC}"
    echo "Make sure you're in the correct directory"
    exit 1
fi

if [ ! -f "rdkit_property_calculator.py" ]; then
    echo -e "${RED}Error: rdkit_property_calculator.py not found${NC}"
    echo "Make sure you're in the correct directory"
    exit 1
fi

# Check if we have the required Python libraries
python -c "import supabase, tqdm, dotenv" 2>/dev/null
if [ $? -ne 0 ]; then
    echo -e "${YELLOW}Installing required Python packages...${NC}"
    pip install python-supabase tqdm python-dotenv
fi

# Make sure the script is executable
chmod +x rdkit_property_calculator.py

# Run the calculator based on command-line arguments
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    python rdkit_property_calculator.py --help
    exit 0
fi

# If no arguments, show usage
if [ $# -eq 0 ]; then
    print_header "CryoProtect RDKit Property Analysis"
    echo "This script calculates molecular properties that predict cryoprotectant effectiveness."
    echo
    echo "Options:"
    echo "  --known     Calculate properties for known cryoprotectants only"
    echo "  --sample N  Calculate for a random sample of N molecules"
    echo "  --limit N   Limit to N molecules"
    echo "  --help      Show detailed help"
    echo
    echo "Example usage:"
    echo "  ./run_rdkit_analysis.sh --known"
    echo "  ./run_rdkit_analysis.sh --sample 50"
    echo
    echo "For more information about the molecular properties and their significance,"
    echo "see RDKit_CRYOPROTECTANT_GUIDE.md"
    exit 0
fi

# Log to both console and file
print_header "Starting RDKit Property Analysis"
echo "Parameters: $@"
echo "Log file: rdkit_analysis.log"
echo

# Run the Python script with the provided arguments
python rdkit_property_calculator.py "$@" 2>&1 | tee -a rdkit_analysis.log

# Check if it was successful
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    print_header "Analysis Completed Successfully"
    echo -e "${GREEN}Results have been stored in the database.${NC}"
    echo "See rdkit_analysis.log for details"
    echo
    echo "Next steps:"
    echo "1. Query the database to find top cryoprotectant candidates"
    echo "2. Use the 'cryoprotectant_score' property to rank molecules"
    echo "3. Examine high-scoring molecules for experimental validation"
else
    print_header "Analysis Failed"
    echo -e "${RED}The RDKit analysis encountered errors.${NC}"
    echo "Check rdkit_analysis.log for details"
fi