#!/bin/bash
# Run ChEMBL Import Verification
set -e

# Define directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPORTS_DIR="${SCRIPT_DIR}/reports"

# Create reports directory if it doesn't exist
mkdir -p "${REPORTS_DIR}"

# Define colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}Starting ChEMBL Import Verification...${NC}"
echo "This script will generate a comprehensive report on the integrity of ChEMBL data imports."

# Check for Python and required dependencies
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Python 3 is required but not found. Please install Python 3 and try again.${NC}"
    exit 1
fi

# Check for virtual environment
if [ -d "venv" ] || [ -d ".venv" ]; then
    # Activate virtual environment if it exists
    if [ -f "venv/bin/activate" ]; then
        echo "Activating virtual environment..."
        source venv/bin/activate
    elif [ -f ".venv/bin/activate" ]; then
        echo "Activating virtual environment..."
        source .venv/bin/activate
    fi
fi

# Install required dependencies if needed
echo -e "${YELLOW}Checking dependencies...${NC}"
python3 -c "
try:
    import psycopg2
except ImportError:
    print('Installing psycopg2...')
    import subprocess
    subprocess.call(['pip', 'install', 'psycopg2-binary'])

try:
    import pandas
except ImportError:
    print('Installing pandas (recommended for statistical analysis)...')
    import subprocess
    subprocess.call(['pip', 'install', 'pandas'])

try:
    from rdkit import Chem
except ImportError:
    print('Note: RDKit not available. Structure validation will be skipped.')
    print('To enable structure validation, install RDKit: pip install rdkit')
"

# Run verification script
echo -e "${GREEN}Running verification script...${NC}"
python3 "${SCRIPT_DIR}/verify_chembl_import_integrity.py" --config=config.py --output="${REPORTS_DIR}"

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Verification completed successfully!${NC}"
    
    # Find the latest report
    LATEST_REPORT=$(ls -t "${REPORTS_DIR}"/chembl_verification_report_*.md 2>/dev/null | head -n 1)
    
    if [ -n "${LATEST_REPORT}" ]; then
        echo -e "${GREEN}Report generated: ${LATEST_REPORT}${NC}"
        echo ""
        echo -e "${YELLOW}Report Summary:${NC}"
        echo "========================================"
        # Extract the General Statistics section
        grep -A 10 "^## General Statistics" "${LATEST_REPORT}"
        echo "========================================"
        echo -e "${YELLOW}Recommendations:${NC}"
        echo "========================================"
        grep -A 20 "^## Recommendations" "${LATEST_REPORT}" | grep -v "^## Specific Issues"
        echo "========================================"
    else
        echo -e "${YELLOW}No report was generated. Check logs for errors.${NC}"
    fi
else
    echo -e "${RED}Verification failed. See error messages above.${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}Complete verification process finished.${NC}"
echo "For detailed results, see the generated reports in the ${REPORTS_DIR} directory."