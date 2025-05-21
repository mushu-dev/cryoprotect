#!/bin/bash
# Run ChEMBL Verification Phase 5
# This script orchestrates the execution of Phase 5 ChEMBL verification components

set -e  # Exit on any error

# Define colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Directory for logs
mkdir -p logs

# Log file
LOG_FILE="logs/chembl_verification_phase5_$(date +%Y%m%d_%H%M%S).log"

# Function to log messages
log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo -e "$msg" | tee -a "$LOG_FILE"
}

log "${GREEN}Starting ChEMBL Verification Phase 5${NC}"

# Step 1: Run the property fixer
log "${YELLOW}Step 1: Fixing missing properties in ChEMBL molecules...${NC}"
python fix_missing_chembl_properties.py --verbose
if [ $? -eq 0 ]; then
    log "${GREEN}✓ Successfully fixed missing properties${NC}"
else
    log "${RED}✗ Failed to fix missing properties${NC}"
    exit 1
fi

# Step 2: Resolve missing PubChem cross-references
log "${YELLOW}Step 2: Resolving missing PubChem cross-references...${NC}"
python resolve_missing_pubchem_references.py --verbose
if [ $? -eq 0 ]; then
    log "${GREEN}✓ Successfully resolved missing PubChem cross-references${NC}"
else
    log "${RED}✗ Failed to resolve missing PubChem cross-references${NC}"
    exit 1
fi

# Step 3: Verify ChEMBL import comprehensively
log "${YELLOW}Step 3: Running comprehensive ChEMBL verification...${NC}"
python comprehensive_chembl_verification.py --generate_visualizations
if [ $? -eq 0 ]; then
    log "${GREEN}✓ Successfully completed comprehensive verification${NC}"
else
    log "${RED}✗ Failed to complete comprehensive verification${NC}"
    exit 1
fi

# Step 4: Generate verification report
log "${YELLOW}Step 4: Generating verification report...${NC}"
python generate_chembl_verification_report.py --visualizations
if [ $? -eq 0 ]; then
    log "${GREEN}✓ Successfully generated verification report${NC}"
else
    log "${RED}✗ Failed to generate verification report${NC}"
    exit 1
fi

# Make the report more accessible
LATEST_REPORT=$(find reports -name "chembl_verification_report_*.html" | sort -r | head -n 1)
if [ -n "$LATEST_REPORT" ]; then
    cp "$LATEST_REPORT" "reports/latest_chembl_verification_report.html"
    log "${GREEN}Latest report copied to reports/latest_chembl_verification_report.html${NC}"
fi

log "${GREEN}ChEMBL Verification Phase 5 completed successfully!${NC}"
log "The verification report is available in the reports directory."
log "To view the report, open the HTML file in a web browser."

# Make script executable
chmod +x run_chembl_verification_phase5.sh