#!/bin/bash
# Script to run ChEMBL import data fixes

# Set up the environment
export PYTHONPATH=$PYTHONPATH:$(pwd)

# Create reports directory if it doesn't exist
mkdir -p reports

# Get the current timestamp for reports
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Set batch size and limit (for demonstration purposes)
BATCH_SIZE=10
LIMIT=20

# Print header
echo "======================================================="
echo "            Running ChEMBL Import Data Fixes           "
echo "======================================================="
echo "Current timestamp: $(date)"
echo "Batch size: $BATCH_SIZE"
echo "Limit: $LIMIT"
echo "-------------------------------------------------------"

# Step 1: Run property fixes
echo "Step 1: Running property fixes..."
python fix_chembl_import_data.py --properties-only --batch-size $BATCH_SIZE --limit $LIMIT --report-file "reports/property_fixes_$TIMESTAMP.json"
echo "Property fixes completed."
echo "-------------------------------------------------------"

# Step 2: Run cross-reference fixes
echo "Step 2: Running cross-reference fixes..."
python fix_chembl_import_data.py --cross-refs-only --batch-size $BATCH_SIZE --limit $LIMIT --report-file "reports/cross_ref_fixes_$TIMESTAMP.json"
echo "Cross-reference fixes completed."
echo "-------------------------------------------------------"

# Step 3: Run full verification to check progress
echo "Step 3: Running verification check..."
python verify_chembl_direct_fixed.py
echo "Verification completed."
echo "-------------------------------------------------------"

echo "All fixes completed!"
echo "Check the 'reports' directory for detailed reports."
echo "======================================================="