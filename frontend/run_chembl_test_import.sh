#!/bin/bash
# Run a test import of ChEMBL data and validate the results

set -e

echo "Starting ChEMBL test import and validation..."

# Set up environment
source .env 2>/dev/null || true
export LOG_LEVEL=${LOG_LEVEL:-INFO}

# Create output directory
TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
OUTPUT_DIR="chembl_test_import_${TIMESTAMP}"
mkdir -p "$OUTPUT_DIR"

echo "Running test import with limited scope (10 compounds)..."
python unified_chembl_import.py --limit 10 --batch-size 5 --checkpoint-interval 5 --output-dir "$OUTPUT_DIR" --test-mode

echo "Validating imported data..."
python test_chembl_import_validation.py --output "${OUTPUT_DIR}/validation_report.json"

echo "Generating summary report..."
echo "# ChEMBL Test Import Results" > "${OUTPUT_DIR}/summary.md"
echo "" >> "${OUTPUT_DIR}/summary.md"
echo "**Date:** $(date "+%Y-%m-%d %H:%M:%S")" >> "${OUTPUT_DIR}/summary.md"
echo "" >> "${OUTPUT_DIR}/summary.md"
echo "## Import Statistics" >> "${OUTPUT_DIR}/summary.md"
grep "Import statistics:" "$OUTPUT_DIR/chembl_import.log" | tail -1 | sed 's/^.*statistics: //' >> "${OUTPUT_DIR}/summary.md"
echo "" >> "${OUTPUT_DIR}/summary.md"
echo "## Validation Summary" >> "${OUTPUT_DIR}/summary.md"
grep -A5 "^## Summary Statistics" "${OUTPUT_DIR}/validation_report.md" >> "${OUTPUT_DIR}/summary.md"
echo "" >> "${OUTPUT_DIR}/summary.md"
echo "## Next Steps" >> "${OUTPUT_DIR}/summary.md"
echo "1. Review validation results in detail" >> "${OUTPUT_DIR}/summary.md"
echo "2. Address any critical issues identified" >> "${OUTPUT_DIR}/summary.md"
echo "3. Scale up import to full dataset if test successful" >> "${OUTPUT_DIR}/summary.md"
echo "" >> "${OUTPUT_DIR}/summary.md"
echo "For detailed validation results, see the [validation report](validation_report.md)" >> "${OUTPUT_DIR}/summary.md"

echo "Test import and validation complete. Results saved to ${OUTPUT_DIR}/"
echo "Summary report: ${OUTPUT_DIR}/summary.md"