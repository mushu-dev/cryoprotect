#!/bin/bash
# run_database_integrity_check_fixed.sh
#
# Runs the fixed database integrity verification script and generates a
# comprehensive report of all issues found.
#
# Usage:
#   ./run_database_integrity_check_fixed.sh [--verbose] [--report=filename.json]

set -e

echo "CryoProtect Database Integrity Verification (Fixed Version)"
echo "=========================================================="

# Parse command line arguments
VERBOSE=""
REPORT_FILE="database_integrity_report_$(date +%Y%m%d_%H%M%S).json"

for arg in "$@"; do
  case $arg in
    --verbose)
      VERBOSE="--verbose"
      shift
      ;;
    --report=*)
      REPORT_FILE="${arg#*=}"
      shift
      ;;
    *)
      # Unknown option
      ;;
  esac
done

echo "Starting verification process..."
echo "Report will be saved to: $REPORT_FILE"

# Check Python environment
if command -v conda &> /dev/null; then
  # Conda environment
  echo "Activating conda environment..."
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate cryoprotect || echo "Warning: Failed to activate conda environment"
fi

# Run the verification script
python verify_database_integrity_fixed.py --output="$REPORT_FILE" $VERBOSE

# Check the exit status
if [ $? -eq 0 ]; then
  echo -e "\n✅ Database integrity verification PASSED"
else
  echo -e "\n❌ Database integrity verification FAILED - See report for details"
fi

# Count issues by severity
if [ -f "$REPORT_FILE" ]; then
  ERROR_COUNT=$(grep -o '"severity": "error"' "$REPORT_FILE" | wc -l)
  WARNING_COUNT=$(grep -o '"severity": "warning"' "$REPORT_FILE" | wc -l)
  
  echo -e "\nSummary:"
  echo "- Error issues: $ERROR_COUNT"
  echo "- Warning issues: $WARNING_COUNT"
  echo -e "\nFull report saved to: $REPORT_FILE"
  
  # Create a markdown summary report
  MD_REPORT="${REPORT_FILE%.json}.md"
  
  echo "# Database Integrity Verification Report" > "$MD_REPORT"
  echo "Generated: $(date)" >> "$MD_REPORT"
  echo "" >> "$MD_REPORT"
  echo "## Summary" >> "$MD_REPORT"
  echo "- Error issues: $ERROR_COUNT" >> "$MD_REPORT"
  echo "- Warning issues: $WARNING_COUNT" >> "$MD_REPORT"
  echo "" >> "$MD_REPORT"
  
  # Extract table counts
  echo "## Data Counts" >> "$MD_REPORT"
  echo "| Table | Row Count |" >> "$MD_REPORT"
  echo "| ----- | --------- |" >> "$MD_REPORT"
  
  # Use jq to extract data counts if available
  if command -v jq &> /dev/null; then
    jq -r '.data_counts | to_entries[] | "| \(.key) | \(.value) |"' "$REPORT_FILE" >> "$MD_REPORT" 2>/dev/null || 
      echo "| Failed to extract table counts | |" >> "$MD_REPORT"
  else
    echo "| Install jq for detailed table counts | |" >> "$MD_REPORT"
  fi
  
  echo "" >> "$MD_REPORT"
  echo "## Issues" >> "$MD_REPORT"
  
  # Use jq to extract issues if available
  if command -v jq &> /dev/null; then
    # Error issues first
    jq -r '.issues[] | select(.severity == "error") | "### ⛔ \(.title)\n\n\(.description)\n\n---\n"' "$REPORT_FILE" >> "$MD_REPORT" 2>/dev/null
    
    # Then warning issues
    jq -r '.issues[] | select(.severity == "warning") | "### ⚠️ \(.title)\n\n\(.description)\n\n---\n"' "$REPORT_FILE" >> "$MD_REPORT" 2>/dev/null
  else
    echo "Install jq for detailed issue listing" >> "$MD_REPORT"
  fi
  
  echo "Markdown summary also saved to: $MD_REPORT"
  
  # Generate HTML report if the generator script exists
  if [ -f "generate_database_integrity_report.py" ]; then
    HTML_REPORT="${REPORT_FILE%.json}.html"
    echo "Generating HTML report..."
    python generate_database_integrity_report.py --input="$REPORT_FILE" --output="$HTML_REPORT"
    echo "HTML report saved to: $HTML_REPORT"
  fi
fi

echo -e "\nVerification process completed"