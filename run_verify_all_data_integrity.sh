#!/bin/bash
# run_verify_all_data_integrity.sh
#
# Master script that runs a complete database verification check
# - Runs enhanced data integrity verification
# - Generates JSON report
# - Creates visual HTML report
# - Saves reports to timestamped files
# - Can be scheduled with cron
#
# Usage:
#   ./run_verify_all_data_integrity.sh [--notify] [--schedule]

set -e

# Parse command line arguments
NOTIFY=false
SCHEDULE=false

for arg in "$@"; do
  case $arg in
    --notify)
      NOTIFY=true
      shift
      ;;
    --schedule)
      SCHEDULE=true
      shift
      ;;
    *)
      # Unknown option
      ;;
  esac
done

# Set file paths with timestamps
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
REPORTS_DIR="reports"
JSON_REPORT="${REPORTS_DIR}/db_integrity_${TIMESTAMP}.json"
HTML_REPORT="${REPORTS_DIR}/db_integrity_${TIMESTAMP}.html"
MD_REPORT="${REPORTS_DIR}/db_integrity_${TIMESTAMP}.md"
LOG_FILE="${REPORTS_DIR}/db_integrity_${TIMESTAMP}.log"

# Create reports directory if it doesn't exist
mkdir -p "$REPORTS_DIR"

# Function to log messages
log() {
  echo "$(date +"%Y-%m-%d %H:%M:%S") - $1" | tee -a "$LOG_FILE"
}

# Set up error handler
error_handler() {
  log "ERROR: Verification process failed at line $1"
  if [ "$NOTIFY" = true ]; then
    echo "Database integrity verification failed. Check log: $LOG_FILE" | mail -s "ALERT: Database Integrity Check Failed" admin@example.com
  fi
  exit 1
}

# Configure error handling
trap 'error_handler $LINENO' ERR

# Banner
log "===================================================="
log "CryoProtect Database Integrity Verification"
log "===================================================="
log "Started at: $(date)"
log "Report files:"
log "- JSON: $JSON_REPORT"
log "- HTML: $HTML_REPORT"
log "- Markdown: $MD_REPORT"
log "- Log: $LOG_FILE"
log "===================================================="

# Activate environment if needed
if command -v conda &> /dev/null; then
  # Conda environment
  log "Activating conda environment..."
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate cryoprotect || log "Warning: Failed to activate conda environment"
fi

# Run data verification
log "Starting database integrity verification..."
python verify_database_integrity_enhanced.py --output="$JSON_REPORT" --verbose

# Check verification status
VERIFICATION_STATUS=$?
if [ $VERIFICATION_STATUS -eq 0 ]; then
  log "✅ Database integrity verification PASSED"
else
  log "❌ Database integrity verification FAILED - See report for details"
fi

# Generate HTML report
log "Generating visual HTML report..."
python generate_database_integrity_report.py --input="$JSON_REPORT" --output="$HTML_REPORT"

# Count issues by severity
if [ -f "$JSON_REPORT" ]; then
  if command -v jq &> /dev/null; then
    ERROR_COUNT=$(jq '.issues | map(select(.severity == "error")) | length' "$JSON_REPORT")
    WARNING_COUNT=$(jq '.issues | map(select(.severity == "warning")) | length' "$JSON_REPORT")
    TOTAL_ISSUES=$((ERROR_COUNT + WARNING_COUNT))
    
    log "Issues summary:"
    log "- Error issues: $ERROR_COUNT"
    log "- Warning issues: $WARNING_COUNT"
    log "- Total issues: $TOTAL_ISSUES"
  else
    log "Install jq for detailed issue counting"
  fi
  
  # Create a markdown summary report
  log "Generating markdown summary report..."
  
  echo "# Database Integrity Verification Report" > "$MD_REPORT"
  echo "Generated: $(date)" >> "$MD_REPORT"
  echo "" >> "$MD_REPORT"
  echo "## Summary" >> "$MD_REPORT"
  
  if command -v jq &> /dev/null; then
    echo "- Status: $(jq -r '.status' "$JSON_REPORT")" >> "$MD_REPORT"
    echo "- Tables checked: $(jq -r '.summary.tables_checked' "$JSON_REPORT")" >> "$MD_REPORT"
    echo "- Rows checked: $(jq -r '.summary.rows_checked' "$JSON_REPORT")" >> "$MD_REPORT"
    echo "- Foreign key constraints: $(jq -r '.summary.foreign_keys_checked' "$JSON_REPORT")" >> "$MD_REPORT"
    echo "- Error issues: $ERROR_COUNT" >> "$MD_REPORT"
    echo "- Warning issues: $WARNING_COUNT" >> "$MD_REPORT"
    echo "- Execution time: $(jq -r '.execution_time_seconds // "unknown"' "$JSON_REPORT") seconds" >> "$MD_REPORT"
  else
    echo "- Status: Manual verification required (install jq for automated summary)" >> "$MD_REPORT"
  fi
  
  echo "" >> "$MD_REPORT"
  echo "## Major Issues" >> "$MD_REPORT"
  
  if command -v jq &> /dev/null; then
    # Extract error issues first
    jq -r '.issues | map(select(.severity == "error")) | .[] | "### ⛔ \(.title)\n\n\(.description)\n\n---\n"' "$JSON_REPORT" >> "$MD_REPORT" 2>/dev/null || 
      echo "Failed to extract error issues details" >> "$MD_REPORT"
      
    # Then extract warning issues (limit to 5)
    jq -r '.issues | map(select(.severity == "warning")) | .[0:5] | .[] | "### ⚠️ \(.title)\n\n\(.description)\n\n---\n"' "$JSON_REPORT" >> "$MD_REPORT" 2>/dev/null || 
      echo "Failed to extract warning issues details" >> "$MD_REPORT"
      
    # Add a note if there are more warnings
    if [ "$WARNING_COUNT" -gt 5 ]; then
      REMAINING=$((WARNING_COUNT - 5))
      echo "**Note**: $REMAINING more warnings not shown. See full report for details." >> "$MD_REPORT"
    fi
  else
    echo "Install jq for detailed issue listing" >> "$MD_REPORT"
  fi
  
  log "Markdown summary saved to: $MD_REPORT"
fi

# Create symlinks to latest reports
ln -sf "$JSON_REPORT" "${REPORTS_DIR}/db_integrity_latest.json"
ln -sf "$HTML_REPORT" "${REPORTS_DIR}/db_integrity_latest.html"
ln -sf "$MD_REPORT" "${REPORTS_DIR}/db_integrity_latest.md"
ln -sf "$LOG_FILE" "${REPORTS_DIR}/db_integrity_latest.log"

log "Created symlinks to latest reports"

# Send notification if requested
if [ "$NOTIFY" = true ]; then
  if [ $VERIFICATION_STATUS -eq 0 ]; then
    echo "Database integrity verification completed successfully. See report: $HTML_REPORT" | mail -s "Database Integrity Check Passed" admin@example.com
  else
    echo "Database integrity verification found issues. See report: $HTML_REPORT" | mail -s "ALERT: Database Integrity Issues Found" admin@example.com
  fi
  log "Notification email sent"
fi

# Add to crontab if requested
if [ "$SCHEDULE" = true ]; then
  SCRIPT_PATH=$(readlink -f "$0")
  (crontab -l 2>/dev/null || echo "") | grep -v "$SCRIPT_PATH" | { cat; echo "0 2 * * * $SCRIPT_PATH --notify"; } | crontab -
  log "Added to crontab - Will run daily at 2:00 AM with notifications"
fi

log "===================================================="
log "Verification process completed at: $(date)"
log "Exit status: $VERIFICATION_STATUS"
log "===================================================="

exit $VERIFICATION_STATUS