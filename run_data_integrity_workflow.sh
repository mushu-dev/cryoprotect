#!/bin/bash
# run_data_integrity_workflow.sh
#
# Master script for the complete database integrity verification 
# and repair workflow. This script runs verification, generates
# reports, and optionally applies fixes.
#
# Usage:
#   ./run_data_integrity_workflow.sh [--auto-fix] [--backup] [--notify=email@example.com]

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Parse command line arguments
AUTO_FIX=false
BACKUP=false
NOTIFY=""
VERBOSE=""

for arg in "$@"; do
  case $arg in
    --auto-fix)
      AUTO_FIX=true
      shift
      ;;
    --backup)
      BACKUP=true
      shift
      ;;
    --notify=*)
      NOTIFY="${arg#*=}"
      shift
      ;;
    --verbose)
      VERBOSE="--verbose"
      shift
      ;;
    *)
      # Unknown option
      ;;
  esac
done

# Print header
echo -e "${BLUE}=================================================="
echo "CryoProtect Database Integrity Workflow"
echo -e "==================================================${NC}"
echo "Started at: $(date)"
echo ""

# Create reports directory if it doesn't exist
mkdir -p reports

# Set file paths with timestamps
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
JSON_REPORT="reports/db_integrity_${TIMESTAMP}.json"
HTML_REPORT="reports/db_integrity_${TIMESTAMP}.html"
MD_REPORT="reports/db_integrity_${TIMESTAMP}.md"
LOG_FILE="reports/db_integrity_${TIMESTAMP}.log"

# Function to notify
send_notification() {
  if [ -n "$NOTIFY" ]; then
    if command -v mail &> /dev/null; then
      echo -e "$1" | mail -s "Database Integrity Report: $2" "$NOTIFY"
      echo "Notification sent to $NOTIFY"
    else
      echo "Warning: 'mail' command not found. Cannot send notification."
    fi
  fi
}

# Function to log
log_message() {
  echo "$(date +"%Y-%m-%d %H:%M:%S") - $1" | tee -a "$LOG_FILE"
}

# Stage 1: Run verification
echo -e "${BLUE}[ Stage 1/3 ]${NC} Running database integrity verification..."
log_message "Starting database integrity verification"

./run_database_integrity_check_fixed.sh $VERBOSE --report="$JSON_REPORT"
VERIFICATION_STATUS=$?

if [ $VERIFICATION_STATUS -eq 0 ]; then
  echo -e "${GREEN}✅ Verification passed${NC}"
  log_message "Verification passed"
  
  if [ -n "$NOTIFY" ]; then
    send_notification "Database integrity verification passed. No issues found." "Successful"
  fi
  
  # Create symlinks to latest reports
  ln -sf "$JSON_REPORT" "reports/db_integrity_latest.json"
  ln -sf "$HTML_REPORT" "reports/db_integrity_latest.html"
  ln -sf "$MD_REPORT" "reports/db_integrity_latest.md"
  ln -sf "$LOG_FILE" "reports/db_integrity_latest.log"
  
  echo -e "${GREEN}Verification complete. All reports saved in 'reports/' directory.${NC}"
  exit 0
else
  echo -e "${YELLOW}⚠️ Verification found issues${NC}"
  log_message "Verification found issues"
  
  # Count issues by severity
  ERROR_COUNT=$(grep -o '"severity": "error"' "$JSON_REPORT" | wc -l)
  WARNING_COUNT=$(grep -o '"severity": "warning"' "$JSON_REPORT" | wc -l)
  
  echo -e "${YELLOW}Found $ERROR_COUNT errors and $WARNING_COUNT warnings${NC}"
  log_message "Found $ERROR_COUNT errors and $WARNING_COUNT warnings"
fi

# Stage 2: Generate reports
echo -e "\n${BLUE}[ Stage 2/3 ]${NC} Generating detailed reports..."
log_message "Generating detailed reports"

if [ -f "generate_database_integrity_report.py" ]; then
  python generate_database_integrity_report.py --input="$JSON_REPORT" --output="$HTML_REPORT"
  echo -e "${GREEN}✅ HTML report generated:${NC} $HTML_REPORT"
  log_message "HTML report generated: $HTML_REPORT"
else
  echo -e "${YELLOW}⚠️ HTML report generator not found${NC}"
  log_message "HTML report generator not found"
fi

# Create symlinks to latest reports
ln -sf "$JSON_REPORT" "reports/db_integrity_latest.json"
ln -sf "$HTML_REPORT" "reports/db_integrity_latest.html"
ln -sf "$MD_REPORT" "reports/db_integrity_latest.md"
ln -sf "$LOG_FILE" "reports/db_integrity_latest.log"

# Stage 3: Fix issues if needed and requested
if [ $VERIFICATION_STATUS -ne 0 ] && [ "$AUTO_FIX" = true ]; then
  echo -e "\n${BLUE}[ Stage 3/3 ]${NC} Fixing database integrity issues..."
  log_message "Applying database fixes"
  
  FIX_OPTS=""
  if [ "$BACKUP" = true ]; then
    FIX_OPTS="--backup"
  fi
  
  ./fix_database_integrity.sh $FIX_OPTS --apply
  FIX_STATUS=$?
  
  if [ $FIX_STATUS -eq 0 ]; then
    echo -e "${GREEN}✅ Database fixes applied successfully${NC}"
    log_message "Database fixes applied successfully"
    
    # Run verification again to check results
    echo -e "\n${BLUE}[ Verification ]${NC} Verifying fixes..."
    log_message "Verifying fixes"
    
    POST_FIX_REPORT="reports/db_integrity_post_fix_${TIMESTAMP}.json"
    ./run_database_integrity_check_fixed.sh --report="$POST_FIX_REPORT"
    POST_VERIFICATION_STATUS=$?
    
    if [ $POST_VERIFICATION_STATUS -eq 0 ]; then
      echo -e "${GREEN}✅ Post-fix verification passed${NC}"
      log_message "Post-fix verification passed"
      
      if [ -n "$NOTIFY" ]; then
        send_notification "Database integrity issues fixed successfully." "Fixed Successfully"
      fi
    else
      echo -e "${RED}❌ Some issues remain after fixes${NC}"
      log_message "Some issues remain after fixes"
      
      if [ -n "$NOTIFY" ]; then
        send_notification "Database integrity issues partially fixed. Manual review needed." "Partial Fix"
      fi
    fi
  else
    echo -e "${RED}❌ Error applying database fixes${NC}"
    log_message "Error applying database fixes"
    
    if [ -n "$NOTIFY" ]; then
      send_notification "Error applying database integrity fixes. Manual intervention required." "Fix Failed"
    fi
  fi
else
  if [ $VERIFICATION_STATUS -ne 0 ]; then
    echo -e "\n${BLUE}[ Stage 3/3 ]${NC} Skipping fix stage (--auto-fix not provided)"
    log_message "Skipping fix stage (--auto-fix not provided)"
    
    if [ -n "$NOTIFY" ]; then
      send_notification "Database integrity issues found. Run fix script to resolve." "Issues Found"
    fi
    
    echo -e "${YELLOW}To fix issues, run:${NC}"
    echo -e "  ./fix_database_integrity.sh --backup --apply"
  fi
fi

echo -e "\n${BLUE}=================================================="
echo "Workflow completed at: $(date)"
echo -e "==================================================${NC}"
echo "Log file: $LOG_FILE"
echo "All reports saved in 'reports/' directory"

exit $VERIFICATION_STATUS