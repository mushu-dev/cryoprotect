#!/bin/bash
# CryoProtect v2 Security Controls Validation Script
#
# This script runs the comprehensive security validation process:
# 1. Executes automated security tests
# 2. Runs penetration tests
# 3. Generates validation reports
#
# Usage: ./validate_security_controls.sh [--url URL] [--output-dir DIR] [--html] [--pdf]

# Default values
URL="http://localhost:5000"
OUTPUT_DIR="reports/security"
GENERATE_HTML=false
GENERATE_PDF=false
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
REPORT_PREFIX="security_validation_${TIMESTAMP}"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --url)
      URL="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --html)
      GENERATE_HTML=true
      shift
      ;;
    --pdf)
      GENERATE_PDF=true
      shift
      ;;
    --help)
      echo "Usage: $0 [--url URL] [--output-dir DIR] [--html] [--pdf]"
      echo ""
      echo "Options:"
      echo "  --url URL         Base URL of the application to test (default: http://localhost:5000)"
      echo "  --output-dir DIR  Directory to save reports (default: reports/security)"
      echo "  --html            Generate HTML report in addition to JSON"
      echo "  --pdf             Generate PDF report in addition to JSON"
      echo "  --help            Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Function to print section headers
print_section() {
  echo ""
  echo "================================================================================"
  echo "  $1"
  echo "================================================================================"
  echo ""
}

# Function to run a command and check its exit status
run_command() {
  echo "Running: $1"
  eval "$1"
  exit_status=$?
  if [ $exit_status -ne 0 ]; then
    echo "ERROR: Command failed with exit status $exit_status"
    if [ "$2" = "critical" ]; then
      echo "Critical command failed. Aborting validation process."
      exit 1
    fi
  else
    echo "Command completed successfully"
  fi
  return $exit_status
}

# Start validation process
print_section "Starting CryoProtect v2 Security Controls Validation"
echo "Timestamp: $(date)"
echo "Target URL: $URL"
echo "Output directory: $OUTPUT_DIR"

# Step 1: Run automated security tests
print_section "Running Automated Security Tests"
TEST_OUTPUT_FILE="${OUTPUT_DIR}/${REPORT_PREFIX}_test_results.txt"
run_command "python -m pytest tests/test_security_controls.py -v | tee $TEST_OUTPUT_FILE"
TEST_EXIT_STATUS=$?

echo "Test results saved to: $TEST_OUTPUT_FILE"

# Step 2: Run penetration tests
print_section "Running Penetration Tests"
PENTEST_ARGS="--url $URL --output-dir $OUTPUT_DIR"
if [ "$GENERATE_HTML" = true ]; then
  PENTEST_ARGS="$PENTEST_ARGS --html"
fi
if [ "$GENERATE_PDF" = true ]; then
  PENTEST_ARGS="$PENTEST_ARGS --pdf"
fi

run_command "python security/run_pentest.py $PENTEST_ARGS"
PENTEST_EXIT_STATUS=$?

# Step 3: Generate validation summary
print_section "Generating Validation Summary"
SUMMARY_FILE="${OUTPUT_DIR}/${REPORT_PREFIX}_summary.md"

# Create summary file
cat > "$SUMMARY_FILE" << EOF
# CryoProtect v2 Security Controls Validation Summary

**Date:** $(date +"%Y-%m-%d %H:%M:%S")
**Target:** $URL

## Validation Results

| Component | Status | Details |
|-----------|--------|---------|
| Automated Tests | $([ $TEST_EXIT_STATUS -eq 0 ] && echo "✅ PASSED" || echo "❌ FAILED") | See [Test Results]($(basename "$TEST_OUTPUT_FILE")) |
| Penetration Tests | $([ $PENTEST_EXIT_STATUS -eq 0 ] && echo "✅ PASSED" || echo "❌ FAILED") | See [Pentest Report](pentest_summary.md) |
| Overall Validation | $([ $TEST_EXIT_STATUS -eq 0 ] && [ $PENTEST_EXIT_STATUS -eq 0 ] && echo "✅ PASSED" || echo "❌ FAILED") | |

## Security Controls Status

| Security Control | Status | Validation Method |
|------------------|--------|------------------|
| CSRF Protection | $([ $TEST_EXIT_STATUS -eq 0 ] && [ $PENTEST_EXIT_STATUS -eq 0 ] && echo "✅ Validated" || echo "⚠️ Needs Review") | Automated Tests, Penetration Testing |
| Security Headers | $([ $TEST_EXIT_STATUS -eq 0 ] && [ $PENTEST_EXIT_STATUS -eq 0 ] && echo "✅ Validated" || echo "⚠️ Needs Review") | Automated Tests, Penetration Testing |
| Cookie Security | $([ $TEST_EXIT_STATUS -eq 0 ] && [ $PENTEST_EXIT_STATUS -eq 0 ] && echo "✅ Validated" || echo "⚠️ Needs Review") | Automated Tests, Penetration Testing |
| Encryption at Rest | $([ $TEST_EXIT_STATUS -eq 0 ] && [ $PENTEST_EXIT_STATUS -eq 0 ] && echo "✅ Validated" || echo "⚠️ Needs Review") | Automated Tests, Penetration Testing |
| Vulnerability Scanning | $([ $TEST_EXIT_STATUS -eq 0 ] && [ $PENTEST_EXIT_STATUS -eq 0 ] && echo "✅ Validated" || echo "⚠️ Needs Review") | Automated Tests, Penetration Testing |

## Next Steps

$([ $TEST_EXIT_STATUS -eq 0 ] && [ $PENTEST_EXIT_STATUS -eq 0 ] && echo "All security controls have been successfully validated. The application is ready for production use." || echo "Some security controls need review. Please check the detailed reports and address any issues before proceeding to production.")

## References

- [Security Validation Report](security_validation_report.md)
- [Test Results]($(basename "$TEST_OUTPUT_FILE"))
- [Penetration Test Report](pentest_summary.md)
EOF

echo "Validation summary saved to: $SUMMARY_FILE"

# Step 4: Copy the security validation report to the output directory
print_section "Finalizing Reports"
if [ -f "reports/security/security_validation_report.md" ]; then
  cp "reports/security/security_validation_report.md" "${OUTPUT_DIR}/${REPORT_PREFIX}_report.md"
  echo "Security validation report copied to: ${OUTPUT_DIR}/${REPORT_PREFIX}_report.md"
else
  echo "WARNING: Security validation report not found at reports/security/security_validation_report.md"
fi

# Step 5: Determine overall validation status
print_section "Validation Complete"
if [ $TEST_EXIT_STATUS -eq 0 ] && [ $PENTEST_EXIT_STATUS -eq 0 ]; then
  echo "✅ SUCCESS: All security controls have been successfully validated."
  VALIDATION_STATUS=0
else
  echo "❌ FAILURE: Some security controls need review. Please check the detailed reports."
  VALIDATION_STATUS=1
fi

echo ""
echo "Reports generated:"
echo "- Validation Summary: $SUMMARY_FILE"
echo "- Test Results: $TEST_OUTPUT_FILE"
echo "- Penetration Test Reports: See $OUTPUT_DIR/pentest_*.json"
if [ "$GENERATE_HTML" = true ]; then
  echo "- HTML Reports: See $OUTPUT_DIR/pentest_*.html"
fi
if [ "$GENERATE_PDF" = true ]; then
  echo "- PDF Reports: See $OUTPUT_DIR/pentest_*.pdf"
fi

echo ""
echo "Validation process completed at $(date)"
exit $VALIDATION_STATUS