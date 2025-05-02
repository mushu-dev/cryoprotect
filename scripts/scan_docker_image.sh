#!/bin/bash
# scan_docker_image.sh - Security scanning for Docker images
# 
# This script scans Docker images for vulnerabilities using Trivy
# and can fail the build if critical vulnerabilities are found.
#
# Usage: ./scan_docker_image.sh [options] IMAGE_NAME
#
# Options:
#   --format FORMAT       Output format (table, json, sarif, cyclonedx) [default: table]
#   --output FILE         Output file name [default: trivy-results.{format}]
#   --severity SEVERITY   Severities to scan for (comma-separated) [default: CRITICAL,HIGH,MEDIUM]
#   --exit-code CODE      Exit code when vulnerabilities are found [default: 0]
#   --fail-on SEVERITY    Fail on specific severity (CRITICAL, HIGH, MEDIUM, LOW, UNKNOWN) [default: CRITICAL]
#   --help                Show this help message and exit

set -e

# Default values
FORMAT="table"
OUTPUT=""
SEVERITY="CRITICAL,HIGH,MEDIUM"
EXIT_CODE=0
FAIL_ON="CRITICAL"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --format)
      FORMAT="$2"
      shift 2
      ;;
    --output)
      OUTPUT="$2"
      shift 2
      ;;
    --severity)
      SEVERITY="$2"
      shift 2
      ;;
    --exit-code)
      EXIT_CODE="$2"
      shift 2
      ;;
    --fail-on)
      FAIL_ON="$2"
      shift 2
      ;;
    --help)
      echo "Usage: ./scan_docker_image.sh [options] IMAGE_NAME"
      echo ""
      echo "Options:"
      echo "  --format FORMAT       Output format (table, json, sarif, cyclonedx) [default: table]"
      echo "  --output FILE         Output file name [default: trivy-results.{format}]"
      echo "  --severity SEVERITY   Severities to scan for (comma-separated) [default: CRITICAL,HIGH,MEDIUM]"
      echo "  --exit-code CODE      Exit code when vulnerabilities are found [default: 0]"
      echo "  --fail-on SEVERITY    Fail on specific severity (CRITICAL, HIGH, MEDIUM, LOW, UNKNOWN) [default: CRITICAL]"
      echo "  --help                Show this help message and exit"
      exit 0
      ;;
    *)
      IMAGE_NAME="$1"
      shift
      ;;
  esac
done

# Check if image name is provided
if [ -z "$IMAGE_NAME" ]; then
  echo "Error: IMAGE_NAME is required"
  echo "Run './scan_docker_image.sh --help' for usage information"
  exit 1
fi

# Set default output file if not provided
if [ -z "$OUTPUT" ]; then
  OUTPUT="trivy-results.${FORMAT}"
fi

echo "Scanning Docker image: $IMAGE_NAME"
echo "Scanning for vulnerabilities with severity: $SEVERITY"
echo "Output format: $FORMAT"
echo "Output file: $OUTPUT"

# Check if Trivy is installed
if ! command -v trivy &> /dev/null; then
  echo "Error: Trivy is not installed"
  echo "Please install Trivy: https://aquasecurity.github.io/trivy/latest/getting-started/installation/"
  exit 1
fi

# Run Trivy scan
trivy image --format "$FORMAT" --output "$OUTPUT" --severity "$SEVERITY" "$IMAGE_NAME"

# Check for vulnerabilities based on FAIL_ON severity
if [ "$FAIL_ON" != "NONE" ]; then
  if [ "$FORMAT" = "json" ]; then
    # For JSON format, use jq if available, otherwise grep
    if command -v jq &> /dev/null; then
      VULN_COUNT=$(jq -r '.Results[] | select(.Vulnerabilities != null) | .Vulnerabilities[] | select(.Severity == "'"$FAIL_ON"'") | .VulnerabilityID' "$OUTPUT" | wc -l)
    else
      VULN_COUNT=$(grep -c "\"Severity\": \"$FAIL_ON\"" "$OUTPUT" || true)
    fi
  else
    # For other formats, use grep
    VULN_COUNT=$(grep -c "$FAIL_ON" "$OUTPUT" || true)
  fi

  if [ "$VULN_COUNT" -gt 0 ]; then
    echo "Found $VULN_COUNT $FAIL_ON severity vulnerabilities!"
    if [ "$EXIT_CODE" -ne 0 ]; then
      echo "Failing build as requested (exit code $EXIT_CODE)"
      exit "$EXIT_CODE"
    else
      echo "Warning: Vulnerabilities found but continuing build (exit code 0)"
    fi
  else
    echo "No $FAIL_ON severity vulnerabilities found."
  fi
fi

echo "Scan completed successfully."