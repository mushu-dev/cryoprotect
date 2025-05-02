#!/bin/bash
# scan_dependency_check.sh - OWASP Dependency-Check scanner
# 
# This script scans project dependencies for known vulnerabilities using OWASP Dependency-Check.
# It can be run as a standalone script or integrated into CI/CD pipelines.
#
# Usage: ./scan_dependency_check.sh [options]
#
# Options:
#   --path PATH           Path to scan (default: current directory)
#   --format FORMAT       Output format (HTML, XML, CSV, JSON, JUNIT, SARIF) (default: HTML)
#   --output DIR          Output directory (default: ./dependency-check-reports)
#   --name NAME           Project name (default: CryoProtect)
#   --exit-on-critical    Exit with code 1 if critical vulnerabilities found (CVSS >= 7.0)
#   --nvd-api-key KEY     NVD API key for better performance
#   --help                Show this help message and exit

set -e

# Default values
SCAN_PATH="."
FORMAT="HTML"
OUTPUT_DIR="./dependency-check-reports"
PROJECT_NAME="CryoProtect"
EXIT_ON_CRITICAL=false
NVD_API_KEY=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --path)
      SCAN_PATH="$2"
      shift 2
      ;;
    --format)
      FORMAT="$2"
      shift 2
      ;;
    --output)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --name)
      PROJECT_NAME="$2"
      shift 2
      ;;
    --exit-on-critical)
      EXIT_ON_CRITICAL=true
      shift
      ;;
    --nvd-api-key)
      NVD_API_KEY="$2"
      shift 2
      ;;
    --help)
      echo "Usage: ./scan_dependency_check.sh [options]"
      echo ""
      echo "Options:"
      echo "  --path PATH           Path to scan (default: current directory)"
      echo "  --format FORMAT       Output format (HTML, XML, CSV, JSON, JUNIT, SARIF) (default: HTML)"
      echo "  --output DIR          Output directory (default: ./dependency-check-reports)"
      echo "  --name NAME           Project name (default: CryoProtect)"
      echo "  --exit-on-critical    Exit with code 1 if critical vulnerabilities found (CVSS >= 7.0)"
      echo "  --nvd-api-key KEY     NVD API key for better performance"
      echo "  --help                Show this help message and exit"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Run './scan_dependency_check.sh --help' for usage information"
      exit 1
      ;;
  esac
done

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Generate timestamp for reports
TIMESTAMP=$(date -u +"%Y%m%d%H%M%S")
REPORT_NAME="${PROJECT_NAME}-${TIMESTAMP}"

# Check if Docker is available
if command -v docker &> /dev/null; then
  echo "Using Docker to run OWASP Dependency-Check..."
  
  # Prepare Docker command
  DOCKER_CMD="docker run --rm"
  DOCKER_CMD+=" -v \"$(pwd)/${SCAN_PATH}:/src\""
  DOCKER_CMD+=" -v \"$(pwd)/${OUTPUT_DIR}:/report\""
  DOCKER_CMD+=" -v \"$(pwd)/dependency-check-data:/usr/share/dependency-check/data\""
  
  # Add NVD API key if provided
  if [ -n "$NVD_API_KEY" ]; then
    DOCKER_CMD+=" -e \"NVD_API_KEY=${NVD_API_KEY}\""
  fi
  
  # Complete the command
  DOCKER_CMD+=" owasp/dependency-check:latest"
  DOCKER_CMD+=" --scan /src"
  DOCKER_CMD+=" --format $FORMAT"
  DOCKER_CMD+=" --project \"$PROJECT_NAME\""
  DOCKER_CMD+=" --out /report"
  DOCKER_CMD+=" --enableExperimental"
  
  # Run the scan
  echo "Running OWASP Dependency-Check scan on $SCAN_PATH..."
  eval $DOCKER_CMD
  
  # Rename the report files with timestamp
  if [ -f "$OUTPUT_DIR/dependency-check-report.html" ]; then
    mv "$OUTPUT_DIR/dependency-check-report.html" "$OUTPUT_DIR/${REPORT_NAME}.html"
  fi
  
  if [ -f "$OUTPUT_DIR/dependency-check-report.json" ]; then
    mv "$OUTPUT_DIR/dependency-check-report.json" "$OUTPUT_DIR/${REPORT_NAME}.json"
  fi
  
  if [ -f "$OUTPUT_DIR/dependency-check-report.xml" ]; then
    mv "$OUTPUT_DIR/dependency-check-report.xml" "$OUTPUT_DIR/${REPORT_NAME}.xml"
  fi
  
  if [ -f "$OUTPUT_DIR/dependency-check-report.csv" ]; then
    mv "$OUTPUT_DIR/dependency-check-report.csv" "$OUTPUT_DIR/${REPORT_NAME}.csv"
  fi
  
  if [ -f "$OUTPUT_DIR/dependency-check-junit.xml" ]; then
    mv "$OUTPUT_DIR/dependency-check-junit.xml" "$OUTPUT_DIR/${REPORT_NAME}-junit.xml"
  fi
  
  if [ -f "$OUTPUT_DIR/dependency-check-sarif.json" ]; then
    mv "$OUTPUT_DIR/dependency-check-sarif.json" "$OUTPUT_DIR/${REPORT_NAME}-sarif.json"
  fi
  
  echo "Scan completed. Reports saved to $OUTPUT_DIR"
  
  # Check for critical vulnerabilities if requested
  if [ "$EXIT_ON_CRITICAL" = true ] && [ -f "$OUTPUT_DIR/${REPORT_NAME}.json" ]; then
    echo "Checking for critical vulnerabilities (CVSS >= 7.0)..."
    
    # Use jq if available, otherwise grep
    if command -v jq &> /dev/null; then
      CRITICAL_COUNT=$(jq -r '.dependencies[] | select(.vulnerabilities != null) | .vulnerabilities[] | select(.cvssv3 != null and .cvssv3.baseScore >= 7.0) | .name' "$OUTPUT_DIR/${REPORT_NAME}.json" | wc -l)
    else
      # This is a rough approximation without jq
      CRITICAL_COUNT=$(grep -c "\"baseScore\": [7-9]" "$OUTPUT_DIR/${REPORT_NAME}.json" || true)
    fi
    
    if [ "$CRITICAL_COUNT" -gt 0 ]; then
      echo "Found $CRITICAL_COUNT critical vulnerabilities (CVSS >= 7.0)!"
      echo "See the report for details: $OUTPUT_DIR/${REPORT_NAME}.html"
      exit 1
    else
      echo "No critical vulnerabilities found."
    fi
  fi
  
else
  # Check if dependency-check script is installed
  if command -v dependency-check.sh &> /dev/null; then
    DEPENDENCY_CHECK_CMD="dependency-check.sh"
  elif [ -f "./dependency-check/bin/dependency-check.sh" ]; then
    DEPENDENCY_CHECK_CMD="./dependency-check/bin/dependency-check.sh"
  else
    echo "Error: OWASP Dependency-Check not found and Docker is not available."
    echo "Please install OWASP Dependency-Check or Docker."
    echo "Installation instructions: https://jeremylong.github.io/DependencyCheck/dependency-check-cli/index.html"
    exit 1
  fi
  
  # Prepare command
  CMD="$DEPENDENCY_CHECK_CMD"
  CMD+=" --scan \"$SCAN_PATH\""
  CMD+=" --format $FORMAT"
  CMD+=" --project \"$PROJECT_NAME\""
  CMD+=" --out \"$OUTPUT_DIR\""
  CMD+=" --enableExperimental"
  
  # Add NVD API key if provided
  if [ -n "$NVD_API_KEY" ]; then
    CMD+=" --nvdApiKey \"$NVD_API_KEY\""
  fi
  
  # Run the scan
  echo "Running OWASP Dependency-Check scan on $SCAN_PATH..."
  eval $CMD
  
  # Rename the report files with timestamp
  if [ -f "$OUTPUT_DIR/dependency-check-report.html" ]; then
    mv "$OUTPUT_DIR/dependency-check-report.html" "$OUTPUT_DIR/${REPORT_NAME}.html"
  fi
  
  if [ -f "$OUTPUT_DIR/dependency-check-report.json" ]; then
    mv "$OUTPUT_DIR/dependency-check-report.json" "$OUTPUT_DIR/${REPORT_NAME}.json"
  fi
  
  if [ -f "$OUTPUT_DIR/dependency-check-report.xml" ]; then
    mv "$OUTPUT_DIR/dependency-check-report.xml" "$OUTPUT_DIR/${REPORT_NAME}.xml"
  fi
  
  if [ -f "$OUTPUT_DIR/dependency-check-report.csv" ]; then
    mv "$OUTPUT_DIR/dependency-check-report.csv" "$OUTPUT_DIR/${REPORT_NAME}.csv"
  fi
  
  if [ -f "$OUTPUT_DIR/dependency-check-junit.xml" ]; then
    mv "$OUTPUT_DIR/dependency-check-junit.xml" "$OUTPUT_DIR/${REPORT_NAME}-junit.xml"
  fi
  
  if [ -f "$OUTPUT_DIR/dependency-check-sarif.json" ]; then
    mv "$OUTPUT_DIR/dependency-check-sarif.json" "$OUTPUT_DIR/${REPORT_NAME}-sarif.json"
  fi
  
  echo "Scan completed. Reports saved to $OUTPUT_DIR"
  
  # Check for critical vulnerabilities if requested
  if [ "$EXIT_ON_CRITICAL" = true ] && [ -f "$OUTPUT_DIR/${REPORT_NAME}.json" ]; then
    echo "Checking for critical vulnerabilities (CVSS >= 7.0)..."
    
    # Use jq if available, otherwise grep
    if command -v jq &> /dev/null; then
      CRITICAL_COUNT=$(jq -r '.dependencies[] | select(.vulnerabilities != null) | .vulnerabilities[] | select(.cvssv3 != null and .cvssv3.baseScore >= 7.0) | .name' "$OUTPUT_DIR/${REPORT_NAME}.json" | wc -l)
    else
      # This is a rough approximation without jq
      CRITICAL_COUNT=$(grep -c "\"baseScore\": [7-9]" "$OUTPUT_DIR/${REPORT_NAME}.json" || true)
    fi
    
    if [ "$CRITICAL_COUNT" -gt 0 ]; then
      echo "Found $CRITICAL_COUNT critical vulnerabilities (CVSS >= 7.0)!"
      echo "See the report for details: $OUTPUT_DIR/${REPORT_NAME}.html"
      exit 1
    else
      echo "No critical vulnerabilities found."
    fi
  fi
fi

# Generate summary file
SUMMARY_FILE="$OUTPUT_DIR/${REPORT_NAME}-summary.json"

echo "{" > "$SUMMARY_FILE"
echo "  \"scanner\": \"owasp-dependency-check\"," >> "$SUMMARY_FILE"
echo "  \"timestamp\": \"$(date -u +"%Y-%m-%dT%H:%M:%SZ")\"," >> "$SUMMARY_FILE"
echo "  \"project\": \"$PROJECT_NAME\"," >> "$SUMMARY_FILE"
echo "  \"scan_path\": \"$SCAN_PATH\"," >> "$SUMMARY_FILE"
echo "  \"report_files\": {" >> "$SUMMARY_FILE"

# Add report files to summary
if [ -f "$OUTPUT_DIR/${REPORT_NAME}.html" ]; then
  echo "    \"html\": \"${REPORT_NAME}.html\"," >> "$SUMMARY_FILE"
fi

if [ -f "$OUTPUT_DIR/${REPORT_NAME}.json" ]; then
  echo "    \"json\": \"${REPORT_NAME}.json\"," >> "$SUMMARY_FILE"
fi

if [ -f "$OUTPUT_DIR/${REPORT_NAME}.xml" ]; then
  echo "    \"xml\": \"${REPORT_NAME}.xml\"," >> "$SUMMARY_FILE"
fi

if [ -f "$OUTPUT_DIR/${REPORT_NAME}.csv" ]; then
  echo "    \"csv\": \"${REPORT_NAME}.csv\"," >> "$SUMMARY_FILE"
fi

if [ -f "$OUTPUT_DIR/${REPORT_NAME}-junit.xml" ]; then
  echo "    \"junit\": \"${REPORT_NAME}-junit.xml\"," >> "$SUMMARY_FILE"
fi

if [ -f "$OUTPUT_DIR/${REPORT_NAME}-sarif.json" ]; then
  echo "    \"sarif\": \"${REPORT_NAME}-sarif.json\"" >> "$SUMMARY_FILE"
else
  # Remove trailing comma from last item
  sed -i 's/,$//' "$SUMMARY_FILE"
fi

echo "  }" >> "$SUMMARY_FILE"
echo "}" >> "$SUMMARY_FILE"

echo "Summary saved to $SUMMARY_FILE"

exit 0