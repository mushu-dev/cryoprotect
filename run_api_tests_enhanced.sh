#!/bin/bash
#
# Enhanced API Endpoint Testing Script
# This script runs the enhanced API endpoint tests with configurable options

# Default settings
BASE_URL="http://localhost:5000"
ENVIRONMENT="local"
PARALLEL=false
TIMEOUT=30
RETRIES=3
VERBOSE=false
GENERATE_REPORT=true
MOCK_SERVER=false  # Use mock server instead of real server

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --url=*)
      BASE_URL="${1#*=}"
      shift
      ;;
    --env=*)
      ENVIRONMENT="${1#*=}"
      shift
      ;;
    --parallel)
      PARALLEL=true
      shift
      ;;
    --timeout=*)
      TIMEOUT="${1#*=}"
      shift
      ;;
    --retries=*)
      RETRIES="${1#*=}"
      shift
      ;;
    --verbose)
      VERBOSE=true
      shift
      ;;
    --no-report)
      GENERATE_REPORT=false
      shift
      ;;
    --mock)
      MOCK_SERVER=true
      shift
      ;;
    *)
      echo "Unknown parameter: $1"
      echo "Usage: $0 [--url=URL] [--env=local|staging|production] [--parallel] [--timeout=SECONDS] [--retries=COUNT] [--verbose] [--no-report] [--mock]"
      exit 1
      ;;
  esac
done

# Check if we're using the mock server
if [ "$MOCK_SERVER" = true ]; then
  echo "Using mock server for testing..."
else
  # Ensure the target server is running
  echo "Checking if API server is running at $BASE_URL..."
  curl -s --head --max-time 5 "$BASE_URL/health" >/dev/null
  if [ $? -ne 0 ]; then
    echo "Warning: API server at $BASE_URL is not responding. Tests may fail."
  else
    echo "API server at $BASE_URL is responding."
  fi
fi

# Build command arguments
ARGS=""
ARGS="$ARGS --base-url $BASE_URL"
ARGS="$ARGS --environment $ENVIRONMENT"
ARGS="$ARGS --timeout $TIMEOUT"
ARGS="$ARGS --retries $RETRIES"

if [ "$PARALLEL" = true ]; then
  ARGS="$ARGS --parallel"
fi

if [ "$VERBOSE" = true ]; then
  ARGS="$ARGS --verbose"
fi

if [ "$MOCK_SERVER" = true ]; then
  ARGS="$ARGS --mock"
fi

# Create timestamp for report file
TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
RESULTS_PATH="memory-bank/verification_results_$TIMESTAMP.json"
ARGS="$ARGS --output $RESULTS_PATH"

if [ "$GENERATE_REPORT" = true ]; then
  REPORT_PATH="memory-bank/API_VERIFICATION_REPORT_$TIMESTAMP.md"
  ARGS="$ARGS --report $REPORT_PATH"
fi

# Ensure memory-bank directory exists
mkdir -p memory-bank

# Run the tests
echo "Starting API endpoint tests with the following settings:"
echo "- Base URL: $BASE_URL"
echo "- Environment: $ENVIRONMENT"
echo "- Parallel: $PARALLEL"
echo "- Timeout: $TIMEOUT seconds"
echo "- Retries: $RETRIES"
echo "- Verbose: $VERBOSE"
echo "- Results: $RESULTS_PATH"
if [ "$GENERATE_REPORT" = true ]; then
  echo "- Report: $REPORT_PATH"
fi

echo ""
python3 test_all_api_endpoints_enhanced.py $ARGS

# Check exit code
if [ $? -eq 0 ]; then
  echo "API testing completed successfully."
  
  if [ "$GENERATE_REPORT" = true ]; then
    echo "Detailed report available at: $REPORT_PATH"
  fi
  
  exit 0
else
  echo "API testing completed with errors."
  exit 1
fi