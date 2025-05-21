#!/bin/bash
# run_connection_pool_stress_test.sh
# Script to run comprehensive connection pool stress tests

set -e

# Ensure necessary directories exist
mkdir -p logs
mkdir -p reports

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== Connection Pool Stress Test Runner ===${NC}"
echo "This script will run comprehensive stress tests on the database connection pool"
echo "to help optimize performance and reliability."
echo ""

# Check if .env file exists and load it
if [ -f .env ]; then
    echo -e "${YELLOW}Loading environment variables from .env file${NC}"
    export $(grep -v '^#' .env | xargs)
fi

# Check for required environment variables
if [ -z "$SUPABASE_DB_HOST" ] || [ -z "$SUPABASE_DB_USER" ] || [ -z "$SUPABASE_DB_PASSWORD" ]; then
    echo -e "${RED}ERROR: Missing required environment variables${NC}"
    echo "Please ensure the following variables are set either in your environment or .env file:"
    echo "  SUPABASE_DB_HOST"
    echo "  SUPABASE_DB_USER"
    echo "  SUPABASE_DB_PASSWORD"
    exit 1
fi

# Parse command line arguments
FULL_SUITE=0
WORKERS=20
DURATION=30
SIMULATE_ERRORS=0
USE_TRANSACTIONS=1
OPTIMIZED=1
MIN_CONN=3
MAX_CONN=20

# Parse arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --full-suite)
            FULL_SUITE=1
            shift
            ;;
        --workers)
            WORKERS="$2"
            shift
            shift
            ;;
        --duration)
            DURATION="$2"
            shift
            shift
            ;;
        --simulate-errors)
            SIMULATE_ERRORS=1
            shift
            ;;
        --no-transactions)
            USE_TRANSACTIONS=0
            shift
            ;;
        --standard)
            OPTIMIZED=0
            shift
            ;;
        --min-connections)
            MIN_CONN="$2"
            shift
            shift
            ;;
        --max-connections)
            MAX_CONN="$2"
            shift
            shift
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --full-suite         Run full test suite with multiple configurations"
            echo "  --workers N          Number of worker threads for load test (default: 20)"
            echo "  --duration N         Duration of test in seconds (default: 30)"
            echo "  --simulate-errors    Simulate random errors during testing"
            echo "  --no-transactions    Don't use transactions for queries"
            echo "  --standard           Use standard connection pool instead of optimized"
            echo "  --min-connections N  Minimum connections in pool (default: 3)"
            echo "  --max-connections N  Maximum connections in pool (default: 20)"
            echo "  --help               Show this help message"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $key${NC}"
            echo "Use --help to see available options"
            exit 1
            ;;
    esac
done

# Print test configuration
echo -e "${BLUE}Test Configuration:${NC}"
if [ $FULL_SUITE -eq 1 ]; then
    echo "- Running full test suite with multiple configurations"
else
    echo "- Workers: $WORKERS"
    echo "- Duration: $DURATION seconds"
    echo "- Pool Type: $([ $OPTIMIZED -eq 1 ] && echo "Optimized" || echo "Standard")"
    echo "- Min Connections: $MIN_CONN"
    echo "- Max Connections: $MAX_CONN"
    echo "- Simulate Errors: $([ $SIMULATE_ERRORS -eq 1 ] && echo "Yes" || echo "No")"
    echo "- Use Transactions: $([ $USE_TRANSACTIONS -eq 1 ] && echo "Yes" || echo "No")"
fi

echo ""
echo -e "${YELLOW}Starting stress tests...${NC}"
echo "This may take some time. Detailed logs will be saved to the logs directory."
echo ""

# Build command
CMD="./stress_test_connection_pool.py"

if [ $FULL_SUITE -eq 1 ]; then
    CMD="$CMD --full-suite"
else
    CMD="$CMD --workers $WORKERS --duration $DURATION"
    [ $OPTIMIZED -eq 0 ] && CMD="$CMD --standard"
    [ $SIMULATE_ERRORS -eq 1 ] && CMD="$CMD --simulate-errors"
    [ $USE_TRANSACTIONS -eq 0 ] && CMD="$CMD --no-transactions"
    CMD="$CMD --min-connections $MIN_CONN --max-connections $MAX_CONN"
fi

# Run the stress test
echo "Running command: $CMD"
echo ""
$CMD

# Check if the test was successful
if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}Stress tests completed successfully!${NC}"
    echo "Detailed reports have been saved to the reports directory."
    echo ""
    echo -e "${BLUE}Next Steps:${NC}"
    echo "1. Review the generated reports in the reports directory"
    echo "2. Update your connection pool configuration based on the recommendations"
    echo "3. Use ./update_connection_pool_config.py to apply the recommended settings"
    echo "   Example: ./update_connection_pool_config.py --medium"
else
    echo ""
    echo -e "${RED}Stress tests failed!${NC}"
    echo "Please check the logs for details on what went wrong."
    exit 1
fi