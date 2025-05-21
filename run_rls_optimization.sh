#!/bin/bash

# Run RLS optimization script
# This script applies the optimized RLS policies to the database

# Set script to exit immediately if a command exits with a non-zero status
set -e

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$SCRIPT_DIR"

# Check if log directory exists
if [ ! -d "logs" ]; then
    mkdir -p logs
fi

# Function to display usage information
show_usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --dry-run               Show what would be done without making changes"
    echo "  --skip-functions        Skip creating security definer functions"
    echo "  --skip-indexes          Skip creating performance indexes"
    echo "  --skip-views            Skip creating materialized views"
    echo "  --skip-policies         Skip creating RLS policies"
    echo "  --verify                Only verify that optimizations have been applied"
    echo "  --performance-test      Run performance tests after applying optimizations"
    echo "  --help                  Display this help message"
}

# Parse command line arguments
DRY_RUN=0
SKIP_FUNCTIONS=0
SKIP_INDEXES=0
SKIP_VIEWS=0
SKIP_POLICIES=0
VERIFY_ONLY=0
PERFORMANCE_TEST=0

while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --skip-functions)
            SKIP_FUNCTIONS=1
            shift
            ;;
        --skip-indexes)
            SKIP_INDEXES=1
            shift
            ;;
        --skip-views)
            SKIP_VIEWS=1
            shift
            ;;
        --skip-policies)
            SKIP_POLICIES=1
            shift
            ;;
        --verify)
            VERIFY_ONLY=1
            shift
            ;;
        --performance-test)
            PERFORMANCE_TEST=1
            shift
            ;;
        --help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Print header
echo "=========================================================="
echo "   CryoProtect RLS Optimization"
echo "=========================================================="
echo

# Log the command execution
echo "Execution started at: $(date)"
if [ $DRY_RUN -eq 1 ]; then
    echo "Mode: Dry run (no changes will be made)"
elif [ $VERIFY_ONLY -eq 1 ]; then
    echo "Mode: Verify only"
else
    echo "Mode: Apply and verify"
fi

# Configure command line arguments for the Python script
ARGS=""

if [ $DRY_RUN -eq 1 ]; then
    ARGS="$ARGS --dry-run"
fi

if [ $SKIP_FUNCTIONS -eq 1 ]; then
    ARGS="$ARGS --skip-functions"
fi

if [ $SKIP_INDEXES -eq 1 ]; then
    ARGS="$ARGS --skip-indexes"
fi

if [ $SKIP_VIEWS -eq 1 ]; then
    ARGS="$ARGS --skip-views"
fi

if [ $SKIP_POLICIES -eq 1 ]; then
    ARGS="$ARGS --skip-policies"
fi

if [ $VERIFY_ONLY -eq 1 ]; then
    ARGS="$ARGS --verify"
fi

if [ $PERFORMANCE_TEST -eq 1 ]; then
    ARGS="$ARGS --performance-test"
fi

echo "Executing: python apply_rls_optimization.py $ARGS"
echo

# Create virtual environment if needed and activate it
if [ ! -d "./quick_env" ]; then
    echo "Creating Python virtual environment..."
    python3 -m venv quick_env
    source quick_env/bin/activate
    pip install --upgrade pip
    pip install psycopg2-binary
else
    source quick_env/bin/activate
fi

# Execute the main script
python apply_rls_optimization.py $ARGS 2>&1 | tee logs/rls_optimization_$(date +%Y%m%d_%H%M%S).log

# Capture the exit status
STATUS=${PIPESTATUS[0]}

# Run performance tests if optimization was successful
if [ $STATUS -eq 0 ] && [ $PERFORMANCE_TEST -eq 1 ]; then
    echo "Running performance tests..."
    python test_rls_optimization.py --benchmark-only 2>&1 | tee logs/rls_performance_$(date +%Y%m%d_%H%M%S).log

    # Update status only if the test fails
    TEST_STATUS=${PIPESTATUS[0]}
    if [ $TEST_STATUS -ne 0 ]; then
        STATUS=$TEST_STATUS
    fi

    # Find the latest performance report
    LATEST_REPORT=$(ls -t reports/rls_performance_report_*.md 2>/dev/null | head -1)
    if [ -n "$LATEST_REPORT" ]; then
        echo "Performance report generated: $LATEST_REPORT"
        echo "Preview:"
        head -10 "$LATEST_REPORT"
        echo "..."
    fi
fi

# Get the return status
STATUS=$?

# Print footer
echo
echo "=========================================================="
if [ $STATUS -eq 0 ]; then
    echo "   RLS Optimization Completed Successfully!"
else
    echo "   RLS Optimization Encountered Issues!"
fi
echo "   See logs folder for details"
echo "=========================================================="

# Deactivate virtual environment
deactivate

exit $STATUS