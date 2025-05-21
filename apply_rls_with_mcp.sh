#!/bin/bash

# Script to apply RLS optimizations using Supabase MCP
# This script uses the Supabase MCP to apply the RLS optimizations

# Set script to exit immediately if a command exits with a non-zero status
set -e

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$SCRIPT_DIR"

# Function to display usage information
show_usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --security-definer-only Only apply the security definer functions"
    echo "  --optimized-only        Only apply the optimized RLS policies"
    echo "  --project-id ID         Supabase project ID (required)"
    echo "  --help                  Display this help message"
}

# Parse command line arguments
SECURITY_DEFINER_ONLY=0
OPTIMIZED_ONLY=0
PROJECT_ID=""

while [ $# -gt 0 ]; do
    case "$1" in
        --security-definer-only)
            SECURITY_DEFINER_ONLY=1
            shift
            ;;
        --optimized-only)
            OPTIMIZED_ONLY=1
            shift
            ;;
        --project-id)
            PROJECT_ID="$2"
            shift 2
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

# Check if project ID is provided
if [ -z "$PROJECT_ID" ]; then
    echo "Error: Supabase project ID is required"
    show_usage
    exit 1
fi

# Directory for migrations
MIGRATIONS_DIR="$SCRIPT_DIR/migrations"

# Check if migrations directory exists
if [ ! -d "$MIGRATIONS_DIR" ]; then
    echo "Error: Migrations directory not found: $MIGRATIONS_DIR"
    exit 1
fi

# Paths to migration files
SECURITY_DEFINER_FILE="$MIGRATIONS_DIR/019_rls_security_definer_functions.sql"
OPTIMIZED_POLICIES_FILE="$MIGRATIONS_DIR/020_optimized_rls_policies.sql"

# Check if migration files exist
if [ $OPTIMIZED_ONLY -eq 0 ] && [ ! -f "$SECURITY_DEFINER_FILE" ]; then
    echo "Error: Security definer functions migration file not found: $SECURITY_DEFINER_FILE"
    exit 1
fi

if [ $SECURITY_DEFINER_ONLY -eq 0 ] && [ ! -f "$OPTIMIZED_POLICIES_FILE" ]; then
    echo "Error: Optimized RLS policies migration file not found: $OPTIMIZED_POLICIES_FILE"
    exit 1
fi

# Ensure reports directory exists
mkdir -p reports

# Print header
echo "=========================================================="
echo "   CryoProtect RLS Optimization - Using Supabase MCP"
echo "=========================================================="
echo

echo "Execution started at: $(date)"
if [ $SECURITY_DEFINER_ONLY -eq 1 ]; then
    echo "Target: Security definer functions only"
elif [ $OPTIMIZED_ONLY -eq 1 ]; then
    echo "Target: Optimized RLS policies only"
else
    echo "Target: Full RLS optimization"
fi

echo "Supabase Project ID: $PROJECT_ID"
echo

# Configure arguments
ARGS="--project-id $PROJECT_ID"

if [ $SECURITY_DEFINER_ONLY -eq 1 ]; then
    ARGS="$ARGS --security-definer-only"
fi

if [ $OPTIMIZED_ONLY -eq 1 ]; then
    ARGS="$ARGS --optimized-only"
fi

echo "Running RLS optimization script using MCP..."
echo

# Execute the main script
python apply_rls_with_mcp.py $ARGS

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
echo "   See reports directory for details"
echo "=========================================================="

exit $STATUS