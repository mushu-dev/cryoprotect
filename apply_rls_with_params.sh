#!/bin/bash

# Script to apply RLS optimizations with connection parameters
# This script allows you to specify database connection details

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
    echo "  --host HOST             Database host (default: localhost)"
    echo "  --port PORT             Database port (default: 5432)"
    echo "  --dbname DBNAME         Database name (default: postgres)"
    echo "  --user USER             Database user (default: postgres)"
    echo "  --password PASSWORD     Database password"
    echo "  --sslmode MODE          SSL mode (disable, allow, prefer, require, verify-ca, verify-full)"
    echo "                          Default: require"
    echo "  --help                  Display this help message"
}

# Parse command line arguments
SECURITY_DEFINER_ONLY=0
OPTIMIZED_ONLY=0
DB_HOST="localhost"
DB_PORT="5432"
DB_NAME="postgres"
DB_USER="postgres"
DB_PASSWORD=""
DB_SSLMODE="require"

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
        --host)
            DB_HOST="$2"
            shift 2
            ;;
        --port)
            DB_PORT="$2"
            shift 2
            ;;
        --dbname)
            DB_NAME="$2"
            shift 2
            ;;
        --user)
            DB_USER="$2"
            shift 2
            ;;
        --password)
            DB_PASSWORD="$2"
            shift 2
            ;;
        --sslmode)
            DB_SSLMODE="$2"
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

# Directory for migrations
MIGRATIONS_DIR="$SCRIPT_DIR/migrations"

# Check if migrations directory exists
if [ ! -d "$MIGRATIONS_DIR" ]; then
    echo "Error: Migrations directory not found: $MIGRATIONS_DIR"
    exit 1
fi

# Paths to migration files
SECURITY_DEFINER_FILE="$MIGRATIONS_DIR/simplified_rls_functions.sql"
OPTIMIZED_POLICIES_FILE="$MIGRATIONS_DIR/simplified_rls_functions.sql"

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
echo "   CryoProtect RLS Optimization - With Parameters"
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

echo "Database: $DB_USER@$DB_HOST:$DB_PORT/$DB_NAME"
echo

# Check if we should prompt for password
if [ -z "$DB_PASSWORD" ]; then
    # Prompt for database password if not provided
    read -s -p "Enter database password: " DB_PASSWORD
    echo
fi

# Configure arguments
ARGS=""

if [ $SECURITY_DEFINER_ONLY -eq 1 ]; then
    ARGS="$ARGS --security-definer-only"
fi

if [ $OPTIMIZED_ONLY -eq 1 ]; then
    ARGS="$ARGS --optimized-only"
fi

# Add connection parameters
ARGS="$ARGS --host $DB_HOST --port $DB_PORT --dbname $DB_NAME --user $DB_USER --sslmode $DB_SSLMODE"

# Add password if provided
if [ -n "$DB_PASSWORD" ]; then
    ARGS="$ARGS --password $DB_PASSWORD"
fi

echo "Running RLS optimization script with parameters..."
echo

# Create Python virtual environment if it doesn't exist
if [ ! -d "quick_env" ]; then
    echo "Creating Python virtual environment..."
    python3 -m venv quick_env
    source quick_env/bin/activate
    pip install --upgrade pip
    pip install psycopg2-binary
else
    source quick_env/bin/activate
fi

# Execute the main script
python apply_rls_with_params.py $ARGS

# Get the return status
STATUS=$?

# Deactivate virtual environment
deactivate

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