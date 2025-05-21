#!/bin/bash
# Run the RLS complex query optimization script with various options

# Set default database parameters
DB_HOST=${DB_HOST:-"localhost"}
DB_PORT=${DB_PORT:-5432}
DB_NAME=${DB_NAME:-"cryoprotect_db"}
DB_USER=${DB_USER:-"postgres"}
DB_PASSWORD=${DB_PASSWORD:-"password"}
PROJECT_ID=${PROJECT_ID:-""}

# Parse command line arguments
POSITIONAL=()
MODE="direct"
VERIFY_ONLY=false
TEST_PERFORMANCE=false
SKIP_FUNCS=false
SKIP_INDEXES=false
SKIP_VIEWS=false
SKIP_POLICIES=false
SKIP_CQ=false
DRY_RUN=false
FORCE=false

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        --supabase)
        MODE="supabase"
        shift
        ;;
        --mcp)
        MODE="mcp"
        shift
        ;;
        --project-id)
        PROJECT_ID="$2"
        shift
        shift
        ;;
        --db-host)
        DB_HOST="$2"
        shift
        shift
        ;;
        --db-port)
        DB_PORT="$2"
        shift
        shift
        ;;
        --db-name)
        DB_NAME="$2"
        shift
        shift
        ;;
        --db-user)
        DB_USER="$2"
        shift
        shift
        ;;
        --db-password)
        DB_PASSWORD="$2"
        shift
        shift
        ;;
        --verify-only)
        VERIFY_ONLY=true
        shift
        ;;
        --test-performance)
        TEST_PERFORMANCE=true
        shift
        ;;
        --skip-functions)
        SKIP_FUNCS=true
        shift
        ;;
        --skip-indexes)
        SKIP_INDEXES=true
        shift
        ;;
        --skip-views)
        SKIP_VIEWS=true
        shift
        ;;
        --skip-policies)
        SKIP_POLICIES=true
        shift
        ;;
        --skip-complex-queries)
        SKIP_CQ=true
        shift
        ;;
        --dry-run)
        DRY_RUN=true
        shift
        ;;
        --force)
        FORCE=true
        shift
        ;;
        --help)
        echo "RLS Complex Query Optimization Script"
        echo ""
        echo "Usage: $0 [options]"
        echo ""
        echo "Connection Options:"
        echo "  --supabase             Use Supabase client (default: direct)"
        echo "  --mcp                  Use Supabase MCP"
        echo "  --project-id ID        Supabase project ID (required for MCP)"
        echo ""
        echo "Database Connection Options (for direct mode):"
        echo "  --db-host HOST         Database host (default: localhost)"
        echo "  --db-port PORT         Database port (default: 5432)"
        echo "  --db-name NAME         Database name (default: cryoprotect_db)"
        echo "  --db-user USER         Database user (default: postgres)"
        echo "  --db-password PASS     Database password (default: password)"
        echo ""
        echo "Optimization Options:"
        echo "  --verify-only          Only verify optimizations without applying them"
        echo "  --test-performance     Test query performance after optimization"
        echo "  --skip-functions       Skip creating security definer functions"
        echo "  --skip-indexes         Skip creating performance indexes"
        echo "  --skip-views           Skip creating materialized views"
        echo "  --skip-policies        Skip creating RLS policies"
        echo "  --skip-complex-queries Skip complex query optimizations"
        echo "  --dry-run              Show what would be done without making changes"
        echo "  --force                Apply optimizations even if they already exist"
        echo ""
        echo "Examples:"
        echo "  $0 --direct --db-host localhost --db-name cryoprotect_db --test-performance"
        echo "  $0 --supabase --test-performance"
        echo "  $0 --mcp --project-id abc123 --verify-only"
        exit 0
        ;;
        *)
        POSITIONAL+=("$1")
        shift
        ;;
    esac
done
set -- "${POSITIONAL[@]}"

# Check if project ID is provided for MCP mode
if [ "$MODE" == "mcp" ] && [ -z "$PROJECT_ID" ]; then
    echo "Error: Project ID is required for MCP mode"
    echo "Use --project-id to specify the project ID"
    exit 1
fi

# Build command line arguments
ARGS=""

# Add connection mode
if [ "$MODE" == "supabase" ]; then
    ARGS="$ARGS --supabase"
elif [ "$MODE" == "mcp" ]; then
    ARGS="$ARGS --mcp --project-id $PROJECT_ID"
else
    ARGS="$ARGS --direct --db-host $DB_HOST --db-port $DB_PORT --db-name $DB_NAME --db-user $DB_USER --db-password $DB_PASSWORD"
fi

# Add optimization options
if [ "$VERIFY_ONLY" = true ]; then
    ARGS="$ARGS --verify-only"
fi

if [ "$TEST_PERFORMANCE" = true ]; then
    ARGS="$ARGS --test-performance"
fi

if [ "$SKIP_FUNCS" = true ]; then
    ARGS="$ARGS --skip-functions"
fi

if [ "$SKIP_INDEXES" = true ]; then
    ARGS="$ARGS --skip-indexes"
fi

if [ "$SKIP_VIEWS" = true ]; then
    ARGS="$ARGS --skip-views"
fi

if [ "$SKIP_POLICIES" = true ]; then
    ARGS="$ARGS --skip-policies"
fi

if [ "$SKIP_CQ" = true ]; then
    ARGS="$ARGS --skip-complex-queries"
fi

if [ "$DRY_RUN" = true ]; then
    ARGS="$ARGS --dry-run"
fi

if [ "$FORCE" = true ]; then
    ARGS="$ARGS --force"
fi

# Make script executable if needed
chmod +x optimize_complex_rls_queries.py

# Run the optimization script
echo "Running: python optimize_complex_rls_queries.py $ARGS"
python optimize_complex_rls_queries.py $ARGS

exit $?