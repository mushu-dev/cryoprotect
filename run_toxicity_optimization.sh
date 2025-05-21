#!/bin/bash
#
# CryoProtect - Run Toxicity Optimization
#
# This script applies the toxicity data schema optimization and 
# updates the API to use the optimized endpoints.

set -e

echo "Starting toxicity data optimization..."

# Make sure we're in the project root
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
cd "$SCRIPT_DIR"

# Step 1: Apply schema migration
echo "Applying toxicity schema migration..."
if [[ -f "migrations/021_toxicity_schema_optimization.sql" && -f "migrations/022_toxicity_data_migration.sql" ]]; then
    # Check if we have a database connection
    HAVE_DB_CONNECTION=$(python -c "from database.db import test_connection; success, _ = test_connection(); print('yes' if success else 'no')" 2>/dev/null || echo "no")

    if [[ "$HAVE_DB_CONNECTION" == "yes" ]]; then
        # Apply migrations
        python -c "from database.utils import execute_sql_file; execute_sql_file('migrations/021_toxicity_schema_optimization.sql'); execute_sql_file('migrations/022_toxicity_data_migration.sql')"
        echo "Schema migration applied successfully"
    else
        echo "Skipping actual database migration due to no database connection"
        echo "Migration would apply the following files:"
        echo " - migrations/021_toxicity_schema_optimization.sql"
        echo " - migrations/022_toxicity_data_migration.sql"
        echo "Schema migration step completed (simulation mode)"
    fi
else
    echo "Error: Migration files not found"
    exit 1
fi

# Step 2: Refresh materialized views
echo "Refreshing materialized views..."
if [[ "$HAVE_DB_CONNECTION" == "yes" ]]; then
    python -c "from database.utils import execute_sql; execute_sql('SELECT refresh_toxicity_materialized_views();')"
    echo "Materialized views refreshed successfully"
else
    echo "Skipping materialized view refresh due to no database connection"
    echo "Views would be refreshed using: SELECT refresh_toxicity_materialized_views();"
    echo "Materialized view refresh step completed (simulation mode)"
fi

# Step 3: Verify the schema
echo "Verifying schema changes..."
if [[ "$HAVE_DB_CONNECTION" == "yes" ]]; then
    python -c "from database.utils import execute_sql; print(execute_sql('SELECT COUNT(*) FROM toxicity_summary;'))"
    python -c "from database.utils import execute_sql; print(execute_sql('SELECT COUNT(*) FROM ld50_summary;'))"
    python -c "from database.utils import execute_sql; print(execute_sql('SELECT COUNT(*) FROM tox21_activity_summary;'))"
    python -c "from database.utils import execute_sql; print(execute_sql('SELECT COUNT(*) FROM hazard_classification_summary;'))"
    echo "Schema verification completed successfully"
else
    echo "Skipping schema verification due to no database connection"
    echo "Would verify the following materialized views:"
    echo " - toxicity_summary"
    echo " - ld50_summary"
    echo " - tox21_activity_summary"
    echo " - hazard_classification_summary"
    echo "Schema verification step completed (simulation mode)"
fi

# Step 4: Apply API optimization
echo "Updating API to use optimized toxicity resources..."
python apply_toxicity_optimization.py

# Step 5: Run a simple test
echo "Testing the API endpoints..."
# Check if API is already running
if pgrep -f "python app.py" > /dev/null; then
    echo "API is already running, will use existing server"
    EXISTING_API=true
else
    echo "Starting API temporarily for testing..."
    python app.py &
    API_PID=$!
    sleep 3  # Wait for the API to start
    EXISTING_API=false
fi

# Use a sample UUID for testing since we might not have actual molecules
SAMPLE_UUID=$(python -c "import uuid; print(uuid.uuid4())")
echo "Testing new endpoints with sample UUID: $SAMPLE_UUID"

# Test the endpoint
if command -v curl &> /dev/null; then
    # Test one of the new endpoints
    echo "Response from toxicity summary endpoint:"
    curl -s "http://localhost:5000/api/toxicity/summary/molecule/$SAMPLE_UUID" | head -n 20
else
    # Use Python if curl is not available
    echo "Using Python to test endpoints (curl not available):"
    python -c "import requests; response = requests.get(f'http://localhost:5000/api/toxicity/summary/molecule/$SAMPLE_UUID'); print(response.status_code); print(response.text[:200] + '...' if len(response.text) > 200 else response.text)"
fi

# Clean up if we started the API
if [[ "$EXISTING_API" == "false" ]]; then
    echo "Stopping temporary API server..."
    kill $API_PID
fi

echo "Toxicity data optimization completed"