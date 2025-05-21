#!/bin/bash
# Apply RLS optimizations using Supabase MCP

# Set default project ID
DEFAULT_PROJECT_ID="tsdlmynydfuypiugmkev"
PROJECT_ID=${1:-$DEFAULT_PROJECT_ID}

# SQL file path
SQL_FILE="migrations/all_complex_rls_optimizations.sql"

# Check if the SQL file exists
if [ ! -f "$SQL_FILE" ]; then
  echo "Error: SQL file not found: $SQL_FILE"
  exit 1
fi

echo "Applying RLS optimizations to project $PROJECT_ID..."

# Split large SQL file into smaller chunks to avoid MCP token limits
mkdir -p migrations/chunks
split -l 300 "$SQL_FILE" migrations/chunks/chunk_

# Apply each chunk sequentially
for chunk in migrations/chunks/chunk_*; do
  echo "Applying chunk: $chunk"
  python -c "
import sys
import os
import json
import time

try:
    from supabase_mcp_tools import execute_sql_on_supabase
    # Read the SQL file
    with open('$chunk', 'r') as f:
        sql = f.read()
    
    # Execute the SQL
    print('Executing SQL chunk...')
    result = execute_sql_on_supabase('$PROJECT_ID', sql)
    print('SQL executed successfully:', result)
except Exception as e:
    print('Error executing SQL:', str(e))
    sys.exit(1)
"
  status=$?
  if [ $status -ne 0 ]; then
    echo "Failed to apply chunk: $chunk"
    exit $status
  fi
  
  # Add a small delay between chunks
  sleep 2
done

# Clean up chunk files
rm -rf migrations/chunks

echo "RLS optimizations applied successfully!"

# Use the Supabase MCP to verify the optimizations
echo "Verifying optimizations..."
python -c "
import sys
import os
import json
import time

try:
    from supabase_mcp_tools import execute_sql_on_supabase
    
    # Define verification queries
    verification_queries = {
        'security_definer_functions': '''
            SELECT proname, prosecdef
            FROM pg_proc p
            JOIN pg_namespace n ON p.pronamespace = n.oid
            WHERE n.nspname = 'public'
            AND proname IN (
                'has_molecule_access', 'has_mixture_access', 'is_team_member',
                'user_has_clearance', 'find_molecules_by_property_range',
                'search_molecules_text', 'filter_accessible_molecules'
            );
        ''',
        
        'performance_indexes': '''
            SELECT indexname, indexdef
            FROM pg_indexes
            WHERE schemaname = 'public'
            AND (
                indexname LIKE 'idx_%' OR
                indexname LIKE '%_molecules_%' OR
                indexname LIKE '%_mixtures_%' OR
                indexname LIKE '%_properties_%'
            )
            LIMIT 10;
        ''',
        
        'materialized_views': '''
            SELECT matviewname
            FROM pg_matviews
            WHERE schemaname = 'public';
        '''
    }
    
    # Execute verification queries
    for name, query in verification_queries.items():
        print(f'\\nVerifying {name}...')
        result = execute_sql_on_supabase('$PROJECT_ID', query)
        if result:
            print(f'✅ {name} verification passed!')
            print(f'Found: {len(result)} items')
        else:
            print(f'❌ {name} verification failed!')
    
    # Test optimized queries
    print('\\nTesting optimized get_molecules_with_properties function...')
    test_query = 'SELECT * FROM get_molecules_with_properties(5, 0);'
    try:
        result = execute_sql_on_supabase('$PROJECT_ID', test_query)
        if result:
            print('✅ Query test passed!')
            print(f'Found: {len(result)} molecules')
        else:
            print('❌ Query test returned no results')
    except Exception as e:
        print('❌ Query test failed:', str(e))
    
except Exception as e:
    print('Error during verification:', str(e))
    sys.exit(1)
"

echo "Verification completed!"
echo "RLS optimization process finished."