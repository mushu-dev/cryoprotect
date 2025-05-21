# Practical Database Approach for Fedora Migration

This document outlines a pragmatic approach for handling database operations in the CryoProtect project on Fedora Linux.

## Current Status

We have successfully tested direct Supabase connectivity in our simplified app:
- REST API connection is working properly
- Environment variables are correctly loaded
- Basic queries are functional

## Recommended Approach

Instead of creating a new abstraction layer or connection mechanism, we should build on what's already working:

1. **Keep Using `simplified_app.py` Pattern**
   - The REST API approach is simple and effective
   - Already tested and working in Fedora environment
   - Minimal dependencies and complexity

2. **Implement Migrations Using Native Tools**
   - Use the existing migration scripts directly
   - Apply migrations sequentially with proper verification
   - Document the migration process for future reference

3. **Document Database Operations**
   - Create a simple reference guide for common database operations
   - Include examples for different query types
   - Keep complexity low and focus on what's needed

## Database Schema Management

For schema management, we'll use a phased approach:

1. **Phase 1: Apply Base Schema**
   - Run initial migrations (001-005) to set up core tables
   - Verify schema using REST API and visual inspection
   - Test basic operations (insert, select, update)

2. **Phase 2: Apply Security Policies**
   - Add RLS policies (006, 018, 019)
   - Implement service role authentication
   - Test access control

3. **Phase 3: Scientific Data**
   - Seed scientific data (007, 008)
   - Verify data integrity
   - Test data retrieval and querying

4. **Phase 4: Additional Features**
   - Add protocol storage (009)
   - Implement performance indexes (010)
   - Add RBAC schema (011)
   - Add toxicity schema (012)
   - Add lab verification schema (013)
   - Apply minor schema updates (015-017)

## Migration Script

We'll create a simple shell script to manage migrations:

```bash
#!/bin/bash
# apply_migrations.sh
# Simple script to apply migrations sequentially

# Load environment variables
source .env

# Set default parameters
START_INDEX=1
END_INDEX=19
DRY_RUN=false
VERIFY=true

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --start=*)
      START_INDEX="${1#*=}"
      shift
      ;;
    --end=*)
      END_INDEX="${1#*=}"
      shift
      ;;
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    --no-verify)
      VERIFY=false
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Function to apply a migration
apply_migration() {
  local index=$1
  local padded_index=$(printf "%03d" $index)
  local migration_file=$(find migrations -name "${padded_index}_*.sql" | sort | head -n 1)
  
  if [ -z "$migration_file" ]; then
    echo "No migration found with index $padded_index"
    return 1
  fi
  
  echo "Applying migration: $migration_file"
  
  if [ "$DRY_RUN" = true ]; then
    echo "DRY RUN: Would apply $migration_file"
  else
    # Create request to apply migration
    curl -X POST \
      "${SUPABASE_URL}/rest/v1/rpc/exec_sql" \
      -H "apikey: ${SUPABASE_SERVICE_KEY}" \
      -H "Authorization: Bearer ${SUPABASE_SERVICE_KEY}" \
      -H "Content-Type: application/json" \
      -d "{\"query\": \"$(cat $migration_file | tr -d '\n' | sed 's/"/\\"/g')\"}"
    
    echo ""
    echo "Migration $padded_index applied."
  fi
}

# Function to verify a migration
verify_migration() {
  local index=$1
  echo "Verifying migration $index..."
  
  # Add verification logic here based on migration index
  # For example, check if tables exist, count records, etc.
}

# Apply migrations in sequence
for i in $(seq $START_INDEX $END_INDEX); do
  apply_migration $i
  
  if [ $? -ne 0 ]; then
    echo "Failed to apply migration $i. Stopping."
    exit 1
  fi
  
  if [ "$VERIFY" = true ]; then
    verify_migration $i
  fi
done

echo "Migration process completed."
```

## Common Database Operations

Here are simplified examples for common operations:

### REST API Operations

```python
# Get data from a table
def get_data(table_name, filters=None, limit=100):
    url = f"{SUPABASE_URL}/rest/v1/{table_name}"
    params = {'limit': limit}
    
    if filters:
        params.update(filters)
    
    headers = {
        'apikey': SUPABASE_KEY,
        'Authorization': f'Bearer {SUPABASE_KEY}'
    }
    
    response = requests.get(url, headers=headers, params=params)
    response.raise_for_status()
    return response.json()

# Insert data
def insert_data(table_name, data):
    url = f"{SUPABASE_URL}/rest/v1/{table_name}"
    
    headers = {
        'apikey': SUPABASE_KEY,
        'Authorization': f'Bearer {SUPABASE_KEY}',
        'Prefer': 'return=representation'
    }
    
    response = requests.post(url, headers=headers, json=data)
    response.raise_for_status()
    return response.json()

# Update data
def update_data(table_name, id, data):
    url = f"{SUPABASE_URL}/rest/v1/{table_name}?id=eq.{id}"
    
    headers = {
        'apikey': SUPABASE_KEY,
        'Authorization': f'Bearer {SUPABASE_KEY}',
        'Prefer': 'return=representation'
    }
    
    response = requests.patch(url, headers=headers, json=data)
    response.raise_for_status()
    return response.json()

# Delete data
def delete_data(table_name, id):
    url = f"{SUPABASE_URL}/rest/v1/{table_name}?id=eq.{id}"
    
    headers = {
        'apikey': SUPABASE_KEY,
        'Authorization': f'Bearer {SUPABASE_KEY}'
    }
    
    response = requests.delete(url, headers=headers)
    response.raise_for_status()
    return True
```

### SQL Operations (When Needed)

For complex operations that require SQL, we can use the `exec_sql` RPC function:

```python
# Execute SQL query
def execute_sql(query):
    url = f"{SUPABASE_URL}/rest/v1/rpc/exec_sql"
    
    headers = {
        'apikey': SUPABASE_SERVICE_KEY,
        'Authorization': f'Bearer {SUPABASE_SERVICE_KEY}',
        'Content-Type': 'application/json'
    }
    
    data = {'query': query}
    
    response = requests.post(url, headers=headers, json=data)
    response.raise_for_status()
    return response.json()
```

## Benefits of This Approach

1. **Simplicity**: Building on what's already working rather than adding complexity
2. **Reliability**: Using well-tested REST API patterns that are standard and stable
3. **Maintainability**: Less code to manage and maintain
4. **Compatibility**: Works across different environments
5. **Focus**: Allows us to focus on actual functionality rather than infrastructure

## Next Steps

1. Create the migration script and test it with individual migrations
2. Verify each migration result to ensure database integrity
3. Document the schema after migrations are applied
4. Implement basic data access functions for common operations
5. Build core API functionality on these simple foundations