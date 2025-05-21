#!/bin/bash
# CryoProtect Fedora - Database Migration Runner Script
# This script runs database migrations directly using the Supabase MCP tool

echo "CryoProtect Fedora - Database Migration Runner"
echo "=============================================="
echo

# Ensure variables are set
if [ -z "$SUPABASE_PROJECT_ID" ] && [ -f .env ]; then
  export $(grep -v '^#' .env | xargs)
fi

if [ -z "$SUPABASE_PROJECT_ID" ]; then
  echo "Error: SUPABASE_PROJECT_ID not found in environment or .env file"
  echo "Please set this variable before running migrations"
  exit 1
fi

# Check if the migrations directory exists
if [ ! -d "migrations" ]; then
  echo "Error: migrations directory not found"
  exit 1
fi

# Parse command line arguments
DRY_RUN=0
MIGRATION_ID=""
FORCE=0
LIST_ONLY=0
VERIFY=0

function print_usage {
  echo "Usage: ./run_fedora_migration.sh [OPTIONS]"
  echo "Options:"
  echo "  --dry-run         Show SQL that would be executed without applying changes"
  echo "  --migration=ID    Run a specific migration by file ID (e.g., --migration=001 for 001_initial_schema.sql)"
  echo "  --force           Skip confirmation prompt and apply migrations"
  echo "  --list            List available migrations without applying them" 
  echo "  --verify          Verify database structure after migrations"
  echo "  --help            Show this help message"
}

for arg in "$@"; do
  case $arg in
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --migration=*)
      MIGRATION_ID="${arg#*=}"
      shift
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --list)
      LIST_ONLY=1
      shift
      ;;
    --verify)
      VERIFY=1
      shift
      ;;
    --help)
      print_usage
      exit 0
      ;;
    *)
      echo "Unknown option: $arg"
      print_usage
      exit 1
      ;;
  esac
done

# List all available migrations
function list_migrations {
  echo "Available migrations:"
  echo "---------------------"
  ls -1 migrations/ | grep "^[0-9].*\.sql$" | sort -n | while read -r file; do
    local id=$(echo "$file" | grep -o "^[0-9]\+")
    local desc=$(echo "$file" | sed -E 's/^[0-9]+_(.*)\.sql/\1/' | tr '_' ' ')
    echo "[$id] $desc"
  done
}

if [ $LIST_ONLY -eq 1 ]; then
  list_migrations
  exit 0
fi

# Get the list of migration files to run
if [ -n "$MIGRATION_ID" ]; then
  MIGRATION_FILES=$(ls -1 migrations/ | grep "^${MIGRATION_ID}.*\.sql$" | sort -n)
  if [ -z "$MIGRATION_FILES" ]; then
    echo "Error: No migration found with ID $MIGRATION_ID"
    echo
    list_migrations
    exit 1
  fi
else
  MIGRATION_FILES=$(ls -1 migrations/ | grep "^[0-9].*\.sql$" | sort -n)
fi

# Show the migrations that will be applied
echo "The following migrations will be applied:"
echo "----------------------------------------"
for file in $MIGRATION_FILES; do
  echo "- $file"
done
echo

# If dry run, print SQL without executing
if [ $DRY_RUN -eq 1 ]; then
  echo "DRY RUN MODE - SQL will not be executed"
  echo "---------------------------------------"
  for file in $MIGRATION_FILES; do
    echo "🔍 Previewing migration: $file"
    echo
    cat "migrations/$file" | sed 's/^/    /'
    echo
    echo "---------------------------------------"
  done
  exit 0
fi

# Confirm before proceeding
if [ $FORCE -eq 0 ]; then
  read -p "Are you sure you want to apply these migrations? (y/n): " CONFIRM
  if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
    echo "Migration cancelled."
    exit 0
  fi
fi

# Apply the migrations
echo
echo "Applying migrations..."
echo

SUCCESS_COUNT=0

for file in $MIGRATION_FILES; do
  echo "📦 Applying migration: $file"
  
  # Generate a unique name for the migration based on the file name
  MIGRATION_NAME=$(echo "$file" | sed -E 's/^[0-9]+_(.*)\.sql/\1/')
  
  # Apply the migration
  echo "python -c \"import os; from supabase_mcp_tools import apply_migration; apply_migration('$SUPABASE_PROJECT_ID', '$MIGRATION_NAME', open('migrations/$file').read())\"" | bash
  
  if [ $? -eq 0 ]; then
    echo "✅ Successfully applied: $file"
    SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
  else
    echo "❌ Failed to apply: $file"
    echo "Migration process aborted."
    exit 1
  fi
  
  echo
done

# Verify database structure
if [ $VERIFY -eq 1 ]; then
  echo "Verifying database structure..."
  echo
  
  # Run verification of tables and schema
  python verify_database_structure.py
  
  if [ $? -eq 0 ]; then
    echo "✅ Database verification successful!"
  else
    echo "⚠️ Database verification completed with warnings."
  fi
  
  echo
fi

echo "Migration completed successfully."
echo "$SUCCESS_COUNT migration(s) applied."
echo
echo "Next steps:"
echo "1. Verify the database schema with the verification script"
echo "2. Run tests to ensure the application works with the new schema"
echo
echo "For more information, see DATABASE_SCHEMA_SUMMARY.md"