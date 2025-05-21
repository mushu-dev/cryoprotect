#!/bin/bash
# fix_database_integrity.sh
#
# Runs the database integrity fix SQL script to repair issues
# found by the verification tool.
#
# Usage:
#   ./fix_database_integrity.sh [--backup] [--apply]

set -e

echo "CryoProtect Database Integrity Fix Tool"
echo "======================================="

# Parse command line arguments
BACKUP=false
APPLY=false

for arg in "$@"; do
  case $arg in
    --backup)
      BACKUP=true
      shift
      ;;
    --apply)
      APPLY=true
      shift
      ;;
    *)
      # Unknown option
      ;;
  esac
done

# Check if environment variables file exists
if [ ! -f .env ]; then
  echo "Error: .env file not found."
  echo "Make sure you are in the correct directory or create an .env file with database credentials."
  exit 1
fi

# Load environment variables
export $(grep -v '^#' .env | xargs)

# Function to check if database credentials are set
check_credentials() {
  if [ -z "$SUPABASE_DB_HOST" ] || [ -z "$SUPABASE_DB_PORT" ] || 
     [ -z "$SUPABASE_DB_NAME" ] || [ -z "$SUPABASE_DB_USER" ] || 
     [ -z "$SUPABASE_DB_PASSWORD" ]; then
    echo "Error: Missing database credentials in .env file."
    echo "Required variables: SUPABASE_DB_HOST, SUPABASE_DB_PORT, SUPABASE_DB_NAME, SUPABASE_DB_USER, SUPABASE_DB_PASSWORD"
    exit 1
  fi
}

# Function to create database backup
create_backup() {
  echo "Creating database backup before applying fixes..."
  BACKUP_FILE="database_backup_$(date +%Y%m%d_%H%M%S).sql"
  
  PGPASSWORD=$SUPABASE_DB_PASSWORD pg_dump -h $SUPABASE_DB_HOST -p $SUPABASE_DB_PORT \
    -U $SUPABASE_DB_USER -d $SUPABASE_DB_NAME -F p > $BACKUP_FILE
  
  if [ $? -eq 0 ]; then
    echo "✅ Backup created successfully: $BACKUP_FILE"
    # Compress the backup
    gzip $BACKUP_FILE
    echo "✅ Backup compressed: ${BACKUP_FILE}.gz"
  else
    echo "❌ Error creating backup. Aborting."
    exit 1
  fi
}

# Function to apply fixes
apply_fixes() {
  echo "Applying database integrity fixes..."
  
  PGPASSWORD=$SUPABASE_DB_PASSWORD psql -h $SUPABASE_DB_HOST -p $SUPABASE_DB_PORT \
    -U $SUPABASE_DB_USER -d $SUPABASE_DB_NAME -f fix_database_integrity_issues.sql
  
  if [ $? -eq 0 ]; then
    echo "✅ Database fixes applied successfully."
  else
    echo "❌ Error applying fixes."
    exit 1
  fi
}

# Check credentials
check_credentials

# Print connection info
echo "Database: $SUPABASE_DB_NAME at $SUPABASE_DB_HOST:$SUPABASE_DB_PORT"

# Create backup if requested
if [ "$BACKUP" = true ]; then
  create_backup
fi

# Show preview or apply fixes
if [ "$APPLY" = true ]; then
  # Double-check with user before applying
  if [ "$BACKUP" = false ]; then
    read -p "WARNING: No backup requested. Apply fixes without backup? (y/N): " confirm
    if [[ "$confirm" != [yY] ]]; then
      echo "Operation cancelled."
      exit 0
    fi
  fi
  
  apply_fixes
  
  echo -e "\nTo verify the fixes, run:"
  echo "./run_database_integrity_check_fixed.sh"
else
  echo -e "\nPreview mode: SQL fixes were not applied."
  echo "To apply fixes, run the script with --apply flag:"
  echo "  ./fix_database_integrity.sh --backup --apply"
  echo -e "\nTo view the SQL fixes that would be applied:"
  echo "  cat fix_database_integrity_issues.sql"
fi