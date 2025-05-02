#!/bin/bash
set -e

# Database Migration Script for Blue/Green Deployment
# This script handles database migrations as part of the deployment process

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BACKUP_DIR="$PROJECT_DIR/database_backups"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --environment)
      ENVIRONMENT="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    --help)
      echo "Usage: $0 [options]"
      echo "Options:"
      echo "  --environment ENV    Specify the environment to run migrations from (blue or green)"
      echo "  --dry-run            Run migrations in dry-run mode without applying changes"
      echo "  --help               Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Check if environment is provided
if [ -z "$ENVIRONMENT" ]; then
  echo "Error: Environment is required. Use --environment to specify the environment (blue or green)."
  exit 1
fi

# Validate environment
if [ "$ENVIRONMENT" != "blue" ] && [ "$ENVIRONMENT" != "green" ]; then
  echo "Error: Invalid environment. Must be 'blue' or 'green'."
  exit 1
fi

# Create backup directory if it doesn't exist
mkdir -p "$BACKUP_DIR"

# Check if the environment is running
if ! docker ps --filter "name=cryoprotect-$ENVIRONMENT" --quiet | grep -q .; then
  echo "Error: The $ENVIRONMENT environment is not running. Cannot run migrations."
  exit 1
fi

# Create database backup before migration
echo "Creating database backup before migration..."
BACKUP_FILE="$BACKUP_DIR/pre_migration_backup_${TIMESTAMP}.sql"

docker exec "$(docker ps --filter "name=cryoprotect-$ENVIRONMENT" --quiet)" \
  /opt/conda/envs/cryoprotect/bin/python -c "
import os
from database.utils.connection import create_connection
conn = create_connection()
result = conn.rpc('create_database_backup', {'filename': 'pre_migration_backup_${TIMESTAMP}'}).execute()
print(f'Backup result: {result.data if hasattr(result, \"data\") else \"No data\"}')
" || {
  echo "Warning: Failed to create database backup. Proceeding with caution..."
}

# Run database migrations
echo "Running database migrations from $ENVIRONMENT environment..."
DRY_RUN_ARG=""
if [ "$DRY_RUN" = true ]; then
  DRY_RUN_ARG="--dry-run"
  echo "Running in dry-run mode. No changes will be applied."
fi

docker exec "$(docker ps --filter "name=cryoprotect-$ENVIRONMENT" --quiet)" \
  /opt/conda/envs/cryoprotect/bin/python -m database.migrations apply $DRY_RUN_ARG || {
  echo "Error: Database migration failed."
  
  # Ask if we should restore from backup
  if [ -z "$DRY_RUN" ]; then
    read -p "Do you want to restore the database from backup? (y/n): " RESTORE
    if [ "$RESTORE" = "y" ] || [ "$RESTORE" = "Y" ]; then
      echo "Restoring database from backup..."
      docker exec "$(docker ps --filter "name=cryoprotect-$ENVIRONMENT" --quiet)" \
        /opt/conda/envs/cryoprotect/bin/python -c "
import os
from database.utils.connection import create_connection
conn = create_connection()
result = conn.rpc('restore_database_backup', {'filename': 'pre_migration_backup_${TIMESTAMP}'}).execute()
print(f'Restore result: {result.data if hasattr(result, \"data\") else \"No data\"}')
"
      if [ $? -eq 0 ]; then
        echo "Database restored successfully."
      else
        echo "Error: Failed to restore database. Manual intervention required."
        exit 1
      fi
    else
      echo "Database not restored. Manual intervention may be required."
      exit 1
    fi
  fi
  
  exit 1
}

if [ "$DRY_RUN" = true ]; then
  echo "Dry run completed. No changes were applied."
else
  echo "Database migrations completed successfully."
fi