#!/bin/bash
# apply_migrations.sh
# Simple script to apply migrations sequentially to Supabase

# Load environment variables
if [ -f .env ]; then
  export $(grep -v '^#' .env | xargs)
else
  echo "Error: .env file not found"
  exit 1
fi

# Validate environment variables
if [ -z "$SUPABASE_URL" ] || [ -z "$SUPABASE_SERVICE_KEY" ]; then
  echo "Error: SUPABASE_URL and SUPABASE_SERVICE_KEY must be set in .env file"
  exit 1
fi

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
    --help)
      echo "Usage: ./apply_migrations.sh [options]"
      echo "Options:"
      echo "  --start=N      Start with migration N (default: 1)"
      echo "  --end=N        End with migration N (default: 19)"
      echo "  --dry-run      Show SQL without executing"
      echo "  --no-verify    Skip verification steps"
      echo "  --help         Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

echo "CryoProtect Database Migration Tool"
echo "=================================="
echo "Applying migrations from $START_INDEX to $END_INDEX"
echo "Dry run: $DRY_RUN"
echo "Verify: $VERIFY"
echo ""

# Create log directory if it doesn't exist
LOG_DIR="migration_logs"
mkdir -p $LOG_DIR

# Create timestamp for log files
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
SUMMARY_LOG="$LOG_DIR/migration_summary_$TIMESTAMP.log"

# Initialize summary log
echo "Migration Summary ($TIMESTAMP)" > $SUMMARY_LOG
echo "===============================" >> $SUMMARY_LOG
echo "" >> $SUMMARY_LOG

# Function to apply a migration
apply_migration() {
  local index=$1
  local padded_index=$(printf "%03d" $index)
  local migration_file=$(find migrations -name "${padded_index}_*.sql" | sort | head -n 1)
  
  if [ -z "$migration_file" ]; then
    echo "No migration found with index $padded_index"
    echo "Migration $padded_index: Not found" >> $SUMMARY_LOG
    return 1
  fi
  
  # Get migration name for logging
  local migration_name=$(basename $migration_file)
  local detail_log="$LOG_DIR/migration_${padded_index}_${TIMESTAMP}.log"
  
  echo "Applying migration: $migration_name"
  echo "Migration $padded_index: $migration_name" >> $SUMMARY_LOG
  
  # Log migration SQL
  echo "SQL for migration $padded_index:" > $detail_log
  cat $migration_file >> $detail_log
  echo "" >> $detail_log
  
  if [ "$DRY_RUN" = true ]; then
    echo "DRY RUN: Would apply $migration_name"
    echo "DRY RUN: Would apply migration" >> $SUMMARY_LOG
    return 0
  fi
  
  # Read migration file and prepare SQL
  MIGRATION_SQL=$(cat $migration_file | tr -d '\n' | sed 's/"/\\"/g')
  
  # Create request to apply migration
  RESPONSE=$(curl -s -X POST \
    "${SUPABASE_URL}/rest/v1/rpc/exec_sql" \
    -H "apikey: ${SUPABASE_SERVICE_KEY}" \
    -H "Authorization: Bearer ${SUPABASE_SERVICE_KEY}" \
    -H "Content-Type: application/json" \
    -d "{\"query\": \"${MIGRATION_SQL}\"}")
  
  # Check for errors in response
  if [[ $RESPONSE == *"error"* ]]; then
    echo "Error applying migration $padded_index:"
    echo $RESPONSE | jq 2>/dev/null || echo $RESPONSE
    echo "Migration $padded_index: FAILED" >> $SUMMARY_LOG
    echo "Error: $RESPONSE" >> $detail_log
    return 1
  fi
  
  echo "Migration $padded_index applied successfully."
  echo "Migration $padded_index: SUCCESS" >> $SUMMARY_LOG
  echo "Response: $RESPONSE" >> $detail_log
  return 0
}

# Function to verify a migration
verify_migration() {
  local index=$1
  local padded_index=$(printf "%03d" $index)
  local detail_log="$LOG_DIR/verification_${padded_index}_${TIMESTAMP}.log"
  
  echo "Verifying migration $padded_index..."
  echo "Verification $padded_index: STARTED" >> $SUMMARY_LOG
  
  # Create verification log
  echo "Verification for migration $padded_index:" > $detail_log
  
  # Switch based on migration index to perform appropriate verification
  case $padded_index in
    "001")
      # Verify initial schema
      echo "Verifying initial schema..." | tee -a $detail_log
      
      # Check if users table exists
      RESPONSE=$(curl -s -X GET \
        "${SUPABASE_URL}/rest/v1/users?limit=1" \
        -H "apikey: ${SUPABASE_KEY}" \
        -H "Authorization: Bearer ${SUPABASE_KEY}")
      
      if [[ $RESPONSE != *"error"* ]]; then
        echo "✓ Users table exists" | tee -a $detail_log
        echo "Verification $padded_index: SUCCESS" >> $SUMMARY_LOG
      else
        echo "✗ Users table not found" | tee -a $detail_log
        echo "Verification $padded_index: FAILED" >> $SUMMARY_LOG
        return 1
      fi
      ;;
      
    "002")
      # Verify projects schema
      echo "Verifying projects schema..." | tee -a $detail_log
      
      # Check if projects table exists
      RESPONSE=$(curl -s -X GET \
        "${SUPABASE_URL}/rest/v1/projects?limit=1" \
        -H "apikey: ${SUPABASE_KEY}" \
        -H "Authorization: Bearer ${SUPABASE_KEY}")
      
      if [[ $RESPONSE != *"error"* ]]; then
        echo "✓ Projects table exists" | tee -a $detail_log
        echo "Verification $padded_index: SUCCESS" >> $SUMMARY_LOG
      else
        echo "✗ Projects table not found" | tee -a $detail_log
        echo "Verification $padded_index: FAILED" >> $SUMMARY_LOG
        return 1
      fi
      ;;
      
    "003")
      # Verify teams schema
      echo "Verifying teams schema..." | tee -a $detail_log
      
      # Check if teams table exists
      RESPONSE=$(curl -s -X GET \
        "${SUPABASE_URL}/rest/v1/teams?limit=1" \
        -H "apikey: ${SUPABASE_KEY}" \
        -H "Authorization: Bearer ${SUPABASE_KEY}")
      
      if [[ $RESPONSE != *"error"* ]]; then
        echo "✓ Teams table exists" | tee -a $detail_log
        echo "Verification $padded_index: SUCCESS" >> $SUMMARY_LOG
      else
        echo "✗ Teams table not found" | tee -a $detail_log
        echo "Verification $padded_index: FAILED" >> $SUMMARY_LOG
        return 1
      fi
      ;;
      
    "005")
      # Verify scientific schema
      echo "Verifying scientific schema..." | tee -a $detail_log
      
      # Check if molecules table exists
      RESPONSE=$(curl -s -X GET \
        "${SUPABASE_URL}/rest/v1/molecules?limit=1" \
        -H "apikey: ${SUPABASE_KEY}" \
        -H "Authorization: Bearer ${SUPABASE_KEY}")
      
      if [[ $RESPONSE != *"error"* ]]; then
        echo "✓ Molecules table exists" | tee -a $detail_log
        echo "Verification $padded_index: SUCCESS" >> $SUMMARY_LOG
      else
        echo "✗ Molecules table not found" | tee -a $detail_log
        echo "Verification $padded_index: FAILED" >> $SUMMARY_LOG
        return 1
      fi
      ;;
      
    "007"|"008")
      # Verify seed data
      echo "Verifying seed data..." | tee -a $detail_log
      
      # Check if molecules table has data
      RESPONSE=$(curl -s -X GET \
        "${SUPABASE_URL}/rest/v1/molecules?limit=1" \
        -H "apikey: ${SUPABASE_KEY}" \
        -H "Authorization: Bearer ${SUPABASE_KEY}")
      
      if [[ $RESPONSE != "[]" ]]; then
        echo "✓ Molecules table has data" | tee -a $detail_log
        echo "Verification $padded_index: SUCCESS" >> $SUMMARY_LOG
      else
        echo "✗ Molecules table is empty" | tee -a $detail_log
        echo "Verification $padded_index: FAILED" >> $SUMMARY_LOG
        return 1
      fi
      ;;
      
    *)
      # Generic verification
      echo "No specific verification for migration $padded_index" | tee -a $detail_log
      echo "Verification $padded_index: SKIPPED" >> $SUMMARY_LOG
      ;;
  esac
  
  return 0
}

# Apply migrations in sequence
FAILED=false

for i in $(seq $START_INDEX $END_INDEX); do
  apply_migration $i
  
  if [ $? -ne 0 ]; then
    echo "Failed to apply migration $i. Stopping."
    FAILED=true
    break
  fi
  
  if [ "$VERIFY" = true ]; then
    verify_migration $i
    
    if [ $? -ne 0 ]; then
      echo "Failed to verify migration $i. Stopping."
      FAILED=true
      break
    fi
  fi
  
  # Sleep briefly between migrations
  sleep 1
done

# Print summary
echo ""
echo "Migration Summary"
echo "================="

if [ "$FAILED" = true ]; then
  echo "Migration process FAILED. See logs for details."
  echo "Migration process FAILED" >> $SUMMARY_LOG
else
  echo "Migration process COMPLETED successfully."
  echo "Migration process COMPLETED successfully" >> $SUMMARY_LOG
fi

echo "Log files:"
echo "- Summary: $SUMMARY_LOG"
echo "- Details: $LOG_DIR/migration_*_$TIMESTAMP.log"
echo "- Verification: $LOG_DIR/verification_*_$TIMESTAMP.log"

# Exit with appropriate status code
if [ "$FAILED" = true ]; then
  exit 1
else
  exit 0
fi