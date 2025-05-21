#!/bin/bash
# run_convex_migrations.sh
# Script to run Convex migrations with unified output format

# Load environment variables
if [ -f .env ]; then
  export $(grep -v '^#' .env | xargs)
else
  echo "Error: .env file not found"
  exit 1
fi

# Validate environment variables
if [ -z "$CONVEX_URL" ]; then
  echo "Error: CONVEX_URL must be set in .env file"
  exit 1
fi

# Set default parameters
ACTION="status"
TARGET_VERSION=""
DRY_RUN=false
FORMAT="text"
MIGRATION_NAME=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    status|apply|rollback|create)
      ACTION=$1
      shift
      ;;
    --target=*)
      TARGET_VERSION="${1#*=}"
      shift
      ;;
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    --format=*)
      FORMAT="${1#*=}"
      shift
      ;;
    --name=*)
      MIGRATION_NAME="${1#*=}"
      shift
      ;;
    --help)
      echo "Usage: ./run_convex_migrations.sh [action] [options]"
      echo ""
      echo "Actions:"
      echo "  status      Show migration status (default)"
      echo "  apply       Apply pending migrations"
      echo "  rollback    Roll back applied migrations"
      echo "  create      Create a new migration"
      echo ""
      echo "Options:"
      echo "  --target=N  Target migration version"
      echo "  --dry-run   Show changes without executing"
      echo "  --format=F  Output format (text or json)"
      echo "  --name=X    Name for new migration (required for create)"
      echo "  --help      Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

echo "CryoProtect Convex Migration Tool"
echo "================================="
echo "Action: $ACTION"
if [ ! -z "$TARGET_VERSION" ]; then
  echo "Target version: $TARGET_VERSION"
fi
if [ "$DRY_RUN" = true ]; then
  echo "Dry run: YES"
fi
echo ""

# Check if Node.js is installed
if ! command -v node &> /dev/null; then
  echo "Error: Node.js is required to run Convex migrations"
  exit 1
fi

# Check if ts-node is installed
if ! command -v ts-node &> /dev/null; then
  echo "Installing ts-node..."
  npm install -g ts-node typescript @types/node
fi

# Create log directory if it doesn't exist
LOG_DIR="migration_logs"
mkdir -p $LOG_DIR

# Create timestamp for log files
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Run the appropriate command
case $ACTION in
  status)
    LOG_FILE="$LOG_DIR/convex_status_$TIMESTAMP.log"
    echo "Getting migration status..."
    
    CMD="node --loader ts-node/esm ./convex/migrations/cli.ts status"
    if [ "$FORMAT" = "json" ]; then
      CMD="$CMD --format json"
    fi
    
    echo "$ $CMD" | tee -a $LOG_FILE
    eval $CMD | tee -a $LOG_FILE
    ;;
    
  apply)
    LOG_FILE="$LOG_DIR/convex_apply_$TIMESTAMP.log"
    echo "Applying migrations..."
    
    CMD="node --loader ts-node/esm ./convex/migrations/cli.ts apply"
    if [ ! -z "$TARGET_VERSION" ]; then
      CMD="$CMD --target $TARGET_VERSION"
    fi
    if [ "$DRY_RUN" = true ]; then
      CMD="$CMD --dry-run"
    fi
    
    echo "$ $CMD" | tee -a $LOG_FILE
    eval $CMD | tee -a $LOG_FILE
    ;;
    
  rollback)
    LOG_FILE="$LOG_DIR/convex_rollback_$TIMESTAMP.log"
    echo "Rolling back migrations..."
    
    CMD="node --loader ts-node/esm ./convex/migrations/cli.ts rollback"
    if [ ! -z "$TARGET_VERSION" ]; then
      CMD="$CMD --target $TARGET_VERSION"
    fi
    if [ "$DRY_RUN" = true ]; then
      CMD="$CMD --dry-run"
    fi
    
    echo "$ $CMD" | tee -a $LOG_FILE
    eval $CMD | tee -a $LOG_FILE
    ;;
    
  create)
    LOG_FILE="$LOG_DIR/convex_create_$TIMESTAMP.log"
    
    if [ -z "$MIGRATION_NAME" ]; then
      echo "Error: --name is required for create action"
      echo "Usage: ./run_convex_migrations.sh create --name=migration_name"
      exit 1
    fi
    
    echo "Creating migration '$MIGRATION_NAME'..."
    
    CMD="node --loader ts-node/esm ./convex/migrations/cli.ts create \"$MIGRATION_NAME\""
    
    echo "$ $CMD" | tee -a $LOG_FILE
    eval $CMD | tee -a $LOG_FILE
    ;;
esac

echo ""
echo "Log saved to: $LOG_FILE"