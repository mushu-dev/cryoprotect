#!/bin/bash
#
# ChEMBL Import Runner Script
#
# This script runs the optimized ChEMBL import process with proper environment setup,
# logging, and error handling.
#

# Ensure script exits on error
set -e

# Configuration
LOG_DIR="logs"
LOG_FILE="$LOG_DIR/chembl_import_$(date +%Y%m%d_%H%M%S).log"
CHECKPOINT_DIR="checkpoints"
CHECKPOINT_FILE="$CHECKPOINT_DIR/chembl_import_checkpoint.json"

# Default parameters (can be overridden with command line arguments)
LIMIT=10000
BATCH_SIZE=100
DRY_RUN=false
MAX_WORKERS=4

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --limit)
      LIMIT="$2"
      shift 2
      ;;
    --batch-size)
      BATCH_SIZE="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    --max-workers)
      MAX_WORKERS="$2"
      shift 2
      ;;
    --checkpoint-file)
      CHECKPOINT_FILE="$2"
      shift 2
      ;;
    --search-terms)
      SEARCH_TERMS="$2"
      shift 2
      ;;
    --help)
      echo "Usage: $0 [options]"
      echo
      echo "Options:"
      echo "  --limit N             Maximum number of compounds to import (default: 10000)"
      echo "  --batch-size N        Batch size for database operations (default: 100)"
      echo "  --dry-run             Don't actually insert data, just simulate"
      echo "  --max-workers N       Maximum number of worker threads/processes (default: 4)"
      echo "  --checkpoint-file F   Path to checkpoint file for resumable imports"
      echo "  --search-terms S      Comma-separated list of search terms"
      echo "  --help                Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Create checkpoint directory if it doesn't exist
mkdir -p "$CHECKPOINT_DIR"

# Print configuration
echo "=== ChEMBL Import Configuration ==="
echo "Limit: $LIMIT"
echo "Batch Size: $BATCH_SIZE"
echo "Dry Run: $DRY_RUN"
echo "Max Workers: $MAX_WORKERS"
echo "Checkpoint File: $CHECKPOINT_FILE"
echo "Search Terms: $SEARCH_TERMS"
echo "Log File: $LOG_FILE"
echo "====================================="

# Build command
CMD="python chembl_direct_import.py --limit $LIMIT --batch-size $BATCH_SIZE --max-workers $MAX_WORKERS --checkpoint-file $CHECKPOINT_FILE"

if [ "$DRY_RUN" = true ]; then
  CMD="$CMD --dry-run"
fi

if [ ! -z "$SEARCH_TERMS" ]; then
  CMD="$CMD --search-terms \"$SEARCH_TERMS\""
fi

# Start import process
echo "Starting ChEMBL import process at $(date)"
echo "Command: $CMD"
echo "Logging to $LOG_FILE"

# Run the import process and log output
$CMD 2>&1 | tee "$LOG_FILE"

# Check exit status
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo "ChEMBL import completed successfully at $(date)"
  exit 0
else
  echo "ChEMBL import failed at $(date)"
  exit 1
fi