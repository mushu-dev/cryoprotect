#!/bin/bash
# Script to run the full ChEMBL import with integrated fixes

# Set the environment
export PYTHONPATH=./:$PYTHONPATH

# Create timestamp for logging
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="logs/chembl_full_import_${TIMESTAMP}.log"

# Create logs directory if it doesn't exist
mkdir -p logs
mkdir -p reports
mkdir -p checkpoints

# Function to display help
function show_help {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --limit NUMBER            Maximum number of compounds to import (default: 1000)"
    echo "  --batch-size NUMBER       Batch size for import operations (default: 50)"
    echo "  --fix-batch-size NUMBER   Batch size for fix operations (default: 50)" 
    echo "  --fix-limit NUMBER        Maximum number of molecules to fix (default: all)"
    echo "  --dry-run                 Don't actually insert data, just simulate"
    echo "  --skip-fixes              Skip applying fixes after import"
    echo "  --properties-only         Only apply property fixes"
    echo "  --cross-refs-only         Only apply cross-reference fixes"
    echo "  --help                    Show this help message and exit"
    echo ""
    echo "Examples:"
    echo "  $0 --limit 100 --batch-size 10"
    echo "  $0 --dry-run"
    echo "  $0 --skip-fixes"
    exit 0
}

# Process command line arguments
LIMIT=1000
BATCH_SIZE=50
FIX_BATCH_SIZE=50
FIX_LIMIT=""
DRY_RUN=""
APPLY_FIXES="--apply-fixes"
PROPERTIES_ONLY=""
CROSS_REFS_ONLY=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --limit)
            LIMIT="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --fix-batch-size)
            FIX_BATCH_SIZE="$2"
            shift 2
            ;;
        --fix-limit)
            FIX_LIMIT="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN="--dry-run"
            shift
            ;;
        --skip-fixes)
            APPLY_FIXES=""
            shift
            ;;
        --properties-only)
            PROPERTIES_ONLY="--properties-only"
            shift
            ;;
        --cross-refs-only)
            CROSS_REFS_ONLY="--cross-refs-only"
            shift
            ;;
        --help)
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# Build the fix limit argument if provided
if [ -n "$FIX_LIMIT" ]; then
    FIX_LIMIT_ARG="--fix-limit $FIX_LIMIT"
else
    FIX_LIMIT_ARG=""
fi

# Display the configuration
echo "Running ChEMBL full import with the following configuration:"
echo "- Limit: $LIMIT compounds"
echo "- Batch size: $BATCH_SIZE"
echo "- Fix batch size: $FIX_BATCH_SIZE"
if [ -n "$FIX_LIMIT" ]; then
    echo "- Fix limit: $FIX_LIMIT molecules"
else
    echo "- Fix limit: all molecules"
fi
if [ -n "$DRY_RUN" ]; then
    echo "- Dry run: yes (no data will be inserted)"
else
    echo "- Dry run: no"
fi
if [ -n "$APPLY_FIXES" ]; then
    echo "- Apply fixes: yes"
    if [ -n "$PROPERTIES_ONLY" ]; then
        echo "  - Properties only: yes"
    elif [ -n "$CROSS_REFS_ONLY" ]; then
        echo "  - Cross-references only: yes"
    else
        echo "  - Both properties and cross-references"
    fi
else
    echo "- Apply fixes: no"
fi
echo "- Log file: $LOG_FILE"
echo ""

# Confirm with the user
read -p "Do you want to continue with these settings? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Import cancelled."
    exit 1
fi

# Run the import
echo "Starting ChEMBL import at $(date)"
echo "This may take a while. You can monitor the progress in $LOG_FILE"

python ChEMBL_Integrated_Import.py --limit $LIMIT --batch-size $BATCH_SIZE --fix-batch-size $FIX_BATCH_SIZE $FIX_LIMIT_ARG $DRY_RUN $APPLY_FIXES 2>&1 | tee -a "$LOG_FILE"

# If fixes were skipped but we want to run them separately
if [ -z "$APPLY_FIXES" ] && [ -z "$DRY_RUN" ]; then
    echo ""
    read -p "Do you want to run fixes separately now? (y/n) " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Running fixes separately at $(date)"
        
        FIX_ARGS="--batch-size $FIX_BATCH_SIZE"
        if [ -n "$FIX_LIMIT" ]; then
            FIX_ARGS="$FIX_ARGS --limit $FIX_LIMIT"
        fi
        if [ -n "$PROPERTIES_ONLY" ]; then
            FIX_ARGS="$FIX_ARGS $PROPERTIES_ONLY"
        fi
        if [ -n "$CROSS_REFS_ONLY" ]; then
            FIX_ARGS="$FIX_ARGS $CROSS_REFS_ONLY"
        fi
        
        python integrated_chembl_import_fix.py $FIX_ARGS 2>&1 | tee -a "$LOG_FILE"
    fi
fi

echo ""
echo "ChEMBL import process completed at $(date)"
echo "See $LOG_FILE for details"

# Run verification if not in dry run mode
if [ -z "$DRY_RUN" ]; then
    echo ""
    read -p "Do you want to run verification on the imported data? (y/n) " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Running verification at $(date)"
        python verify_chembl_import_data.py 2>&1 | tee -a "$LOG_FILE"
    fi
fi