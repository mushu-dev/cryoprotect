#\!/bin/bash
# Script to run the unified ChEMBL import

# Set the environment
export PYTHONPATH=./:$PYTHONPATH

# Create timestamp for logging
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="logs/unified_chembl_import_${TIMESTAMP}.log"

# Create logs directory if it doesn't exist
mkdir -p logs
mkdir -p reports
mkdir -p checkpoints

# Function to display help
function show_help {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --limit NUMBER            Maximum number of compounds to import (default: 5000)"
    echo "  --batch-size NUMBER       Batch size for import operations (default: 50)"
    echo "  --dry-run                 Don't actually insert data, just simulate"
    echo "  --verify-only             Only verify existing data, don't import"
    echo "  --verify-limit NUMBER     Maximum number of molecules to verify (default: all)"
    echo "  --help                    Show this help message and exit"
    echo ""
    echo "Examples:"
    echo "  $0 --limit 100 --batch-size 10"
    echo "  $0 --dry-run"
    echo "  $0 --verify-only"
    exit 0
}

# Process command line arguments
LIMIT=5000
BATCH_SIZE=50
DRY_RUN=""
VERIFY_ONLY=""
VERIFY_LIMIT=""

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
        --verify-limit)
            VERIFY_LIMIT="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN="--dry-run"
            shift
            ;;
        --verify-only)
            VERIFY_ONLY="--verify-only"
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

# Build the verify limit argument if provided
if [ -n "$VERIFY_LIMIT" ]; then
    VERIFY_LIMIT_ARG="--verify-limit $VERIFY_LIMIT"
else
    VERIFY_LIMIT_ARG=""
fi

# Display the configuration
echo "Running Unified ChEMBL Import with the following configuration:"
echo "- Limit: $LIMIT compounds"
echo "- Batch size: $BATCH_SIZE"
if [ -n "$DRY_RUN" ]; then
    echo "- Dry run: yes (no data will be inserted)"
else
    echo "- Dry run: no"
fi
if [ -n "$VERIFY_ONLY" ]; then
    echo "- Verify only: yes (no data will be imported)"
    if [ -n "$VERIFY_LIMIT" ]; then
        echo "- Verify limit: $VERIFY_LIMIT molecules"
    else
        echo "- Verify limit: all molecules"
    fi
else
    echo "- Verify only: no"
fi
echo "- Log file: $LOG_FILE"
echo ""

# Skip confirmation for automated tests if AUTO_CONFIRM is set
if [[ -z "$AUTO_CONFIRM" ]]; then
    # Confirm with the user
    read -p "Do you want to continue with these settings? (y/n) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Import cancelled."
        exit 1
    fi
else
    echo "Auto-confirm enabled, proceeding without confirmation"
fi

# Run the import
echo "Starting Unified ChEMBL Import at $(date)"
echo "This may take a while. You can monitor the progress in $LOG_FILE"

python unified_chembl_import.py --limit $LIMIT --batch-size $BATCH_SIZE $DRY_RUN $VERIFY_ONLY $VERIFY_LIMIT_ARG 2>&1  < /dev/null |  tee -a "$LOG_FILE"

echo ""
echo "ChEMBL import process completed at $(date)"
echo "See $LOG_FILE for details"
