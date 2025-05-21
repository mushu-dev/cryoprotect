#!/bin/bash
# Database maintenance script for CryoProtect
# This script runs maintenance tasks on the database

# Set working directory to the script directory
cd "$(dirname "$0")"

# Create log directory if it doesn't exist
mkdir -p logs/maintenance

# Generate timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="logs/maintenance/maintenance_${TIMESTAMP}.log"

# Function to log messages
log() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] $1" | tee -a "$LOG_FILE"
}

log "Starting database maintenance..."

# Check if this is a full maintenance run
# Full maintenance runs on weekends, standard runs on weekdays
if [[ $(date +"%u") -ge 6 ]]; then
    FULL_MAINTENANCE="--full"
    log "Running FULL maintenance (weekend schedule)"
else
    FULL_MAINTENANCE=""
    log "Running standard maintenance (weekday schedule)"
fi

# Run the maintenance script
log "Executing maintenance tasks..."
python3 -m database.maintenance $FULL_MAINTENANCE --report --report-file="logs/maintenance/report_${TIMESTAMP}.txt" 2>&1 | tee -a "$LOG_FILE"
RESULT=$?

if [ $RESULT -eq 0 ]; then
    log "Maintenance completed successfully"
else
    log "Maintenance failed with exit code $RESULT"
    
    # Send an alert (example: you could add email notification here)
    echo "Database maintenance failed. Check the logs at $LOG_FILE" > logs/maintenance/alert_${TIMESTAMP}.txt
fi

# Also clean up old log files (keep last 30 days)
find logs/maintenance -name "*.log" -type f -mtime +30 -delete
find logs/maintenance -name "*.txt" -type f -mtime +30 -delete

log "Maintenance script completed"
exit $RESULT