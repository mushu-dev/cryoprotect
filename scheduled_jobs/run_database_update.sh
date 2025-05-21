#!/bin/bash

# Run Database Update Job
# This script is designed to be executed by a scheduler (cron)

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
LOG_DIR="${SCRIPT_DIR}/logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/db_update_run_${TIMESTAMP}.log"

# Ensure log directory exists
mkdir -p "${LOG_DIR}"

# Log start time
echo "Starting database update job at $(date)" | tee -a "${LOG_FILE}"

# Activate virtual environment if needed
if [ -d "${SCRIPT_DIR}/../venv" ]; then
    echo "Activating virtual environment" | tee -a "${LOG_FILE}"
    source "${SCRIPT_DIR}/../venv/bin/activate"
fi

# Run the database update script
echo "Running database update script" | tee -a "${LOG_FILE}"
python3 "${SCRIPT_DIR}/database_update.py" 2>&1 | tee -a "${LOG_FILE}"

# Check exit status
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo "Database update completed successfully at $(date)" | tee -a "${LOG_FILE}"
    exit 0
else
    echo "Database update failed at $(date)" | tee -a "${LOG_FILE}"
    
    # Send notification of failure
    # In a production environment, this would send an email or Slack notification
    echo "Sending failure notification" | tee -a "${LOG_FILE}"
    
    exit 1
fi