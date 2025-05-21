#!/bin/bash
# Setup cron jobs for database maintenance and cache processing

# Set working directory to the script directory
cd "$(dirname "$0")"

# Make sure the scripts are executable
chmod +x run_database_maintenance.sh
chmod +x run_cache_processor.sh

# Get the absolute path of the scripts
SCRIPT_DIR=$(pwd)
MAINTENANCE_SCRIPT="${SCRIPT_DIR}/run_database_maintenance.sh"
CACHE_PROCESSOR_SCRIPT="${SCRIPT_DIR}/run_cache_processor.sh"

# Create the cron entries
# Run database maintenance daily at 3:00 AM
MAINTENANCE_CRON="0 3 * * * ${MAINTENANCE_SCRIPT} >> ${SCRIPT_DIR}/logs/cron_maintenance.log 2>&1"

# Check if the cron jobs already exist
EXISTING_CRONS=$(crontab -l 2>/dev/null || echo "")

if echo "$EXISTING_CRONS" | grep -q "run_database_maintenance.sh"; then
    echo "Database maintenance cron job already exists."
else
    # Add the cron jobs to the current user's crontab
    (echo "$EXISTING_CRONS"; echo "$MAINTENANCE_CRON") | crontab -
    echo "Database maintenance cron job has been added."
fi

# Create the systemd service file for cache processor
echo "Creating systemd service for cache processor..."

SYSTEMD_SERVICE_FILE="${SCRIPT_DIR}/cryoprotect-cache-processor.service"
cat > "$SYSTEMD_SERVICE_FILE" << EOF
[Unit]
Description=CryoProtect Cache Processor Service
After=network.target postgresql.service redis.service

[Service]
Type=simple
User=$(whoami)
WorkingDirectory=${SCRIPT_DIR}
ExecStart=${SCRIPT_DIR}/run_cache_processor.sh
ExecStop=${SCRIPT_DIR}/stop_cache_processor.sh
Restart=on-failure
RestartSec=10
StandardOutput=append:${SCRIPT_DIR}/logs/cache_processor_service.log
StandardError=append:${SCRIPT_DIR}/logs/cache_processor_service.log

[Install]
WantedBy=multi-user.target
EOF

echo "To install the systemd service, run the following commands as root:"
echo "  sudo cp \"$SYSTEMD_SERVICE_FILE\" /etc/systemd/system/"
echo "  sudo systemctl daemon-reload"
echo "  sudo systemctl enable cryoprotect-cache-processor.service"
echo "  sudo systemctl start cryoprotect-cache-processor.service"

echo "Setup completed!"