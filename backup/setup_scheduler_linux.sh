#!/bin/bash
# Setup script for automated backup scheduling on Linux/macOS using cron

# Get the absolute path to the project directory
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BACKUP_SCRIPT="$PROJECT_DIR/backup/backup_manager.py"
CONFIG_FILE="$PROJECT_DIR/backup/backup_config.json"
LOG_FILE="$PROJECT_DIR/logs/backup.log"

# Ensure the logs directory exists
mkdir -p "$PROJECT_DIR/logs"

# Create the cron job command
CRON_CMD="cd $PROJECT_DIR && python -m backup.backup_manager backup --config $CONFIG_FILE >> $LOG_FILE 2>&1"

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Config file not found. Creating from template..."
    if [ -f "$PROJECT_DIR/backup/backup_config.json.template" ]; then
        cp "$PROJECT_DIR/backup/backup_config.json.template" "$CONFIG_FILE"
        echo "Created config file from template. Please edit $CONFIG_FILE to customize your backup settings."
    else
        echo "Error: Template config file not found."
        exit 1
    fi
fi

# Get the backup schedule from the config file
DAILY_TIME=$(grep -o '"daily": *"[^"]*"' "$CONFIG_FILE" | cut -d'"' -f4)
WEEKLY_DAY=$(grep -o '"weekly": *"[^"]*"' "$CONFIG_FILE" | cut -d'"' -f4)
MONTHLY_DAY=$(grep -o '"monthly": *[0-9]*' "$CONFIG_FILE" | grep -o '[0-9]*')

# Parse the daily time
HOUR=$(echo $DAILY_TIME | cut -d':' -f1)
MINUTE=$(echo $DAILY_TIME | cut -d':' -f2)

# Map day of week to cron format (0-6, where 0 is Sunday)
case $WEEKLY_DAY in
    "sunday") DOW=0 ;;
    "monday") DOW=1 ;;
    "tuesday") DOW=2 ;;
    "wednesday") DOW=3 ;;
    "thursday") DOW=4 ;;
    "friday") DOW=5 ;;
    "saturday") DOW=6 ;;
    *) DOW=0 ;; # Default to Sunday
esac

# Create the cron entries
DAILY_CRON="$MINUTE $HOUR * * * $CRON_CMD --type daily"
WEEKLY_CRON="$MINUTE $HOUR * * $DOW $CRON_CMD --type weekly"
MONTHLY_CRON="$MINUTE $HOUR $MONTHLY_DAY * * $CRON_CMD --type monthly"

# Create a temporary file with the current crontab
crontab -l > /tmp/crontab.tmp 2>/dev/null || echo "" > /tmp/crontab.tmp

# Check if the cron jobs already exist
if grep -q "$BACKUP_SCRIPT" /tmp/crontab.tmp; then
    echo "Backup cron jobs already exist. Updating..."
    # Remove existing backup cron jobs
    grep -v "$BACKUP_SCRIPT" /tmp/crontab.tmp > /tmp/crontab.new
    mv /tmp/crontab.new /tmp/crontab.tmp
fi

# Add the new cron jobs
echo "# CryoProtect v2 Backup System" >> /tmp/crontab.tmp
echo "$DAILY_CRON" >> /tmp/crontab.tmp
echo "$WEEKLY_CRON" >> /tmp/crontab.tmp
echo "$MONTHLY_CRON" >> /tmp/crontab.tmp

# Install the new crontab
crontab /tmp/crontab.tmp
rm /tmp/crontab.tmp

echo "Backup scheduler setup complete."
echo "Daily backup: $MINUTE $HOUR * * * (every day at $DAILY_TIME)"
echo "Weekly backup: $MINUTE $HOUR * * $DOW (every $WEEKLY_DAY at $DAILY_TIME)"
echo "Monthly backup: $MINUTE $HOUR $MONTHLY_DAY * * (every month on day $MONTHLY_DAY at $DAILY_TIME)"
echo "Logs will be written to: $LOG_FILE"