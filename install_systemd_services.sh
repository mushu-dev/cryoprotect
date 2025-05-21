#!/bin/bash
# install_systemd_services.sh
#
# This script installs the CryoProtect systemd services for automatic startup.
#
# Usage:
#   ./install_systemd_services.sh [--user|--system]

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default to user mode
USER_MODE=true

# Parse command line arguments
for arg in "$@"; do
    case $arg in
        --user)
            USER_MODE=true
            shift
            ;;
        --system)
            USER_MODE=false
            shift
            ;;
        *)
            # Unknown option
            shift
            ;;
    esac
done

echo -e "${BLUE}======================================"
echo "CryoProtect Systemd Services Installation"
echo -e "======================================${NC}"

# Get current directory
CURRENT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$CURRENT_DIR"

# Service files to install
SERVICE_FILES=(
    "cryoprotect-postgres.service"
    "cryoprotect-app.service"
    "cryoprotect-db-integrity-check.service"
    "cryoprotect-db-integrity-check.timer"
    "cryoprotect-postgres-health.service"
    "cryoprotect-postgres-health.timer"
)

# Create directory for user systemd services if needed
if [ "$USER_MODE" = true ]; then
    mkdir -p ~/.config/systemd/user
    echo -e "${BLUE}Installing services in user mode${NC}"
    SYSTEMD_DIR=~/.config/systemd/user
    SYSTEMCTL_CMD="systemctl --user"
else
    echo -e "${BLUE}Installing services in system mode${NC}"
    SYSTEMD_DIR=/etc/systemd/system
    SYSTEMCTL_CMD="sudo systemctl"
    
    # Check if running as root when using system mode
    if [ "$EUID" -ne 0 ]; then
        echo -e "${YELLOW}Warning: Not running as root. Using sudo for system-level operations.${NC}"
    fi
fi

# Copy service files to systemd directory
echo "Copying service files to $SYSTEMD_DIR..."
for service_file in "${SERVICE_FILES[@]}"; do
    if [ -f "$service_file" ]; then
        echo "Installing $service_file..."
        if [ "$USER_MODE" = true ]; then
            cp "$service_file" "$SYSTEMD_DIR/"
        else
            sudo cp "$service_file" "$SYSTEMD_DIR/"
        fi
    else
        echo -e "${YELLOW}Warning: Service file $service_file not found!${NC}"
    fi
done

# Make script files executable
echo "Making scripts executable..."
chmod +x run_optimized_postgres.sh
chmod +x check_postgres_health.py
if [ -f run_verify_all_data_integrity.sh ]; then
    chmod +x run_verify_all_data_integrity.sh
fi

# Reload systemd configuration
echo "Reloading systemd configuration..."
if [ "$USER_MODE" = true ]; then
    systemctl --user daemon-reload
else
    sudo systemctl daemon-reload
fi

# Enable and start services
echo "Enabling and starting services..."

# Skip PostgreSQL service setup - using existing database
echo "Skipping PostgreSQL service setup - using existing database..."

# Main application service
echo "Enabling main application service..."
$SYSTEMCTL_CMD enable cryoprotect-app.service
echo "Starting main application service..."
$SYSTEMCTL_CMD start cryoprotect-app.service

# Integrity check timer
echo "Enabling integrity check timer..."
$SYSTEMCTL_CMD enable cryoprotect-db-integrity-check.timer
echo "Starting integrity check timer..."
$SYSTEMCTL_CMD start cryoprotect-db-integrity-check.timer

# PostgreSQL health check timer
echo "Enabling PostgreSQL health check timer..."
$SYSTEMCTL_CMD enable cryoprotect-postgres-health.timer
echo "Starting PostgreSQL health check timer..."
$SYSTEMCTL_CMD start cryoprotect-postgres-health.timer

# Show status of services
echo -e "\n${BLUE}Service Status:${NC}"
echo -e "${GREEN}PostgreSQL Service:${NC}"
$SYSTEMCTL_CMD status cryoprotect-postgres.service --no-pager

echo -e "\n${GREEN}Main Application Service:${NC}"
$SYSTEMCTL_CMD status cryoprotect-app.service --no-pager

echo -e "\n${GREEN}Database Integrity Check Timer:${NC}"
$SYSTEMCTL_CMD status cryoprotect-db-integrity-check.timer --no-pager

echo -e "\n${GREEN}PostgreSQL Health Check Timer:${NC}"
$SYSTEMCTL_CMD status cryoprotect-postgres-health.timer --no-pager

echo -e "\n${BLUE}======================================"
echo "Installation Complete!"
echo -e "======================================${NC}"

echo -e "${GREEN}✓${NC} Installed systemd services in ${USER_MODE:+user mode}${USER_MODE:+system mode}"
echo -e "${GREEN}✓${NC} Services will start automatically on system boot"

echo -e "\n${BLUE}Usage:${NC}"
echo "- Check service status: $SYSTEMCTL_CMD status cryoprotect-app.service"
echo "- Start a service: $SYSTEMCTL_CMD start cryoprotect-app.service"
echo "- Stop a service: $SYSTEMCTL_CMD stop cryoprotect-app.service"
echo "- Restart a service: $SYSTEMCTL_CMD restart cryoprotect-app.service"
echo "- View logs: $SYSTEMCTL_CMD ${USER_MODE:+--user }journalctl -u cryoprotect-app.service"

# Update database connection environment variables
echo -e "${BLUE}Setting environment variables for database connection${NC}"
if [ "$USER_MODE" = true ]; then
  mkdir -p ~/.config/environment.d/
  echo "SUPABASE_DB_PORT=5433" > ~/.config/environment.d/cryoprotect.conf
  echo -e "${GREEN}✓${NC} Added environment variables to ~/.config/environment.d/cryoprotect.conf"
else
  sudo bash -c 'echo "SUPABASE_DB_PORT=5433" > /etc/environment.d/cryoprotect.conf'
  echo -e "${GREEN}✓${NC} Added environment variables to /etc/environment.d/cryoprotect.conf"
fi