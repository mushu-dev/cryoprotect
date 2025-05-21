#!/bin/bash
# Setup SELinux contexts for CryoProtect on Fedora

BLUE='\033[0;34m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to print section headers
print_section() {
    echo -e "\n${BLUE}================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}================================================${NC}\n"
}

# Function to print success messages
print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

# Function to print warning messages
print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

# Function to print error messages
print_error() {
    echo -e "${RED}✗ $1${NC}"
}

# Check if script is run as root
if [ "$EUID" -ne 0 ]; then
    print_error "This script must be run as root (with sudo)"
    exit 1
fi

# Create required directories if they don't exist
print_section "Creating directories if needed"
mkdir -p logs
mkdir -p backup/data
mkdir -p database/verification/reports

# Set SELinux contexts for directories
print_section "Setting SELinux contexts for directories"

# Main data directories
echo "Setting context for logs directory..."
chcon -Rt container_file_t ./logs
print_success "Logs directory context set"

echo "Setting context for backup directory..."
chcon -Rt container_file_t ./backup
print_success "Backup directory context set"

echo "Setting context for data directory..."
chcon -Rt container_file_t ./data
print_success "Data directory context set"

echo "Setting context for database directory..."
chcon -Rt container_file_t ./database
print_success "Database directory context set"

# Check if contexts were applied correctly
print_section "Verifying SELinux contexts"

if [ "$(ls -lZ ./logs | grep container_file_t)" ]; then
    print_success "Logs directory context verified"
else
    print_error "Failed to set logs directory context"
fi

if [ "$(ls -lZ ./backup | grep container_file_t)" ]; then
    print_success "Backup directory context verified"
else
    print_error "Failed to set backup directory context"
fi

if [ "$(ls -lZ ./data | grep container_file_t)" ]; then
    print_success "Data directory context verified"
else
    print_error "Failed to set data directory context"
fi

if [ "$(ls -lZ ./database | grep container_file_t)" ]; then
    print_success "Database directory context verified"
else
    print_error "Failed to set database directory context"
fi

print_section "SELinux Context Setup Complete"
echo "Your directories now have the correct SELinux contexts for Podman."
echo "To run containers with these volumes, use the ':Z' suffix in your volume mounts."
echo "Example: -v ./logs:/app/logs:Z"