#!/bin/bash
# CryoProtect Fedora Setup Script
# A comprehensive setup script for running CryoProtect on Fedora Linux

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Function to print section headers
print_section() {
    echo
    echo -e "${BLUE}=================================${NC}"
    echo -e "${BLUE}  $1${NC}"
    echo -e "${BLUE}=================================${NC}"
    echo
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

# Function to confirm actions
confirm() {
    read -p "$(echo -e $CYAN"$1 [y/N]: "$NC)" choice
    case "$choice" in
        y|Y|yes|Yes|YES ) return 0;;
        * ) return 1;;
    esac
}

# Welcome message
clear
echo -e "${BOLD}${BLUE}"
echo "  _____ _____          _____  _____   ____  _______ ______ _____ _______ "
echo " / ____|  __ \        |  __ \|  __ \ / __ \|__   __|  ____|  __ \__   __|"
echo "| |    | |__) |_   _  | |__) | |__) | |  | |  | |  | |__  | |  | | | |   "
echo "| |    |  _  /| | | | |  ___/|  _  /| |  | |  | |  |  __| | |  | | | |   "
echo "| |____| | \ \| |_| | | |    | | \ \| |__| |  | |  | |____| |__| | | |   "
echo " \_____|_|  \_\\__, | |_|    |_|  \_\\____/   |_|  |______|_____/  |_|   "
echo "                __/ |                                                     "
echo "               |___/          Fedora Linux Setup Guide                   "
echo -e "${NC}"
echo
echo "This script will guide you through setting up CryoProtect on Fedora Linux."
echo "It will help you configure SELinux, Podman, and test connectivity to Supabase."
echo
echo -e "${YELLOW}Note: Some commands may require sudo privileges.${NC}"
echo

# Check if running as root and warn if so
if [ "$EUID" -eq 0 ]; then
    print_warning "This script is running as root. It's recommended to run it as a regular user."
    if ! confirm "Continue anyway?"; then
        echo "Exiting..."
        exit 1
    fi
fi

# Main menu
while true; do
    print_section "MAIN MENU"
    echo "Select an option to proceed:"
    echo "1. Run environment readiness check"
    echo "2. Configure SELinux for application directories"
    echo "3. Create and configure IPv4-only Podman network"
    echo "4. Update Supabase credentials in .env file"
    echo "5. Test Supabase connectivity"
    echo "6. Test database connection"
    echo "7. Start application with Podman"
    echo "8. Exit"
    echo
    read -p "$(echo -e $CYAN"Enter your choice [1-8]: "$NC)" choice
    
    case $choice in
        1)
            print_section "ENVIRONMENT READINESS CHECK"
            if confirm "Run the environment readiness test?"; then
                echo "Running test_podman_readiness.sh..."
                ./test_podman_readiness.sh
            fi
            ;;
        2)
            print_section "SELINUX CONFIGURATION"
            if confirm "Configure SELinux contexts for application directories?"; then
                echo "This will set the proper SELinux contexts for CryoProtect directories."
                echo "Since this requires sudo privileges, you'll be prompted for your password."
                
                if confirm "Proceed with setting SELinux contexts?"; then
                    echo "Creating sudo-enabled setup_selinux_contexts.sh script..."
                    
                    # Create temporary script
                    cat > /tmp/setup_selinux_contexts.sh << 'EOF'
#!/bin/bash
# Setup SELinux contexts for CryoProtect directories

# Create required directories if they don't exist
mkdir -p logs
mkdir -p backup/data
mkdir -p data

# Set SELinux contexts
echo "Setting context for logs directory..."
chcon -Rt container_file_t ./logs

echo "Setting context for backup directory..."
chcon -Rt container_file_t ./backup

echo "Setting context for data directory..."
chcon -Rt container_file_t ./data

echo "Setting context for database directory..."
chcon -Rt container_file_t ./database 2>/dev/null || true

echo "Done setting SELinux contexts."
EOF
                    
                    chmod +x /tmp/setup_selinux_contexts.sh
                    echo "Running sudo /tmp/setup_selinux_contexts.sh..."
                    sudo /tmp/setup_selinux_contexts.sh
                    rm /tmp/setup_selinux_contexts.sh
                    
                    print_success "SELinux contexts configured successfully."
                fi
            fi
            ;;
        3)
            print_section "PODMAN NETWORK CONFIGURATION"
            if confirm "Create and configure IPv4-only Podman network?"; then
                echo "Checking existing Podman networks..."
                podman network ls
                
                # Check if cryoprotect-net already exists
                if podman network ls | grep -q "cryoprotect-net"; then
                    print_warning "A network named 'cryoprotect-net' already exists."
                    if confirm "Remove existing network and create a new one?"; then
                        podman network rm cryoprotect-net
                    else
                        continue
                    fi
                fi
                
                echo "Creating IPv4-only network 'cryoprotect-net'..."
                podman network create --ipv6=false cryoprotect-net
                
                # Update podman-compose.minimal.yml to use this network
                if [ -f "podman-compose.minimal.yml" ]; then
                    echo "Updating podman-compose.minimal.yml to use the cryoprotect-net network..."
                    sed -i 's/networks:\n  app-network:\n    driver: bridge/networks:\n  app-network:\n    external: true\n    name: cryoprotect-net/' podman-compose.minimal.yml
                    print_success "Podman network configuration updated."
                else
                    print_warning "podman-compose.minimal.yml not found. Network created but Compose file not updated."
                fi
            fi
            ;;
        4)
            print_section "UPDATE SUPABASE CREDENTIALS"
            if confirm "Update Supabase credentials in .env file?"; then
                echo "Running scripts/update_env.sh..."
                cd scripts && ./update_env.sh
                cd ..
            fi
            ;;
        5)
            print_section "TEST SUPABASE CONNECTIVITY"
            if confirm "Test connectivity to Supabase?"; then
                echo "This will test if your containers can successfully connect to Supabase."
                echo "Running scripts/test_supabase_connectivity.sh..."
                cd scripts && ./test_supabase_connectivity.sh
                cd ..
            fi
            ;;
        6)
            print_section "TEST DATABASE CONNECTION"
            if confirm "Test connection to Supabase PostgreSQL database?"; then
                echo "This will test direct connection to your Supabase PostgreSQL database."
                echo
                echo -e "${YELLOW}Note: You'll need to provide your database password.${NC}"
                echo
                
                if confirm "Continue with database connection test?"; then
                    python3 scripts/test_db_connection.py
                fi
            fi
            ;;
        7)
            print_section "START APPLICATION"
            if confirm "Start the CryoProtect application using Podman?"; then
                echo "This will start the application using the minimal podman-compose configuration."
                echo "The application will be available at http://localhost:5000"
                echo
                
                if confirm "Start application now?"; then
                    ./quickstart_podman.sh
                fi
            fi
            ;;
        8)
            echo "Exiting..."
            exit 0
            ;;
        *)
            print_error "Invalid option. Please select a number between 1 and 8."
            ;;
    esac
    
    echo
    read -p "Press Enter to continue..."
done