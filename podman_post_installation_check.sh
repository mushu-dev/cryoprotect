#!/bin/bash
#
# podman_post_installation_check.sh - Verify Podman installation and setup for CryoProtect
#
# This script performs post-installation verification after migrating to Podman,
# checking container functionality and network connectivity.

set -e

# Text formatting
BOLD="\033[1m"
GREEN="\033[32m"
YELLOW="\033[33m"
RED="\033[31m"
RESET="\033[0m"

# Function to log messages
log() {
    echo -e "[$(date +"%Y-%m-%d %H:%M:%S")] $1"
}

# Function to log success messages
log_success() {
    echo -e "[$(date +"%Y-%m-%d %H:%M:%S")] ${GREEN}${BOLD}$1${RESET}"
}

# Function to log warning messages
log_warning() {
    echo -e "[$(date +"%Y-%m-%d %H:%M:%S")] ${YELLOW}${BOLD}$1${RESET}"
}

# Function to log error messages
log_error() {
    echo -e "[$(date +"%Y-%m-%d %H:%M:%S")] ${RED}${BOLD}$1${RESET}"
}

# Function to check IPv4 Supabase connectivity
check_ipv4_connectivity() {
    local host=$(echo $SUPABASE_URL | sed -e 's|^[^/]*//||' -e 's|/.*$||')
    log "Testing IPv4 connectivity to Supabase host: $host"

    if ping -c 3 -4 $host > /dev/null 2>&1; then
        log_success "IPv4 connectivity to $host confirmed"
        return 0
    else
        log_error "Cannot reach $host via IPv4"
        return 1
    fi
}

# Function to verify container functionality
verify_container() {
    log "Verifying container setup..."
    
    # Check for running containers
    if podman ps | grep -q cryoprotect; then
        log_success "CryoProtect container is running"
    else
        log_warning "CryoProtect container is not running. Starting minimal setup..."
        
        # Try to start the minimal container
        if [ -f podman-compose.minimal.yml ]; then
            podman-compose -f podman-compose.minimal.yml up -d
            sleep 5
            
            if podman ps | grep -q cryoprotect; then
                log_success "Successfully started CryoProtect container"
            else
                log_error "Failed to start CryoProtect container"
                return 1
            fi
        else
            log_error "podman-compose.minimal.yml not found"
            return 1
        fi
    fi
    
    # Check container health
    local container_id=$(podman ps | grep cryoprotect | awk '{print $1}')
    if [ -n "$container_id" ]; then
        if podman healthcheck run $container_id &> /dev/null; then
            log_success "Container health check passed"
        else
            log_warning "Container health check failed"
        fi
    fi
    
    return 0
}

# Function to check SELinux contexts
check_selinux_contexts() {
    log "Checking SELinux contexts for volume mounts..."
    
    if ! command -v getenforce &> /dev/null; then
        log_warning "SELinux management tools not found"
        return 0
    fi
    
    if [ "$(getenforce)" == "Disabled" ]; then
        log_warning "SELinux is disabled on this system"
        return 0
    fi
    
    # Check context of volume mount directories
    local dirs=("./logs" "./backup/data" "./redis-data")
    local has_errors=0
    
    for dir in "${dirs[@]}"; do
        if [ -d "$dir" ]; then
            local context=$(ls -Zd $dir | awk '{print $1}')
            if [[ "$context" =~ container_file_t ]]; then
                log_success "Directory $dir has correct SELinux context: $context"
            else
                log_warning "Directory $dir has incorrect SELinux context: $context"
                log "  Run: chcon -Rt container_file_t $dir"
                has_errors=1
            fi
        else
            log_warning "Directory $dir does not exist"
        fi
    done
    
    return $has_errors
}

# Function to check network configuration
check_network_config() {
    log "Checking container network configuration..."
    
    # Check if container network exists
    if podman network ls | grep -q app-network; then
        log_success "Container network 'app-network' exists"
    else
        log_warning "Container network 'app-network' does not exist"
        log "  This will be created automatically when using podman-compose"
    fi
    
    # Check container connectivity
    local container_id=$(podman ps | grep cryoprotect | awk '{print $1}')
    if [ -n "$container_id" ]; then
        # Check container IP address
        local container_ip=$(podman inspect -f '{{.NetworkSettings.IPAddress}}' $container_id)
        log "Container IP address: $container_ip"
        
        # Check port mapping
        local port_mapping=$(podman port $container_id | grep 5000)
        if [ -n "$port_mapping" ]; then
            log_success "Port mapping configured correctly: $port_mapping"
        else
            log_warning "Port 5000 not mapped correctly"
        fi
        
        # Check localhost connectivity
        if curl -s http://localhost:5000/health &> /dev/null; then
            log_success "Application is accessible via localhost:5000"
        else
            log_warning "Application is not accessible via localhost:5000"
        fi
    else
        log_warning "No running CryoProtect container found"
    fi
    
    return 0
}

# Function to verify Podman version
check_podman_version() {
    log "Checking Podman version..."
    
    local podman_version=$(podman --version | awk '{print $3}')
    local major_version=$(echo $podman_version | cut -d. -f1)
    local minor_version=$(echo $podman_version | cut -d. -f2)
    
    log "Detected Podman version: $podman_version"
    
    # Check if version is recent enough (Podman 3.0.0+)
    if [ "$major_version" -ge 3 ]; then
        log_success "Podman version $podman_version is suitable for production use"
    else
        log_warning "Podman version $podman_version might be outdated"
        log "  Consider upgrading to Podman 3.0.0 or newer"
    fi
    
    return 0
}

# Function to verify resource usage
check_resource_usage() {
    log "Checking container resource usage..."
    
    local container_id=$(podman ps | grep cryoprotect | awk '{print $1}')
    if [ -n "$container_id" ]; then
        # Get container stats
        local stats=$(podman stats --no-stream $container_id)
        log "Container resource usage:"
        echo "$stats"
        
        # Extract CPU and memory usage
        local cpu_usage=$(echo "$stats" | tail -n 1 | awk '{print $3}')
        local mem_usage=$(echo "$stats" | tail -n 1 | awk '{print $5" "$6}')
        
        log "CPU Usage: $cpu_usage"
        log "Memory Usage: $mem_usage"
    else
        log_warning "No running CryoProtect container found"
    fi
    
    return 0
}

# Function to check application logs for errors
check_logs() {
    log "Checking application logs for errors..."
    
    local container_id=$(podman ps | grep cryoprotect | awk '{print $1}')
    if [ -n "$container_id" ]; then
        # Check for errors in logs
        local error_count=$(podman logs $container_id 2>&1 | grep -i "error\|exception\|fail" | wc -l)
        
        if [ "$error_count" -eq 0 ]; then
            log_success "No obvious errors found in container logs"
        else
            log_warning "Found $error_count potential error messages in logs"
            log "  Run 'podman logs $container_id' to view the full logs"
            
            # Show the first few errors
            log "Sample error messages:"
            podman logs $container_id 2>&1 | grep -i "error\|exception\|fail" | head -5
        fi
    else
        log_warning "No running CryoProtect container found"
    fi
    
    return 0
}

# Main function to run checks
main() {
    echo -e "${BOLD}CryoProtect Podman Post-Installation Verification${RESET}"
    echo "This script will verify your Podman installation and container setup."
    echo "---------------------------------------------------------------------"
    
    # Check if .env file exists
    if [ ! -f ./.env ]; then
        log_warning ".env file not found. Some checks may fail."
        log "Create a .env file with at least SUPABASE_URL and SUPABASE_KEY."
    else
        # Load environment variables
        set -a
        source ./.env
        set +a
    fi
    
    # Run all verification steps
    check_podman_version
    echo ""
    verify_container
    echo ""
    check_selinux_contexts
    echo ""
    check_network_config
    echo ""
    check_resource_usage
    echo ""
    
    # If SUPABASE_URL is defined, check connectivity
    if [ -n "$SUPABASE_URL" ]; then
        check_ipv4_connectivity
        echo ""
    fi
    
    check_logs
    echo ""
    
    # Final summary
    echo "---------------------------------------------------------------------"
    echo -e "${BOLD}Post-Installation Verification Complete${RESET}"
    echo "Please review any warnings or errors flagged above."
    echo ""
    echo "Next steps:"
    echo "1. Run integration tests to ensure all functionality works with Podman"
    echo "2. Check application performance under load"
    echo "3. Consider creating a production deployment with proper secrets"
    echo ""
    echo "For more information, see the PODMAN_DEPLOYMENT_GUIDE.md documentation."
}

# Run the main function
main