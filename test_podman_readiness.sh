#!/bin/bash
#
# test_podman_readiness.sh - Test if the system is ready for Podman migration
#
# This script performs comprehensive checks to determine if a system
# is ready for migrating from Docker to Podman, with special focus
# on Fedora-specific requirements.

set -e

# Text formatting
BOLD="\033[1m"
GREEN="\033[32m"
YELLOW="\033[33m"
RED="\033[31m"
RESET="\033[0m"

# Counter for passed and failed tests
PASSED=0
FAILED=0
WARNINGS=0

# Function to log test results
log_test() {
    local status=$1
    local message=$2
    
    if [ "$status" == "PASS" ]; then
        echo -e "[${GREEN}${BOLD}PASS${RESET}] $message"
        PASSED=$((PASSED+1))
    elif [ "$status" == "WARN" ]; then
        echo -e "[${YELLOW}${BOLD}WARN${RESET}] $message"
        WARNINGS=$((WARNINGS+1))
    else
        echo -e "[${RED}${BOLD}FAIL${RESET}] $message"
        FAILED=$((FAILED+1))
    fi
}

# Function to check if a command exists
command_exists() {
    command -v "$1" &> /dev/null
}

# Function to check if a file exists
file_exists() {
    [ -f "$1" ]
}

# Function to check if a directory exists
dir_exists() {
    [ -d "$1" ]
}

echo -e "${BOLD}CryoProtect Podman Migration Readiness Test${RESET}"
echo "This script will check if your system is ready for migrating from Docker to Podman."
echo "---------------------------------------------------------------------"

# Check if running on Fedora
if [ -f /etc/fedora-release ]; then
    FEDORA_VERSION=$(cat /etc/fedora-release | grep -oP '(?<=release )[0-9]+')
    log_test "PASS" "Running on Fedora $FEDORA_VERSION"
    
    # Check if Fedora version is 35 or higher (recommended for Podman)
    if [ "$FEDORA_VERSION" -ge 35 ]; then
        log_test "PASS" "Fedora version is $FEDORA_VERSION (recommended: >=35)"
    else
        log_test "WARN" "Fedora version is $FEDORA_VERSION (recommended: >=35)"
    fi
else
    log_test "WARN" "Not running on Fedora. This script is optimized for Fedora systems."
fi

# Check if Podman is installed
if command_exists podman; then
    PODMAN_VERSION=$(podman --version | awk '{print $3}')
    log_test "PASS" "Podman is installed (version $PODMAN_VERSION)"
else
    log_test "FAIL" "Podman is not installed. Run: sudo dnf install -y podman"
fi

# Check if podman-compose is installed
if command_exists podman-compose; then
    PODMAN_COMPOSE_VERSION=$(podman-compose --version | awk '{print $3}')
    log_test "PASS" "podman-compose is installed (version $PODMAN_COMPOSE_VERSION)"
else
    log_test "FAIL" "podman-compose is not installed. Run: sudo dnf install -y podman-compose"
fi

# Check if Docker is installed (not required, but good to know)
if command_exists docker; then
    DOCKER_VERSION=$(docker --version | awk '{print $3}' | tr -d ',')
    log_test "WARN" "Docker is installed (version $DOCKER_VERSION). Consider removing after migration."
else
    log_test "PASS" "Docker is not installed. Clean Podman-only setup possible."
fi

# Check if SELinux is enabled (important for Podman on Fedora)
if command_exists getenforce; then
    SELINUX_STATUS=$(getenforce)
    if [ "$SELINUX_STATUS" == "Enforcing" ]; then
        log_test "PASS" "SELinux is enabled and enforcing (good for Podman security)"
    elif [ "$SELINUX_STATUS" == "Permissive" ]; then
        log_test "WARN" "SELinux is in permissive mode. Consider enabling enforcing mode."
    else
        log_test "WARN" "SELinux is disabled. Podman works better with SELinux enabled."
    fi
else
    log_test "WARN" "Cannot determine SELinux status."
fi

# Check for container-selinux package
if rpm -q container-selinux &> /dev/null; then
    log_test "PASS" "container-selinux package is installed"
else
    log_test "FAIL" "container-selinux package is not installed. Run: sudo dnf install container-selinux"
fi

# Check if project structure exists
if [ -f "app.py" ]; then
    log_test "PASS" "CryoProtect project files found"
else
    log_test "WARN" "CryoProtect project files not found. Are you in the project directory?"
fi

# Check for Docker files to convert
if file_exists "docker-compose.yml"; then
    log_test "PASS" "docker-compose.yml found (will be converted to podman-compose.yml)"
else
    log_test "WARN" "docker-compose.yml not found. Manual configuration may be needed."
fi

if file_exists "Dockerfile"; then
    log_test "PASS" "Dockerfile found (compatible with Podman)"
else
    log_test "WARN" "Dockerfile not found. Container build may not be possible."
fi

# Check for firewalld (common on Fedora)
if command_exists firewall-cmd; then
    log_test "PASS" "firewalld is installed (may need configuration for container networking)"
else
    log_test "WARN" "firewalld not found. This is unusual for Fedora."
fi

# Check usernamespace support (important for rootless podman)
if [ -f /proc/sys/kernel/unprivileged_userns_clone ] && [ "$(cat /proc/sys/kernel/unprivileged_userns_clone)" -eq 1 ]; then
    log_test "PASS" "Unprivileged user namespaces are enabled (good for rootless Podman)"
else
    log_test "WARN" "Unprivileged user namespaces may not be enabled. Rootless Podman might have issues."
fi

# Check system resources
MEMORY=$(free -m | awk '/^Mem:/{print $2}')
if [ "$MEMORY" -ge 2048 ]; then
    log_test "PASS" "System has sufficient memory ($MEMORY MB)"
else
    log_test "WARN" "System memory may be insufficient ($MEMORY MB, recommended: >=2048 MB)"
fi

DISK_SPACE=$(df -h . | awk 'NR==2 {print $4}')
log_test "PASS" "Available disk space: $DISK_SPACE"

# Test basic Podman functionality if installed
if command_exists podman; then
    if podman info &> /dev/null; then
        log_test "PASS" "Podman is functioning properly"
    else
        log_test "FAIL" "Podman is installed but not functioning properly"
    fi
fi

# Summary
echo "---------------------------------------------------------------------"
echo -e "${BOLD}Test Summary:${RESET}"
echo -e "${GREEN}Passed: $PASSED${RESET}"
echo -e "${YELLOW}Warnings: $WARNINGS${RESET}"
echo -e "${RED}Failed: $FAILED${RESET}"

if [ $FAILED -eq 0 ]; then
    if [ $WARNINGS -eq 0 ]; then
        echo -e "\n${GREEN}${BOLD}All tests passed! Your system is ready for Podman migration.${RESET}"
        echo "You can now run the migrate_to_podman.sh script."
    else
        echo -e "\n${YELLOW}${BOLD}Your system can proceed with Podman migration, but there are some warnings to consider.${RESET}"
        echo "Review the warnings above before running migrate_to_podman.sh."
    fi
else
    echo -e "\n${RED}${BOLD}Your system has some issues that need to be addressed before migration.${RESET}"
    echo "Fix the failed tests before running migrate_to_podman.sh."
fi

exit $FAILED