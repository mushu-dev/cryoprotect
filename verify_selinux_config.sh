#!/bin/bash
# verify_selinux_config.sh
#
# This script verifies the SELinux configuration for the CryoProtect application.
# It checks contexts, booleans, ports, and policies to ensure proper security.
#
# Usage:
#   ./verify_selinux_config.sh [--generate-report]

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}======================================"
echo "CryoProtect SELinux Configuration Verification"
echo -e "======================================${NC}"

# Parse arguments
GENERATE_REPORT=false

for arg in "$@"; do
    case $arg in
        --generate-report)
            GENERATE_REPORT=true
            shift
            ;;
        *)
            # Unknown option
            shift
            ;;
    esac
done

# Function to check a condition and report status
check_condition() {
    local condition=$1
    local success_msg=$2
    local failure_msg=$3
    
    if eval "$condition"; then
        echo -e "${GREEN}✓${NC} $success_msg"
        return 0
    else
        echo -e "${RED}✗${NC} $failure_msg"
        return 1
    fi
}

# Initialize arrays to store test results
passed_tests=()
failed_tests=()

# ===== 1. Check SELinux Status =====
echo -e "\n${BLUE}[1/5]${NC} Checking SELinux status..."

# Check if SELinux is enabled
if check_condition "sestatus | grep -q 'SELinux status: *enabled'"; then
    passed_tests+=("SELinux is enabled")
else
    failed_tests+=("SELinux is not enabled")
fi

# Check if SELinux is in enforcing mode
if check_condition "sestatus | grep -q 'Current mode: *enforcing'"; then
    passed_tests+=("SELinux is in enforcing mode")
else
    failed_tests+=("SELinux is not in enforcing mode")
fi

# ===== 2. Check SELinux Booleans =====
echo -e "\n${BLUE}[2/5]${NC} Checking SELinux booleans..."

# Check container-related booleans
if check_condition "getsebool container_manage_cgroup | grep -q 'on'"; then
    passed_tests+=("container_manage_cgroup boolean is enabled")
else
    failed_tests+=("container_manage_cgroup boolean is not enabled")
fi

if check_condition "getsebool container_use_devices | grep -q 'on'"; then
    passed_tests+=("container_use_devices boolean is enabled")
else
    failed_tests+=("container_use_devices boolean is not enabled")
fi

# ===== 3. Check Directory Contexts =====
echo -e "\n${BLUE}[3/5]${NC} Checking directory contexts..."

# Application directory
APP_DIR="/home/mushu/Projects/CryoProtect"

# Function to check directory context
check_dir_context() {
    local dir=$1
    local expected_context=$2
    local dir_exists=false
    
    if [ -d "$dir" ]; then
        dir_exists=true
        if ls -ldZ "$dir" | grep -q "$expected_context"; then
            passed_tests+=("Directory $dir has correct context: $expected_context")
            return 0
        else
            actual_context=$(ls -ldZ "$dir" | awk '{print $4}')
            failed_tests+=("Directory $dir has wrong context: $actual_context (expected: $expected_context)")
            return 1
        fi
    else
        echo -e "${YELLOW}⚠${NC} Directory $dir does not exist, skipping context check"
        return 0
    fi
}

# Check contexts of key directories
check_dir_context "$APP_DIR/database" "container_file_t"
check_dir_context "$APP_DIR/backups" "container_file_t"
check_dir_context "$APP_DIR/reports" "user_home_t"

# ===== 4. Check Port Contexts =====
echo -e "\n${BLUE}[4/5]${NC} Checking port contexts..."

# Check if PostgreSQL port is properly labeled
if check_condition "semanage port -l | grep -qE '5432.*postgresql_port_t'"; then
    passed_tests+=("PostgreSQL port 5432 has correct context: postgresql_port_t")
else
    if check_condition "semanage port -l | grep -qE '5432.*http_port_t'"; then
        passed_tests+=("PostgreSQL port 5432 has alternative context: http_port_t")
    else
        failed_tests+=("PostgreSQL port 5432 is not properly labeled")
    fi
fi

# Check if application port is properly labeled
if check_condition "semanage port -l | grep -qE '8000.*http_port_t'"; then
    passed_tests+=("Application port 8000 has correct context: http_port_t")
else
    failed_tests+=("Application port 8000 is not properly labeled")
fi

# ===== 5. Check SELinux Policies =====
echo -e "\n${BLUE}[5/5]${NC} Checking SELinux policies..."

# Check for custom policy modules
if check_condition "semodule -l | grep -qE 'cryoprotect'"; then
    passed_tests+=("Custom CryoProtect SELinux policy is installed")
else
    failed_tests+=("No custom CryoProtect SELinux policy found")
fi

# Check for SELinux denials related to the application
echo "Checking for SELinux denials in the last 24 hours..."
DENIALS=$(ausearch -m AVC -ts yesterday 2>/dev/null | grep -i -E "cryoprotect|flask|python|postgres|container" | grep -i "denied" || echo "")

if [ -z "$DENIALS" ]; then
    echo -e "${GREEN}No SELinux denials found for CryoProtect.${NC}"
    passed_tests+=("No SELinux denials found")
else
    echo -e "${RED}SELinux denials found for CryoProtect!${NC}"
    echo "$DENIALS" | head -n 3
    if [ $(echo "$DENIALS" | wc -l) -gt 3 ]; then
        echo -e "${YELLOW}(showing first 3 of $(echo "$DENIALS" | wc -l) denials)${NC}"
    fi
    failed_tests+=("SELinux denials found: $(echo "$DENIALS" | wc -l) denials")
fi

# ===== Generate Report =====
if [ "$GENERATE_REPORT" = true ]; then
    echo -e "\n${BLUE}Generating SELinux verification report...${NC}"
    
    REPORT_FILE="$APP_DIR/reports/selinux_verification_$(date +%Y%m%d_%H%M%S).md"
    mkdir -p "$APP_DIR/reports"
    
    # Create report
    cat > "$REPORT_FILE" << EOL
# CryoProtect SELinux Verification Report

Generated: $(date)

## Summary

* **Passed Tests:** ${#passed_tests[@]}
* **Failed Tests:** ${#failed_tests[@]}
* **Overall Status:** $([ ${#failed_tests[@]} -eq 0 ] && echo "✅ PASS" || echo "❌ FAIL")

## System Information

\`\`\`
$(sestatus)
\`\`\`

## Passed Tests

$(for test in "${passed_tests[@]}"; do echo "* ✅ $test"; done)

## Failed Tests

$(if [ ${#failed_tests[@]} -eq 0 ]; then echo "* None"; else for test in "${failed_tests[@]}"; do echo "* ❌ $test"; done; fi)

## SELinux Booleans

\`\`\`
$(getsebool -a | grep container)
\`\`\`

## Directory Contexts

\`\`\`
$(ls -ldZ "$APP_DIR/database" "$APP_DIR/backups" "$APP_DIR/reports" 2>/dev/null || echo "Directories not found")
\`\`\`

## Port Contexts

\`\`\`
$(semanage port -l | grep -E '5432|5433|8000')
\`\`\`

## Active SELinux Policies

\`\`\`
$(semodule -l | grep -i "cryoprotect")
\`\`\`

## Recommendations

$(if [ ${#failed_tests[@]} -eq 0 ]; then
    echo "* SELinux configuration is correct, no changes needed."
else
    echo "* Run the SELinux configuration script: \`sudo ./setup_selinux.sh\`"
    echo "* For Podman-specific settings: \`sudo ./setup_podman_selinux.sh --apply\`"
    if echo "${failed_tests[@]}" | grep -q "denials"; then
        echo "* Check SELinux denials: \`sudo ./monitor_selinux_denials.sh\`"
    fi
fi)

EOL
    
    echo -e "${GREEN}Report generated:${NC} $REPORT_FILE"
    
    # Create a symlink to the latest report
    ln -sf "$REPORT_FILE" "$APP_DIR/reports/selinux_verification_latest.md"
fi

# ===== Summary =====
echo -e "\n${BLUE}======================================"
echo "SELinux Verification Summary"
echo -e "======================================${NC}"
echo -e "Passed Tests: ${GREEN}${#passed_tests[@]}${NC}"
echo -e "Failed Tests: ${RED}${#failed_tests[@]}${NC}"

if [ ${#failed_tests[@]} -eq 0 ]; then
    echo -e "\n${GREEN}✅ All SELinux verification tests passed!${NC}"
    echo -e "The CryoProtect application has proper SELinux configuration."
else
    echo -e "\n${RED}❌ Some SELinux verification tests failed.${NC}"
    echo -e "To fix the issues, run the SELinux configuration scripts:"
    echo "  sudo ./setup_selinux.sh"
    echo "  sudo ./setup_podman_selinux.sh --apply"
fi

# Exit with status based on whether all tests passed
if [ ${#failed_tests[@]} -eq 0 ]; then
    exit 0
else
    exit 1
fi