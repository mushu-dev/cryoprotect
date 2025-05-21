#!/bin/bash
# setup_selinux.sh
#
# This script configures SELinux for the CryoProtect application.
# It sets up appropriate contexts, booleans, and ports for secure operation.
#
# Usage:
#   ./setup_selinux.sh [--full-setup]

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if running as root
if [ "$(id -u)" -ne 0 ]; then
    echo -e "${RED}Error: This script must be run as root. Use 'sudo ./setup_selinux.sh'${NC}"
    exit 1
fi

echo -e "${BLUE}======================================"
echo "CryoProtect SELinux Configuration"
echo -e "======================================${NC}"
echo "Configuring SELinux for secure operation..."

# Parse arguments
FULL_SETUP=false

for arg in "$@"; do
    case $arg in
        --full-setup)
            FULL_SETUP=true
            shift
            ;;
        *)
            # Unknown option
            shift
            ;;
    esac
done

# Check SELinux status
SELINUX_STATUS=$(getenforce)
if [ "$SELINUX_STATUS" != "Enforcing" ]; then
    echo -e "${YELLOW}Warning: SELinux is not in enforcing mode. Current mode: $SELINUX_STATUS${NC}"
    echo "Recommended: Set SELinux to enforcing mode for better security"
    echo "You can change this in /etc/selinux/config and reboot"
fi

# ===== 1. Container SELinux Settings =====

echo -e "\n${BLUE}[1/6]${NC} Configuring container SELinux settings..."

# Set container-related booleans
echo "Setting container SELinux booleans..."
setsebool -P container_manage_cgroup on
setsebool -P container_use_devices on

# Check status of settings
echo -e "${GREEN}Container SELinux booleans set:${NC}"
getsebool container_manage_cgroup container_use_devices

# ===== 2. Database Directory Contexts =====

echo -e "\n${BLUE}[2/6]${NC} Setting up database directory contexts..."

# Define application directories
APP_DIR="/home/mushu/Projects/CryoProtect"
DATABASE_DIR="$APP_DIR/database"
BACKUPS_DIR="$APP_DIR/backups"
REPORTS_DIR="$APP_DIR/reports"

# Set SELinux contexts for database directories
echo "Setting contexts for database directories..."
semanage fcontext -a -t container_file_t "$DATABASE_DIR(/.*)?"
semanage fcontext -a -t container_file_t "$BACKUPS_DIR(/.*)?"
restorecon -Rv "$DATABASE_DIR" "$BACKUPS_DIR"

# Apply special context for PostgreSQL data directory if it exists
POSTGRES_DATA_DIR="$APP_DIR/postgres_data"
if [ -d "$POSTGRES_DATA_DIR" ]; then
    echo "Setting contexts for PostgreSQL data directory..."
    semanage fcontext -a -t container_file_t "$POSTGRES_DATA_DIR(/.*)?"
    restorecon -Rv "$POSTGRES_DATA_DIR"
fi

echo -e "${GREEN}Database directory contexts set${NC}"

# ===== 3. Port Configurations =====

echo -e "\n${BLUE}[3/6]${NC} Configuring SELinux port settings..."

# Configure ports used by the application
echo "Setting up port contexts..."
semanage port -a -t http_port_t -p tcp 5432 || semanage port -m -t http_port_t -p tcp 5432
semanage port -a -t http_port_t -p tcp 5433 || semanage port -m -t http_port_t -p tcp 5433
semanage port -a -t http_port_t -p tcp 8000 || semanage port -m -t http_port_t -p tcp 8000

echo -e "${GREEN}Port contexts set${NC}"

# ===== 4. Volume Mount Contexts =====

echo -e "\n${BLUE}[4/6]${NC} Setting up volume mount contexts..."

# Special context for systemd service files
SYSTEMD_DIR="/etc/systemd/system"
if [ -f "$SYSTEMD_DIR/cryoprotect.service" ]; then
    echo "Setting contexts for systemd service files..."
    semanage fcontext -a -t systemd_unit_file_t "$SYSTEMD_DIR/cryoprotect*.service"
    semanage fcontext -a -t systemd_unit_file_t "$SYSTEMD_DIR/cryoprotect*.timer"
    restorecon -v "$SYSTEMD_DIR/cryoprotect"*.service "$SYSTEMD_DIR/cryoprotect"*.timer 2>/dev/null || true
fi

echo -e "${GREEN}Volume mount contexts set${NC}"

# ===== 5. Create Custom SELinux Policy =====

if [ "$FULL_SETUP" = true ]; then
    echo -e "\n${BLUE}[5/6]${NC} Creating custom SELinux policy module..."

    # Create temporary policy file
    POLICY_FILE="cryoprotect_custom"
    cat > ${POLICY_FILE}.te << 'EOL'
module cryoprotect_custom 1.0;

require {
    type container_t;
    type http_port_t;
    type postgresql_port_t;
    type container_file_t;
    type user_home_t;
    class tcp_socket name_connect;
    class file { read write getattr open };
    class dir { search read write add_name remove_name };
}

# Allow containers to connect to PostgreSQL port
allow container_t postgresql_port_t:tcp_socket name_connect;

# Allow containers to read/write container_file_t files and directories
allow container_t container_file_t:file { read write getattr open };
allow container_t container_file_t:dir { search read write add_name remove_name };

# Allow containers to read user home files (for configuration)
allow container_t user_home_t:file { read getattr open };
allow container_t user_home_t:dir search;
EOL

    # Compile and install the policy module
    echo "Compiling and installing SELinux policy module..."
    checkmodule -M -m -o ${POLICY_FILE}.mod ${POLICY_FILE}.te
    semodule_package -o ${POLICY_FILE}.pp -m ${POLICY_FILE}.mod
    semodule -i ${POLICY_FILE}.pp

    # Clean up temporary files
    rm -f ${POLICY_FILE}.te ${POLICY_FILE}.mod ${POLICY_FILE}.pp

    echo -e "${GREEN}Custom SELinux policy module installed${NC}"
else
    echo -e "\n${BLUE}[5/6]${NC} Skipping custom SELinux policy creation (use --full-setup to enable)"
fi

# ===== 6. Setup SELinux Auditing =====

echo -e "\n${BLUE}[6/6]${NC} Setting up SELinux auditing..."

# Create a script to monitor SELinux denials
cat > "$APP_DIR/monitor_selinux_denials.sh" << 'EOL'
#!/bin/bash
# monitor_selinux_denials.sh
#
# This script monitors and reports SELinux denials related to CryoProtect.
# It can also generate policy recommendations to fix issues.

# Colors
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

REPORT_DIR="reports/selinux"
mkdir -p "$REPORT_DIR"

echo -e "${BLUE}CryoProtect SELinux Denial Monitor${NC}"
echo "Checking for SELinux denials..."

# Check recent denials
DENIALS=$(ausearch -m AVC -ts recent 2>/dev/null | grep -i "denied" | grep -iE "cryoprotect|flask|python|postgres|container" || echo "")

if [ -z "$DENIALS" ]; then
    echo -e "${GREEN}No SELinux denials found related to CryoProtect.${NC}"
else
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    REPORT_FILE="$REPORT_DIR/selinux_denials_$TIMESTAMP.txt"
    POLICY_FILE="$REPORT_DIR/selinux_policy_$TIMESTAMP.te"
    
    echo -e "${YELLOW}SELinux denials found!${NC}"
    
    # Save denials to report file
    echo "SELinux Denials Report - $(date)" > "$REPORT_FILE"
    echo "=================================" >> "$REPORT_FILE"
    echo "$DENIALS" >> "$REPORT_FILE"
    
    # Generate policy recommendations
    echo -e "\n${BLUE}Generating policy recommendations...${NC}"
    ausearch -m AVC -ts recent 2>/dev/null | audit2allow -R > "$POLICY_FILE"
    
    echo -e "\nDenials saved to: $REPORT_FILE"
    echo -e "Policy recommendations saved to: $POLICY_FILE"
    
    # Show instructions
    echo -e "\n${YELLOW}To fix these denials, you can create a custom policy module:${NC}"
    echo "sudo ausearch -m AVC -ts recent | audit2allow -M cryoprotect_fix"
    echo "sudo semodule -i cryoprotect_fix.pp"
fi

# Check for permissive domains
PERMISSIVE=$(semodule -l | grep permissive)
if [ -n "$PERMISSIVE" ]; then
    echo -e "\n${YELLOW}Warning: The following domains are in permissive mode:${NC}"
    echo "$PERMISSIVE"
fi

echo -e "\n${BLUE}SELinux Context Information:${NC}"
echo -e "${YELLOW}Container Processes:${NC}"
ps -eZ | grep -E "container|postgres|flask|python" | grep -v grep || echo "No matching processes found"

echo -e "\n${YELLOW}Important Directories:${NC}"
ls -ldZ database backups reports 2>/dev/null || echo "Directories not found"
EOL

# Make the monitoring script executable
chmod +x "$APP_DIR/monitor_selinux_denials.sh"

echo -e "${GREEN}SELinux monitoring script created: $APP_DIR/monitor_selinux_denials.sh${NC}"

# ===== Final Configuration =====

echo -e "\n${BLUE}Running final SELinux configuration checks...${NC}"

# Check for any remaining issues
echo "Checking for SELinux issues..."
ISSUES=$(ausearch -m AVC -ts today 2>/dev/null | grep -i "denied" | grep -i "cryoprotect" || echo "")

if [ -n "$ISSUES" ]; then
    echo -e "${YELLOW}Warning: SELinux denials detected for CryoProtect.${NC}"
    echo "Run the monitoring script for details: ./monitor_selinux_denials.sh"
else
    echo -e "${GREEN}No SELinux issues detected for CryoProtect.${NC}"
fi

# Create documentation file
cat > "$APP_DIR/SELINUX_CONFIGURATION_GUIDE.md" << 'EOL'
# CryoProtect SELinux Configuration Guide

This document describes the SELinux configuration applied to secure the CryoProtect application.

## Current Configuration

The SELinux configuration for CryoProtect includes:

1. **Container Security Settings**
   - `container_manage_cgroup` - Allows containers to manage cgroup files
   - `container_use_devices` - Allows containers to use devices

2. **Directory Contexts**
   - Database directories: `container_file_t`
   - Backup directories: `container_file_t`
   - PostgreSQL data directory: `container_file_t`

3. **Port Configurations**
   - PostgreSQL ports (5432, 5433): `http_port_t`
   - Application port (8000): `http_port_t`

4. **Custom Policy**
   - Custom policy module for CryoProtect-specific rules
   - Allows containers to connect to PostgreSQL
   - Allows containers to read/write data directories

## Managing SELinux for CryoProtect

### Checking Status

To check the status of SELinux settings:

```bash
# Check SELinux status
sestatus

# Check booleans related to containers
getsebool -a | grep container

# Check port assignments
semanage port -l | grep -E "5432|5433|8000"

# Check file contexts
ls -ldZ database backups postgres_data
```

### Monitoring Denials

A monitoring script is provided to check for SELinux denials:

```bash
./monitor_selinux_denials.sh
```

This script will report any denials and suggest policy fixes.

### Troubleshooting Common Issues

1. **Container Cannot Connect to PostgreSQL**
   - Check if the PostgreSQL port is properly labeled:
     ```bash
     semanage port -l | grep 5432
     ```
   - If missing, add it:
     ```bash
     sudo semanage port -a -t postgresql_port_t -p tcp 5432
     ```

2. **Container Cannot Access Files**
   - Check file contexts:
     ```bash
     ls -lZ /path/to/file
     ```
   - Fix context if needed:
     ```bash
     sudo semanage fcontext -a -t container_file_t "/path/to/file(/.*)?"
     sudo restorecon -Rv /path/to/file
     ```

3. **Generating Policy from Denials**
   - Use the following to create a policy from denials:
     ```bash
     sudo ausearch -m AVC -ts recent | audit2allow -M mycustommodule
     sudo semodule -i mycustommodule.pp
     ```

## Best Practices

1. Keep SELinux in enforcing mode
2. Use volumes with `:Z` suffix for dedicated container volumes
3. Use `:z` suffix for shared container volumes
4. Regularly run the monitoring script to check for issues
5. After making application changes, check for new denials
6. Update the SELinux policy when adding new features

## Resources

- [Fedora SELinux Guide](https://docs.fedoraproject.org/en-US/quick-docs/selinux/)
- [Container SELinux Contexts](https://www.redhat.com/sysadmin/container-selinux-types)
- [SELinux Troubleshooting](https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/7/html/selinux_users_and_administrators_guide/chap-security-enhanced_linux-troubleshooting)
EOL

echo -e "${GREEN}SELinux configuration guide created: $APP_DIR/SELINUX_CONFIGURATION_GUIDE.md${NC}"

# ===== Summary =====

echo -e "\n${BLUE}======================================"
echo "SELinux Configuration Summary"
echo -e "======================================${NC}"
echo -e "${GREEN}✓${NC} Container SELinux booleans configured"
echo -e "${GREEN}✓${NC} Database directory contexts set"
echo -e "${GREEN}✓${NC} Port contexts configured"
echo -e "${GREEN}✓${NC} Volume mount contexts set"
if [ "$FULL_SETUP" = true ]; then
    echo -e "${GREEN}✓${NC} Custom SELinux policy module installed"
else
    echo -e "${YELLOW}⚠${NC} Custom SELinux policy module not installed (use --full-setup)"
fi
echo -e "${GREEN}✓${NC} SELinux monitoring tools created"
echo -e "${GREEN}✓${NC} SELinux configuration guide created"

echo -e "\n${BLUE}Next Steps:${NC}"
echo "1. Run the monitoring script to check for SELinux denials:"
echo "   sudo ./monitor_selinux_denials.sh"
echo "2. Restart your application to apply all SELinux changes"
echo "3. Review the SELinux configuration guide for more information:"
echo "   less SELINUX_CONFIGURATION_GUIDE.md"

echo -e "\n${GREEN}SELinux configuration completed successfully.${NC}"