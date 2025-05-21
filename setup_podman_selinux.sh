#!/bin/bash
# setup_podman_selinux.sh
#
# This script configures SELinux specifically for Podman containers in the CryoProtect application.
# It addresses common SELinux issues with Podman and ensures correct volume labeling.
#
# Usage:
#   ./setup_podman_selinux.sh [--apply]

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if running as root
if [ "$(id -u)" -ne 0 ]; then
    echo -e "${RED}Error: This script must be run as root. Use 'sudo ./setup_podman_selinux.sh'${NC}"
    exit 1
fi

echo -e "${BLUE}======================================"
echo "CryoProtect Podman SELinux Configuration"
echo -e "======================================${NC}"
echo "Configuring SELinux for Podman containers..."

# Parse arguments
APPLY=false

for arg in "$@"; do
    case $arg in
        --apply)
            APPLY=true
            shift
            ;;
        *)
            # Unknown option
            shift
            ;;
    esac
done

# Application directory
APP_DIR="/home/mushu/Projects/CryoProtect"

# ===== 1. Podman SELinux Booleans =====

echo -e "\n${BLUE}[1/5]${NC} Configuring Podman SELinux booleans..."

echo "Current container boolean settings:"
getsebool -a | grep container | sort

if [ "$APPLY" = true ]; then
    echo "Setting Podman SELinux booleans..."
    
    # Essential booleans for Podman
    setsebool -P container_manage_cgroup on
    setsebool -P container_use_devices on
    
    # Additional useful booleans
    setsebool -P container_user_exec_content on
    
    echo -e "${GREEN}Podman SELinux booleans set.${NC}"
else
    echo -e "${YELLOW}Preview: Would set these booleans:${NC}"
    echo "container_manage_cgroup -> on"
    echo "container_use_devices -> on"
    echo "container_user_exec_content -> on"
    echo -e "${YELLOW}Use --apply to make these changes${NC}"
fi

# ===== 2. Podman Volume Management =====

echo -e "\n${BLUE}[2/5]${NC} Setting up Podman volume SELinux management..."

# Identify volume directories
VOLUME_DIRS=(
    "$APP_DIR/database"
    "$APP_DIR/backups"
    "$APP_DIR/reports"
)

if [ "$APPLY" = true ]; then
    for dir in "${VOLUME_DIRS[@]}"; do
        if [ -d "$dir" ]; then
            echo "Setting container_file_t context for $dir..."
            semanage fcontext -a -t container_file_t "$dir(/.*)?"
            restorecon -Rv "$dir"
        else
            mkdir -p "$dir"
            echo "Created and set container_file_t context for $dir..."
            semanage fcontext -a -t container_file_t "$dir(/.*)?"
            restorecon -Rv "$dir"
        fi
    done
    
    echo -e "${GREEN}Podman volume contexts set.${NC}"
else
    echo -e "${YELLOW}Preview: Would set container_file_t on these directories:${NC}"
    for dir in "${VOLUME_DIRS[@]}"; do
        echo "- $dir"
    done
    echo -e "${YELLOW}Use --apply to make these changes${NC}"
fi

# ===== 3. Create Podman-Specific Helper Scripts =====

echo -e "\n${BLUE}[3/5]${NC} Creating Podman SELinux helper scripts..."

# Create script for running PostgreSQL container with correct SELinux context
cat > "$APP_DIR/run_postgres_container.sh" << 'EOL'
#!/bin/bash
# run_postgres_container.sh
#
# Runs PostgreSQL in a Podman container with proper SELinux contexts.
#
# Usage: 
#   ./run_postgres_container.sh

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Starting PostgreSQL container with SELinux contexts...${NC}"

# Database directory
DB_DIR="./database/postgres_data"
mkdir -p "$DB_DIR"

# Check if container is already running
CONTAINER_ID=$(podman ps -q --filter name=cryoprotect-postgres)
if [ -n "$CONTAINER_ID" ]; then
    echo -e "${YELLOW}PostgreSQL container is already running with ID: $CONTAINER_ID${NC}"
    exit 0
fi

# Set context labels on volume
sudo chcon -Rt container_file_t "$DB_DIR"

# Run PostgreSQL container with proper volume mapping
podman run -d \
    --name cryoprotect-postgres \
    --security-opt label=type:container_t \
    -e POSTGRES_PASSWORD=postgres \
    -e POSTGRES_USER=postgres \
    -e POSTGRES_DB=postgres \
    -p 5432:5432 \
    -v "$DB_DIR:/var/lib/postgresql/data:Z" \
    postgres:13

echo -e "${GREEN}PostgreSQL container started with correct SELinux contexts.${NC}"
EOL

# Create script for running Flask app container with correct SELinux context
cat > "$APP_DIR/run_flask_container.sh" << 'EOL'
#!/bin/bash
# run_flask_container.sh
#
# Runs Flask app in a Podman container with proper SELinux contexts.
#
# Usage: 
#   ./run_flask_container.sh

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Starting Flask container with SELinux contexts...${NC}"

# App directories
APP_DIR="."
REPORTS_DIR="./reports"
mkdir -p "$REPORTS_DIR"

# Check if container is already running
CONTAINER_ID=$(podman ps -q --filter name=cryoprotect-app)
if [ -n "$CONTAINER_ID" ]; then
    echo -e "${YELLOW}Flask container is already running with ID: $CONTAINER_ID${NC}"
    exit 0
fi

# Set context labels on volume
sudo chcon -Rt container_file_t "$REPORTS_DIR"

# Build the container if not exists
if ! podman image exists cryoprotect-app:latest; then
    echo -e "${YELLOW}Building cryoprotect-app image...${NC}"
    podman build -t cryoprotect-app:latest .
fi

# Run Flask container with proper volume mapping
podman run -d \
    --name cryoprotect-app \
    --security-opt label=type:container_t \
    -p 8000:8000 \
    -v "$APP_DIR:/app:Z" \
    -v "$REPORTS_DIR:/app/reports:Z" \
    --network=host \
    cryoprotect-app:latest

echo -e "${GREEN}Flask container started with correct SELinux contexts.${NC}"
EOL

# Make the scripts executable
chmod +x "$APP_DIR/run_postgres_container.sh"
chmod +x "$APP_DIR/run_flask_container.sh"

echo -e "${GREEN}Podman helper scripts created in $APP_DIR${NC}"

# ===== 4. Create Podman Compose File with SELinux Labels =====

echo -e "\n${BLUE}[4/5]${NC} Creating Podman Compose file with SELinux configuration..."

cat > "$APP_DIR/podman-compose.selinux.yml" << 'EOL'
version: '3'

services:
  postgres:
    image: postgres:13
    container_name: cryoprotect-postgres
    restart: always
    environment:
      POSTGRES_USER: postgres
      POSTGRES_PASSWORD: postgres
      POSTGRES_DB: postgres
    ports:
      - "5432:5432"
    volumes:
      - ./database/postgres_data:/var/lib/postgresql/data:Z
    security_opt:
      - label=type:container_t
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U postgres"]
      interval: 10s
      timeout: 5s
      retries: 5

  pgadmin:
    image: dpage/pgadmin4
    container_name: cryoprotect-pgadmin
    restart: always
    environment:
      PGADMIN_DEFAULT_EMAIL: admin@example.com
      PGADMIN_DEFAULT_PASSWORD: admin
    ports:
      - "5050:80"
    volumes:
      - ./database/pgadmin_data:/var/lib/pgadmin:Z
    depends_on:
      - postgres
    security_opt:
      - label=type:container_t

  cryoprotect:
    build: .
    container_name: cryoprotect-app
    restart: always
    ports:
      - "8000:8000"
    volumes:
      - .:/app:Z
      - ./reports:/app/reports:Z
    depends_on:
      - postgres
    environment:
      - DATABASE_URL=postgresql://postgres:postgres@postgres:5432/postgres
      - FLASK_ENV=production
    security_opt:
      - label=type:container_t
EOL

echo -e "${GREEN}Podman Compose file with SELinux labels created:${NC} $APP_DIR/podman-compose.selinux.yml"

# ===== 5. Create SELinux Policy for Podman =====

echo -e "\n${BLUE}[5/5]${NC} Creating custom SELinux policy for Podman..."

# Create temporary policy file
POLICY_FILE="cryoprotect_podman_policy"
cat > ${POLICY_FILE}.te << 'EOL'
module cryoprotect_podman_policy 1.0;

require {
    type container_t;
    type container_file_t;
    type postgresql_port_t;
    type http_port_t;
    type user_home_t;
    class tcp_socket name_connect;
    class file { open read write getattr append create unlink };
    class dir { add_name read remove_name search write };
    class unix_stream_socket connectto;
    class process { setrlimit };
}

# Allow containers to connect to PostgreSQL
allow container_t postgresql_port_t:tcp_socket name_connect;

# Allow containers to connect to web ports
allow container_t http_port_t:tcp_socket name_connect;

# Allow containers to work with container_file_t
allow container_t container_file_t:file { open read write getattr append create unlink };
allow container_t container_file_t:dir { add_name read remove_name search write };

# Allow containers to work with user home dirs
allow container_t user_home_t:file { read getattr open };
allow container_t user_home_t:dir search;

# Allow containers to set resource limits
allow container_t self:process setrlimit;

# Allow container to connect to container
allow container_t self:unix_stream_socket connectto;
EOL

if [ "$APPLY" = true ]; then
    echo "Compiling and installing SELinux policy for Podman..."
    
    # Compile and install the policy module
    checkmodule -M -m -o ${POLICY_FILE}.mod ${POLICY_FILE}.te
    semodule_package -o ${POLICY_FILE}.pp -m ${POLICY_FILE}.mod
    semodule -i ${POLICY_FILE}.pp
    
    # Preserve the policy file
    cp ${POLICY_FILE}.te "$APP_DIR/${POLICY_FILE}.te"
    
    # Clean up temporary files
    rm -f ${POLICY_FILE}.mod ${POLICY_FILE}.pp
    
    echo -e "${GREEN}SELinux policy for Podman installed.${NC}"
else
    # Just save the policy file
    mv ${POLICY_FILE}.te "$APP_DIR/${POLICY_FILE}.te"
    
    echo -e "${YELLOW}Preview: Created SELinux policy file:${NC} $APP_DIR/${POLICY_FILE}.te"
    echo -e "${YELLOW}Use --apply to compile and install this policy${NC}"
fi

# ===== Create Documentation =====

echo -e "\n${BLUE}Creating SELinux for Podman documentation...${NC}"

cat > "$APP_DIR/PODMAN_SELINUX_GUIDE.md" << 'EOL'
# Podman SELinux Configuration Guide for CryoProtect

This guide explains how to use Podman with proper SELinux configuration for the CryoProtect application.

## Understanding SELinux & Podman

Podman with SELinux provides strong security isolation between containers and the host system. Key concepts:

1. **Container Process Labels**: Determine what the container process can do
2. **Volume Labels**: Control access to mounted volumes
3. **Port Labels**: Control network access

## Using Volume Mounts with SELinux

When mounting volumes in Podman, use these options:

- `:Z` - For dedicated volumes (used by a single container)
- `:z` - For shared volumes (used by multiple containers)

Example:
```bash
# Mount with dedicated label
podman run -v /host/dir:/container/dir:Z ...

# Mount with shared label
podman run -v /host/dir:/container/dir:z ...
```

## Podman Compose with SELinux

We provide a SELinux-aware Podman Compose file: `podman-compose.selinux.yml`

To use it:
```bash
podman-compose -f podman-compose.selinux.yml up -d
```

## Helper Scripts

We provide scripts to run containers with proper SELinux contexts:

1. **run_postgres_container.sh** - Runs PostgreSQL with correct contexts
2. **run_flask_container.sh** - Runs the Flask app with correct contexts

## Troubleshooting SELinux with Podman

### Common Issues and Solutions

1. **"Permission denied" accessing volumes**
   - Check the SELinux context of the volume directory:
     ```bash
     ls -ldZ /path/to/volume
     ```
   - Apply the container_file_t context:
     ```bash
     sudo chcon -Rt container_file_t /path/to/volume
     ```
   - Or use `:Z` in your volume mount

2. **Container can't connect to PostgreSQL**
   - Ensure the PostgreSQL port has the correct label:
     ```bash
     sudo semanage port -l | grep 5432
     ```
   - Add the label if missing:
     ```bash
     sudo semanage port -a -t postgresql_port_t -p tcp 5432
     ```
   - Apply the custom SELinux policy:
     ```bash
     sudo semodule -i cryoprotect_podman_policy.pp
     ```

3. **Network connectivity issues**
   - Allow container networking:
     ```bash
     sudo setsebool -P container_connect_any on
     ```

### Viewing SELinux Denials

To see if SELinux is blocking Podman:
```bash
sudo ausearch -m AVC -ts recent | grep -i container
```

### Creating Custom Rules from Denials

If you see denials, you can create custom rules:
```bash
sudo ausearch -m AVC -ts recent | audit2allow -M my_podman_policy
sudo semodule -i my_podman_policy.pp
```

## Best Practices

1. Always use `:Z` or `:z` for volume mounts
2. Run containers with `--security-opt label=type:container_t`
3. Apply appropriate context to host directories before mounting
4. Use the provided helper scripts for consistent SELinux contexts
5. Check for SELinux denials after making changes to container configurations

## Resources

- [Podman and SELinux](https://www.redhat.com/sysadmin/podman-security-selinux)
- [Container SELinux Policies](https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/8/html/security_hardening/using-selinux_security-hardening#container-selinux-policies_using-selinux)
- [Udica: Custom SELinux policies for containers](https://github.com/containers/udica)
EOL

echo -e "${GREEN}Podman SELinux guide created:${NC} $APP_DIR/PODMAN_SELINUX_GUIDE.md"

# ===== Summary =====

echo -e "\n${BLUE}======================================"
echo "Podman SELinux Configuration Summary"
echo -e "======================================${NC}"

if [ "$APPLY" = true ]; then
    echo -e "${GREEN}✓${NC} Podman SELinux booleans set"
    echo -e "${GREEN}✓${NC} Podman volume contexts configured"
    echo -e "${GREEN}✓${NC} Podman SELinux policy installed"
else
    echo -e "${YELLOW}⚠${NC} Changes previewed but not applied (use --apply to make changes)"
fi

echo -e "${GREEN}✓${NC} Podman helper scripts created"
echo -e "${GREEN}✓${NC} Podman Compose file with SELinux labels created"
echo -e "${GREEN}✓${NC} Podman SELinux documentation created"

echo -e "\n${BLUE}Next Steps:${NC}"
echo "1. Apply the SELinux configuration (if not done):"
echo "   sudo ./setup_podman_selinux.sh --apply"
echo "2. Use the provided helper scripts to start containers:"
echo "   ./run_postgres_container.sh"
echo "   ./run_flask_container.sh"
echo "3. Or use Podman Compose with SELinux configuration:"
echo "   podman-compose -f podman-compose.selinux.yml up -d"
echo "4. Read the Podman SELinux guide for more information:"
echo "   less PODMAN_SELINUX_GUIDE.md"

echo -e "\n${GREEN}Podman SELinux configuration completed.${NC}"