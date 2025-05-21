#!/bin/bash
# optimize_postgresql_fedora.sh
#
# This script optimizes PostgreSQL for the CryoProtect application on Fedora.
# It configures performance settings, security, and system integration.
#
# Usage:
#   ./optimize_postgresql_fedora.sh [--apply] [--container] [--local]

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default options
APPLY=false
CONTAINER_MODE=false
LOCAL_MODE=false

# Parse command line arguments
for arg in "$@"; do
    case $arg in
        --apply)
            APPLY=true
            shift
            ;;
        --container)
            CONTAINER_MODE=true
            shift
            ;;
        --local)
            LOCAL_MODE=true
            shift
            ;;
        *)
            # Unknown option
            shift
            ;;
    esac
done

# Function to check if running as root
check_root() {
    if [ "$APPLY" = true ] && [ "$EUID" -ne 0 ] && [ "$CONTAINER_MODE" = false ]; then
        echo -e "${RED}Error: This script must be run as root when using --apply.${NC}"
        echo -e "Please run with 'sudo $0 --apply' or use '--container' for container mode."
        exit 1
    fi
}

# Check environment
check_environment() {
    if [ "$CONTAINER_MODE" = true ] && [ "$LOCAL_MODE" = true ]; then
        echo -e "${RED}Error: Cannot use both --container and --local modes.${NC}"
        exit 1
    fi
    
    # Default to container mode if neither is specified
    if [ "$CONTAINER_MODE" = false ] && [ "$LOCAL_MODE" = false ]; then
        CONTAINER_MODE=true
        echo -e "${YELLOW}No mode specified, defaulting to container mode.${NC}"
    fi
}

# Print header
echo -e "${BLUE}======================================"
echo "CryoProtect PostgreSQL Optimization for Fedora"
echo -e "======================================${NC}"

# Check root access and environment
check_root
check_environment

# Set paths based on mode
if [ "$CONTAINER_MODE" = true ]; then
    echo -e "${BLUE}Running in container mode${NC}"
    PG_CONF_PATH="/var/lib/postgresql/data/postgresql.conf"
    PG_HBA_PATH="/var/lib/postgresql/data/pg_hba.conf"
    PG_DATA_DIR="/var/lib/postgresql/data"
elif [ "$LOCAL_MODE" = true ]; then
    echo -e "${BLUE}Running in local mode${NC}"
    PG_CONF_PATH="/var/lib/pgsql/data/postgresql.conf"
    PG_HBA_PATH="/var/lib/pgsql/data/pg_hba.conf"
    PG_DATA_DIR="/var/lib/pgsql/data"
fi

# Create the directory for PostgreSQL custom configuration
CONF_DIR="./postgres_config"
mkdir -p "$CONF_DIR"

# ===== 1. Generate PostgreSQL Configuration =====

echo -e "\n${BLUE}[1/5]${NC} Generating PostgreSQL configuration..."

# Determine system memory
TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
TOTAL_MEM_MB=$((TOTAL_MEM_KB / 1024))
TOTAL_MEM_GB=$((TOTAL_MEM_MB / 1024))

# Calculate settings based on system memory
if [ "$TOTAL_MEM_GB" -ge 32 ]; then
    # High memory (32GB+)
    SHARED_BUFFERS="8GB"
    EFFECTIVE_CACHE_SIZE="24GB"
    WORK_MEM="64MB"
    MAINTENANCE_WORK_MEM="2GB"
    MAX_CONNECTIONS=200
    WAL_BUFFERS="16MB"
    MAX_WORKER_PROCESSES=8
    MAX_PARALLEL_WORKERS_PER_GATHER=4
    MAX_PARALLEL_WORKERS=8
    RANDOM_PAGE_COST=1.1
elif [ "$TOTAL_MEM_GB" -ge 16 ]; then
    # Medium memory (16GB)
    SHARED_BUFFERS="4GB"
    EFFECTIVE_CACHE_SIZE="12GB"
    WORK_MEM="32MB"
    MAINTENANCE_WORK_MEM="1GB"
    MAX_CONNECTIONS=100
    WAL_BUFFERS="16MB"
    MAX_WORKER_PROCESSES=4
    MAX_PARALLEL_WORKERS_PER_GATHER=2
    MAX_PARALLEL_WORKERS=4
    RANDOM_PAGE_COST=1.1
elif [ "$TOTAL_MEM_GB" -ge 8 ]; then
    # Low memory (8GB)
    SHARED_BUFFERS="2GB"
    EFFECTIVE_CACHE_SIZE="6GB"
    WORK_MEM="16MB"
    MAINTENANCE_WORK_MEM="512MB"
    MAX_CONNECTIONS=50
    WAL_BUFFERS="16MB"
    MAX_WORKER_PROCESSES=2
    MAX_PARALLEL_WORKERS_PER_GATHER=1
    MAX_PARALLEL_WORKERS=2
    RANDOM_PAGE_COST=1.1
else
    # Minimal memory (<8GB)
    SHARED_BUFFERS="1GB"
    EFFECTIVE_CACHE_SIZE="3GB"
    WORK_MEM="8MB"
    MAINTENANCE_WORK_MEM="256MB"
    MAX_CONNECTIONS=30
    WAL_BUFFERS="8MB"
    MAX_WORKER_PROCESSES=1
    MAX_PARALLEL_WORKERS_PER_GATHER=0
    MAX_PARALLEL_WORKERS=1
    RANDOM_PAGE_COST=2.0
fi

# Create optimized PostgreSQL configuration
cat > "$CONF_DIR/postgresql.optimized.conf" << EOL
# CryoProtect PostgreSQL Optimized Configuration for Fedora
# Automatically generated for a system with ${TOTAL_MEM_GB}GB of RAM

# Connection Settings
listen_addresses = '*'
max_connections = ${MAX_CONNECTIONS}
superuser_reserved_connections = 3

# Memory Settings
shared_buffers = ${SHARED_BUFFERS}
effective_cache_size = ${EFFECTIVE_CACHE_SIZE}
work_mem = ${WORK_MEM}
maintenance_work_mem = ${MAINTENANCE_WORK_MEM}
huge_pages = try   # Use huge pages if available (improves performance)
temp_buffers = 16MB
max_prepared_transactions = 0 # Disable if not needed for better performance

# Parallelism Settings
max_worker_processes = ${MAX_WORKER_PROCESSES}
max_parallel_workers_per_gather = ${MAX_PARALLEL_WORKERS_PER_GATHER}
max_parallel_workers = ${MAX_PARALLEL_WORKERS}
max_parallel_maintenance_workers = ${MAX_PARALLEL_WORKERS_PER_GATHER}
dynamic_shared_memory_type = posix  # Works better with SELinux on Fedora

# Write Ahead Log (WAL) Settings
wal_level = replica
max_wal_size = 1GB
min_wal_size = 80MB
checkpoint_completion_target = 0.9
wal_buffers = ${WAL_BUFFERS}
synchronous_commit = off  # Improves performance but with small durability trade-off

# Query Optimization
random_page_cost = ${RANDOM_PAGE_COST}  # Lower for SSDs
effective_io_concurrency = 200  # Higher for SSDs/NVMe
default_statistics_target = 100  # Balanced for OLTP with some analysis

# Fedora/SELinux Specific Settings
unix_socket_directories = '/var/run/postgresql, /tmp'
unix_socket_permissions = 0777
bonjour = off   # Disable bonjour for security

# Logging
log_destination = 'stderr'
logging_collector = on
log_directory = 'log'
log_filename = 'postgresql-%a.log'
log_truncate_on_rotation = on
log_rotation_age = 1d
log_rotation_size = 0
log_line_prefix = '%m [%p] '
log_timezone = 'UTC'

# OLTP optimization 
default_transaction_isolation = 'read committed'
geqo = on
jit = on

# Autovacuum settings - adjusted for OLTP
autovacuum = on
autovacuum_max_workers = ${MAX_WORKER_PROCESSES}
autovacuum_naptime = 10min
autovacuum_vacuum_threshold = 50
autovacuum_analyze_threshold = 50
autovacuum_vacuum_scale_factor = 0.05
autovacuum_analyze_scale_factor = 0.025
autovacuum_vacuum_cost_delay = 20ms
autovacuum_vacuum_cost_limit = 2000

# Fedora/RHEL specific settings
data_directory = '${PG_DATA_DIR}'
timezone = 'UTC'
lc_messages = 'en_US.UTF-8'
lc_monetary = 'en_US.UTF-8'
lc_numeric = 'en_US.UTF-8'
lc_time = 'en_US.UTF-8'

# Security settings for CryoProtect
password_encryption = scram-sha-256
ssl = on
ssl_prefer_server_ciphers = on
ssl_ecdh_curve = 'prime256v1'
EOL

# Create optimized pg_hba.conf
cat > "$CONF_DIR/pg_hba.optimized.conf" << EOL
# CryoProtect pg_hba.conf Optimized Configuration
# TYPE  DATABASE        USER            ADDRESS                 METHOD

# "local" is for Unix domain socket connections only
local   all             all                                     peer
# IPv4 local connections:
host    all             all             127.0.0.1/32            scram-sha-256
# IPv6 local connections:
host    all             all             ::1/128                 scram-sha-256
# Allow connections from Docker/Podman containers
host    all             all             172.16.0.0/12           scram-sha-256
host    all             all             192.168.0.0/16          scram-sha-256
host    all             all             10.0.0.0/8              scram-sha-256
# Allow connections from the app server
host    all             postgres        all                     scram-sha-256
EOL

echo -e "${GREEN}Optimized PostgreSQL configuration generated in $CONF_DIR/${NC}"

# ===== 2. Create Podman Container Configuration =====

echo -e "\n${BLUE}[2/5]${NC} Creating Podman container configuration..."

# Create podman-compose file for PostgreSQL
cat > "$CONF_DIR/postgres-podman-compose.yml" << EOL
version: '3.8'

services:
  postgres:
    image: postgres:14-alpine
    container_name: cryoprotect-postgres
    restart: always
    environment:
      POSTGRES_USER: postgres
      POSTGRES_PASSWORD: postgres
      POSTGRES_DB: postgres
      # Performance optimization flags for container
      PGSSLMODE: prefer
      PGCLIENTENCODING: UTF8
      WORK_MEM: ${WORK_MEM}
      MAINTENANCE_WORK_MEM: ${MAINTENANCE_WORK_MEM}
      PGDATESTYLE: "iso, ymd"
      PGTZ: "UTC"
    command: 
      - "postgres"
      - "-c"
      - "config_file=/etc/postgresql/postgresql.conf"
    ports:
      - "5432:5432"
    volumes:
      - ./postgres_data:/var/lib/postgresql/data:Z
      - ./postgres_config/postgresql.optimized.conf:/etc/postgresql/postgresql.conf:Z
      - ./postgres_config/pg_hba.optimized.conf:/etc/postgresql/pg_hba.conf:Z
    # Map container user to postgres UID and GID
    user: "postgres:postgres"
    security_opt:
      - label=type:container_file_t
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U postgres"]
      interval: 10s
      timeout: 5s
      retries: 5
    networks:
      - postgres-network
    # Apply resource limitations based on system memory
    deploy:
      resources:
        limits:
          cpus: '2.0'
          memory: ${SHARED_BUFFERS}
        reservations:
          cpus: '0.5'
          memory: 512M

networks:
  postgres-network:
    driver: bridge
EOL

# Create helper script for Podman
cat > "./run_optimized_postgres.sh" << EOL
#!/bin/bash
# run_optimized_postgres.sh
#
# Runs PostgreSQL with optimized settings for CryoProtect on Fedora.
#
# Usage:
#   ./run_optimized_postgres.sh [--start|--stop|--restart|--status]

set -e

# Directory setup
POSTGRES_DATA_DIR="./postgres_data"
POSTGRES_CONFIG_DIR="./postgres_config"
mkdir -p "\$POSTGRES_DATA_DIR"

# Default command
COMMAND=\${1:-"--start"}

case "\$COMMAND" in
  --start)
    echo "Starting optimized PostgreSQL container..."
    podman-compose -f "\$POSTGRES_CONFIG_DIR/postgres-podman-compose.yml" up -d
    echo "PostgreSQL container started"
    ;;
  --stop)
    echo "Stopping PostgreSQL container..."
    podman-compose -f "\$POSTGRES_CONFIG_DIR/postgres-podman-compose.yml" down
    echo "PostgreSQL container stopped"
    ;;
  --restart)
    echo "Restarting PostgreSQL container..."
    podman-compose -f "\$POSTGRES_CONFIG_DIR/postgres-podman-compose.yml" restart
    echo "PostgreSQL container restarted"
    ;;
  --status)
    echo "PostgreSQL container status:"
    podman ps -f name=cryoprotect-postgres
    
    # Check if container is running
    if podman ps -f name=cryoprotect-postgres --format '{{.Names}}' | grep -q 'cryoprotect-postgres'; then
      echo "Container is running"
      echo "PostgreSQL connection info:"
      echo "  Host: localhost"
      echo "  Port: 5432"
      echo "  User: postgres"
      echo "  Database: postgres"
      
      # Check if PostgreSQL is responsive
      if podman exec cryoprotect-postgres pg_isready -U postgres > /dev/null 2>&1; then
        echo "PostgreSQL server is accepting connections"
      else
        echo "PostgreSQL server is not responding"
      fi
    else
      echo "Container is not running"
    fi
    ;;
  *)
    echo "Unknown command: \$COMMAND"
    echo "Usage: \$0 [--start|--stop|--restart|--status]"
    exit 1
    ;;
esac
EOL

chmod +x "./run_optimized_postgres.sh"
echo -e "${GREEN}Created helper script: run_optimized_postgres.sh${NC}"

# ===== 3. Create systemd service for PostgreSQL =====

echo -e "\n${BLUE}[3/5]${NC} Creating systemd service..."

# Create systemd service file
cat > "$CONF_DIR/cryoprotect-postgres.service" << EOL
[Unit]
Description=CryoProtect PostgreSQL Container
After=network.target
Documentation=https://github.com/yourusername/cryoprotect

[Service]
Type=forking
User=mushu
WorkingDirectory=/home/mushu/Projects/CryoProtect
ExecStart=/home/mushu/Projects/CryoProtect/run_optimized_postgres.sh --start
ExecStop=/home/mushu/Projects/CryoProtect/run_optimized_postgres.sh --stop
ExecReload=/home/mushu/Projects/CryoProtect/run_optimized_postgres.sh --restart
Restart=on-failure
RestartSec=30s
TimeoutStartSec=180s
TimeoutStopSec=120s

[Install]
WantedBy=multi-user.target
EOL

echo -e "${GREEN}Created systemd service file: $CONF_DIR/cryoprotect-postgres.service${NC}"

# ===== 4. Create PostgreSQL Health Check Tool =====

echo -e "\n${BLUE}[4/5]${NC} Creating PostgreSQL health check tool..."

# Create a health check script
cat > "./check_postgres_health.py" << EOL
#!/usr/bin/env python3
"""
PostgreSQL Health Check Tool for CryoProtect.

This script checks the health and performance of the PostgreSQL database
and provides optimization recommendations based on the results.
"""

import os
import sys
import psycopg2
import psycopg2.extras
import json
import subprocess
from datetime import datetime
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Database connection parameters
DB_HOST = os.getenv("SUPABASE_DB_HOST", "localhost")
DB_PORT = os.getenv("SUPABASE_DB_PORT", "5432")
DB_NAME = os.getenv("SUPABASE_DB_NAME", "postgres")
DB_USER = os.getenv("SUPABASE_DB_USER", "postgres")
DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD", "postgres")

class PostgresHealthCheck:
    """PostgreSQL health check tool."""
    
    def __init__(self):
        """Initialize the health check tool."""
        self.conn = None
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "connection": {"status": "unknown"},
            "version": "unknown",
            "settings": {},
            "statistics": {},
            "issues": [],
            "recommendations": []
        }
        
    def connect(self):
        """Connect to the PostgreSQL database."""
        try:
            self.conn = psycopg2.connect(
                host=DB_HOST,
                port=DB_PORT,
                dbname=DB_NAME,
                user=DB_USER,
                password=DB_PASSWORD
            )
            self.results["connection"]["status"] = "ok"
            self.results["connection"]["details"] = {
                "host": DB_HOST,
                "port": DB_PORT,
                "dbname": DB_NAME,
                "user": DB_USER
            }
            return True
        except Exception as e:
            self.results["connection"]["status"] = "failed"
            self.results["connection"]["error"] = str(e)
            return False
            
    def check_version(self):
        """Check the PostgreSQL version."""
        try:
            cursor = self.conn.cursor()
            cursor.execute("SELECT version()")
            version = cursor.fetchone()[0]
            self.results["version"] = version
            cursor.close()
        except Exception as e:
            self.results["issues"].append({
                "severity": "error",
                "message": f"Failed to get PostgreSQL version: {e}"
            })
            
    def check_settings(self):
        """Check important PostgreSQL settings."""
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
            
            # Get the most important settings
            settings_query = """
                SELECT name, setting, unit, context, short_desc
                FROM pg_settings
                WHERE name IN (
                    'shared_buffers', 'work_mem', 'maintenance_work_mem',
                    'effective_cache_size', 'max_connections',
                    'random_page_cost', 'effective_io_concurrency',
                    'max_worker_processes', 'max_parallel_workers',
                    'max_parallel_workers_per_gather', 'autovacuum',
                    'default_statistics_target', 'synchronous_commit',
                    'wal_buffers', 'checkpoint_completion_target',
                    'max_wal_size', 'min_wal_size'
                )
                ORDER BY name;
            """
            cursor.execute(settings_query)
            settings = cursor.fetchall()
            
            for setting in settings:
                self.results["settings"][setting["name"]] = {
                    "value": setting["setting"],
                    "unit": setting["unit"],
                    "context": setting["context"],
                    "description": setting["short_desc"]
                }
                
            cursor.close()
        except Exception as e:
            self.results["issues"].append({
                "severity": "error",
                "message": f"Failed to check PostgreSQL settings: {e}"
            })
            
    def check_statistics(self):
        """Check database statistics."""
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
            
            # Database size
            cursor.execute("SELECT pg_size_pretty(pg_database_size(current_database())) as size")
            db_size = cursor.fetchone()
            self.results["statistics"]["database_size"] = db_size["size"]
            
            # Table counts
            cursor.execute("SELECT COUNT(*) as count FROM information_schema.tables WHERE table_schema = 'public'")
            table_count = cursor.fetchone()
            self.results["statistics"]["table_count"] = table_count["count"]
            
            # Connection statistics
            cursor.execute("""
                SELECT 
                    count(*) AS active_connections,
                    sum(CASE WHEN state = 'active' THEN 1 ELSE 0 END) AS active_queries,
                    sum(CASE WHEN state = 'idle' THEN 1 ELSE 0 END) AS idle_connections
                FROM pg_stat_activity
            """)
            connections = cursor.fetchone()
            self.results["statistics"]["connections"] = connections
            
            # Cache hit ratio
            cursor.execute("""
                SELECT 
                    sum(heap_blks_read) as heap_read,
                    sum(heap_blks_hit) as heap_hit,
                    sum(heap_blks_hit) / (sum(heap_blks_hit) + sum(heap_blks_read)) as ratio
                FROM pg_statio_user_tables
            """)
            cache = cursor.fetchone()
            self.results["statistics"]["cache_hit_ratio"] = cache["ratio"] if cache["ratio"] is not None else 0
            
            # Get largest tables
            cursor.execute("""
                SELECT
                    tablename AS table_name,
                    pg_size_pretty(pg_total_relation_size('"' || tablename || '"')) AS size
                FROM pg_tables
                WHERE schemaname = 'public'
                ORDER BY pg_total_relation_size('"' || tablename || '"') DESC
                LIMIT 5
            """)
            largest_tables = cursor.fetchall()
            self.results["statistics"]["largest_tables"] = largest_tables
            
            cursor.close()
        except Exception as e:
            self.results["issues"].append({
                "severity": "error",
                "message": f"Failed to check database statistics: {e}"
            })
            
    def analyze_performance(self):
        """Analyze database performance and suggest optimizations."""
        try:
            # Check cache hit ratio
            if "cache_hit_ratio" in self.results["statistics"]:
                ratio = self.results["statistics"]["cache_hit_ratio"]
                if ratio < 0.95:
                    self.results["issues"].append({
                        "severity": "warning",
                        "message": f"Cache hit ratio is low: {ratio:.2%}",
                        "recommendation": "Consider increasing shared_buffers"
                    })
                    self.results["recommendations"].append(
                        "Increase shared_buffers to improve cache hit ratio"
                    )
                    
            # Check connections
            if "connections" in self.results["statistics"]:
                conn_data = self.results["statistics"]["connections"]
                max_conn = int(self.results["settings"].get("max_connections", {}).get("value", 100))
                conn_percent = conn_data["active_connections"] / max_conn
                
                if conn_percent > 0.7:
                    self.results["issues"].append({
                        "severity": "warning",
                        "message": f"High connection usage: {conn_percent:.2%} of max_connections",
                        "recommendation": "Consider implementing connection pooling"
                    })
                    self.results["recommendations"].append(
                        "Implement connection pooling to reduce the number of active connections"
                    )
                    
            # Check for missing indexes on foreign keys
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
            cursor.execute("""
                SELECT
                    conrelid::regclass AS table_from,
                    a.attname AS column,
                    confrelid::regclass AS table_to,
                    af.attname AS ftable_column,
                    (SELECT COUNT(1) FROM pg_index i WHERE indrelid=conrelid AND (a.attnum = ANY(i.indkey))) AS has_index
                FROM
                    pg_constraint c
                    JOIN pg_attribute a ON a.attnum = ANY(c.conkey) AND a.attrelid = c.conrelid
                    JOIN pg_attribute af ON af.attnum = ANY(c.confkey) AND af.attrelid = c.confrelid
                WHERE
                    c.contype = 'f'
                    AND (SELECT COUNT(1) FROM pg_index i WHERE indrelid=conrelid AND (a.attnum = ANY(i.indkey))) = 0
                ORDER BY
                    conrelid::regclass::text,
                    a.attname;
            """)
            missing_indexes = cursor.fetchall()
            
            if missing_indexes:
                self.results["issues"].append({
                    "severity": "warning",
                    "message": f"Found {len(missing_indexes)} foreign keys without indexes",
                    "details": [f"{row['table_from']}.{row['column']} -> {row['table_to']}.{row['ftable_column']}" 
                               for row in missing_indexes[:5]]
                })
                self.results["recommendations"].append(
                    "Add indexes to foreign key columns to improve join performance"
                )
            
            cursor.close()
            
        except Exception as e:
            self.results["issues"].append({
                "severity": "error",
                "message": f"Failed to analyze performance: {e}"
            })
            
    def run_checks(self):
        """Run all health checks."""
        if not self.connect():
            return
            
        self.check_version()
        self.check_settings()
        self.check_statistics()
        self.analyze_performance()
        
        # Close connection
        if self.conn:
            self.conn.close()
            
    def print_report(self):
        """Print the health check report to the console."""
        print("=" * 70)
        print("PostgreSQL Health Check Report")
        print("=" * 70)
        print(f"Timestamp: {self.results['timestamp']}")
        print(f"Connection: {self.results['connection']['status']}")
        print(f"Version: {self.results['version']}")
        print()
        
        print("Database Statistics:")
        print(f"  Database Size: {self.results['statistics'].get('database_size', 'Unknown')}")
        print(f"  Table Count: {self.results['statistics'].get('table_count', 'Unknown')}")
        if "connections" in self.results["statistics"]:
            conn = self.results["statistics"]["connections"]
            print(f"  Active Connections: {conn.get('active_connections', 'Unknown')}")
            print(f"  Active Queries: {conn.get('active_queries', 'Unknown')}")
        if "cache_hit_ratio" in self.results["statistics"]:
            print(f"  Cache Hit Ratio: {self.results['statistics']['cache_hit_ratio']:.2%}")
        print()
        
        print("Key Settings:")
        for name, details in self.results["settings"].items():
            value = details["value"]
            unit = details["unit"] if details["unit"] else ""
            print(f"  {name}: {value}{unit}")
        print()
        
        if self.results["issues"]:
            print("Issues Found:")
            for issue in self.results["issues"]:
                severity = issue["severity"].upper()
                message = issue["message"]
                print(f"  [{severity}] {message}")
                if "details" in issue:
                    for detail in issue["details"]:
                        print(f"    - {detail}")
            print()
            
        if self.results["recommendations"]:
            print("Recommendations:")
            for i, rec in enumerate(self.results["recommendations"], 1):
                print(f"  {i}. {rec}")
            print()
        
        print("=" * 70)
        
    def save_report(self, filename):
        """Save the health check report to a file."""
        with open(filename, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"Report saved to {filename}")

def main():
    """Main function."""
    health_check = PostgresHealthCheck()
    health_check.run_checks()
    health_check.print_report()
    
    # Save report to file
    report_dir = "reports"
    os.makedirs(report_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{report_dir}/postgres_health_check_{timestamp}.json"
    health_check.save_report(filename)
    
    # Create symlink to latest report
    latest_link = f"{report_dir}/postgres_health_latest.json"
    if os.path.exists(latest_link):
        os.remove(latest_link)
    os.symlink(filename, latest_link)
    
    # Return exit code based on issues
    critical_issues = sum(1 for issue in health_check.results["issues"] if issue["severity"] == "error")
    return 1 if critical_issues > 0 else 0

if __name__ == "__main__":
    sys.exit(main())
EOL

chmod +x "./check_postgres_health.py"
echo -e "${GREEN}Created PostgreSQL health check tool: check_postgres_health.py${NC}"

# ===== 5. Create Documentation =====

echo -e "\n${BLUE}[5/5]${NC} Creating documentation..."

# Create README file
cat > "./POSTGRESQL_FEDORA_GUIDE.md" << EOL
# PostgreSQL Optimization for CryoProtect on Fedora

This guide describes the optimized PostgreSQL configuration for running CryoProtect on Fedora Linux.

## Configuration Overview

The configuration has been optimized based on your system specifications:
- Total Memory: ${TOTAL_MEM_GB}GB
- Shared Buffers: ${SHARED_BUFFERS}
- Effective Cache Size: ${EFFECTIVE_CACHE_SIZE}
- Work Memory: ${WORK_MEM}
- Maintenance Work Memory: ${MAINTENANCE_WORK_MEM}

## Getting Started

### Option 1: Run with Podman (Recommended)

1. Start the optimized PostgreSQL container:
   \`\`\`bash
   ./run_optimized_postgres.sh --start
   \`\`\`

2. Check the status of the container:
   \`\`\`bash
   ./run_optimized_postgres.sh --status
   \`\`\`

3. Stop the container:
   \`\`\`bash
   ./run_optimized_postgres.sh --stop
   \`\`\`

### Option 2: Install as a Systemd Service

1. Install the systemd service:
   \`\`\`bash
   sudo cp ${CONF_DIR}/cryoprotect-postgres.service /etc/systemd/system/
   sudo systemctl daemon-reload
   sudo systemctl enable cryoprotect-postgres.service
   sudo systemctl start cryoprotect-postgres.service
   \`\`\`

2. Check service status:
   \`\`\`bash
   sudo systemctl status cryoprotect-postgres.service
   \`\`\`

## Health Checks and Monitoring

Run the health check tool to analyze your PostgreSQL performance:

\`\`\`bash
./check_postgres_health.py
\`\`\`

This will generate a report with statistics and recommendations.

## Configuration Files

- **PostgreSQL Configuration**: \`${CONF_DIR}/postgresql.optimized.conf\`
- **Client Authentication**: \`${CONF_DIR}/pg_hba.optimized.conf\`
- **Podman Compose**: \`${CONF_DIR}/postgres-podman-compose.yml\`

## Performance Optimization Details

### Memory Settings

- **shared_buffers**: ${SHARED_BUFFERS}
  Determines how much memory is dedicated to PostgreSQL for data caching.

- **effective_cache_size**: ${EFFECTIVE_CACHE_SIZE}
  Helps PostgreSQL estimate how much memory is available for disk caching.

- **work_mem**: ${WORK_MEM}
  Memory used for sort operations and hash tables per operation.

- **maintenance_work_mem**: ${MAINTENANCE_WORK_MEM}
  Memory used for maintenance operations like VACUUM and CREATE INDEX.

### Parallelism Settings

- **max_worker_processes**: ${MAX_WORKER_PROCESSES}
- **max_parallel_workers_per_gather**: ${MAX_PARALLEL_WORKERS_PER_GATHER}
- **max_parallel_workers**: ${MAX_PARALLEL_WORKERS}

### Workload Optimization

The configuration is optimized for OLTP (Online Transaction Processing) with some analytical workloads, 
which is the typical pattern for CryoProtect application.

## Fedora-Specific Considerations

### SELinux

If you encounter SELinux-related issues, you might need to:

1. Check the SELinux context of the data directory:
   \`\`\`bash
   ls -ldZ ./postgres_data
   \`\`\`

2. Set the correct SELinux context:
   \`\`\`bash
   sudo chcon -Rt container_file_t ./postgres_data
   \`\`\`

### Firewall

If connecting remotely, ensure the PostgreSQL port is open:

\`\`\`bash
sudo firewall-cmd --permanent --add-port=5432/tcp
sudo firewall-cmd --reload
\`\`\`

## Troubleshooting

If you encounter issues:

1. Check PostgreSQL logs:
   \`\`\`bash
   podman logs cryoprotect-postgres
   \`\`\`

2. Run the health check tool:
   \`\`\`bash
   ./check_postgres_health.py
   \`\`\`

3. Verify the container is running:
   \`\`\`bash
   podman ps -a
   \`\`\`

4. Check for SELinux denials:
   \`\`\`bash
   sudo ausearch -m avc --start recent
   \`\`\`
EOL

echo -e "${GREEN}Created documentation file: POSTGRESQL_FEDORA_GUIDE.md${NC}"

# ===== Apply Settings (if requested) =====

if [ "$APPLY" = true ]; then
    echo -e "\n${BLUE}Applying PostgreSQL optimizations...${NC}"
    
    if [ "$CONTAINER_MODE" = true ]; then
        # Create data directory if it doesn't exist
        if [ ! -d "./postgres_data" ]; then
            echo "Creating PostgreSQL data directory..."
            mkdir -p "./postgres_data"
            chmod 700 "./postgres_data"
        fi
        
        # Copy configuration files to proper directory
        echo "Applying container configuration..."
        
        # Run container if requested
        echo "Starting PostgreSQL container with optimized settings..."
        ./run_optimized_postgres.sh --start
        
        echo -e "${GREEN}PostgreSQL container started with optimized settings${NC}"
    elif [ "$LOCAL_MODE" = true ]; then
        # Copy configuration files to system location
        echo "Applying system configuration..."
        if [ -f "$PG_CONF_PATH" ]; then
            echo "Backing up existing PostgreSQL configuration..."
            cp "$PG_CONF_PATH" "${PG_CONF_PATH}.backup.$(date +%Y%m%d%H%M%S)"
        fi
        
        if [ -f "$PG_HBA_PATH" ]; then
            echo "Backing up existing HBA configuration..."
            cp "$PG_HBA_PATH" "${PG_HBA_PATH}.backup.$(date +%Y%m%d%H%M%S)"
        fi
        
        echo "Copying optimized configuration files..."
        cp "$CONF_DIR/postgresql.optimized.conf" "$PG_CONF_PATH"
        cp "$CONF_DIR/pg_hba.optimized.conf" "$PG_HBA_PATH"
        
        echo "Setting correct permissions..."
        chown postgres:postgres "$PG_CONF_PATH" "$PG_HBA_PATH"
        chmod 600 "$PG_CONF_PATH" "$PG_HBA_PATH"
        
        echo "Restarting PostgreSQL service..."
        systemctl restart postgresql
        
        echo -e "${GREEN}PostgreSQL optimized configuration applied to system${NC}"
    fi
else
    echo -e "\n${YELLOW}Skipping application of settings (use --apply to apply)${NC}"
fi

# ===== Summary =====

echo -e "\n${BLUE}======================================"
echo "PostgreSQL Optimization Summary"
echo -e "======================================${NC}"
echo -e "${GREEN}✓${NC} Generated optimized PostgreSQL configuration"
echo -e "${GREEN}✓${NC} Created Podman container configuration"
echo -e "${GREEN}✓${NC} Created systemd service file"
echo -e "${GREEN}✓${NC} Created PostgreSQL health check tool"
echo -e "${GREEN}✓${NC} Created documentation"

if [ "$APPLY" = true ]; then
    if [ "$CONTAINER_MODE" = true ]; then
        echo -e "${GREEN}✓${NC} Started PostgreSQL container with optimized settings"
    elif [ "$LOCAL_MODE" = true ]; then
        echo -e "${GREEN}✓${NC} Applied optimized settings to system PostgreSQL"
    fi
fi

echo -e "\n${BLUE}Next Steps:${NC}"
echo "1. Review the configuration in $CONF_DIR"
echo "2. Check PostgreSQL health with ./check_postgres_health.py"
echo "3. Read the documentation in POSTGRESQL_FEDORA_GUIDE.md"

echo -e "\n${GREEN}PostgreSQL optimization for Fedora completed successfully.${NC}"