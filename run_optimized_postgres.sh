#!/bin/bash
# run_optimized_postgres.sh
#
# Runs PostgreSQL with optimized settings for CryoProtect on Fedora.
#
# Usage:
#   ./run_optimized_postgres.sh [--start|--stop|--restart|--status]

set -e

# Get the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Directory setup
POSTGRES_DATA_DIR="${SCRIPT_DIR}/postgres_data"
POSTGRES_CONFIG_DIR="${SCRIPT_DIR}/postgres_config"
mkdir -p "$POSTGRES_DATA_DIR"

# Default command
COMMAND=${1:-"--start"}

case "$COMMAND" in
  --start)
    echo "Starting optimized PostgreSQL container..."
    # Create network if it doesn't exist
    podman network exists postgres-network || podman network create postgres-network

    # Run PostgreSQL container directly
    podman run -d --name cryoprotect-postgres \
      --restart always \
      -e POSTGRES_USER=postgres \
      -e POSTGRES_PASSWORD=postgres \
      -e POSTGRES_DB=postgres \
      -e PGSSLMODE=prefer \
      -e PGCLIENTENCODING=UTF8 \
      -e PGDATESTYLE="iso, ymd" \
      -e PGTZ=UTC \
      -p 5433:5432 \
      -v "${POSTGRES_DATA_DIR}:/var/lib/postgresql/data" \
      -v "${POSTGRES_CONFIG_DIR}/postgresql.optimized.conf:/etc/postgresql/postgresql.conf" \
      -v "${POSTGRES_CONFIG_DIR}/pg_hba.optimized.conf:/etc/postgresql/pg_hba.conf" \
      --network postgres-network \
      docker.io/library/postgres:14-alpine \
      postgres -c config_file=/etc/postgresql/postgresql.conf

    echo "PostgreSQL container started"
    ;;
  --stop)
    echo "Stopping PostgreSQL container..."
    podman stop cryoprotect-postgres && podman rm -f cryoprotect-postgres
    echo "PostgreSQL container stopped"
    ;;
  --restart)
    echo "Restarting PostgreSQL container..."
    podman restart cryoprotect-postgres
    echo "PostgreSQL container restarted"
    ;;
  --status)
    echo "PostgreSQL container status:"
    podman ps -f name=cryoprotect-postgres

    # Check if container is running
    if podman ps | grep -q 'cryoprotect-postgres'; then
      echo "Container is running"
      echo "PostgreSQL connection info:"
      echo "  Host: localhost"
      echo "  Port: 5433"
      echo "  User: postgres"
      echo "  Database: postgres"
      
      # Check if PostgreSQL is responsive (with a retry)
      for i in {1..5}; do
        if podman exec cryoprotect-postgres pg_isready -U postgres > /dev/null 2>&1; then
          echo "PostgreSQL server is accepting connections"
          break
        else
          echo "Waiting for PostgreSQL server to start (attempt $i/5)..."
          sleep 2
        fi

        if [ $i -eq 5 ]; then
          echo "PostgreSQL server is not responding after multiple attempts"
        fi
      done
    else
      echo "Container is not running"
    fi
    ;;
  *)
    echo "Unknown command: $COMMAND"
    echo "Usage: $0 [--start|--stop|--restart|--status]"
    exit 1
    ;;
esac
