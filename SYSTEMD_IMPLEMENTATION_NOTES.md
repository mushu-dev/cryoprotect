# CryoProtect Systemd Implementation Notes

## Overview

This document provides important notes about the systemd service implementation for CryoProtect on Fedora Linux.

## Database Configuration

During testing, we encountered permission issues with running PostgreSQL in a Podman container. These issues are related to SELinux constraints in the Fedora environment. The following alternatives are available:

1. **Use Existing PostgreSQL Service**: If PostgreSQL is already running on the system, you can configure CryoProtect to use it directly. This is the recommended approach for development environments.

2. **Run PostgreSQL with SELinux in Permissive Mode**: For production environments, you can run PostgreSQL in a container with SELinux in permissive mode:
   ```bash
   sudo setenforce 0
   ./run_optimized_postgres.sh --start
   sudo setenforce 1
   ```

3. **Create Custom SELinux Policy**: For the most secure approach, create a custom SELinux policy for the CryoProtect containers:
   ```bash
   sudo ausearch -m avc --start recent | audit2allow -M cryoprotect
   sudo semodule -i cryoprotect.pp
   ```

## Environment Variables

When running with the existing PostgreSQL instance, make sure to set the correct environment variables:

```bash
export SUPABASE_DB_HOST=localhost
export SUPABASE_DB_PORT=5432  # Use default PostgreSQL port
export SUPABASE_DB_NAME=postgres
export SUPABASE_DB_USER=postgres
export SUPABASE_DB_PASSWORD=postgres
```

## Systemd Service Configuration

The provided systemd services are configured as follows:

1. **cryoprotect-app.service**: Main Flask application
   - Depends on PostgreSQL service
   - Automatically restarts on failure
   - Uses conda environment for Python dependencies

2. **cryoprotect-db-integrity-check.service**: Database integrity checker
   - Runs daily via timer (cryoprotect-db-integrity-check.timer)
   - Generates reports in the reports/ directory

3. **cryoprotect-postgres-health.service**: Database health monitor
   - Runs hourly via timer (cryoprotect-postgres-health.timer)
   - Monitors PostgreSQL performance and generates recommendations

## Installation Options

The install_systemd_services.sh script provides two installation modes:

1. **User Mode**: Services run under the current user and start on user login
   ```bash
   ./install_systemd_services.sh --user
   ```

2. **System Mode**: Services run system-wide and start on system boot (requires sudo)
   ```bash
   sudo ./install_systemd_services.sh --system
   ```

## Troubleshooting

If you encounter issues with the systemd services:

1. Check service status:
   ```bash
   systemctl --user status cryoprotect-app.service
   ```

2. View logs:
   ```bash
   journalctl --user-unit=cryoprotect-app.service -f
   ```

3. Verify PostgreSQL connectivity:
   ```bash
   psql -h localhost -p 5432 -U postgres -d postgres -c "\l"
   ```

4. Check SELinux denials:
   ```bash
   sudo ausearch -m avc --start recent
   ```

## Future Improvements

For future versions, consider implementing:

1. Socket activation for on-demand service startup
2. Resource limits for improved container isolation
3. Backup service with automatic rotation
4. Monitoring integration with Prometheus/Grafana
5. Log aggregation with Elasticsearch/Kibana