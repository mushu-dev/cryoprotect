# PostgreSQL Optimization for CryoProtect on Fedora

This guide describes the optimized PostgreSQL configuration for running CryoProtect on Fedora Linux.

## Configuration Overview

The configuration has been optimized based on your system specifications:
- Total Memory: 15GB
- Shared Buffers: 2GB
- Effective Cache Size: 6GB
- Work Memory: 16MB
- Maintenance Work Memory: 512MB

## Getting Started

### Option 1: Run with Podman (Recommended)

1. Start the optimized PostgreSQL container:
   ```bash
   ./run_optimized_postgres.sh --start
   ```

2. Check the status of the container:
   ```bash
   ./run_optimized_postgres.sh --status
   ```

3. Stop the container:
   ```bash
   ./run_optimized_postgres.sh --stop
   ```

### Option 2: Install as a Systemd Service

1. Install the systemd service:
   ```bash
   sudo cp ./postgres_config/cryoprotect-postgres.service /etc/systemd/system/
   sudo systemctl daemon-reload
   sudo systemctl enable cryoprotect-postgres.service
   sudo systemctl start cryoprotect-postgres.service
   ```

2. Check service status:
   ```bash
   sudo systemctl status cryoprotect-postgres.service
   ```

## Health Checks and Monitoring

Run the health check tool to analyze your PostgreSQL performance:

```bash
./check_postgres_health.py
```

This will generate a report with statistics and recommendations.

## Configuration Files

- **PostgreSQL Configuration**: `./postgres_config/postgresql.optimized.conf`
- **Client Authentication**: `./postgres_config/pg_hba.optimized.conf`
- **Podman Compose**: `./postgres_config/postgres-podman-compose.yml`

## Performance Optimization Details

### Memory Settings

- **shared_buffers**: 2GB
  Determines how much memory is dedicated to PostgreSQL for data caching.

- **effective_cache_size**: 6GB
  Helps PostgreSQL estimate how much memory is available for disk caching.

- **work_mem**: 16MB
  Memory used for sort operations and hash tables per operation.

- **maintenance_work_mem**: 512MB
  Memory used for maintenance operations like VACUUM and CREATE INDEX.

### Parallelism Settings

- **max_worker_processes**: 2
- **max_parallel_workers_per_gather**: 1
- **max_parallel_workers**: 2

### Workload Optimization

The configuration is optimized for OLTP (Online Transaction Processing) with some analytical workloads, 
which is the typical pattern for CryoProtect application.

## Fedora-Specific Considerations

### SELinux

If you encounter SELinux-related issues, you might need to:

1. Check the SELinux context of the data directory:
   ```bash
   ls -ldZ ./postgres_data
   ```

2. Set the correct SELinux context:
   ```bash
   sudo chcon -Rt container_file_t ./postgres_data
   ```

### Firewall

If connecting remotely, ensure the PostgreSQL port is open:

```bash
sudo firewall-cmd --permanent --add-port=5432/tcp
sudo firewall-cmd --reload
```

## Troubleshooting

If you encounter issues:

1. Check PostgreSQL logs:
   ```bash
   podman logs cryoprotect-postgres
   ```

2. Run the health check tool:
   ```bash
   ./check_postgres_health.py
   ```

3. Verify the container is running:
   ```bash
   podman ps -a
   ```

4. Check for SELinux denials:
   ```bash
   sudo ausearch -m avc --start recent
   ```
