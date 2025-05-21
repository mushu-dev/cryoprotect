# Database Maintenance Guide

This document describes the database maintenance system implemented for the CryoProtect application.

## Overview

The database maintenance system is designed to keep the PostgreSQL database running efficiently with minimal manual intervention. It includes:

1. Automatic vacuuming of tables with many dead tuples
2. Regular analysis of tables to update statistics for the query planner
3. Reindexing of important tables to maintain index efficiency
4. Monitoring and cancellation of long-running queries
5. Cleaning of old cache invalidation events, sessions, and logs
6. Identification of database issues like unused indexes or missing foreign key indexes

## Maintenance Tasks

The system performs two types of maintenance runs:

### Standard Maintenance (Daily)

Runs every weekday and includes:
- Vacuum tables with high dead tuple percentages
- Analyze tables to update statistics
- Cancel queries running longer than 10 minutes
- Clean up old invalidation events, sessions, and audit logs
- Identify potential database issues

### Full Maintenance (Weekly)

Runs on weekends and includes all standard maintenance tasks plus:
- Reindex important tables
- More aggressive vacuuming
- Full database statistics analysis

## Scheduling

The maintenance system is scheduled using cron:

- **Standard Maintenance**: Runs daily at 3:00 AM on weekdays
- **Full Maintenance**: Runs at 3:00 AM on weekends

## Manual Execution

You can manually run maintenance tasks using:

```bash
# Standard maintenance
./run_database_maintenance.sh

# Full maintenance
./run_database_maintenance.sh --full
```

## Maintenance Scripts

The system consists of several files:

- `database/maintenance.py` - Core maintenance functionality and database utilities
- `run_database_maintenance.sh` - Shell script that runs maintenance tasks
- `setup_maintenance_cron.sh` - Script to set up the cron job

## Cache Processor Service

The cache processor is designed to run as a systemd service:

```bash
# Start the cache processor service
sudo systemctl start cryoprotect-cache-processor.service

# Check status
sudo systemctl status cryoprotect-cache-processor.service

# Stop the service
sudo systemctl stop cryoprotect-cache-processor.service
```

## Logs and Reports

Each maintenance run generates:

1. A log file in `logs/maintenance/maintenance_TIMESTAMP.log`
2. A detailed report in `logs/maintenance/report_TIMESTAMP.txt`

Logs older than 30 days are automatically deleted.

## Monitoring Maintenance

To check if maintenance is running properly:

1. Check the log files in `logs/maintenance/`
2. Look for any alert files in `logs/maintenance/alert_*.txt`
3. Run `systemctl status cryoprotect-cache-processor.service` to verify the cache processor is running

## Database Issues Detection

The maintenance system can detect several common database issues:

- **Unused Indexes**: Indexes with very few scans that consume space
- **Duplicate Indexes**: Multiple indexes covering the same columns
- **Foreign Keys Without Indexes**: Foreign keys that should have indexes for performance
- **Table Bloat**: Tables with excessive dead space
- **Long-Running Queries**: Queries that might be problematic

## Best Practices

1. **Review Maintenance Reports**: Regularly check the maintenance reports for issues
2. **Monitor Disk Space**: Ensure adequate disk space is available for database growth
3. **Address Detected Issues**: Fix issues like unused indexes or foreign keys without indexes
4. **Adjust Maintenance Schedule**: Modify the cron schedule if needed for your server's peak hours
5. **Backup Database**: Always ensure you have proper backups before major maintenance

## Maintenance Parameters

Several parameters can be adjusted in `database/maintenance.py`:

- `DEFAULT_VACUUM_THRESHOLD`: Vacuum tables with dead tuples above this percentage (default: 20%)
- `DEFAULT_ANALYZE_THRESHOLD`: Analyze tables with modified tuples above this percentage (default: 10%)
- `DEFAULT_MAX_RUNTIME`: Maximum runtime for maintenance in seconds (default: 1 hour)
- `DEFAULT_QUERY_TIMEOUT`: Timeout for long-running queries in seconds (default: 10 minutes)

## Troubleshooting

### Maintenance Fails to Run

1. Check if PostgreSQL is running: `systemctl status postgresql`
2. Verify cron jobs are properly set up: `crontab -l`
3. Check logs for specific errors

### Cache Processor Issues

1. Check if Redis is running: `systemctl status redis`
2. Verify the service is running: `systemctl status cryoprotect-cache-processor.service`
3. Check logs for specific errors

### Database Performance Issues

If performance issues persist after maintenance:

1. Check for long-running queries: `python3 -m database.maintenance --report`
2. Consider adjusting vacuum and analyze thresholds
3. Review PostgreSQL configuration for potential improvements