# CryoProtect v2 Database Maintenance Guide

## Overview

This document provides a comprehensive guide for maintaining the CryoProtect v2 database. It covers routine maintenance tasks, troubleshooting common issues, and best practices to ensure the database remains stable, performant, and secure.

## Routine Maintenance Tasks

### 1. Database Integrity Checks (Weekly)

- **Purpose**: Ensure data consistency and identify potential corruption.
- **Tasks**:
  - Run `VACUUM ANALYZE` to reclaim storage and update query planner statistics.
  - Check for orphaned records using SQL queries (e.g., identify `mixture_components` without corresponding `mixtures` or `molecules`).
  - Verify foreign key constraints are intact and not violated.
- **Scripts**:
  - `batch_scripts/run_weekly_maintenance.sh` (example script to automate tasks)
  - `supabase_database_audit.py` (script for database auditing and integrity checks)
- **Verification**:
  - Review logs for any errors or warnings during maintenance tasks.
  - Check database statistics for anomalies (e.g., table sizes, index usage).

### 2. Performance Monitoring (Monthly)

- **Purpose**: Track database performance and identify potential bottlenecks.
- **Tasks**:
  - Monitor query performance using Supabase dashboard or pgAdmin.
  - Identify slow queries using query logs or performance monitoring tools.
  - Analyze query plans for inefficient execution strategies.
  - Check index usage and identify missing or redundant indexes.
- **Tools**:
  - Supabase Performance Dashboard
  - pgAdmin Query Analyzer
  - `test_database_performance.py` (script for performance benchmarking)
- **Metrics to Monitor**:
  - Average query execution time
  - P95 and P99 response times
  - Query throughput
  - CPU and memory usage
  - Disk I/O
- **Actions based on Monitoring**:
  - Optimize slow queries by rewriting or adding indexes.
  - Review and adjust database configuration parameters.
  - Consider database scaling if performance degrades consistently.

### 3. Security Audits (Quarterly)

- **Purpose**: Ensure database security and compliance with security policies.
- **Tasks**:
  - Review Row Level Security (RLS) policies for all tables.
  - Verify user roles and permissions are correctly configured.
  - Audit authentication logs for suspicious activity.
  - Check for any security vulnerabilities in database extensions or configurations.
  - Review service role usage and ensure it's only used for necessary tasks.
- **Tools**:
  - Supabase Security Dashboard
  - `supabase_database_audit.py` (script for security policy auditing)
- **Verification**:
  - Document audit findings and remediation actions.
  - Update security policies and procedures as needed.

### 4. Backup and Recovery Testing (Monthly)

- **Purpose**: Ensure backups are reliable and recovery procedures are effective.
- **Tasks**:
  - Perform regular database backups (daily or more frequently for critical data).
  - Test backup restoration to a staging environment.
  - Verify data integrity after restoration.
  - Document backup and recovery procedures.
- **Scripts**:
  - `create_production_backup.sh` (script for creating database backups)
  - `batch_scripts/test_backup_restore.sh` (example script for testing restore)
- **Verification**:
  - Confirm successful backup creation and restoration.
  - Compare restored data with backup data to ensure consistency.
  - Document any issues encountered and resolutions.

### 5. Schema Evolution Management

- **Purpose**: Manage database schema changes in a controlled and documented manner.
- **Tasks**:
  - Use migration scripts for all schema changes (DDL operations).
  - Review and test migration scripts in a development environment.
  - Apply migrations to staging and production environments in a phased approach.
  - Document all schema changes and their impact.
  - Maintain backward compatibility when possible.
- **Tools**:
  - Supabase Migrations CLI
  - `migrations/` directory (for storing migration scripts)
- **Best Practices**:
  - Version control migration scripts.
  - Use descriptive names for migration files.
  - Test migrations in non-production environments first.
  - Implement rollback procedures for migrations.

## Troubleshooting Common Issues

### 1. Slow Query Performance

- **Symptoms**: API endpoints are slow, application response times are high, database CPU usage is high.
- **Possible Causes**:
  - Missing indexes on frequently queried columns.
  - Inefficient query design (e.g., complex joins, full table scans).
  - Outdated query planner statistics.
  - Database server resource limitations.
- **Troubleshooting Steps**:
  1. Identify slow queries using query logs or performance monitoring tools.
  2. Analyze query execution plans using `EXPLAIN ANALYZE`.
  3. Add missing indexes on relevant columns.
  4. Rewrite inefficient queries to optimize performance.
  5. Run `VACUUM ANALYZE` to update query planner statistics.
  6. Check database server resource usage (CPU, memory, disk I/O).
  7. Consider database scaling or read replicas for high load scenarios.
- **Tools**:
  - Supabase Performance Dashboard
  - pgAdmin Query Analyzer
  - `test_database_performance.py` (script for performance benchmarking)

### 2. RLS Policy Issues

- **Symptoms**: Users cannot access data they should be able to, or vice versa; unexpected authorization errors.
- **Possible Causes**:
  - Incorrectly configured RLS policies.
  - Missing RLS policies on new tables.
  - Conflicts between RLS policies.
  - Performance issues with complex RLS policies.
- **Troubleshooting Steps**:
  1. Review RLS policies for the affected tables using Supabase dashboard or SQL queries.
  2. Test RLS policies using different user roles and scenarios.
  3. Simplify complex RLS policies if possible.
  4. Ensure indexes are in place to support RLS policy evaluation (e.g., on `created_by` columns).
  5. Check for policy conflicts or overlapping policies.
- **Tools**:
  - Supabase Security Dashboard
  - `supabase_database_audit.py` (script for security policy auditing)

### 3. Database Connection Errors

- **Symptoms**: Application cannot connect to the database; API endpoints return database connection errors.
- **Possible Causes**:
  - Incorrect database connection parameters (URL, credentials).
  - Database server is down or unreachable.
  - Network connectivity issues.
  - Database connection limits exceeded.
- **Troubleshooting Steps**:
  1. Verify database connection parameters in application configuration.
  2. Check database server status using Supabase dashboard or monitoring tools.
  3. Test network connectivity to the database server (ping, traceroute).
  4. Review database connection limits and increase if necessary.
  5. Check database logs for connection errors or server issues.
- **Tools**:
  - Supabase Dashboard (connection status)
  - Network diagnostic tools (ping, traceroute)

### 4. Data Integrity Issues

- **Symptoms**: Data inconsistencies, orphaned records, foreign key violations, incorrect data values.
- **Possible Causes**:
  - Application bugs causing data corruption.
  - Incomplete or incorrect data migration scripts.
  - Manual data manipulation errors.
  - Database integrity constraint violations.
- **Troubleshooting Steps**:
  1. Run database integrity checks (see Routine Maintenance Tasks).
  2. Review application code for potential data corruption issues.
  3. Examine data migration scripts for errors.
  4. Check database logs for constraint violations or errors.
  5. Implement data validation and sanitization in the application.
  6. Restore from backup if data corruption is severe.
- **Scripts**:
  - `supabase_database_audit.py` (script for database auditing and integrity checks)
  - `remediate_database_integrity.py` (script for fixing database integrity issues)

### 5. Backup and Restore Failures

- **Symptoms**: Backups fail to complete; database restoration fails or results in data loss.
- **Possible Causes**:
  - Insufficient storage space for backups.
  - Backup process errors (e.g., permissions, network issues).
  - Corrupted backup files.
  - Incorrect restoration procedures.
- **Troubleshooting Steps**:
  1. Verify sufficient storage space for backups.
  2. Check backup process logs for errors.
  3. Test backup integrity by attempting restoration in a staging environment.
  4. Review and correct backup and restore procedures.
  5. Implement backup monitoring and alerting.
- **Scripts**:
  - `create_production_backup.sh` (script for creating database backups)
  - `batch_scripts/test_backup_restore.sh` (example script for testing restore)

## Best Practices for Database Maintenance

### 1. Automate Routine Tasks

- Automate routine maintenance tasks like integrity checks, performance monitoring, and backups using scripts and scheduled jobs.
- Use tools like cron jobs or Supabase Functions to schedule tasks.
- Implement monitoring and alerting for automated tasks to detect failures.

### 2. Version Control Database Changes

- Use migration scripts for all schema changes and store them in version control (e.g., Git).
- Review and test migration scripts before applying them to production.
- Implement rollback procedures for migrations.

### 3. Monitor Database Performance Continuously

- Set up continuous performance monitoring using Supabase dashboard or external tools.
- Establish performance baselines and alerts for critical metrics.
- Regularly review performance data and identify trends.

### 4. Implement Comprehensive Security Policies

- Enforce Row Level Security (RLS) on all tables to control data access.
- Regularly audit RLS policies and user permissions.
- Follow security best practices for database configuration and access management.

### 5. Test Backup and Recovery Procedures Regularly

- Test backup and restore procedures monthly to ensure they are effective.
- Document backup and recovery procedures and keep them up to date.
- Store backups in a secure and redundant location.

### 6. Document Database Schema and Maintenance Procedures

- Maintain comprehensive documentation of the database schema (as provided in `DATABASE_SCHEMA_DOCUMENTATION.md`).
- Document all maintenance procedures, troubleshooting steps, and best practices in this guide.
- Keep documentation up to date with any changes to the database or maintenance processes.

### 7. Use Staging Environment for Testing

- Always test database changes, migrations, and maintenance tasks in a staging environment before applying them to production.
- Use a staging environment that is as similar to production as possible.
- Verify changes thoroughly in staging before deploying to production.

### 8. Train Database Administrators and Developers

- Ensure database administrators and developers are trained on database maintenance procedures, security best practices, and troubleshooting techniques.
- Provide access to documentation and training resources.
- Conduct regular knowledge sharing sessions to improve team expertise.

By following this maintenance guide and implementing these best practices, you can ensure the CryoProtect v2 database remains a reliable and performant foundation for the application. Regular maintenance, proactive monitoring, and adherence to security policies are crucial for long-term database health and stability.