# Audit Logging System

This document describes the comprehensive audit logging system implemented in the CryoProtect database.

## Overview

The audit logging system provides:

1. Automatic tracking of all data changes (INSERT, UPDATE, DELETE)
2. Application event logging for custom events and operations
3. User activity monitoring and analytics
4. Complete data change history for regulatory compliance
5. Performance-optimized design with proper indexing

## Database Structure

### Schemas

The system uses a dedicated `audit` schema to separate audit data from the application data.

### Core Tables

#### `audit.audit_log`

Tracks all data changes:

```sql
CREATE TABLE audit.audit_log (
    id BIGSERIAL PRIMARY KEY,
    table_name TEXT NOT NULL,
    table_schema TEXT NOT NULL DEFAULT 'public',
    user_id UUID,
    action TEXT NOT NULL CHECK (action IN ('INSERT', 'UPDATE', 'DELETE', 'TRUNCATE')),
    row_id TEXT,
    old_data JSONB,
    new_data JSONB,
    changed_fields TEXT[],
    ip_address TEXT,
    user_agent TEXT,
    app_context JSONB,
    recorded_at TIMESTAMP WITH TIME ZONE NOT NULL DEFAULT NOW()
);
```

#### `audit.application_log`

Tracks application events and operations:

```sql
CREATE TABLE audit.application_log (
    id BIGSERIAL PRIMARY KEY,
    event_type TEXT NOT NULL,
    event_severity TEXT NOT NULL CHECK (event_severity IN ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL')),
    event_source TEXT,
    user_id UUID,
    message TEXT NOT NULL,
    details JSONB,
    ip_address TEXT,
    user_agent TEXT,
    recorded_at TIMESTAMP WITH TIME ZONE NOT NULL DEFAULT NOW()
);
```

## Key Components

### Trigger-Based Auditing

The system uses PostgreSQL triggers to automatically capture changes:

1. Each audited table has a trigger attached
2. The trigger fires on INSERT, UPDATE, and DELETE actions
3. For each change, a record is stored in the `audit_log` table

### Helper Functions

#### `audit.audit_trigger_func()`

Core trigger function that captures change data.

#### `audit.enable_table_auditing()`

```sql
SELECT audit.enable_table_auditing(
    target_table := 'table_name',
    target_schema := 'public',
    pk_column := 'id',
    excluded_cols := ARRAY['updated_at']
);
```

Enables auditing for a specific table with these parameters:
- `target_table`: The table to audit
- `target_schema`: The schema containing the table (default: 'public')
- `pk_column`: The primary key column used for row identification
- `excluded_cols`: Columns to exclude from auditing (e.g., timestamps)
- `included_cols`: Only audit these columns (if specified)
- `log_row_data`: Whether to log the actual row data (default: true)
- `log_changed_only`: Whether to log only if fields changed (default: true)

#### `audit.enable_schema_auditing()`

```sql
SELECT audit.enable_schema_auditing(
    target_schema := 'public',
    excluded_tables := ARRAY['large_table', 'temp_table'],
    excluded_cols := ARRAY['created_at', 'updated_at']
);
```

Enables auditing for all tables in a schema.

#### `audit.get_record_history()`

```sql
SELECT * FROM audit.get_record_history(
    target_table := 'molecules',
    record_id := '123e4567-e89b-12d3-a456-426614174000'
);
```

Retrieves the complete audit history for a specific record.

#### `audit.log_application_event()`

```sql
SELECT audit.log_application_event(
    event_type := 'USER_LOGIN',
    event_severity := 'INFO',
    message := 'User logged in successfully',
    event_source := 'AUTH_SYSTEM',
    details := '{"ip": "192.168.1.1", "method": "password"}'::jsonb
);
```

Logs an application event.

#### `audit.cleanup_old_logs()`

```sql
SELECT audit.cleanup_old_logs(retention_days := 90);
```

Cleans up old audit logs based on a retention policy.

### Monitoring Views

#### `audit.recent_activity`

Shows recent audit logs with user information.

#### `audit.activity_summary`

Summarizes audit activity by table and action type.

#### `audit.user_activity`

Summarizes user activity by action type.

#### `audit.recent_app_logs`

Shows recent application logs with user information.

#### `audit.app_log_summary`

Summarizes application logs by event type and severity.

## Security

The audit system is secured with Row Level Security (RLS):

- Only admins and service roles can access audit logs
- All functions use `SECURITY DEFINER` to ensure consistent privilege evaluation
- The audit schema is separate from application data

## Usage Examples

### Viewing Audit History for a Record

```sql
-- Get the complete history for a molecule
SELECT * FROM audit.get_record_history('molecules', record_id := 'your-molecule-id');
```

### Logging Application Events

```sql
-- Log a user action
SELECT audit.log_application_event(
    'EXPORT_DATA',
    'INFO',
    'User exported molecule data',
    'EXPORT_SYSTEM',
    jsonb_build_object(
        'format', 'CSV',
        'record_count', 42,
        'filters', 'molecular_weight > 100'
    )
);

-- Log an error
SELECT audit.log_application_event(
    'API_ERROR',
    'ERROR',
    'Failed to process API request',
    'API_GATEWAY',
    jsonb_build_object(
        'endpoint', '/api/molecules',
        'method', 'POST',
        'error_code', 500,
        'error_message', 'Database connection failed'
    )
);
```

### Finding Changed Records

```sql
-- Find all molecules updated in the last 24 hours
SELECT 
    a.table_name,
    a.row_id,
    a.changed_fields,
    a.recorded_at,
    u.email as modified_by
FROM 
    audit.audit_log a
JOIN
    auth.users u ON a.user_id = u.id
WHERE 
    a.table_name = 'molecules'
    AND a.action = 'UPDATE'
    AND a.recorded_at > NOW() - INTERVAL '24 hours'
ORDER BY
    a.recorded_at DESC;
```

### Getting User Activity

```sql
-- Find what a specific user has been doing
SELECT 
    a.table_name,
    a.action,
    a.row_id,
    a.changed_fields,
    a.recorded_at
FROM 
    audit.audit_log a
JOIN
    auth.users u ON a.user_id = u.id
WHERE 
    u.email = 'user@example.com'
ORDER BY
    a.recorded_at DESC
LIMIT 100;
```

## Performance Considerations

1. **Indexes**: The audit tables are properly indexed for performance
2. **Selective Auditing**: You can exclude unnecessary columns
3. **Log Cleanup**: Regular cleanup of old logs is recommended
4. **Changed Fields Only**: By default, only records actions that actually change data

## Best Practices

1. **Sensitive Data**: Consider excluding sensitive fields from audit logs
2. **Retention Policy**: Implement a retention policy using `audit.cleanup_old_logs()`
3. **Application Logging**: Use `audit.log_application_event()` for important application events
4. **Regular Reviews**: Periodically review audit logs for security concerns

## Adding Auditing to New Tables

When creating new tables that should be audited:

```sql
-- Create the table
CREATE TABLE new_table (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    name TEXT NOT NULL,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

-- Enable auditing for the table
SELECT audit.enable_table_auditing(
    'new_table',
    excluded_cols := ARRAY['updated_at']
);
```

## Future Enhancements

Potential future enhancements include:

1. Admin UI for audit log review
2. Scheduled audit reports
3. Anomaly detection based on audit patterns
4. Enhanced visualization of record history
5. Integration with external SIEM systems