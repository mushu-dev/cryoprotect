-- Migration: 031_implement_audit_logging.sql
-- Purpose: Implement comprehensive audit logging system
-- This migration adds audit schema, tables, and triggers to track changes

-- Step 1: Create audit schema if it doesn't exist
CREATE SCHEMA IF NOT EXISTS audit;

-- Step 2: Create the audit_log table to store all changes
CREATE TABLE IF NOT EXISTS audit.audit_log (
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

-- Add indexes to the audit_log table
CREATE INDEX IF NOT EXISTS idx_audit_log_table_name ON audit.audit_log(table_name);
CREATE INDEX IF NOT EXISTS idx_audit_log_user_id ON audit.audit_log(user_id);
CREATE INDEX IF NOT EXISTS idx_audit_log_action ON audit.audit_log(action);
CREATE INDEX IF NOT EXISTS idx_audit_log_recorded_at ON audit.audit_log(recorded_at);
CREATE INDEX IF NOT EXISTS idx_audit_log_row_id ON audit.audit_log(row_id);

-- Enable RLS on the audit_log table (only admins and service roles can access)
ALTER TABLE audit.audit_log ENABLE ROW LEVEL SECURITY;

DROP POLICY IF EXISTS audit_log_admin_policy ON audit.audit_log;
CREATE POLICY audit_log_admin_policy ON audit.audit_log
    FOR SELECT
    USING (
        auth.is_admin()
        OR auth.is_service_role_cached()
    );

-- Step 3: Create a function to extract request context information
CREATE OR REPLACE FUNCTION audit.get_request_context() RETURNS JSONB AS $$
DECLARE
    context JSONB;
BEGIN
    -- Initialize with empty context
    context := '{}'::JSONB;
    
    -- Try to get IP address if the current_setting exists
    BEGIN
        context := context || jsonb_build_object('ip_address', current_setting('request.headers')::JSONB->>'x-forwarded-for');
    EXCEPTION WHEN OTHERS THEN
        -- Ignore error if setting doesn't exist
    END;
    
    -- Try to get user agent if the current_setting exists
    BEGIN
        context := context || jsonb_build_object('user_agent', current_setting('request.headers')::JSONB->>'user-agent');
    EXCEPTION WHEN OTHERS THEN
        -- Ignore error if setting doesn't exist
    END;
    
    -- Get the jwt claims if available
    BEGIN
        context := context || jsonb_build_object('jwt_claims', current_setting('request.jwt.claims', true)::JSONB);
    EXCEPTION WHEN OTHERS THEN
        -- Ignore error if setting doesn't exist
    END;
    
    RETURN context;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Step 4: Create a function to find the difference between two records
CREATE OR REPLACE FUNCTION audit.changed_fields(old_data JSONB, new_data JSONB) RETURNS TEXT[] AS $$
DECLARE
    result TEXT[] := '{}'::TEXT[];
    key TEXT;
BEGIN
    -- For each key in new_data
    FOR key IN SELECT jsonb_object_keys(new_data)
    LOOP
        -- If the key doesn't exist in old_data or values are different
        IF (old_data->key) IS DISTINCT FROM (new_data->key) THEN
            result := array_append(result, key);
        END IF;
    END LOOP;
    
    -- For each key in old_data that might have been removed
    FOR key IN SELECT jsonb_object_keys(old_data)
    LOOP
        -- If the key doesn't exist in new_data and wasn't already added
        IF (new_data->key) IS NULL AND NOT (key = ANY(result)) THEN
            result := array_append(result, key);
        END IF;
    END LOOP;
    
    RETURN result;
END;
$$ LANGUAGE plpgsql IMMUTABLE;

-- Step 5: Create the core audit trigger function
CREATE OR REPLACE FUNCTION audit.audit_trigger_func() RETURNS TRIGGER AS $$
DECLARE
    excluded_cols TEXT[] := ARRAY[]::TEXT[];
    included_cols TEXT[] := ARRAY[]::TEXT[];
    log_row_data BOOLEAN := TRUE;
    log_changed_only BOOLEAN := TRUE;
    row_id TEXT;
    old_data JSONB := '{}'::JSONB;
    new_data JSONB := '{}'::JSONB;
    changed_fields TEXT[] := '{}'::TEXT[];
    app_context JSONB;
    log_action TEXT;
    row_to_json JSONB;
BEGIN
    -- Set the log action based on the operation
    log_action := TG_OP;
    
    -- Get the context information
    app_context := audit.get_request_context();

    -- Set the primary key for the record
    IF TG_ARGV[0] IS NOT NULL THEN
        EXECUTE 'SELECT ($1).' || TG_ARGV[0]::TEXT
        INTO row_id
        USING CASE WHEN TG_OP = 'DELETE' THEN OLD ELSE NEW END;
    ELSE
        row_id := NULL;
    END IF;
    
    -- Determine which columns to include/exclude
    IF TG_ARGV[1] IS NOT NULL THEN
        excluded_cols := TG_ARGV[1]::TEXT[];
    END IF;
    
    IF TG_ARGV[2] IS NOT NULL THEN
        included_cols := TG_ARGV[2]::TEXT[];
    END IF;
    
    -- Only log data if specified
    IF TG_ARGV[3] IS NOT NULL THEN
        log_row_data := TG_ARGV[3]::BOOLEAN;
    END IF;
    
    -- Only log changed fields if specified
    IF TG_ARGV[4] IS NOT NULL THEN
        log_changed_only := TG_ARGV[4]::BOOLEAN;
    END IF;
    
    -- For INSERT, only log new data
    IF (TG_OP = 'INSERT') THEN
        row_to_json := to_jsonb(NEW);
        
        IF log_row_data THEN
            -- Include only the desired columns
            IF array_length(included_cols, 1) > 0 THEN
                SELECT jsonb_object_agg(key, value)
                INTO new_data
                FROM jsonb_each(row_to_json)
                WHERE key = ANY(included_cols);
            -- Exclude specified columns
            ELSIF array_length(excluded_cols, 1) > 0 THEN
                SELECT jsonb_object_agg(key, value)
                INTO new_data
                FROM jsonb_each(row_to_json)
                WHERE key <> ALL(excluded_cols);
            -- Include all columns
            ELSE
                new_data := row_to_json;
            END IF;
        END IF;
        
        changed_fields := array_agg(key) FROM jsonb_object_keys(new_data) key;
        
    -- For DELETE, only log old data
    ELSIF (TG_OP = 'DELETE') THEN
        row_to_json := to_jsonb(OLD);
        
        IF log_row_data THEN
            -- Include only the desired columns
            IF array_length(included_cols, 1) > 0 THEN
                SELECT jsonb_object_agg(key, value)
                INTO old_data
                FROM jsonb_each(row_to_json)
                WHERE key = ANY(included_cols);
            -- Exclude specified columns
            ELSIF array_length(excluded_cols, 1) > 0 THEN
                SELECT jsonb_object_agg(key, value)
                INTO old_data
                FROM jsonb_each(row_to_json)
                WHERE key <> ALL(excluded_cols);
            -- Include all columns
            ELSE
                old_data := row_to_json;
            END IF;
        END IF;
        
        changed_fields := array_agg(key) FROM jsonb_object_keys(old_data) key;
        
    -- For UPDATE, log both old and new data
    ELSIF (TG_OP = 'UPDATE') THEN
        -- Convert both old and new records to JSON
        row_to_json := to_jsonb(OLD);
        
        IF log_row_data THEN
            -- Include only the desired columns for the old data
            IF array_length(included_cols, 1) > 0 THEN
                SELECT jsonb_object_agg(key, value)
                INTO old_data
                FROM jsonb_each(row_to_json)
                WHERE key = ANY(included_cols);
            -- Exclude specified columns for the old data
            ELSIF array_length(excluded_cols, 1) > 0 THEN
                SELECT jsonb_object_agg(key, value)
                INTO old_data
                FROM jsonb_each(row_to_json)
                WHERE key <> ALL(excluded_cols);
            -- Include all columns for the old data
            ELSE
                old_data := row_to_json;
            END IF;
            
            -- Convert new record to JSON
            row_to_json := to_jsonb(NEW);
            
            -- Include only the desired columns for the new data
            IF array_length(included_cols, 1) > 0 THEN
                SELECT jsonb_object_agg(key, value)
                INTO new_data
                FROM jsonb_each(row_to_json)
                WHERE key = ANY(included_cols);
            -- Exclude specified columns for the new data
            ELSIF array_length(excluded_cols, 1) > 0 THEN
                SELECT jsonb_object_agg(key, value)
                INTO new_data
                FROM jsonb_each(row_to_json)
                WHERE key <> ALL(excluded_cols);
            -- Include all columns for the new data
            ELSE
                new_data := row_to_json;
            END IF;
        END IF;
        
        -- Calculate the changed fields
        changed_fields := audit.changed_fields(old_data, new_data);
        
        -- If no fields were changed, don't log this update
        IF log_changed_only AND array_length(changed_fields, 1) IS NULL THEN
            RETURN NULL;
        END IF;
    END IF;
    
    -- Insert the audit log entry
    INSERT INTO audit.audit_log (
        table_schema,
        table_name,
        user_id,
        action,
        row_id,
        old_data,
        new_data,
        changed_fields,
        ip_address,
        user_agent,
        app_context
    ) VALUES (
        TG_TABLE_SCHEMA,
        TG_TABLE_NAME,
        (SELECT id FROM auth.users WHERE id = auth.uid()),
        log_action,
        row_id,
        old_data,
        new_data,
        changed_fields,
        (app_context->>'ip_address')::TEXT,
        (app_context->>'user_agent')::TEXT,
        app_context
    );
    
    RETURN NULL;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Step 6: Create a helper function to easily add audit triggers to tables
CREATE OR REPLACE FUNCTION audit.enable_table_auditing(
    target_table TEXT,
    target_schema TEXT DEFAULT 'public',
    pk_column TEXT DEFAULT 'id',
    excluded_cols TEXT[] DEFAULT NULL,
    included_cols TEXT[] DEFAULT NULL,
    log_row_data BOOLEAN DEFAULT TRUE,
    log_changed_only BOOLEAN DEFAULT TRUE
) RETURNS VOID AS $$
DECLARE
    trigger_name TEXT;
    fully_qualified_table TEXT;
BEGIN
    fully_qualified_table := format('%I.%I', target_schema, target_table);
    trigger_name := format('audit_trigger_%s', target_table);
    
    EXECUTE format('DROP TRIGGER IF EXISTS %I ON %s', trigger_name, fully_qualified_table);
    
    EXECUTE format('
        CREATE TRIGGER %I
        AFTER INSERT OR UPDATE OR DELETE ON %s
        FOR EACH ROW EXECUTE FUNCTION audit.audit_trigger_func(%L, %L, %L, %L, %L);
    ', 
    trigger_name, 
    fully_qualified_table, 
    pk_column, 
    excluded_cols, 
    included_cols, 
    log_row_data, 
    log_changed_only);
    
    RAISE NOTICE 'Audit trigger created for %', fully_qualified_table;
END;
$$ LANGUAGE plpgsql;

-- Step 7: Create a function to enable auditing for all tables in a schema
CREATE OR REPLACE FUNCTION audit.enable_schema_auditing(
    target_schema TEXT DEFAULT 'public',
    excluded_tables TEXT[] DEFAULT NULL,
    excluded_cols TEXT[] DEFAULT ARRAY['created_at', 'updated_at']
) RETURNS VOID AS $$
DECLARE
    tables RECORD;
BEGIN
    FOR tables IN 
        SELECT tablename 
        FROM pg_tables 
        WHERE schemaname = target_schema
          AND tablename <> ALL(COALESCE(excluded_tables, ARRAY[]::TEXT[]))
    LOOP
        BEGIN
            PERFORM audit.enable_table_auditing(
                target_table := tables.tablename,
                target_schema := target_schema,
                excluded_cols := excluded_cols
            );
        EXCEPTION WHEN OTHERS THEN
            RAISE WARNING 'Failed to enable auditing for table %: %', tables.tablename, SQLERRM;
        END;
    END LOOP;
END;
$$ LANGUAGE plpgsql;

-- Step 8: Create monitoring views for audit data
CREATE OR REPLACE VIEW audit.recent_activity AS
SELECT
    a.id,
    a.table_schema,
    a.table_name,
    a.action,
    a.row_id,
    a.changed_fields,
    a.recorded_at,
    u.email as user_email,
    a.ip_address
FROM
    audit.audit_log a
LEFT JOIN
    auth.users u ON a.user_id = u.id
ORDER BY
    a.recorded_at DESC
LIMIT 1000;

CREATE OR REPLACE VIEW audit.activity_summary AS
SELECT
    table_schema,
    table_name,
    action,
    COUNT(*) as action_count,
    MIN(recorded_at) as first_action,
    MAX(recorded_at) as last_action
FROM
    audit.audit_log
GROUP BY
    table_schema,
    table_name,
    action
ORDER BY
    table_schema,
    table_name,
    action;

CREATE OR REPLACE VIEW audit.user_activity AS
SELECT
    u.email as user_email,
    a.action,
    COUNT(*) as action_count,
    MIN(a.recorded_at) as first_action,
    MAX(a.recorded_at) as last_action
FROM
    audit.audit_log a
JOIN
    auth.users u ON a.user_id = u.id
GROUP BY
    u.email,
    a.action
ORDER BY
    u.email,
    a.action;

-- Step 9: Create a function to retrieve audit history for a specific record
CREATE OR REPLACE FUNCTION audit.get_record_history(
    target_table TEXT,
    target_schema TEXT DEFAULT 'public',
    record_id TEXT
) RETURNS TABLE (
    audit_id BIGINT,
    action TEXT,
    changed_fields TEXT[],
    old_data JSONB,
    new_data JSONB,
    user_email TEXT,
    recorded_at TIMESTAMPTZ
) AS $$
BEGIN
    RETURN QUERY
    SELECT
        a.id as audit_id,
        a.action,
        a.changed_fields,
        a.old_data,
        a.new_data,
        u.email as user_email,
        a.recorded_at
    FROM
        audit.audit_log a
    LEFT JOIN
        auth.users u ON a.user_id = u.id
    WHERE
        a.table_schema = target_schema
        AND a.table_name = target_table
        AND a.row_id = record_id
    ORDER BY
        a.recorded_at DESC;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Step 10: Enable auditing for core tables

-- Enable auditing for molecules table
SELECT audit.enable_table_auditing('molecules', excluded_cols := ARRAY['updated_at']);

-- Enable auditing for mixtures table
SELECT audit.enable_table_auditing('mixtures', excluded_cols := ARRAY['updated_at']);

-- Enable auditing for molecular_properties table
SELECT audit.enable_table_auditing('molecular_properties', excluded_cols := ARRAY['updated_at']);

-- Enable auditing for teams table
SELECT audit.enable_table_auditing('teams', excluded_cols := ARRAY['updated_at']);

-- Enable auditing for team_members table
SELECT audit.enable_table_auditing('team_members', pk_column := 'team_id');

-- Enable auditing for mixture_components table
SELECT audit.enable_table_auditing('mixture_components', pk_column := 'mixture_id');

-- Step 11: Create a specialized logging table for application events
CREATE TABLE IF NOT EXISTS audit.application_log (
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

-- Add indexes to the application_log table
CREATE INDEX IF NOT EXISTS idx_app_log_event_type ON audit.application_log(event_type);
CREATE INDEX IF NOT EXISTS idx_app_log_event_severity ON audit.application_log(event_severity);
CREATE INDEX IF NOT EXISTS idx_app_log_user_id ON audit.application_log(user_id);
CREATE INDEX IF NOT EXISTS idx_app_log_recorded_at ON audit.application_log(recorded_at);

-- Enable RLS on the application_log table
ALTER TABLE audit.application_log ENABLE ROW LEVEL SECURITY;

DROP POLICY IF EXISTS application_log_admin_policy ON audit.application_log;
CREATE POLICY application_log_admin_policy ON audit.application_log
    FOR SELECT
    USING (
        auth.is_admin()
        OR auth.is_service_role_cached()
    );

-- Create a function to log application events
CREATE OR REPLACE FUNCTION audit.log_application_event(
    event_type TEXT,
    event_severity TEXT,
    message TEXT,
    event_source TEXT DEFAULT NULL,
    details JSONB DEFAULT NULL
) RETURNS BIGINT AS $$
DECLARE
    new_log_id BIGINT;
    context JSONB;
BEGIN
    -- Get request context
    context := audit.get_request_context();
    
    -- Insert the log entry
    INSERT INTO audit.application_log (
        event_type,
        event_severity,
        event_source,
        user_id,
        message,
        details,
        ip_address,
        user_agent
    ) VALUES (
        event_type,
        event_severity,
        event_source,
        (SELECT id FROM auth.users WHERE id = auth.uid()),
        message,
        details,
        (context->>'ip_address')::TEXT,
        (context->>'user_agent')::TEXT
    )
    RETURNING id INTO new_log_id;
    
    RETURN new_log_id;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Step 12: Create monitoring views for application logs
CREATE OR REPLACE VIEW audit.recent_app_logs AS
SELECT
    a.id,
    a.event_type,
    a.event_severity,
    a.event_source,
    a.message,
    a.details,
    a.recorded_at,
    u.email as user_email,
    a.ip_address
FROM
    audit.application_log a
LEFT JOIN
    auth.users u ON a.user_id = u.id
ORDER BY
    a.recorded_at DESC
LIMIT 1000;

CREATE OR REPLACE VIEW audit.app_log_summary AS
SELECT
    event_type,
    event_severity,
    COUNT(*) as event_count,
    MIN(recorded_at) as first_event,
    MAX(recorded_at) as last_event
FROM
    audit.application_log
GROUP BY
    event_type,
    event_severity
ORDER BY
    event_severity,
    event_type;

-- Step 13: Create a clean-up function for old audit logs
CREATE OR REPLACE FUNCTION audit.cleanup_old_logs(
    retention_days INTEGER DEFAULT 90
) RETURNS INTEGER AS $$
DECLARE
    audit_deleted INTEGER;
    app_log_deleted INTEGER;
    cutoff_date TIMESTAMP WITH TIME ZONE;
BEGIN
    cutoff_date := NOW() - (retention_days || ' days')::INTERVAL;
    
    -- Delete old audit logs
    DELETE FROM audit.audit_log
    WHERE recorded_at < cutoff_date
    RETURNING COUNT(*) INTO audit_deleted;
    
    -- Delete old application logs
    DELETE FROM audit.application_log
    WHERE recorded_at < cutoff_date
    RETURNING COUNT(*) INTO app_log_deleted;
    
    -- Log the cleanup
    PERFORM audit.log_application_event(
        'LOG_CLEANUP',
        'INFO',
        format('Cleaned up %s audit logs and %s application logs older than %s days',
               audit_deleted, app_log_deleted, retention_days),
        'SYSTEM',
        jsonb_build_object(
            'audit_deleted', audit_deleted,
            'app_log_deleted', app_log_deleted,
            'retention_days', retention_days,
            'cutoff_date', cutoff_date
        )
    );
    
    RETURN audit_deleted + app_log_deleted;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;