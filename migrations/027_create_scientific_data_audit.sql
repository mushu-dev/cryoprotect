-- Check if table exists
DO $$ 
BEGIN
    -- Create scientific_data_audit table if it doesn't exist
    IF NOT EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'scientific_data_audit') THEN
        CREATE TABLE IF NOT EXISTS public.scientific_data_audit (
            id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
            table_name TEXT NOT NULL,
            record_id TEXT NOT NULL,
            operation TEXT NOT NULL,
            old_value JSONB,
            new_value JSONB,
            user_id UUID,
            timestamp TIMESTAMP WITH TIME ZONE DEFAULT NOW()
        );
        
        -- Create indexes for efficient querying
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_table_name 
        ON public.scientific_data_audit(table_name);
        
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_record_id 
        ON public.scientific_data_audit(record_id);
        
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_timestamp 
        ON public.scientific_data_audit(timestamp);
        
        RAISE NOTICE 'Created scientific_data_audit table with indexes.';
    ELSE
        RAISE NOTICE 'scientific_data_audit table already exists.';
    END IF;
    
    -- Create RLS policies for scientific_data_audit if needed
    -- This can be extended as needed for your specific security requirements
END $$;
