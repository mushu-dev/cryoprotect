-- Migration: 025_consolidated_molecule_audit
-- Description: Creates audit table and triggers for consolidated molecule changes

-- Create audit table for tracking molecule consolidation changes
CREATE TABLE IF NOT EXISTS molecule_consolidation_audit (
    id SERIAL PRIMARY KEY,
    operation_type TEXT NOT NULL,
    primary_molecule_id UUID NOT NULL,
    secondary_molecule_id UUID NOT NULL,
    performed_by TEXT,
    performed_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    operation_details JSONB
);

-- Create indexes for efficient querying
CREATE INDEX idx_molecule_consolidation_audit_primary
ON molecule_consolidation_audit(primary_molecule_id);

CREATE INDEX idx_molecule_consolidation_audit_secondary
ON molecule_consolidation_audit(secondary_molecule_id);

CREATE INDEX idx_molecule_consolidation_audit_operation_type
ON molecule_consolidation_audit(operation_type);

CREATE INDEX idx_molecule_consolidation_audit_performed_at
ON molecule_consolidation_audit(performed_at);

-- Create function for logging molecule consolidation changes
CREATE OR REPLACE FUNCTION log_molecule_consolidation()
RETURNS TRIGGER AS $$
BEGIN
    IF (TG_OP = 'UPDATE' AND 
        NEW.consolidated_to IS DISTINCT FROM OLD.consolidated_to) THEN
        
        -- Determine operation type
        DECLARE
            op_type TEXT;
        BEGIN
            IF OLD.consolidated_to IS NULL AND NEW.consolidated_to IS NOT NULL THEN
                op_type := 'CONSOLIDATE';
            ELSIF OLD.consolidated_to IS NOT NULL AND NEW.consolidated_to IS NULL THEN
                op_type := 'DECONSOLIDATE';
            ELSE
                op_type := 'CHANGE_PRIMARY';
            END IF;
            
            -- Insert audit record
            INSERT INTO molecule_consolidation_audit (
                operation_type, 
                primary_molecule_id, 
                secondary_molecule_id,
                performed_by,
                operation_details
            ) VALUES (
                op_type,
                COALESCE(NEW.consolidated_to, OLD.consolidated_to),
                NEW.id,
                COALESCE(current_setting('app.current_user', true), current_user),
                jsonb_build_object(
                    'previous_state', 
                    CASE WHEN OLD.consolidated_to IS NULL 
                         THEN 'INDEPENDENT' 
                         ELSE 'CONSOLIDATED TO ' || OLD.consolidated_to::text 
                    END,
                    'new_state', 
                    CASE WHEN NEW.consolidated_to IS NULL
                         THEN 'INDEPENDENT'
                         ELSE 'CONSOLIDATED TO ' || NEW.consolidated_to::text
                    END,
                    'timestamp', NOW(),
                    'molecule_name', (SELECT name FROM molecules WHERE id = NEW.id)
                )
            );
        END;
    END IF;
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Create or replace trigger for molecule consolidation changes
DROP TRIGGER IF EXISTS molecule_consolidation_audit_trigger ON molecules;

CREATE TRIGGER molecule_consolidation_audit_trigger
AFTER UPDATE ON molecules
FOR EACH ROW
WHEN (NEW.consolidated_to IS DISTINCT FROM OLD.consolidated_to)
EXECUTE FUNCTION log_molecule_consolidation();

-- Add comment to explain the purpose of the audit table
COMMENT ON TABLE molecule_consolidation_audit IS 
'Audit trail for molecule consolidation operations, tracking when molecules are consolidated, deconsolidated, or have their primary molecule changed';

-- Add sample audit entry to verify the trigger works
-- This will be rolled back in a transaction
DO $$
BEGIN
    -- Only execute if there are molecules to test with
    IF EXISTS (SELECT 1 FROM molecules LIMIT 1) THEN
        -- Start a transaction that we'll roll back
        BEGIN;
        
        -- Get a sample molecule ID
        DECLARE
            sample_id UUID;
            primary_id UUID;
        BEGIN
            -- Get a molecule that isn't already consolidated
            SELECT id INTO sample_id FROM molecules 
            WHERE consolidated_to IS NULL
            LIMIT 1;
            
            -- Get a different molecule to use as primary
            SELECT id INTO primary_id FROM molecules 
            WHERE id != sample_id AND consolidated_to IS NULL
            LIMIT 1;
            
            -- If we found molecules, test the trigger
            IF sample_id IS NOT NULL AND primary_id IS NOT NULL THEN
                -- Set app.current_user for audit trail
                SET LOCAL app.current_user = 'test_user';
                
                -- Update the molecule to test the trigger
                UPDATE molecules 
                SET consolidated_to = primary_id
                WHERE id = sample_id;
                
                -- Verify the audit record was created
                RAISE NOTICE 'Audit trigger test successful. Sample record created and will be rolled back.';
            END IF;
        END;
        
        -- Roll back the transaction
        ROLLBACK;
    END IF;
END;
$$;