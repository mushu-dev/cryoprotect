-- Migration 017: Fix Foreign Key Relationships
-- Phase 2: Essential Maintenance Utilities

-- 1. Fix lab_verifications.verifier to reference auth.users(id)
-- (Assumes all verifier values are valid user UUIDs or will be migrated)
ALTER TABLE public.lab_verifications
    ALTER COLUMN verifier TYPE UUID USING verifier::UUID,
    ALTER COLUMN verifier SET NOT NULL;

ALTER TABLE public.lab_verifications
    ADD CONSTRAINT fk_lab_verifications_verifier
    FOREIGN KEY (verifier) REFERENCES auth.users(id) ON DELETE SET NULL;

-- 2. Add foreign key constraints to shares.item_id and shared_resources.resource_id
-- Note: These are polymorphic references; true FK constraints require table restructuring.
-- Instead, add NOT VALID constraints for the most common types, and document the limitation.

-- Example: If most shares are for molecules, add a partial constraint (PostgreSQL supports this with DEFERRABLE/NOT VALID)
-- Uncomment and adapt as needed:
-- ALTER TABLE public.shares
--     ADD CONSTRAINT fk_shares_item_id_molecules
--     FOREIGN KEY (item_id) REFERENCES public.molecules(id) NOT VALID;

-- For shared_resources, similar approach:
-- ALTER TABLE public.shared_resources
--     ADD CONSTRAINT fk_shared_resources_resource_id_project
--     FOREIGN KEY (resource_id) REFERENCES public.project(id) NOT VALID;

-- 3. Ensure ON DELETE CASCADE/SET NULL is set appropriately
-- (Already present in most tables; review and update as needed)

-- 4. Add comments for future maintainers
COMMENT ON COLUMN public.lab_verifications.verifier IS 'References auth.users(id); was previously VARCHAR, now UUID FK.';

-- 5. (Optional) Data migration for verifier if needed (manual step if existing data is not UUIDs)

-- End of Migration 017