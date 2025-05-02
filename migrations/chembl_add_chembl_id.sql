-- Add chembl_id column to molecules table
ALTER TABLE public.molecules ADD COLUMN IF NOT EXISTS chembl_id VARCHAR(20);

-- Backfill chembl_id from data_source
UPDATE public.molecules
SET chembl_id = SUBSTRING(data_source, 12)
WHERE data_source LIKE 'ChEMBL ID: %' AND chembl_id IS NULL;

-- Create index for fast lookup
CREATE INDEX IF NOT EXISTS idx_molecules_chembl_id ON public.molecules(chembl_id);