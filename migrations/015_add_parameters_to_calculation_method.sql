-- Migration to add 'parameters' column to calculation_method table for toxicity scoring compatibility

ALTER TABLE calculation_method
ADD COLUMN IF NOT EXISTS parameters JSONB;