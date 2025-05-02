-- CryoProtect Analyzer Database Schema Migration - Protocol Storage
-- Compatible with Supabase and PostgreSQL

-- Create protocols table
CREATE TABLE public.protocols (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    mixture_id UUID NOT NULL REFERENCES public.mixtures(id) ON DELETE CASCADE,
    name TEXT NOT NULL,
    description TEXT,
    target_concentration NUMERIC NOT NULL,
    sample_type TEXT NOT NULL,
    starting_temperature NUMERIC NOT NULL,
    target_temperature NUMERIC,
    step_count INTEGER NOT NULL,
    steps JSONB NOT NULL,
    custom_sensitivity JSONB,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id)
);

-- Create index for performance
CREATE INDEX idx_protocols_mixture_id ON public.protocols(mixture_id);

-- Create trigger for updated_at timestamp
CREATE TRIGGER set_timestamp_protocols
BEFORE UPDATE ON public.protocols
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

-- Enable Row Level Security (RLS)
ALTER TABLE public.protocols ENABLE ROW LEVEL SECURITY;

-- Create RLS Policies
-- Protocols: everyone can view, only authenticated users can insert, only creators can update/delete
CREATE POLICY "Protocols are viewable by everyone" 
ON public.protocols FOR SELECT USING (true);

CREATE POLICY "Protocols can be inserted by authenticated users" 
ON public.protocols FOR INSERT WITH CHECK (auth.role() = 'authenticated');

CREATE POLICY "Protocols can be updated by the creator" 
ON public.protocols FOR UPDATE USING (auth.uid() = created_by);

CREATE POLICY "Protocols can be deleted by the creator" 
ON public.protocols FOR DELETE USING (auth.uid() = created_by);

-- Create a view for protocols with mixture information
CREATE OR REPLACE VIEW public.protocol_with_mixture AS
SELECT 
    p.id,
    p.name,
    p.description,
    p.mixture_id,
    m.name AS mixture_name,
    p.target_concentration,
    p.sample_type,
    p.starting_temperature,
    p.target_temperature,
    p.step_count,
    p.steps,
    p.custom_sensitivity,
    p.created_at,
    p.updated_at,
    p.created_by
FROM 
    public.protocols p
JOIN 
    public.mixtures m ON p.mixture_id = m.id;

-- Add comment for documentation
COMMENT ON TABLE public.protocols IS 'Stores cryoprotectant protocols for stepwise addition/removal';