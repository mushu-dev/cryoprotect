-- Migration: Create predictions and experiments tables for CryoProtect v2

-- Table: predictions
CREATE TABLE IF NOT EXISTS public.predictions (
    id SERIAL PRIMARY KEY,
    mixture_id INTEGER NOT NULL REFERENCES public.mixture(id) ON DELETE CASCADE,
    property_type_id INTEGER NOT NULL REFERENCES public.property_type(id) ON DELETE CASCADE,
    value DOUBLE PRECISION NOT NULL,
    units VARCHAR(32) NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- Table: experiments
CREATE TABLE IF NOT EXISTS public.experiments (
    id SERIAL PRIMARY KEY,
    mixture_id INTEGER NOT NULL REFERENCES public.mixture(id) ON DELETE CASCADE,
    experiment_date DATE,
    result JSONB,
    notes TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);