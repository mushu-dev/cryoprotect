-- Create lab_verifications table
CREATE TABLE IF NOT EXISTS public.lab_verifications (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    experiment_id UUID NOT NULL REFERENCES public.experiments(id) ON DELETE CASCADE,
    verification_status VARCHAR(20) NOT NULL CHECK (verification_status IN ('pending', 'verified', 'rejected')),
    verifier VARCHAR(255) NOT NULL,
    equipment_used VARCHAR(255) NOT NULL,
    comments TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- Add RLS policies
ALTER TABLE public.lab_verifications ENABLE ROW LEVEL SECURITY;

CREATE POLICY lab_verifications_select ON public.lab_verifications
    FOR SELECT USING (auth.uid() = verifier OR auth.uid() IN (
        SELECT user_id FROM user_team_membership utm
        JOIN team_roles tr ON utm.role_id = tr.id
        WHERE tr.can_view_verifications = true
    ));

CREATE POLICY lab_verifications_insert ON public.lab_verifications
    FOR INSERT WITH CHECK (auth.uid() = verifier);

CREATE POLICY lab_verifications_update ON public.lab_verifications
    FOR UPDATE USING (auth.uid() = verifier OR auth.uid() IN (
        SELECT user_id FROM user_team_membership utm
        JOIN team_roles tr ON utm.role_id = tr.id
        WHERE tr.can_update_verifications = true
    ));