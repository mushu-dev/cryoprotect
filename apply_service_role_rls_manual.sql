-- CryoProtect Analyzer Service Role RLS Policies Migration
-- Enables service role to insert data into tables with RLS enabled
-- This is necessary for scripts like populate_molecules.py to work properly

-- =========================
-- 1. molecule table
-- =========================
CREATE POLICY "Allow service role inserts on molecule" 
  ON public.molecule 
  FOR INSERT 
  WITH CHECK (auth.role() = 'service_role');
COMMENT ON POLICY "Allow service role inserts on molecule" ON molecule IS
  'Allows service role to insert molecules regardless of project membership.';

-- =========================
-- 2. mixture table
-- =========================
CREATE POLICY "Allow service role inserts on mixture" 
  ON public.mixture 
  FOR INSERT 
  WITH CHECK (auth.role() = 'service_role');
COMMENT ON POLICY "Allow service role inserts on mixture" ON mixture IS
  'Allows service role to insert mixtures regardless of project membership.';

-- =========================
-- 3. mixture_component table
-- =========================
CREATE POLICY "Allow service role inserts on mixture_component" 
  ON public.mixture_component 
  FOR INSERT 
  WITH CHECK (auth.role() = 'service_role');
COMMENT ON POLICY "Allow service role inserts on mixture_component" ON mixture_component IS
  'Allows service role to insert mixture components regardless of project membership.';

-- =========================
-- 4. experiment table
-- =========================
CREATE POLICY "Allow service role inserts on experiment" 
  ON public.experiment 
  FOR INSERT 
  WITH CHECK (auth.role() = 'service_role');
COMMENT ON POLICY "Allow service role inserts on experiment" ON experiment IS
  'Allows service role to insert experiments regardless of project membership.';

-- =========================
-- 5. molecular_property table
-- =========================
CREATE POLICY "Allow service role inserts on molecular_property" 
  ON public.molecular_property 
  FOR INSERT 
  WITH CHECK (auth.role() = 'service_role');
COMMENT ON POLICY "Allow service role inserts on molecular_property" ON molecular_property IS
  'Allows service role to insert molecular properties regardless of project membership.';

-- =========================
-- 6. prediction table
-- =========================
CREATE POLICY "Allow service role inserts on prediction" 
  ON public.prediction 
  FOR INSERT 
  WITH CHECK (auth.role() = 'service_role');
COMMENT ON POLICY "Allow service role inserts on prediction" ON prediction IS
  'Allows service role to insert predictions regardless of project membership.';

-- =========================
-- 7. experiment_property table
-- =========================
CREATE POLICY "Allow service role inserts on experiment_property" 
  ON public.experiment_property 
  FOR INSERT 
  WITH CHECK (auth.role() = 'service_role');
COMMENT ON POLICY "Allow service role inserts on experiment_property" ON experiment_property IS
  'Allows service role to insert experiment properties regardless of project membership.';