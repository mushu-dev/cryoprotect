// Common Types

export type ApiResponse<T> = {
  status: 'success' | 'error';
  data?: T;
  error?: {
    code: string;
    message: string;
    details?: any;
  };
  pagination?: {
    limit: number;
    offset: number;
    total: number;
  };
};

// Molecule Types

export type MoleculeStatus = 'original' | 'primary' | 'duplicate';

export type Molecule = {
  id: string;
  name: string;
  inchikey: string;
  smiles: string;
  molecular_formula: string;
  molecular_weight: number;
  is_consolidated: boolean;
  molecule_status: MoleculeStatus;
  primary_molecule_id?: string;
  primary_molecule_name?: string;
  created_at: string;
  updated_at: string;
  properties?: MoleculeProperty[];
  duplicate_molecules?: { id: string; name: string }[];
  duplicate_count?: number;
};

export type MoleculeProperty = {
  id: string;
  molecule_id: string;
  property_type_id: string;
  property_name: string;
  property_type: string;
  numeric_value?: number;
  text_value?: string;
  unit?: string;
  source: string;
  created_at: string;
  updated_at: string;
};

// Mixture Types

export type MixtureComponent = {
  id: string;
  mixture_id: string;
  molecule_id: string;
  molecule_name: string;
  concentration: number;
  concentration_unit: string;
  role?: string;
};

export type Mixture = {
  id: string;
  name: string;
  description: string;
  created_at: string;
  updated_at: string;
  components: MixtureComponent[];
  properties?: MoleculeProperty[];
  component_count?: number;
};

// User Types

export type User = {
  id: string;
  email: string;
  role: string;
  roles?: string[];
  permissions?: string[];
};

// Consolidated Molecule Types

export type ConsolidatedMolecule = Molecule & {
  primary_note?: string;
  audit_history?: AuditRecord[];
};

export type AuditRecord = {
  id: string;
  table_name: string;
  record_id: string;
  operation: string;
  old_value: string;
  new_value: string;
  user_id: string;
  timestamp: string;
};

// Property Migration

export type PropertyMigration = {
  source_molecule_id: string;
  target_molecule_id: string;
  property_ids?: string[];
  migrated_properties?: {
    source_property_id: string;
    target_property_id: string;
    property_name: string;
  }[];
  skipped_properties?: string[];
  migrated_count: number;
  skipped_count: number;
};

// Dashboard Types

export type DashboardStats = {
  molecules: {
    total: number;
    consolidated: number;
    unique: number;
  };
  mixtures: {
    total: number;
    public: number;
    private: number;
  };
  experiments: {
    total: number;
    successful: number;
    failed: number;
  };
  predictions: {
    total: number;
    high_confidence: number;
    medium_confidence: number;
    low_confidence: number;
  };
};

export type Activity = {
  id: string;
  action: string;
  entity_type: string;
  entity_id: string;
  entity_name: string;
  timestamp: string;
  user_id: string;
  user_email: string;
};

// RDKit Types

export type RDKitProperties = {
  molecular_weight: number;
  logp: number;
  num_atoms: number;
  num_bonds: number;
  num_rings: number;
  tpsa: number;
  h_bond_donors: number;
  h_bond_acceptors: number;
  rotatable_bonds: number;
};

export type VisualizationFormat = 'svg' | 'png' | 'base64';

export type MoleculeVisualization = {
  visualization: string;
  format: VisualizationFormat;
};

export type StructureMatch = {
  id: string;
  name: string;
  smiles: string;
  similarity: number;
};