# Experimental Data Schema

This document outlines the unified schema for the experimental data system.

## Overview

The experimental data schema is designed to support:

- Comprehensive experiment metadata
- Protocol versioning and templating
- Detailed result tracking with uncertainty
- Proper relationships with molecules and mixtures
- Scientific method documentation and reproducibility

## Tables

### Core Tables

#### experiments

The central table for experimental data.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| name | TEXT | Name of the experiment |
| description | TEXT | Description of the experiment |
| experiment_type_id | UUID | Reference to experiment_types |
| tissue_type_id | UUID | Reference to tissue_types |
| protocol_id | UUID | Reference to protocol used |
| preparation_protocol | TEXT | Textual description of protocol |
| temperature | FLOAT | Temperature in Celsius |
| temperature_unit | TEXT | Temperature unit (C, K, F) |
| pressure | FLOAT | Pressure value |
| pressure_unit | TEXT | Pressure unit (kPa, atm) |
| created_by | UUID | User who created the experiment |
| created_at | TIMESTAMP | Creation timestamp |
| updated_at | TIMESTAMP | Last update timestamp |
| data_source | TEXT | Source of the experimental data |
| version | INTEGER | Version number |
| modification_history | JSONB | Detailed modification history |
| uncertainty_model | TEXT | Type of uncertainty model used |
| data_provenance | JSONB | Data provenance information |
| quality_score | FLOAT | Quality score for the data |

#### experiment_types

Types of cryopreservation experiments.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| name | TEXT | Name of the experiment type |
| description | TEXT | Description of the type |
| protocol_details | JSONB | Default protocol parameters |
| created_by | UUID | User who created the type |
| created_at | TIMESTAMP | Creation timestamp |
| updated_at | TIMESTAMP | Last update timestamp |

#### tissue_types

Biological samples used in experiments.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| name | TEXT | Name of the tissue type |
| description | TEXT | Description of the tissue |
| species | TEXT | Species of origin |
| taxonomy_id | INTEGER | Taxonomy identifier |
| created_by | UUID | User who created the type |
| created_at | TIMESTAMP | Creation timestamp |
| updated_at | TIMESTAMP | Last update timestamp |

#### protocols

Stepwise procedures for experiments.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| name | TEXT | Name of the protocol |
| description | TEXT | Description of the protocol |
| template_id | UUID | Reference to template (if derived) |
| steps | JSONB | Ordered protocol steps |
| parameters | JSONB | Protocol parameters |
| version | TEXT | Semantic version number |
| created_by | UUID | User who created the protocol |
| created_at | TIMESTAMP | Creation timestamp |
| updated_at | TIMESTAMP | Last update timestamp |

#### experiment_results

Outcomes and measurements from experiments.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| experiment_id | UUID | Reference to experiments |
| molecule_id | UUID | Reference to molecules (optional) |
| mixture_id | UUID | Reference to mixtures (optional) |
| tissue_type_id | UUID | Reference to tissue_types |
| property_type | TEXT | Type of property measured |
| value | FLOAT | Measured value |
| uncertainty | FLOAT | Measurement uncertainty |
| unit | TEXT | Unit of measurement |
| method_id | UUID | Measurement method used |
| provenance | TEXT | Data provenance information |
| created_by | UUID | User who created the result |
| created_at | TIMESTAMP | Creation timestamp |
| updated_at | TIMESTAMP | Last update timestamp |

### Bridge Tables

#### molecule_experiment_bridge

Links molecules to experiments.

| Column | Type | Description |
|--------|------|-------------|
| molecule_id | UUID | Reference to molecules |
| experiment_id | UUID | Reference to experiments |
| relationship_type | TEXT | Type of relationship |

#### mixture_experiment_bridge

Links mixtures to experiments.

| Column | Type | Description |
|--------|------|-------------|
| mixture_id | UUID | Reference to mixtures |
| experiment_id | UUID | Reference to experiments |
| relationship_type | TEXT | Type of relationship |

## Indexes

The following indexes are created for performance optimization:

```sql
-- Experiment indexes
CREATE INDEX idx_experiments_type ON experiments(experiment_type_id);
CREATE INDEX idx_experiments_tissue ON experiments(tissue_type_id);
CREATE INDEX idx_experiments_created_by ON experiments(created_by);
CREATE INDEX idx_experiments_quality ON experiments(quality_score);

-- Result indexes
CREATE INDEX idx_experiment_results_experiment ON experiment_results(experiment_id);
CREATE INDEX idx_experiment_results_molecule ON experiment_results(molecule_id);
CREATE INDEX idx_experiment_results_mixture ON experiment_results(mixture_id);
CREATE INDEX idx_experiment_results_property ON experiment_results(property_type);

-- Bridge table indexes
CREATE INDEX idx_molecule_experiment_molecule ON molecule_experiment_bridge(molecule_id);
CREATE INDEX idx_molecule_experiment_experiment ON molecule_experiment_bridge(experiment_id);
CREATE INDEX idx_mixture_experiment_mixture ON mixture_experiment_bridge(mixture_id);
CREATE INDEX idx_mixture_experiment_experiment ON mixture_experiment_bridge(experiment_id);
```

## Views

The following views are created for convenient access:

```sql
-- Molecule experiment results view
CREATE OR REPLACE VIEW molecule_experiment_results AS
SELECT 
    er.id,
    m.name as molecule_name,
    m.smiles as molecule_smiles,
    et.name as experiment_type,
    tt.name as tissue_type,
    tt.species,
    e.temperature,
    er.property_type,
    er.value,
    er.uncertainty,
    er.unit,
    er.method_id,
    e.created_at
FROM 
    experiment_results er
JOIN 
    experiments e ON er.experiment_id = e.id
JOIN 
    experiment_types et ON e.experiment_type_id = et.id
JOIN 
    tissue_types tt ON er.tissue_type_id = tt.id
JOIN 
    molecules m ON er.molecule_id = m.id
WHERE 
    er.molecule_id IS NOT NULL;

-- Mixture experiment results view
CREATE OR REPLACE VIEW mixture_experiment_results AS
SELECT 
    er.id,
    mx.name as mixture_name,
    et.name as experiment_type,
    tt.name as tissue_type,
    tt.species,
    e.temperature,
    er.property_type,
    er.value,
    er.uncertainty,
    er.unit,
    er.method_id,
    e.created_at
FROM 
    experiment_results er
JOIN 
    experiments e ON er.experiment_id = e.id
JOIN 
    experiment_types et ON e.experiment_type_id = et.id
JOIN 
    tissue_types tt ON er.tissue_type_id = tt.id
JOIN 
    mixtures mx ON er.mixture_id = mx.id
WHERE 
    er.mixture_id IS NOT NULL;
```

## Schema Evolution

The schema is designed to be extensible for future needs. Planned extensions include:

1. Support for equipment-specific configurations
2. External database cross-references
3. Expanded protocol parameters for specialized protocols
4. Enhanced uncertainty modeling with covariance data

## Data Migration Strategy

Migration from the current schema will follow these principles:

1. Preserve all existing data
2. Enrich data with additional metadata where possible
3. Create bridge table entries to maintain relationships
4. Convert text fields to structured data where appropriate
5. Apply defaults for new required fields

Detailed migration scripts will be provided for each stage of the migration process.