# Unified Experimental Schema Design

## Overview

This document presents the complete schema design for the enhanced experimental data system in CryoProtect. The design unifies and extends the current experimental schema to provide a more robust, scientifically accurate framework for managing cryoprotection experimental data.

## Design Goals

1. **Data Integrity**: Ensure data consistency through properly defined relationships and constraints
2. **Scientific Accuracy**: Support uncertainty measurements, provenance tracking, and physical model validation 
3. **Extensibility**: Design for extensible protocols, property types, and experimental conditions
4. **Compatibility**: Maintain backward compatibility with existing data and APIs
5. **Traceability**: Enable full traceability of experimental data from raw measurements to final results

## Schema Diagram

```
┌───────────────────────┐     ┌─────────────────────────┐     ┌───────────────────────┐
│ experiment_types      │     │ experiments             │     │ experiment_results    │
├───────────────────────┤     ├─────────────────────────┤     ├───────────────────────┤
│ id                    │     │ id                      │     │ id                    │
│ name                  │     │ experiment_type_id      │◄────┤ experiment_id         │
│ description           │     │ title                   │     │ molecule_id           │
│ protocol_details      │     │ description             │     │ mixture_id            │
│ created_by            │     │ protocol_id             │     │ tissue_type_id        │◄─┐
│ created_at            │     │ date_performed          │     │ concentration         │  │
│ updated_at            │     │ temperature             │     │ concentration_unit    │  │
└───────────┬───────────┘     │ cooling_rate            │     │ viability_percentage  │  │
            │                 │ thawing_rate            │     │ recovery_rate         │  │
            │                 │ parameters              │     │ functionality_score   │  │
            └────────────────►│ version                 │     │ result_details        │  │
                              │ provenance              │     │ uncertainty           │  │
                              │ created_by              │     │ provenance            │  │
                              │ created_at              │     │ notes                 │  │
                              │ updated_at              │     │ created_by            │  │
                              └─────────────┬───────────┘     │ created_at            │  │
                                            │                 │ updated_at            │  │
                                            │                 └───────────────────────┘  │
┌───────────────────────┐                  │                                            │
│ protocols             │                  │                                            │
├───────────────────────┤                  │                 ┌───────────────────────┐  │
│ id                    │                  │                 │ tissue_types          │  │
│ name                  │                  │                 ├───────────────────────┤  │
│ description           │                  │                 │ id                    │──┘
│ steps                 │                  │                 │ name                  │
│ parameters            │◄─────────────────┘                 │ description           │
│ version               │                                    │ species               │
│ parent_id             │                                    │ taxonomy_id           │
│ created_by            │                                    │ properties            │
│ created_at            │                                    │ created_by            │
│ updated_at            │                                    │ created_at            │
└───────────────────────┘                                    │ updated_at            │
                                                             └───────────────────────┘

┌───────────────────────┐     ┌─────────────────────────┐     ┌───────────────────────┐
│ equipment_types       │     │ equipment               │     │ experiment_equipment  │
├───────────────────────┤     ├─────────────────────────┤     ├───────────────────────┤
│ id                    │     │ id                      │     │ id                    │
│ name                  │─────►│ equipment_type_id      │◄────┤ experiment_id         │
│ description           │     │ name                    │     │ equipment_id          │
│ manufacturer          │     │ model                   │     │ parameters            │
│ created_by            │     │ description             │     │ created_by            │
│ created_at            │     │ parameters              │     │ created_at            │
│ updated_at            │     │ created_by              │     │ updated_at            │
└───────────────────────┘     │ created_at              │     └───────────────────────┘
                              │ updated_at              │
                              └─────────────────────────┘
                              
┌───────────────────────┐     ┌─────────────────────────┐     ┌───────────────────────┐
│ time_series           │     │ time_series_data        │     │ validation_rules      │
├───────────────────────┤     ├─────────────────────────┤     ├───────────────────────┤
│ id                    │     │ id                      │     │ id                    │
│ experiment_id         │     │ time_series_id          │◄────┤ property_type         │
│ name                  │─────►│ timestamp              │     │ rule_type             │
│ description           │     │ value                   │     │ parameters            │
│ property_type         │     │ uncertainty             │     │ description           │
│ unit                  │     │ provenance              │     │ created_by            │
│ metadata              │     │ created_by              │     │ created_at            │
│ created_by            │     │ created_at              │     │ updated_at            │
│ created_at            │     │ updated_at              │     └───────────────────────┘
│ updated_at            │     └─────────────────────────┘
└───────────────────────┘
```

## Table Definitions

### experiment_types
Categorizes different types of cryopreservation experiments.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| name | VARCHAR(255) | Name of the experiment type |
| description | TEXT | Detailed description |
| protocol_details | JSONB | Default protocol parameters |
| created_by | UUID | User who created the record |
| created_at | TIMESTAMPTZ | Creation timestamp |
| updated_at | TIMESTAMPTZ | Last update timestamp |

### protocols
Template-based protocol definitions with versioning support.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| name | VARCHAR(255) | Protocol name |
| description | TEXT | Detailed description |
| steps | JSONB | Ordered list of protocol steps |
| parameters | JSONB | Protocol parameters |
| version | INTEGER | Version number |
| parent_id | UUID | Reference to previous version |
| created_by | UUID | User who created the record |
| created_at | TIMESTAMPTZ | Creation timestamp |
| updated_at | TIMESTAMPTZ | Last update timestamp |

### tissue_types
Biological samples used in experiments.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| name | VARCHAR(255) | Tissue name |
| description | TEXT | Detailed description |
| species | VARCHAR(255) | Species name |
| taxonomy_id | INTEGER | NCBI taxonomy ID |
| properties | JSONB | Additional properties (e.g., cell density) |
| created_by | UUID | User who created the record |
| created_at | TIMESTAMPTZ | Creation timestamp |
| updated_at | TIMESTAMPTZ | Last update timestamp |

### experiments
Core experiment records.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| experiment_type_id | UUID | Reference to experiment_types |
| title | VARCHAR(255) | Experiment title |
| description | TEXT | Detailed description |
| protocol_id | UUID | Reference to protocols |
| date_performed | DATE | Date when experiment was performed |
| temperature | NUMERIC | Temperature value |
| cooling_rate | NUMERIC | Rate of temperature decrease |
| thawing_rate | NUMERIC | Rate of temperature increase |
| parameters | JSONB | Additional parameters |
| version | INTEGER | Version number |
| provenance | JSONB | Data source information |
| created_by | UUID | User who created the record |
| created_at | TIMESTAMPTZ | Creation timestamp |
| updated_at | TIMESTAMPTZ | Last update timestamp |

### experiment_results
Results from experiments on specific molecules or mixtures.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| experiment_id | UUID | Reference to experiments |
| molecule_id | UUID | Reference to molecules (optional) |
| mixture_id | UUID | Reference to mixtures (optional) |
| tissue_type_id | UUID | Reference to tissue_types |
| concentration | NUMERIC | Concentration value |
| concentration_unit | VARCHAR(50) | Unit of concentration |
| viability_percentage | NUMERIC | Cell survival rate |
| recovery_rate | NUMERIC | Recovery rate after thawing |
| functionality_score | NUMERIC | Functional assessment score |
| result_details | JSONB | Additional result data |
| uncertainty | JSONB | Uncertainty measurements |
| provenance | JSONB | Data source information |
| notes | TEXT | Additional notes |
| created_by | UUID | User who created the record |
| created_at | TIMESTAMPTZ | Creation timestamp |
| updated_at | TIMESTAMPTZ | Last update timestamp |

### equipment_types
Categories of laboratory equipment.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| name | VARCHAR(255) | Equipment type name |
| description | TEXT | Detailed description |
| manufacturer | VARCHAR(255) | Equipment manufacturer |
| created_by | UUID | User who created the record |
| created_at | TIMESTAMPTZ | Creation timestamp |
| updated_at | TIMESTAMPTZ | Last update timestamp |

### equipment
Specific equipment instances.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| equipment_type_id | UUID | Reference to equipment_types |
| name | VARCHAR(255) | Equipment name |
| model | VARCHAR(255) | Model number |
| description | TEXT | Detailed description |
| parameters | JSONB | Equipment parameters |
| created_by | UUID | User who created the record |
| created_at | TIMESTAMPTZ | Creation timestamp |
| updated_at | TIMESTAMPTZ | Last update timestamp |

### experiment_equipment
Links experiments to equipment used.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| experiment_id | UUID | Reference to experiments |
| equipment_id | UUID | Reference to equipment |
| parameters | JSONB | Equipment settings for this experiment |
| created_by | UUID | User who created the record |
| created_at | TIMESTAMPTZ | Creation timestamp |
| updated_at | TIMESTAMPTZ | Last update timestamp |

### time_series
Time-series data definitions.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| experiment_id | UUID | Reference to experiments |
| name | VARCHAR(255) | Time series name |
| description | TEXT | Detailed description |
| property_type | VARCHAR(255) | Type of property measured |
| unit | VARCHAR(50) | Unit of measurement |
| metadata | JSONB | Additional metadata |
| created_by | UUID | User who created the record |
| created_at | TIMESTAMPTZ | Creation timestamp |
| updated_at | TIMESTAMPTZ | Last update timestamp |

### time_series_data
Individual time-series data points.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| time_series_id | UUID | Reference to time_series |
| timestamp | TIMESTAMPTZ | Measurement timestamp |
| value | NUMERIC | Measured value |
| uncertainty | NUMERIC | Measurement uncertainty |
| provenance | JSONB | Data source information |
| created_by | UUID | User who created the record |
| created_at | TIMESTAMPTZ | Creation timestamp |
| updated_at | TIMESTAMPTZ | Last update timestamp |

### validation_rules
Rules for validating experimental data.

| Column | Type | Description |
|--------|------|-------------|
| id | UUID | Primary key |
| property_type | VARCHAR(255) | Type of property to validate |
| rule_type | VARCHAR(50) | Type of validation rule |
| parameters | JSONB | Rule parameters |
| description | TEXT | Detailed description |
| created_by | UUID | User who created the record |
| created_at | TIMESTAMPTZ | Creation timestamp |
| updated_at | TIMESTAMPTZ | Last update timestamp |

## Migration Strategy

The migration strategy will ensure seamless transition from the existing schema to the new unified schema:

1. **Create New Tables**: Add all new tables without dropping existing ones
2. **Copy Data**: Migrate data from existing tables to new tables
3. **Create Views**: Maintain compatibility with existing queries through views
4. **Verify Data**: Ensure data integrity after migration
5. **Update APIs**: Modify API endpoints to use new schema
6. **Test**: Comprehensive testing to ensure functionality
7. **Deprecate Old Tables**: Mark old tables as deprecated with a timeline for removal

## Data Validation

The schema includes built-in validation mechanisms:

1. **Foreign Key Constraints**: Ensure referential integrity
2. **Check Constraints**: Validate data values (e.g., temperature ranges)
3. **Uniqueness Constraints**: Prevent duplicate records
4. **Validation Rules**: Custom validation rules for specific properties
5. **Uncertainty Tracking**: Support for error propagation

## Indexes

Strategic indexes will optimize query performance:

1. **Primary Keys**: UUID indexes on all ID fields
2. **Foreign Keys**: Indexes on all foreign key references
3. **Search Fields**: Indexes on frequently searched fields (name, title)
4. **Range Queries**: Indexes on fields used in range queries (date, temperature)
5. **JSON Indexes**: GIN indexes on JSONB fields for efficient querying

## Views

Views will be created to simplify common queries and maintain compatibility:

1. **experiment_summary**: Aggregate view of experiment results
2. **molecular_experiment_results**: Results for specific molecules
3. **mixture_experiment_results**: Results for specific mixtures
4. **protocol_templates**: Simplified view of protocol templates
5. **experiment_history**: Time-series of all changes to experiments

## Security Considerations

The schema design incorporates proper security measures:

1. **Row-Level Security**: Policies to restrict access to data
2. **Audit Trails**: Tracking of all data modifications
3. **Versioning**: Support for non-destructive updates
4. **User Permissions**: Role-based access control

## Performance Considerations

To ensure optimal performance:

1. **Partitioning**: Time-series data will be partitioned by date
2. **Materialized Views**: For frequently accessed aggregated data
3. **Selective Indexes**: To balance between query speed and write performance
4. **JSON Compression**: For efficient storage of JSONB fields
5. **Connection Pooling**: For efficient database connections

## Implementation Plan

The implementation will proceed in phases:

1. **Phase 1**: Core schema tables and migration scripts
2. **Phase 2**: Advanced features (protocols, time-series)
3. **Phase 3**: Validation and analytics components
4. **Phase 4**: Integration with external systems

Each phase will include thorough testing and validation before proceeding to the next.