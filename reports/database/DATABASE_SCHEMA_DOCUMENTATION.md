# CryoProtect v2 Database Schema Documentation

## Overview

This document provides comprehensive documentation of the CryoProtect v2 database schema following the remediation project. It includes detailed information about tables, relationships, RLS policies, indexes, and application roles.

## Database Architecture

The CryoProtect v2 database is built on Supabase (PostgreSQL) and follows a normalized relational design. The schema implements:

- **Third Normal Form (3NF)** to eliminate data redundancy
- **Row Level Security (RLS)** for fine-grained access control
- **Junction tables** for many-to-many relationships
- **Foreign key constraints** for data integrity
- **Performance indexes** for query optimization

## Table Descriptions and Relationships

### Core Tables

#### 1. molecules

Primary table storing information about cryoprotectant molecules.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| pubchem_cid | INTEGER | PubChem Compound ID | UNIQUE |
| name | VARCHAR(255) | Molecule name | NOT NULL |
| formula | VARCHAR(100) | Chemical formula | |
| smiles | TEXT | SMILES notation | |
| inchikey | VARCHAR(27) | InChI Key | UNIQUE, NOT NULL |
| molecular_weight | NUMERIC | Molecular weight (g/mol) | |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- UNIQUE INDEX on `inchikey`
- UNIQUE INDEX on `pubchem_cid`
- INDEX on `name` (with trigram for search)
- INDEX on `created_by` (for RLS)

**Foreign Keys:**
- `created_by` REFERENCES auth.users(id)

#### 2. property_types

Defines types of properties that can be measured or calculated.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| name | VARCHAR(100) | Property name | UNIQUE, NOT NULL |
| data_type | VARCHAR(50) | Data type (numeric, text, boolean) | NOT NULL |
| description | TEXT | Property description | |
| units | VARCHAR(50) | Standard units | |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- UNIQUE INDEX on `name`
- INDEX on `created_by` (for RLS)

**Foreign Keys:**
- `created_by` REFERENCES auth.users(id)

#### 3. molecular_properties

Stores properties of individual molecules.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| molecule_id | UUID | Reference to molecule | FK, NOT NULL |
| property_type_id | UUID | Reference to property type | FK, NOT NULL |
| numeric_value | NUMERIC | Numeric property value | |
| text_value | TEXT | Text property value | |
| boolean_value | BOOLEAN | Boolean property value | |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- INDEX on `molecule_id`
- INDEX on `property_type_id`
- INDEX on `created_by` (for RLS)
- COMPOSITE INDEX on (molecule_id, property_type_id)

**Foreign Keys:**
- `molecule_id` REFERENCES molecules(id) ON DELETE CASCADE
- `property_type_id` REFERENCES property_types(id)
- `created_by` REFERENCES auth.users(id)

**Constraints:**
- CHECK constraint ensuring only one value type is not null

#### 4. mixtures

Information about mixtures of molecules.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| name | VARCHAR(255) | Mixture name | NOT NULL |
| description | TEXT | Mixture description | |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- INDEX on `name` (with trigram for search)
- INDEX on `created_by` (for RLS)

**Foreign Keys:**
- `created_by` REFERENCES auth.users(id)

#### 5. mixture_components

Junction table linking molecules to mixtures with concentration information.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| mixture_id | UUID | Reference to mixture | FK, NOT NULL |
| molecule_id | UUID | Reference to molecule | FK, NOT NULL |
| concentration | NUMERIC | Concentration value | NOT NULL |
| concentration_unit | VARCHAR(50) | Unit of concentration | NOT NULL |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- INDEX on `mixture_id`
- INDEX on `molecule_id`
- INDEX on `created_by` (for RLS)
- UNIQUE INDEX on (mixture_id, molecule_id)

**Foreign Keys:**
- `mixture_id` REFERENCES mixtures(id) ON DELETE CASCADE
- `molecule_id` REFERENCES molecules(id) ON DELETE CASCADE
- `created_by` REFERENCES auth.users(id)

### Prediction and Experiment Tables

#### 6. calculation_methods

Methods used for property predictions.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| name | VARCHAR(100) | Method name | UNIQUE, NOT NULL |
| description | TEXT | Method description | |
| version | VARCHAR(50) | Version information | |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- UNIQUE INDEX on `name`
- INDEX on `created_by` (for RLS)

**Foreign Keys:**
- `created_by` REFERENCES auth.users(id)

#### 7. predictions

Computational predictions of mixture or molecule properties.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| molecule_id | UUID | Reference to molecule | FK, NULL |
| mixture_id | UUID | Reference to mixture | FK, NULL |
| property_type_id | UUID | Reference to property type | FK, NOT NULL |
| calculation_method_id | UUID | Reference to calculation method | FK, NOT NULL |
| numeric_value | NUMERIC | Numeric prediction value | |
| text_value | TEXT | Text prediction value | |
| boolean_value | BOOLEAN | Boolean prediction value | |
| confidence | NUMERIC | Confidence level (0-1) | |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- INDEX on `molecule_id`
- INDEX on `mixture_id`
- INDEX on `property_type_id`
- INDEX on `calculation_method_id`
- INDEX on `created_by` (for RLS)
- COMPOSITE INDEX on (mixture_id, property_type_id, calculation_method_id)
- COMPOSITE INDEX on (molecule_id, property_type_id, calculation_method_id)

**Foreign Keys:**
- `molecule_id` REFERENCES molecules(id) ON DELETE CASCADE
- `mixture_id` REFERENCES mixtures(id) ON DELETE CASCADE
- `property_type_id` REFERENCES property_types(id)
- `calculation_method_id` REFERENCES calculation_methods(id)
- `created_by` REFERENCES auth.users(id)

**Constraints:**
- CHECK constraint ensuring either molecule_id or mixture_id is not null, but not both

#### 8. experiments

Experimental data for mixtures or molecules.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| mixture_id | UUID | Reference to mixture | FK, NULL |
| molecule_id | UUID | Reference to molecule | FK, NULL |
| property_type_id | UUID | Reference to property type | FK, NOT NULL |
| numeric_value | NUMERIC | Numeric experimental value | |
| text_value | TEXT | Text experimental value | |
| boolean_value | BOOLEAN | Boolean experimental value | |
| experimental_conditions | TEXT | Conditions of the experiment | |
| date_performed | DATE | Date experiment was performed | |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- INDEX on `mixture_id`
- INDEX on `molecule_id`
- INDEX on `property_type_id`
- INDEX on `created_by` (for RLS)
- INDEX on `date_performed`

**Foreign Keys:**
- `mixture_id` REFERENCES mixtures(id) ON DELETE CASCADE
- `molecule_id` REFERENCES molecules(id) ON DELETE CASCADE
- `property_type_id` REFERENCES property_types(id)
- `created_by` REFERENCES auth.users(id)

**Constraints:**
- CHECK constraint ensuring either molecule_id or mixture_id is not null, but not both

### Junction Tables

#### 9. molecule_proteins

Junction table linking molecules to proteins with binding information.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| molecule_id | UUID | Reference to molecule | FK, NOT NULL |
| protein_id | UUID | Reference to protein | FK, NOT NULL |
| binding_affinity | NUMERIC | Binding affinity value | |
| interaction_type | VARCHAR(100) | Type of interaction | |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- INDEX on `molecule_id`
- INDEX on `protein_id`
- INDEX on `created_by` (for RLS)
- UNIQUE INDEX on (molecule_id, protein_id)

**Foreign Keys:**
- `molecule_id` REFERENCES molecules(id) ON DELETE CASCADE
- `protein_id` REFERENCES proteins(id) ON DELETE CASCADE
- `created_by` REFERENCES auth.users(id)

#### 10. molecule_experiments

Junction table linking molecules directly to experiments.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| molecule_id | UUID | Reference to molecule | FK, NOT NULL |
| experiment_id | UUID | Reference to experiment | FK, NOT NULL |
| role | VARCHAR(100) | Role of molecule in experiment | |
| concentration | NUMERIC | Concentration value | |
| concentration_unit | VARCHAR(50) | Unit of concentration | |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- INDEX on `molecule_id`
- INDEX on `experiment_id`
- INDEX on `created_by` (for RLS)
- UNIQUE INDEX on (molecule_id, experiment_id)

**Foreign Keys:**
- `molecule_id` REFERENCES molecules(id) ON DELETE CASCADE
- `experiment_id` REFERENCES experiments(id) ON DELETE CASCADE
- `created_by` REFERENCES auth.users(id)

### Organization Tables

#### 11. projects

Research projects grouping data.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| name | VARCHAR(255) | Project name | NOT NULL |
| description | TEXT | Project description | |
| team_id | UUID | Reference to team | FK |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- INDEX on `team_id`
- INDEX on `created_by` (for RLS)

**Foreign Keys:**
- `team_id` REFERENCES teams(id)
- `created_by` REFERENCES auth.users(id)

#### 12. teams

Research teams.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| name | VARCHAR(255) | Team name | NOT NULL |
| description | TEXT | Team description | |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |
| created_by | UUID | User who created the record | FK to auth.users |

**Indexes:**
- PRIMARY KEY on `id`
- INDEX on `created_by` (for RLS)

**Foreign Keys:**
- `created_by` REFERENCES auth.users(id)

#### 13. user_profile

User accounts and metadata.

| Column | Type | Description | Constraints |
|--------|------|-------------|-------------|
| id | UUID | Primary key | PK, NOT NULL |
| auth_user_id | UUID | Reference to auth.users | FK, NOT NULL |
| name | VARCHAR(255) | User's full name | |
| email | VARCHAR(255) | User's email | |
| team_id | UUID | Reference to team | FK |
| created_at | TIMESTAMPTZ | Creation timestamp | DEFAULT NOW() |
| updated_at | TIMESTAMPTZ | Last update timestamp | DEFAULT NOW() |

**Indexes:**
- PRIMARY KEY on `id`
- UNIQUE INDEX on `auth_user_id`
- INDEX on `team_id`
- INDEX on `email`

**Foreign Keys:**
- `auth_user_id` REFERENCES auth.users(id) ON DELETE CASCADE
- `team_id` REFERENCES teams(id)

## Entity-Relationship Diagram

```
+----------------+      +---------------------+      +------------------+
|   molecules    |<-----| molecular_properties|----->|  property_types  |
+----------------+      +---------------------+      +------------------+
        ^                        ^
        |                        |
        v                        |
+----------------+      +---------------------+      +------------------+
|mixture_components|<----|    mixtures        |      | calculation_methods|
+----------------+      +---------------------+      +------------------+
        ^                        ^                        ^
        |                        |                        |
        v                        v                        |
+----------------+      +---------------------+      +------------------+
|  predictions   |<-----|   experiments       |----->|experiment_properties|
+----------------+      +---------------------+      +------------------+
        ^                        ^
        |                        |
        v                        v
+----------------+      +---------------------+
|molecule_proteins|      |molecule_experiments|
+----------------+      +---------------------+
        ^                        ^
        |                        |
        v                        |
+----------------+      +---------------------+
|   projects     |<-----|      teams          |
+----------------+      +---------------------+
        ^                        ^
        |                        |
        v                        v
+----------------+
| user_profile   |
+----------------+
```

## Row Level Security (RLS) Policies

RLS is enabled on all tables to control data access based on user identity. The following policies are implemented:

### Standard RLS Policies

For tables with a `created_by` column:

1. **Users access their own data**
   - Policy name: `{table_name}_user_access`
   - Operation: SELECT
   - Using expression: `auth.uid() = created_by`
   - Description: Users can view records they created

2. **Users modify their own data**
   - Policy name: `{table_name}_user_modify`
   - Operations: UPDATE, DELETE
   - Using expression: `auth.uid() = created_by`
   - Description: Users can modify or delete records they created

3. **Users insert data**
   - Policy name: `{table_name}_user_insert`
   - Operation: INSERT
   - With check expression: `auth.uid() = created_by`
   - Description: Users can insert records with their user ID

### Special RLS Policies

For tables without a `created_by` column or with special access requirements:

1. **Public read-only access**
   - Policy name: `{table_name}_public_read`
   - Operation: SELECT
   - Using expression: `true`
   - Description: Anyone can view these records (e.g., property_types)

2. **Team-based access**
   - Policy name: `{table_name}_team_access`
   - Operation: SELECT
   - Using expression: `auth.uid() IN (SELECT auth_user_id FROM user_profile WHERE team_id = {table_name}.team_id)`
   - Description: Users can view records associated with their team

3. **Service role access**
   - Policy name: `{table_name}_service_role`
   - Operations: ALL
   - Using expression: `auth.role() = 'service_role'`
   - Description: Service role can perform all operations

## Indexes and Performance Considerations

### Primary Indexes

- Every table has a PRIMARY KEY index on its `id` column

### Foreign Key Indexes

- All foreign key columns have indexes to improve join performance
- Composite indexes are created for commonly joined columns

### Performance Indexes

1. **Text Search Indexes**
   - Trigram indexes on `molecules.name` and `mixtures.name` for efficient text search
   - Full-text search indexes on description fields

2. **RLS Performance Indexes**
   - Indexes on all `created_by` columns to improve RLS policy evaluation
   - Index on `user_profile.auth_user_id` for efficient user lookups

3. **Query Optimization Indexes**
   - Composite indexes for common query patterns
   - Indexes on date fields for time-based queries

### Performance Considerations

1. **Query Patterns**
   - Indexes are optimized for common query patterns
   - Complex joins are minimized through proper normalization

2. **Data Volume**
   - Partitioning strategy for large tables (e.g., molecular_properties)
   - Selective indexing to balance performance and storage

3. **RLS Impact**
   - RLS policies are designed to minimize performance impact
   - Indexes support efficient policy evaluation

## Application Roles and Permissions

### User Roles

1. **Anonymous Users**
   - Role: `anon`
   - Permissions: Limited read-only access to non-sensitive data
   - Tables accessible: property_types, calculation_methods

2. **Authenticated Users**
   - Role: `authenticated`
   - Permissions: 
     - Read access to all public data
     - Create, read, update, delete access to their own data
     - Read access to team data (if part of a team)
   - Tables accessible: All tables with appropriate RLS policies

3. **Service Role**
   - Role: `service_role`
   - Permissions: Full access to all tables
   - Used for: Backend services, data migration, maintenance tasks

### Permission Management

1. **Default Grants**
   ```sql
   -- Grant usage on schema
   GRANT USAGE ON SCHEMA public TO authenticated, anon;
   
   -- Grant select on public tables
   GRANT SELECT ON property_types, calculation_methods TO anon;
   
   -- Grant all on all tables to authenticated users (RLS will restrict)
   GRANT ALL ON ALL TABLES IN SCHEMA public TO authenticated;
   ```

2. **RLS Enforcement**
   - RLS is enabled on all tables
   - Policies restrict access based on user identity
   - Service role bypasses RLS for administrative functions

## Views and Functions

### Views

1. **molecule_with_properties**
   - Combines molecules with their properties in JSON format
   - Simplifies querying of molecules with all their properties

2. **mixture_with_components**
   - Combines mixtures with their components in JSON format
   - Includes molecule details for each component

3. **experiment_with_results**
   - Combines experiments with their results and conditions
   - Simplifies querying of experimental data

### Functions

1. **import_molecule_from_pubchem(cid INTEGER)**
   - Imports molecule data from PubChem by CID
   - Returns the ID of the created or existing molecule

2. **calculate_mixture_score(mixture_id UUID)**
   - Calculates a cryoprotection score for a mixture
   - Based on molecular properties and experimental data

3. **compare_prediction_with_experiment(prediction_id UUID, experiment_id UUID)**
   - Compares predicted values with experimental results
   - Returns difference metrics and statistical analysis

## Maintenance Considerations

1. **Schema Evolution**
   - Use migration scripts for schema changes
   - Maintain backward compatibility when possible
   - Document schema changes thoroughly

2. **Index Maintenance**
   - Regularly analyze and reindex tables
   - Monitor index usage and performance
   - Add or remove indexes based on query patterns

3. **RLS Policy Updates**
   - Update RLS policies when adding new tables
   - Test policy effectiveness with different user roles
   - Document policy changes

4. **Performance Monitoring**
   - Monitor query performance regularly
   - Identify and optimize slow queries
   - Adjust indexes based on performance data