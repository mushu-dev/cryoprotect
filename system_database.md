# CryoProtect Analyzer Database Structure

## Overview

CryoProtect's database is designed following scientific data management best practices, implementing a flexible, relational schema optimized for molecular and experimental data. The database includes comprehensive audit trails, versioning, and fine-grained access control through Row Level Security (RLS).

## Core Database Principles

1. **Scientific Data Integrity**: Precise numeric types, explicit units, proper versioning
2. **Flexible Property Storage**: Support for various property types (numeric, text, boolean)
3. **Graph-based Relationships**: Many-to-many relationships between entities
4. **Multi-tenant Security**: Row-level security policies for collaborative access
5. **Audit and Provenance**: Track data lineage and modifications

## Technology

- **PostgreSQL 15.8+**: Advanced relational database
- **Supabase**: PostgreSQL with built-in Auth, RLS, and API
- **JSON/JSONB**: For flexible property storage and metadata
- **Direct Connection**: Performance-optimized PostgreSQL connection pooling
- **Fedora PostgreSQL**: Integration with Fedora-native PostgreSQL packages

## Database Schema

### Core Scientific Tables

#### `molecules`
Stores molecular structure and identification data.
- `id` (UUID, PK): Unique identifier
- `name` (VARCHAR): Molecule name
- `smiles` (VARCHAR): Canonical SMILES representation
- `inchi` (TEXT): InChI string
- `inchikey` (VARCHAR): InChIKey for indexing
- `formula` (VARCHAR): Molecular formula
- `molecular_weight` (NUMERIC): Calculated molecular weight
- `pubchem_cid` (INTEGER): PubChem Compound ID
- `cid` (INTEGER): Alternative PubChem ID
- `chembl_id` (VARCHAR): ChEMBL database ID
- `created_by` (UUID): User who created the record
- `is_public` (BOOLEAN): Public visibility flag
- `data_source` (VARCHAR): Origin of the data
- `version` (INTEGER): Version number
- `modification_history` (JSONB): Record of changes
- `properties` (JSONB): Cached property values
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `molecular_properties`
Stores properties of molecules with scientific precision.
- `id` (UUID, PK): Unique identifier
- `molecule_id` (UUID, FK): Reference to molecule
- `property_type_id` (UUID, FK): Type of property
- `numeric_value` (NUMERIC): Numeric property value
- `text_value` (TEXT): Text property value
- `boolean_value` (BOOLEAN): Boolean property value
- `unit` (VARCHAR): Unit of measurement
- `property_name` (VARCHAR): Name of property
- `property_type` (VARCHAR): Type of property data
- `source` (VARCHAR): Data source
- `sensitivity_level` (TEXT): Data sensitivity classification
- `created_by` (UUID): Creator reference
- `version` (INTEGER): Version number
- `modification_history` (JSONB): Change history
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `property_types`
Defines the types of properties that can be measured.
- `id` (UUID, PK): Unique identifier
- `name` (VARCHAR): Property name
- `description` (TEXT): Description of property
- `units` (VARCHAR): Default unit of measurement
- `data_type` (VARCHAR): Data type (numeric, text, boolean)
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `mixtures`
Stores information about cryoprotectant mixtures.
- `id` (UUID, PK): Unique identifier
- `name` (VARCHAR): Mixture name
- `description` (TEXT): Description
- `is_public` (BOOLEAN): Public visibility flag
- `project_id` (UUID, FK): Project reference
- `created_by` (UUID): Creator reference
- `properties` (JSONB): Cached property values
- `version` (INTEGER): Version number
- `modification_history` (JSONB): Change history
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `mixture_components`
Defines which molecules are in each mixture and their concentrations.
- `id` (UUID, PK): Unique identifier
- `mixture_id` (UUID, FK): Mixture reference
- `molecule_id` (UUID, FK): Molecule reference
- `concentration` (NUMERIC): Component concentration
- `concentration_unit` (VARCHAR): Unit of concentration
- `created_by` (UUID): Creator reference
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `experiments`
Stores experimental data from laboratory tests.
- `id` (UUID, PK): Unique identifier
- `mixture_id` (UUID, FK): Mixture reference
- `molecule_id` (UUID, FK): Optional direct molecule reference
- `property_type_id` (UUID, FK): Property type reference
- `numeric_value` (NUMERIC): Numeric result
- `text_value` (TEXT): Text result
- `boolean_value` (BOOLEAN): Boolean result
- `experimental_conditions` (TEXT): Conditions description
- `created_by` (UUID): Creator reference
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `experiment_properties`
Additional properties associated with experiments.
- `id` (UUID, PK): Unique identifier
- `experiment_id` (UUID, FK): Experiment reference
- `property_type_id` (UUID, FK): Property type reference
- `numeric_value` (NUMERIC): Numeric value
- `text_value` (TEXT): Text value
- `boolean_value` (BOOLEAN): Boolean value
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `predictions`
Stores computational predictions of mixture/molecule properties.
- `id` (UUID, PK): Unique identifier
- `mixture_id` (UUID, FK): Mixture reference
- `molecule_id` (UUID, FK): Optional direct molecule reference
- `property_type_id` (UUID, FK): Property type reference
- `calculation_method_id` (UUID, FK): Method reference
- `numeric_value` (NUMERIC): Numeric prediction
- `text_value` (TEXT): Text prediction
- `boolean_value` (BOOLEAN): Boolean prediction
- `confidence` (NUMERIC): Confidence level (0-1)
- `created_by` (UUID): Creator reference
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `calculation_methods`
Stores different methods used for property predictions.
- `id` (UUID, PK): Unique identifier
- `name` (VARCHAR): Method name
- `description` (TEXT): Method description
- `version` (VARCHAR): Version information
- `parameters` (JSONB): Method parameters
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

### Collaboration and Organization Tables

#### `user_profile`
Extended user profile information.
- `id` (UUID, PK): Unique identifier
- `auth_user_id` (UUID): Supabase Auth user ID
- `display_name` (VARCHAR): User display name
- `email` (VARCHAR): User email
- `affiliation` (VARCHAR): Organizational affiliation
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `teams`
Research teams for collaborative work.
- `id` (UUID, PK): Unique identifier
- `name` (VARCHAR): Team name
- `description` (TEXT): Team description
- `created_by` (UUID): Creator reference
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `team_members`
Maps users to teams.
- `id` (UUID, PK): Unique identifier
- `team_id` (UUID, FK): Team reference
- `user_id` (UUID, FK): User reference
- `role` (VARCHAR): Role within team
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `projects`
Research projects managed by teams.
- `id` (UUID, PK): Unique identifier
- `name` (VARCHAR): Project name
- `description` (TEXT): Project description
- `team_id` (UUID, FK): Team reference
- `created_by` (UUID): Creator reference
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `shared_resources`
Resources shared between users or teams.
- `id` (UUID, PK): Unique identifier
- `resource_type` (VARCHAR): Type of resource
- `resource_id` (UUID): ID of the shared resource
- `team_id` (UUID, FK): Team reference
- `created_by` (UUID, FK): Creator reference
- `permissions` (JSONB): Permission details
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

### Special-Purpose Tables

#### `toxicity_data_source`
Sources for toxicity information.
- `id` (UUID, PK): Unique identifier
- `name` (VARCHAR): Source name
- `url` (VARCHAR): Source URL
- `description` (TEXT): Source description
- `reliability_score` (NUMERIC): Reliability rating
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `proteins`
Protein targets for molecular interactions.
- `id` (UUID, PK): Unique identifier
- `name` (VARCHAR): Protein name
- `uniprot_id` (VARCHAR): UniProt identifier
- `sequence` (TEXT): Amino acid sequence
- `structure` (JSONB): Structure information
- `created_by` (UUID, FK): Creator reference
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `molecule_proteins`
Relations between molecules and proteins.
- `id` (UUID, PK): Unique identifier
- `molecule_id` (UUID, FK): Molecule reference
- `protein_id` (UUID, FK): Protein reference
- `interaction_type` (VARCHAR): Type of interaction
- `binding_affinity` (NUMERIC): Binding affinity value
- `created_by` (UUID, FK): Creator reference
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

#### `lab_verifications`
Verification of experimental results.
- `id` (UUID, PK): Unique identifier
- `experiment_id` (UUID, FK): Experiment reference
- `verifier` (UUID, FK): User who verified
- `verification_date` (TIMESTAMP): When verified
- `notes` (TEXT): Verification notes
- `status` (VARCHAR): Verification status
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

### Process Tables

#### `chembl_import_batches`
Tracks ChEMBL data import operations.
- `batch_id` (UUID, PK): Unique batch identifier
- `start_time` (TIMESTAMP): Start time of import
- `end_time` (TIMESTAMP): End time of import
- `status` (VARCHAR): Import status
- `molecule_count` (INTEGER): Number of molecules
- `success_count` (INTEGER): Successful imports
- `error_count` (INTEGER): Error count
- `metadata` (JSONB): Additional metadata

#### `chembl_import_logs`
Detailed logs for ChEMBL imports.
- `id` (UUID, PK): Unique identifier
- `batch_id` (UUID, FK): Batch reference
- `molecule_id` (UUID, FK): Molecule reference (if successful)
- `chembl_id` (VARCHAR): ChEMBL ID
- `status` (VARCHAR): Status (success/error)
- `message` (TEXT): Log message
- `timestamp` (TIMESTAMP): Event timestamp
- `details` (JSONB): Additional details

#### `property_calculation_queue`
Queue for asynchronous property calculations.
- `id` (UUID, PK): Unique identifier
- `molecule_id` (UUID, FK): Molecule reference
- `property_type` (VARCHAR): Property to calculate
- `status` (VARCHAR): Job status
- `priority` (INTEGER): Job priority
- `created_at`, `updated_at` (TIMESTAMP): Timestamps

### Database Views

#### `molecule_with_properties`
Provides molecules with their properties in JSON format.
```sql
SELECT 
    m.id, m.name, m.smiles, m.inchikey, m.formula, m.molecular_weight,
    m.pubchem_cid, m.chembl_id, m.created_by, m.is_public,
    m.created_at, m.updated_at,
    jsonb_object_agg(pt.name, 
        CASE 
            WHEN pt.data_type = 'numeric' THEN to_jsonb(mp.numeric_value)
            WHEN pt.data_type = 'text' THEN to_jsonb(mp.text_value)
            WHEN pt.data_type = 'boolean' THEN to_jsonb(mp.boolean_value)
            ELSE NULL
        END
    ) AS properties
FROM 
    public.molecules m
LEFT JOIN 
    public.molecular_properties mp ON m.id = mp.molecule_id
LEFT JOIN 
    public.property_types pt ON mp.property_type_id = pt.id
GROUP BY 
    m.id;
```

#### `mixture_with_components`
Provides mixtures with their components in JSON format.
```sql
SELECT 
    mix.id, mix.name, mix.description, mix.is_public,
    mix.project_id, mix.created_by, mix.created_at, mix.updated_at,
    jsonb_agg(
        jsonb_build_object(
            'molecule_id', mol.id,
            'name', mol.name,
            'concentration', mc.concentration,
            'concentration_unit', mc.concentration_unit
        )
    ) AS components
FROM 
    public.mixtures mix
LEFT JOIN 
    public.mixture_components mc ON mix.id = mc.mixture_id
LEFT JOIN 
    public.molecules mol ON mc.molecule_id = mol.id
GROUP BY 
    mix.id;
```

#### `experiment_with_results`
Provides experiments with their results in JSON format.
```sql
SELECT 
    e.id, e.mixture_id, e.molecule_id, e.property_type_id,
    e.created_by, e.created_at, e.updated_at,
    pt.name AS property_name, pt.data_type,
    CASE 
        WHEN pt.data_type = 'numeric' THEN e.numeric_value
        WHEN pt.data_type = 'text' THEN e.text_value
        WHEN pt.data_type = 'boolean' THEN e.boolean_value
        ELSE NULL
    END AS value,
    pt.units,
    mix.name AS mixture_name,
    mol.name AS molecule_name
FROM 
    public.experiments e
JOIN 
    public.property_types pt ON e.property_type_id = pt.id
LEFT JOIN 
    public.mixtures mix ON e.mixture_id = mix.id
LEFT JOIN 
    public.molecules mol ON e.molecule_id = mol.id;
```

## Database Migrations

The database evolution is managed through a series of migration scripts:

1. `001_initial_schema.sql`: Basic tables (molecules, properties, mixtures)
2. `002_projects_schema.sql`: Project management functionality
3. `003_teams_schema.sql`: Team collaboration schema
4. `004_export_sharing_schema.sql`: Data sharing functionality
5. `005_scientific_schema.sql`: Enhanced scientific data schema
6. `006_rls_policies.sql`: Row Level Security implementation
7. `007_seed_scientific_data.sql`: Initial reference data
8. `009_protocol_storage_schema.sql`: Protocol design schema
9. `010_performance_indexes.sql`: Performance optimization indexes
10. `011_rbac_schema.sql`: Role-based access control
11. `012_toxicity_schema.sql`: Toxicity data schema
12. `013_lab_verification_schema.sql`: Lab verification workflow

Additional migrations handle:
- Foreign key fixes
- RLS policy enhancements
- Performance optimizations
- Schema standardization

## Row Level Security (RLS)

CryoProtect implements fine-grained Row Level Security policies to protect data:

### Basic RLS Pattern
```sql
ALTER TABLE [table_name] ENABLE ROW LEVEL SECURITY;

-- Read access policy
CREATE POLICY "[policy_name]" 
  ON [table_name]
  FOR SELECT
  USING ([access_condition]);

-- Write access policy  
CREATE POLICY "[policy_name]" 
  ON [table_name]
  FOR INSERT/UPDATE/DELETE
  USING ([access_condition]);
```

### Access Control Patterns

1. **Creator-based access**:
   ```sql
   (auth.uid() = created_by)
   ```

2. **Team-based access**:
   ```sql
   EXISTS (
     SELECT 1 FROM team_members
     WHERE team_members.team_id = [table].team_id
       AND team_members.user_id = auth.uid()
   )
   ```

3. **Project-based access**:
   ```sql
   EXISTS (
     SELECT 1 FROM projects
     JOIN team_members ON team_members.team_id = projects.team_id
     WHERE projects.id = [table].project_id
       AND team_members.user_id = auth.uid()
   )
   ```

4. **Public/private distinction**:
   ```sql
   ([table].is_public = true OR auth.uid() = [table].created_by)
   ```

## Database Functions

Key PostgreSQL functions enhance the database functionality:

1. **`trigger_set_timestamp()`**: Updates `updated_at` timestamp on record modification
   
2. **`calculate_mixture_score(mixture_id)`**: Calculates weighted score for mixture components

3. **`compare_prediction_with_experiment(mixture_id, property_id)`**: Compares prediction with experimental results

4. **`import_molecule_from_pubchem(cid, user_id)`**: Imports molecule from PubChem

5. **`get_molecular_similarity(molecule1_id, molecule2_id)`**: Calculates molecular similarity

## Molecular Filtering and Search

### Structure-Based Similarity Search

The database implements several molecular similarity search methods:

1. **Morgan/ECFP Fingerprints**: Extended-Connectivity Fingerprints for structural similarity
   ```sql
   CREATE OR REPLACE FUNCTION morganbv_similarity(smiles1 TEXT, smiles2 TEXT)
   RETURNS FLOAT AS $$
     from rdkit import Chem
     from rdkit.Chem import AllChem
     import rdkit.DataStructs

     # Generate molecules from SMILES
     mol1 = Chem.MolFromSmiles(smiles1)
     mol2 = Chem.MolFromSmiles(smiles2)
     
     if mol1 is None or mol2 is None:
       return 0.0
       
     # Generate Morgan fingerprints (ECFP4)
     fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
     fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
     
     # Calculate Tanimoto similarity
     similarity = rdkit.DataStructs.TanimotoSimilarity(fp1, fp2)
     return similarity
   $$ LANGUAGE plpython3u;
   ```

2. **Topological Fingerprints**: Path-based structural features
   ```sql
   CREATE OR REPLACE FUNCTION topological_similarity(smiles1 TEXT, smiles2 TEXT)
   RETURNS FLOAT AS $$
     from rdkit import Chem
     from rdkit.Chem.Fingerprints import FingerprintMols
     import rdkit.DataStructs

     mol1 = Chem.MolFromSmiles(smiles1)
     mol2 = Chem.MolFromSmiles(smiles2)
     
     if mol1 is None or mol2 is None:
       return 0.0
       
     fp1 = FingerprintMols.FingerprintMol(mol1)
     fp2 = FingerprintMols.FingerprintMol(mol2)
     
     similarity = rdkit.DataStructs.TanimotoSimilarity(fp1, fp2)
     return similarity
   $$ LANGUAGE plpython3u;
   ```

3. **MACCS Keys**: 166 structural keys for focused similarity
   ```sql
   CREATE OR REPLACE FUNCTION maccs_similarity(smiles1 TEXT, smiles2 TEXT)
   RETURNS FLOAT AS $$
     from rdkit import Chem
     from rdkit.Chem import MACCSkeys
     import rdkit.DataStructs

     mol1 = Chem.MolFromSmiles(smiles1)
     mol2 = Chem.MolFromSmiles(smiles2)
     
     if mol1 is None or mol2 is None:
       return 0.0
       
     fp1 = MACCSkeys.GenMACCSKeys(mol1)
     fp2 = MACCSkeys.GenMACCSKeys(mol2)
     
     similarity = rdkit.DataStructs.TanimotoSimilarity(fp1, fp2)
     return similarity
   $$ LANGUAGE plpython3u;
   ```

### Property-Based Molecular Filtering

Molecules are filtered based on scientific properties relevant for cryoprotection:

1. **Physicochemical Properties**:
   - Molecular weight (preferred range: 76-150 Da)
   - LogP (water solubility indicator)
   - Hydrogen bond donors/acceptors
   - Polar surface area
   - Rotatable bonds

2. **Cryoprotection-Specific Properties**:
   - Glass transition temperature
   - Ice formation inhibition
   - Membrane permeability
   - Cytotoxicity profiles
   - Viscosity at low temperatures

3. **Filter Implementation**:
   ```sql
   CREATE OR REPLACE FUNCTION filter_cryoprotectant_candidates()
   RETURNS TABLE(molecule_id UUID) AS $$
   BEGIN
     RETURN QUERY
     SELECT m.id
     FROM molecules m
     JOIN molecular_properties mp_mw ON m.id = mp_mw.molecule_id
     JOIN molecular_properties mp_hba ON m.id = mp_hba.molecule_id
     JOIN molecular_properties mp_hbd ON m.id = mp_hbd.molecule_id
     JOIN molecular_properties mp_logp ON m.id = mp_logp.molecule_id
     WHERE mp_mw.property_name = 'molecular_weight' 
       AND mp_mw.numeric_value BETWEEN 76 AND 150
       AND mp_hba.property_name = 'h_bond_acceptors'
       AND mp_hba.numeric_value >= 3
       AND mp_hbd.property_name = 'h_bond_donors'
       AND mp_hbd.numeric_value >= 2
       AND mp_logp.property_name = 'logp'
       AND mp_logp.numeric_value BETWEEN -3 AND 1;
   END;
   $$ LANGUAGE plpgsql;
   ```

## Performance Considerations

1. **Indexes**:
   - Primary indexes on all ID fields
   - Foreign key indexes for efficient joins
   - Specialized indexes for common query patterns:
     - `idx_molecules_inchikey` for chemical structure lookup
     - `idx_molecular_properties_molecule_id_property_type_id` for property access
     - `idx_mixture_components_mixture_id` for mixture composition

2. **Query Optimization**:
   - Views with pre-joined data for common queries
   - JSON aggregation for denormalized data retrieval
   - Property caching in JSONB columns for frequent access patterns

3. **Connection Pooling**:
   - Min connections: 2
   - Max connections: 10
   - Connection timeout: 30 seconds

## Fedora-Specific Database Configuration

### PostgreSQL Installation
```bash
# Install PostgreSQL 15
sudo dnf install -y postgresql15-server postgresql15-contrib postgresql15-devel

# Initialize the database
sudo /usr/pgsql-15/bin/postgresql-15-setup initdb

# Enable and start PostgreSQL service
sudo systemctl enable postgresql-15
sudo systemctl start postgresql-15
```

### SELinux Configuration
```bash
# Set appropriate SELinux context for PostgreSQL data directory
sudo semanage fcontext -a -t postgresql_db_t "/var/lib/pgsql/15/data(/.*)?)"
sudo restorecon -Rv /var/lib/pgsql/15/data

# Allow network connections if needed
sudo setsebool -P postgresql_can_connect_all 1
```

### IPv4/IPv6 Configuration

PostgreSQL can be configured to listen on specific IP addresses and protocols. For Supabase connectivity, it's often necessary to ensure IPv4 compatibility.

To configure PostgreSQL for IPv4 only:

1. Edit the PostgreSQL configuration file:
```bash
sudo vi /var/lib/pgsql/15/data/postgresql.conf
```

2. Set the listen_addresses parameter to use IPv4 only:
```
listen_addresses = '0.0.0.0'   # listens on all IPv4 addresses
```

3. Configure client authentication in pg_hba.conf:
```bash
sudo vi /var/lib/pgsql/15/data/pg_hba.conf
```

4. Add IPv4-specific entries:
```
# IPv4 local connections
host    cryoprotect     cryoprotect     127.0.0.1/32            md5
# IPv4 remote connections
host    cryoprotect     cryoprotect     0.0.0.0/0               md5
```

5. Restart PostgreSQL:
```bash
sudo systemctl restart postgresql-15
```

If dual-stack (IPv4/IPv6) is required, use these settings:
```
# In postgresql.conf
listen_addresses = '*'   # listens on all IPv4 and IPv6 addresses

# In pg_hba.conf
# IPv4 connections
host    cryoprotect     cryoprotect     127.0.0.1/32            md5
host    cryoprotect     cryoprotect     0.0.0.0/0               md5
# IPv6 connections
host    cryoprotect     cryoprotect     ::1/128                 md5
host    cryoprotect     cryoprotect     ::/0                    md5
```

### User Setup
```bash
# Create database user
sudo -u postgres psql -c "CREATE USER cryoprotect WITH PASSWORD 'your_secure_password';"
sudo -u postgres psql -c "CREATE DATABASE cryoprotect OWNER cryoprotect;"
sudo -u postgres psql -c "GRANT ALL PRIVILEGES ON DATABASE cryoprotect TO cryoprotect;"
```

### Connection Configuration
```bash
# Edit pg_hba.conf to allow password authentication
sudo vi /var/lib/pgsql/15/data/pg_hba.conf

# Add the following line for local access
# TYPE  DATABASE        USER            ADDRESS                 METHOD
local   cryoprotect     cryoprotect                            md5
host    cryoprotect     cryoprotect     127.0.0.1/32            md5
host    cryoprotect     cryoprotect     ::1/128                 md5

# Restart PostgreSQL
sudo systemctl restart postgresql-15
```

## Data Source Integration

The database integrates data from multiple sources:

1. **PubChem**: Public chemical database
   - Direct API integration via `ResilientPubChemClient`
   - Regular import of new cryoprotectant candidates
   - Property extraction and normalization

2. **ChEMBL**: Bioactivity database
   - API integration via `ResilientChEMBLClient`
   - Import of relevant bioactivity data
   - Reference compound integration

3. **Manual Entry**: User-contributed data
   - Web interface for data entry
   - Bulk import capabilities
   - Validation workflows

## Backup and Recovery

Database backup strategy includes:

1. **Scheduled Backups**:
   - Daily backups (retention: 7 days)
   - Weekly backups (retention: 4 weeks)
   - Monthly backups (retention: 6 months)
   - Yearly backups (retention: 2 years)

2. **Backup Verification**:
   - Checksums for integrity verification
   - Test restoration process
   - Automated monitoring

3. **Recovery Procedures**:
   - Point-in-time recovery capability
   - Restoration validation

### Fedora-Specific Backup Configuration

```bash
# Create backup directory with appropriate permissions
sudo mkdir -p /var/backups/cryoprotect
sudo chown postgres:postgres /var/backups/cryoprotect
sudo chmod 750 /var/backups/cryoprotect

# Set up SELinux context for backup directory
sudo semanage fcontext -a -t postgresql_db_t "/var/backups/cryoprotect(/.*)?)"
sudo restorecon -Rv /var/backups/cryoprotect

# Create a systemd timer for automated backups
sudo vi /etc/systemd/system/cryoprotect-backup.service
# [Unit]
# Description=CryoProtect Database Backup
# 
# [Service]
# Type=oneshot
# User=postgres
# ExecStart=/usr/bin/pg_dump -Fc cryoprotect -f /var/backups/cryoprotect/cryoprotect_$(date +\%Y\%m\%d_\%H\%M\%S).dump
# 
# [Install]
# WantedBy=multi-user.target

sudo vi /etc/systemd/system/cryoprotect-backup.timer
# [Unit]
# Description=Run CryoProtect backup daily
# 
# [Timer]
# OnCalendar=*-*-* 02:00:00
# Persistent=true
# 
# [Install]
# WantedBy=timers.target

sudo systemctl daemon-reload
sudo systemctl enable cryoprotect-backup.timer
sudo systemctl start cryoprotect-backup.timer
```

## Implementation Priorities

For simplified implementation, we'll focus on these database priorities:

1. **Core Table Structure**: Start with molecules, properties, and mixtures
2. **Basic RLS Policies**: Implement simple access control first
3. **Essential Indexes**: Create indexes for most common queries
4. **Direct Connection**: Set up reliable connection pooling
5. **Minimal Backup Strategy**: Implement basic backup/restore

## Database Evolution Strategy

1. **Schema Migration**:
   - Versioned migration scripts
   - Forward-only migrations
   - Backward compatibility with views

2. **Data Migration**:
   - Batch processing for large datasets
   - Checkpointing for resumable operations
   - Validation after migration

3. **RLS Evolution**:
   - Automated policy testing
   - Access verification reports

## Known Limitations and Gaps

1. **SELinux Integration**: Need to complete SELinux policy setup
2. **Podman DB Container Support**: Need to adapt Docker Compose for Podman-compatible PostgreSQL containers
3. **Systemd Service Configuration**: Need to complete service configuration for automated processes
4. **Fedora Path Adjustments**: Update backup paths to match Fedora filesystem layout

## Future Database Enhancements

1. **Advanced Analytics**: Materialized views for analytical queries
2. **Full-Text Search**: Specialized indexing for text search capabilities
3. **Molecular Similarity Search**: Chemical fingerprint indexing
4. **Partitioned Tables**: Performance optimization for large tables