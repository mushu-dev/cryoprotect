# CryoProtect Database Population Guide Using MCP Tools

This guide provides instructions for Roo agents to implement the simplified database population approach using Supabase MCP tools. Follow this step-by-step process to ensure a reliable, observable data pipeline.

## 1. Prerequisites

- Familiarity with the core philosophy: **Simplify, Verify, Then Optimize**
- Understanding of the project's verification criteria
- Access to Supabase MCP tools and SQL execution capabilities
- Knowledge of CryoProtect database schema

## 2. Core Workflow Overview

1. Verify database connection and schema
2. Populate reference compounds (9 total)
3. Verify reference compound properties
4. Import PubChem data incrementally
5. Import ChEMBL data incrementally
6. Verify all molecule data meets criteria
7. Optimize performance with indexes

## 3. Using MCP Tools for Database Operations

### 3.1 Supabase Project Access

```
mcp__supabase__list_projects
```

This will give you the project IDs. Select the CryoProtect project ID for subsequent operations.

### 3.2 Executing SQL with MCP

Always prefer `mcp__supabase__execute_sql` over direct database connections for better reliability:

```python
# Example MCP SQL execution
result = mcp__supabase__execute_sql(
    project_id="your_project_id",
    query="SELECT COUNT(*) FROM molecules"
)
```

### 3.3 Schema Verification

Use this SQL query with MCP to verify the required tables exist:

```sql
SELECT 
    table_name,
    EXISTS(SELECT 1 FROM information_schema.tables WHERE table_name = 'molecules') AS molecules_exists,
    EXISTS(SELECT 1 FROM information_schema.tables WHERE table_name = 'molecular_properties') AS properties_exists,
    EXISTS(SELECT 1 FROM information_schema.tables WHERE table_name = 'property_types') AS property_types_exists
FROM 
    information_schema.tables
WHERE 
    table_name IN ('molecules', 'molecular_properties', 'property_types')
GROUP BY 
    table_name;
```

## 4. Step-by-Step Implementation

### 4.1 Reference Compounds Population

1. Use this SQL to insert reference compounds (one at a time for simplicity):

```sql
-- Insert a reference compound
INSERT INTO molecules 
    (name, chembl_id, smiles, inchi, inchikey, formula, data_source)
VALUES 
    ('Glycerol', 'CHEMBL388978', 'C(C(CO)O)O', 
     'InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2', 
     'PEDCQBHIVMGVHV-UHFFFAOYSA-N', 'C3H8O3', 'reference')
ON CONFLICT (chembl_id) DO UPDATE SET
    name = EXCLUDED.name,
    smiles = EXCLUDED.smiles,
    inchi = EXCLUDED.inchi,
    inchikey = EXCLUDED.inchikey,
    formula = EXCLUDED.formula,
    updated_at = NOW()
RETURNING id;
```

2. For each molecule, insert the required properties:

```sql
-- First create/get the property type
INSERT INTO property_types (name, data_type) 
VALUES ('logP', 'numeric')
ON CONFLICT (name) DO UPDATE SET updated_at = NOW()
RETURNING id;

-- Then insert the property value
INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
VALUES 
    ('molecule_id_from_previous_query', 'property_type_id_from_previous_query', -1.76)
ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
    numeric_value = EXCLUDED.numeric_value,
    updated_at = NOW();
```

### 4.2 Repeat for All Reference Compounds

The reference compounds with their complete properties are:

| Name | ChEMBL ID | logP | H-Bond Donors | H-Bond Acceptors |
|------|-----------|------|--------------|------------------|
| Glycerol | CHEMBL388978 | -1.76 | 3 | 3 |
| DMSO | CHEMBL1098659 | -1.35 | 0 | 1 |
| beta-Alanine | CHEMBL66195 | -3.17 | 2 | 3 |
| tert-Butanol | CHEMBL500033 | 0.35 | 1 | 1 |
| Urea | CHEMBL1487 | -2.11 | 2 | 1 |
| Ethylene glycol | CHEMBL6196 | -1.36 | 2 | 2 |
| Propylene glycol | CHEMBL967 | -0.92 | 2 | 2 |
| Trehalose | CHEMBL262548 | -3.77 | 8 | 11 |
| Glycine | CHEMBL6752 | -3.21 | 1 | 3 |

### 4.3 Verify Reference Compounds

Use this SQL query with MCP to verify reference compounds have all required properties:

```sql
SELECT 
    m.name,
    m.chembl_id,
    COUNT(DISTINCT pt.name) AS property_count,
    BOOL_OR(pt.name = 'logP') AS has_logp,
    BOOL_OR(pt.name = 'h_bond_donors') AS has_hbd,
    BOOL_OR(pt.name = 'h_bond_acceptors') AS has_hba
FROM 
    molecules m
LEFT JOIN 
    molecular_properties mp ON m.id = mp.molecule_id
LEFT JOIN 
    property_types pt ON mp.property_type_id = pt.id
WHERE 
    m.chembl_id IN (
        'CHEMBL388978', 'CHEMBL1098659', 'CHEMBL66195', 'CHEMBL500033',
        'CHEMBL1487', 'CHEMBL6196', 'CHEMBL967', 'CHEMBL262548', 'CHEMBL6752'
    )
GROUP BY 
    m.id, m.name, m.chembl_id;
```

### 4.4 Import PubChem Data

Use the PubChem API for small, focused batches:

1. **Fetch data from PubChem API**:
   ```python
   # Example PubChem API request
   url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,IsomericSMILES,InChI,InChIKey,Title/JSON"
   response = requests.get(url)
   data = response.json()["PropertyTable"]["Properties"][0]
   ```

2. **Insert molecule data using MCP SQL**:
   ```python
   # Example MCP SQL execution for PubChem molecule insertion
   sql = """
   INSERT INTO molecules 
       (pubchem_cid, name, smiles, inchi, inchikey, formula, molecular_weight, data_source)
   VALUES 
       (%s, %s, %s, %s, %s, %s, %s, 'PubChem')
   ON CONFLICT (pubchem_cid) DO UPDATE SET
       name = EXCLUDED.name,
       smiles = EXCLUDED.smiles,
       inchi = EXCLUDED.inchi,
       inchikey = EXCLUDED.inchikey,
       formula = EXCLUDED.formula,
       molecular_weight = EXCLUDED.molecular_weight,
       updated_at = NOW()
   RETURNING id;
   """
   
   result = mcp__supabase__execute_sql(
       project_id="your_project_id",
       query=sql,
       params=[data["CID"], data.get("Title", f"CID {data['CID']}"), 
               data["IsomericSMILES"], data["InChI"], data["InChIKey"], 
               data["MolecularFormula"], data.get("MolecularWeight")]
   )
   ```

### 4.5 Property Handling for PubChem Data

Focus on getting the 3 critical properties correctly:

```python
# Example MCP SQL execution for property insertion
properties = {
    "logP": float(data["XLogP"]),
    "h_bond_donors": int(data["HBondDonorCount"]),
    "h_bond_acceptors": int(data["HBondAcceptorCount"])
}

for prop_name, value in properties.items():
    # Get or create property type
    prop_type_sql = """
    INSERT INTO property_types (name, data_type)
    VALUES (%s, 'numeric')
    ON CONFLICT (name) DO UPDATE SET updated_at = NOW()
    RETURNING id;
    """
    
    prop_type_result = mcp__supabase__execute_sql(
        project_id="your_project_id",
        query=prop_type_sql,
        params=[prop_name]
    )
    
    property_type_id = prop_type_result[0]["id"]
    
    # Insert property value
    property_sql = """
    INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
    VALUES (%s, %s, %s)
    ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
        numeric_value = EXCLUDED.numeric_value,
        updated_at = NOW();
    """
    
    mcp__supabase__execute_sql(
        project_id="your_project_id",
        query=property_sql,
        params=[molecule_id, property_type_id, value]
    )
```

### 4.6 ChEMBL Data Integration

Similar to PubChem, but using ChEMBL API:

```python
# Example ChEMBL API request
url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
response = requests.get(url)
data = response.json()

# Extract relevant properties from the data
molecule_props = data.get('molecule_properties', {})
logp = molecule_props.get('alogp')
hbd = molecule_props.get('hbd')
hba = molecule_props.get('hba')
```

### 4.7 Verification Queries

Use these queries to verify progress at different stages:

```sql
-- Verify total molecules and property completeness
SELECT 
    COUNT(*) AS total_molecules,
    COUNT(DISTINCT m.id) FILTER (
        WHERE EXISTS (
            SELECT 1 FROM molecular_properties mp1
            JOIN property_types pt1 ON mp1.property_type_id = pt1.id AND pt1.name = 'logP'
            WHERE mp1.molecule_id = m.id
        ) AND EXISTS (
            SELECT 1 FROM molecular_properties mp2
            JOIN property_types pt2 ON mp2.property_type_id = pt2.id AND pt2.name = 'h_bond_donors'
            WHERE mp2.molecule_id = m.id
        ) AND EXISTS (
            SELECT 1 FROM molecular_properties mp3
            JOIN property_types pt3 ON mp3.property_type_id = pt3.id AND pt3.name = 'h_bond_acceptors'
            WHERE mp3.molecule_id = m.id
        )
    ) AS complete_molecules
FROM 
    molecules m;
```

### 4.8 Performance Optimization

Once data is successfully imported, add indexes:

```sql
-- Add indexes for commonly queried fields
CREATE INDEX IF NOT EXISTS idx_molecules_pubchem_cid ON molecules (pubchem_cid);
CREATE INDEX IF NOT EXISTS idx_molecules_chembl_id ON molecules (chembl_id);
CREATE INDEX IF NOT EXISTS idx_molecules_inchikey ON molecules (inchikey);
CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON molecular_properties (molecule_id);
CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_type_id ON molecular_properties (property_type_id);
CREATE INDEX IF NOT EXISTS idx_property_types_name ON property_types (name);

-- Add compound index for frequent join
CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_property ON molecular_properties (molecule_id, property_type_id);

-- Update statistics
ANALYZE;
```

## 5. Tracking Progress

Create checkpoints to track progress throughout the process:

```python
def save_checkpoint(project_id, position, processed_count, success_count):
    """Save checkpoint data to the database."""
    checkpoint_sql = """
    INSERT INTO import_checkpoints (import_type, position, processed_count, success_count, timestamp)
    VALUES ('pubchem', %s, %s, %s, NOW())
    RETURNING id;
    """
    
    mcp__supabase__execute_sql(
        project_id=project_id,
        query=checkpoint_sql,
        params=[position, processed_count, success_count]
    )
```

## 6. Success Criteria Verification

Before concluding, verify that all requirements are met:

1. **Molecule Count**: At least 5,000 molecules in database
2. **Reference Compounds**: All 9 reference compounds with complete properties
3. **Property Completeness**: At least 90% of molecules have all required properties
4. **Query Performance**: Average query time under 50ms

Use this final verification query:

```sql
WITH 
molecule_count AS (
    SELECT COUNT(*) AS count FROM molecules
),
ref_compound_check AS (
    SELECT 
        COUNT(DISTINCT m.id) = 9 AS all_refs_present,
        COUNT(DISTINCT m.id) FILTER (
            WHERE EXISTS (
                SELECT 1 FROM molecular_properties mp1
                JOIN property_types pt1 ON mp1.property_type_id = pt1.id AND pt1.name = 'logP'
                WHERE mp1.molecule_id = m.id
            ) AND EXISTS (
                SELECT 1 FROM molecular_properties mp2
                JOIN property_types pt2 ON mp2.property_type_id = pt2.id AND pt2.name = 'h_bond_donors'
                WHERE mp2.molecule_id = m.id
            ) AND EXISTS (
                SELECT 1 FROM molecular_properties mp3
                JOIN property_types pt3 ON mp3.property_type_id = pt3.id AND pt3.name = 'h_bond_acceptors'
                WHERE mp3.molecule_id = m.id
            )
        ) = 9 AS all_refs_complete
    FROM 
        molecules m
    WHERE 
        m.chembl_id IN (
            'CHEMBL388978', 'CHEMBL1098659', 'CHEMBL66195', 'CHEMBL500033',
            'CHEMBL1487', 'CHEMBL6196', 'CHEMBL967', 'CHEMBL262548', 'CHEMBL6752'
        )
),
property_completeness AS (
    SELECT 
        COUNT(DISTINCT m.id) AS complete_molecules,
        (SELECT count FROM molecule_count) AS total_molecules,
        (COUNT(DISTINCT m.id)::float / (SELECT count FROM molecule_count)) * 100 AS completeness_pct
    FROM 
        molecules m
    WHERE 
        EXISTS (
            SELECT 1 FROM molecular_properties mp1
            JOIN property_types pt1 ON mp1.property_type_id = pt1.id AND pt1.name = 'logP'
            WHERE mp1.molecule_id = m.id
        ) AND EXISTS (
            SELECT 1 FROM molecular_properties mp2
            JOIN property_types pt2 ON mp2.property_type_id = pt2.id AND pt2.name = 'h_bond_donors'
            WHERE mp2.molecule_id = m.id
        ) AND EXISTS (
            SELECT 1 FROM molecular_properties mp3
            JOIN property_types pt3 ON mp3.property_type_id = pt3.id AND pt3.name = 'h_bond_acceptors'
            WHERE mp3.molecule_id = m.id
        )
)
SELECT 
    (SELECT count FROM molecule_count) AS total_molecules,
    (SELECT count >= 5000 FROM molecule_count) AS molecule_count_met,
    (SELECT all_refs_present FROM ref_compound_check) AS all_refs_present,
    (SELECT all_refs_complete FROM ref_compound_check) AS all_refs_complete,
    (SELECT completeness_pct FROM property_completeness) AS property_completeness_pct,
    (SELECT completeness_pct >= 90 FROM property_completeness) AS property_completeness_met;
```

## 7. Implementation Phases

### Phase 1: Reference Compounds
1. Insert the 9 reference compounds
2. Create property types for logP, h_bond_donors, h_bond_acceptors
3. Insert all properties for each compound
4. Verify all 9 compounds have complete properties

### Phase 2: PubChem Data (Batch Processing)
1. Process in small batches (10-20 compounds)
2. Focus on compounds with complete property data
3. Create checkpoints after each batch
4. Incrementally verify progress

### Phase 3: ChEMBL Data (Focused Approach)
1. Target specific cryoprotectant classes
2. Use the reference compounds as starting points
3. Import only compounds with complete property data
4. Verify ChEMBL data quality

### Phase 4: Performance Optimization
1. Only after data population is successful
2. Add targeted indexes for specific query patterns
3. Run ANALYZE to update statistics
4. Verify query performance improvements