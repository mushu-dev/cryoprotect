# ChEMBL Molecules Schema Inspection Report

**Date:** 2025-04-26  
**Task:** task-imp-chembl-schema-001  
**Author:** Apex Implementer  
**References:**
- `.specs/chembl_molecules_schema_remediation.md`
- `migrations/001_initial_schema.sql`

## 1. Schema Inspection Results

### 1.1 Canonical Schema (from migrations/001_initial_schema.sql)

```sql
CREATE TABLE public.molecules (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    cid INTEGER UNIQUE NOT NULL,  -- PubChem Compound ID
    name TEXT,
    molecular_formula TEXT,
    smiles TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id),
    pubchem_link TEXT GENERATED ALWAYS AS ('https://pubchem.ncbi.nlm.nih.gov/compound/' || cid) STORED
);
```

### 1.2 Current Database Schema

```sql
CREATE TABLE public.molecules (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name CHARACTER VARYING NOT NULL,
    smiles CHARACTER VARYING,
    inchi TEXT,
    inchikey CHARACTER VARYING,
    formula CHARACTER VARYING,
    molecular_weight NUMERIC,
    created_by UUID,
    is_public BOOLEAN DEFAULT false,
    data_source CHARACTER VARYING,
    version INTEGER,
    modification_history JSONB,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);
```

## 2. Schema Drift Analysis

### 2.1 Missing Columns

The following canonical columns are **missing** from the current schema:

| Column Name | Data Type | Constraints | Purpose |
|-------------|-----------|-------------|---------|
| cid | INTEGER | UNIQUE NOT NULL | PubChem Compound ID |
| molecular_formula | TEXT | | Chemical formula |
| pubchem_link | TEXT | GENERATED | Link to PubChem compound page |

### 2.2 Additional Columns

The current schema has the following **additional** columns not in the canonical schema:

| Column Name | Data Type | Purpose |
|-------------|-----------|---------|
| inchi | TEXT | International Chemical Identifier |
| inchikey | CHARACTER VARYING | InChI Key for the molecule |
| formula | CHARACTER VARYING | Chemical formula (similar to molecular_formula) |
| molecular_weight | NUMERIC | Molecular weight |
| is_public | BOOLEAN | Flag for public visibility |
| data_source | CHARACTER VARYING | Source of the molecule data |
| version | INTEGER | Version tracking |
| modification_history | JSONB | History of changes |

### 2.3 Type Differences

| Column Name | Canonical Type | Current Type | Compatible? |
|-------------|----------------|--------------|------------|
| name | TEXT | CHARACTER VARYING NOT NULL | Yes, but with added NOT NULL constraint |
| smiles | TEXT | CHARACTER VARYING | Yes |

### 2.4 Missing Indexes

The canonical schema defines an index on the `cid` column which is missing:

```sql
CREATE INDEX idx_molecules_cid ON public.molecules(cid);
```

## 3. Code Analysis

### 3.1 `cid` vs `pubchem_cid` Usage

- **Canonical Schema:** Uses `cid` for PubChem Compound ID
- **Code References:** Multiple Python files reference `pubchem_cid` instead of `cid`
  - Found in: test_migration_script.py, populate_database_supabase.py, chemical_data/toxicity modules
- **SQL References:** The original migration file (001_initial_schema.sql) uses `cid`

### 3.2 No Alternative Storage

- No other tables or columns in the database appear to be storing PubChem IDs
- No molecule_identifier_mapping table exists that might be mapping molecules to PubChem IDs

## 4. Conclusion

There is significant schema drift between the canonical schema and the current implementation:

1. **Critical Missing Columns:**
   - The `cid` column (or `pubchem_cid`) is completely missing, which is essential for PubChem integration
   - The generated `pubchem_link` column is missing

2. **Code-Schema Mismatch:**
   - The code expects a `pubchem_cid` column
   - The canonical schema defines a `cid` column
   - The actual database has neither column

3. **Schema Evolution:**
   - The current schema has evolved to include additional molecular identifiers (inchi, inchikey)
   - It also includes metadata fields not in the canonical schema

## 5. Recommended Remediation Approach

Based on the analysis, the recommended approach is:

1. **Add `pubchem_cid` Column:**
   - Add a `pubchem_cid` INTEGER column to align with code expectations
   - Make it UNIQUE and allow NULL initially to support gradual data population

2. **Add `cid` as Generated Column:**
   - Add `cid` as a generated column based on `pubchem_cid` to maintain canonical schema compatibility
   - This approach allows both column names to be used in code

3. **Add `pubchem_link` Column:**
   - Add the generated `pubchem_link` column referencing `pubchem_cid`

4. **Create Required Indexes:**
   - Create indexes on both `cid` and `pubchem_cid` columns

This approach balances the need to align with the canonical schema while minimizing code changes and supporting the evolved schema.