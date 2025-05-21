# Consolidated Molecules Schema

This document describes the schema design for the consolidated molecules feature in CryoProtect.

## Table Structure

The `consolidated_molecules` table has the following structure:

```sql
CREATE TABLE public.consolidated_molecules (
    id UUID PRIMARY KEY REFERENCES public.molecules(id),
    name TEXT NOT NULL,
    smiles TEXT,
    inchi TEXT,
    inchikey TEXT,
    formula TEXT,
    molecular_weight DOUBLE PRECISION,
    pubchem_cid TEXT,
    is_public BOOLEAN DEFAULT TRUE,
    data_source TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    is_consolidated BOOLEAN DEFAULT FALSE,
    primary_molecule_id UUID REFERENCES public.molecules(id),
    primary_molecule_name TEXT,
    primary_pubchem_cid TEXT,
    molecule_status TEXT DEFAULT 'original',
    CONSTRAINT valid_status CHECK (molecule_status IN ('original', 'primary', 'duplicate'))
);
```

## Entity Relationship Diagram

```
┌─────────────────────┐               ┌─────────────────────┐
│                     │               │                     │
│      molecules      │               │  scientific_data_   │
│                     │               │       audit         │
└─────────┬───────────┘               └─────────────────────┘
          │                                     ▲
          │                                     │
          │ 1                                   │ records changes
          ▼                                     │
┌─────────────────────┐               ┌─────────────────────┐
│                     │ redirects to  │                     │
│   consolidated_     │◄──────────────┤   dependent tables  │
│     molecules       │               │                     │
└─────────────────────┘               └─────────────────────┘
```

## Field Descriptions

### Core Fields (copied from molecules table)
- `id`: Primary key, identical to the molecule ID in the molecules table
- `name`: Molecule name
- `smiles`: SMILES notation
- `inchi`: InChI identifier
- `inchikey`: InChIKey (used for deduplication)
- `formula`: Molecular formula
- `molecular_weight`: Molecular weight
- `pubchem_cid`: PubChem Compound ID
- `is_public`: Whether the molecule is publicly available
- `data_source`: Source of the molecule data
- `created_at`: When the record was created
- `updated_at`: When the record was last updated

### Consolidation-specific Fields
- `is_consolidated`: Whether this molecule has been involved in consolidation
- `primary_molecule_id`: For duplicates, references the primary molecule ID
- `primary_molecule_name`: Name of the primary molecule (for easier reporting)
- `primary_pubchem_cid`: PubChem CID of the primary molecule (if available)
- `molecule_status`: One of:
  - `'original'`: Not yet involved in consolidation
  - `'primary'`: This is the canonical molecule for an InChIKey
  - `'duplicate'`: This is a duplicate of a primary molecule

## Indexes

```sql
CREATE INDEX idx_consolidated_molecules_inchikey ON public.consolidated_molecules(inchikey);
CREATE INDEX idx_consolidated_molecules_primary_id ON public.consolidated_molecules(primary_molecule_id);
CREATE INDEX idx_consolidated_molecules_status ON public.consolidated_molecules(molecule_status);
```

## Audit Trail

The `scientific_data_audit` table records all consolidation operations:

```sql
CREATE TABLE public.scientific_data_audit (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    table_name TEXT NOT NULL,
    record_id TEXT NOT NULL,
    operation TEXT NOT NULL,
    old_value JSONB,
    new_value JSONB,
    user_id UUID,
    timestamp TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);
```

Each consolidation operation creates an audit entry with:
- `table_name`: 'molecules'
- `record_id`: The duplicate molecule ID
- `operation`: 'consolidated'
- `old_value`: JSON with the duplicate molecule details
- `new_value`: JSON with the primary molecule details
- `user_id`: The user performing the consolidation
- `timestamp`: When the consolidation occurred

## Query Examples

### Find all primary molecules

```sql
SELECT * FROM consolidated_molecules WHERE molecule_status = 'primary';
```

### Find all duplicates of a specific primary molecule

```sql
SELECT * FROM consolidated_molecules 
WHERE primary_molecule_id = '00000000-0000-0000-0000-000000000000'
AND molecule_status = 'duplicate';
```

### Count duplicates per primary molecule

```sql
SELECT primary_molecule_id, primary_molecule_name, COUNT(*) as duplicate_count
FROM consolidated_molecules
WHERE molecule_status = 'duplicate'
GROUP BY primary_molecule_id, primary_molecule_name
ORDER BY duplicate_count DESC;
```