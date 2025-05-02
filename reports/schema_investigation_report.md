# Database Schema Investigation Report

## Task Overview
This report documents the investigation of the actual schema of the 'molecules' table in the CryoProtect database, with specific focus on the existence and structure of a column intended to store JSON properties. This investigation was prompted by the failure of the verification task (task-dbpop-06), which indicated a schema mismatch preventing data insertion/updates.

## Database Schema Analysis

### Molecules Table Structure
The 'molecules' table has the following structure:

| Column Name | Data Type | Nullable |
|-------------|-----------|----------|
| id | uuid | NO |
| name | character varying | NO |
| smiles | character varying | YES |
| inchi | text | YES |
| inchikey | character varying | YES |
| formula | character varying | YES |
| molecular_weight | numeric | YES |
| created_by | uuid | YES |
| is_public | boolean | YES |
| data_source | character varying | YES |
| version | integer | YES |
| modification_history | jsonb | YES |
| created_at | timestamp with time zone | NO |
| updated_at | timestamp with time zone | NO |
| pubchem_cid | integer | YES |
| cid | integer | YES |
| pubchem_link | text | YES |
| molecular_formula | text | YES |
| chembl_id | character varying | YES |

**Key Finding**: The 'molecules' table does NOT have a 'properties' column of type JSONB as expected in the directive file. The only JSONB column is 'modification_history', which serves a different purpose.

### Property Storage Architecture
Instead of storing properties directly in the 'molecules' table, the database uses a normalized approach with separate tables:

1. **property_types**: Defines the types of properties that can be stored
   - Contains definitions for properties like LogP, H-bond donors/acceptors, etc.

2. **molecular_properties**: Stores the actual property values with foreign key relationships
   - Structure includes molecule_id, property_type_id, and various value columns (numeric_value, text_value, boolean_value)

This architecture follows a more normalized database design pattern compared to the denormalized approach (using a JSONB column) that was expected by the import scripts.

## Current Data State

### Molecule Counts
- Total molecules in database: 723
- Molecules with PubChem CIDs: 693
- Molecules with ChEMBL IDs: 10
- Molecules with both PubChem CID and ChEMBL ID (cross-referenced): 0

### Property Data
- Total property records in molecular_properties table: 72
- Molecules with LogP properties: 46
- Molecules with Hydrogen Bond Donor Count: 1
- Molecules with Hydrogen Bond Acceptor Count: 1
- Molecules with PubChem CIDs that have properties: 0

### Reference Compounds
None of the reference compounds specified in the directive (CHEMBL388978, CHEMBL1098659, CHEMBL66195, CHEMBL500033, CHEMBL1487, CHEMBL6196, CHEMBL967, CHEMBL262548, CHEMBL6752) are present in the database.

## Identified Issues

1. **Schema Mismatch**: The expected 'properties' JSONB column in the 'molecules' table doesn't exist. Instead, properties are stored in a separate 'molecular_properties' table.

2. **Missing Reference Compounds**: The reference compounds specified in the directive are not present in the database.

3. **Incomplete ChEMBL Data Import**: Only 10 molecules with ChEMBL IDs are present, far below the expected minimum of 500.

4. **Missing PubChem Property Enhancement**: None of the molecules with PubChem CIDs have properties in the molecular_properties table.

5. **Missing Cross-References**: No molecules have both ChEMBL ID and PubChem CID, indicating that cross-referencing was not completed.

## Conclusion

The verification task (task-dbpop-06) failed primarily due to a fundamental schema mismatch. The import scripts were written assuming a 'properties' JSONB column in the 'molecules' table, but the actual database schema uses a separate 'molecular_properties' table with a different structure.

To resolve this issue, the import scripts need to be modified to:
1. Recognize the actual database schema
2. Insert property data into the 'molecular_properties' table instead of trying to update a non-existent 'properties' column
3. Establish proper relationships between molecules and their properties using the foreign key structure

This schema difference explains all the verification failures reported in task-dbpop-06, including the missing reference compounds, incomplete ChEMBL data, missing property data for PubChem molecules, and lack of cross-references.