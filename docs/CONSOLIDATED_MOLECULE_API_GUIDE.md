# Consolidated Molecule API Guide

## Overview

This guide explains how to work with the CryoProtect API when dealing with consolidated molecules and differentiation groups. These concepts were introduced in Phase 2 of the Molecule Data Quality Enhancement Plan.

## Consolidated Molecule Concepts

### What are Consolidated Molecules?

In the CryoProtect database, "consolidated molecules" are instances where duplicate molecule records have been identified and merged:

- **Primary Molecule**: The main molecule record that remains active
- **Secondary Molecule (Consolidated)**: Duplicate molecules that point to the primary molecule

When you request a consolidated (secondary) molecule through the API, the system will automatically redirect to the primary molecule.

### What are Differentiation Groups?

Differentiation groups are sets of molecules that have similar names but different chemical structures:

- Molecules in the same differentiation group are related but distinct chemical entities
- Each molecule in a group has a differentiation description explaining how it differs from others
- Unlike consolidated molecules, differentiated molecules remain separate records

## API Endpoints

### Consolidated Molecule Endpoints

#### Get a Molecule with Consolidated Handling

```
GET /api/v1/consolidated/molecules/{molecule_id}
```

This endpoint retrieves a molecule with consolidated molecule handling:
- If the requested molecule is a secondary (consolidated) molecule, it returns the primary molecule
- The response includes information about the consolidation
- Example response field: `is_consolidated`, `consolidated_to`, `consolidated_molecules`

#### Batch Operations with Consolidated Molecules

```
POST /api/v1/consolidated/batch
```

This endpoint handles batch operations with consolidated molecules:
- Any consolidated molecules in the request are automatically redirected to their primaries
- The response includes information about which molecules were redirected
- Example request body:
  ```json
  {
    "molecule_ids": ["123e4567-e89b-12d3-a456-426614174000", "123e4567-e89b-12d3-a456-426614174001"]
  }
  ```

#### Get Primary Molecule

```
GET /api/v1/molecules/{molecule_id}/primary
```

This endpoint retrieves the primary molecule for a given molecule ID:
- If the molecule is already a primary, it returns itself
- If the molecule is consolidated, it returns its primary
- Includes detailed consolidation information

#### List All Consolidated Molecules

```
GET /api/v1/consolidated
```

This endpoint returns all consolidated molecule relationships:
- Lists all primary molecules and their secondary (consolidated) molecules
- Includes count information

### Differentiation Group Endpoints

#### List All Differentiation Groups

```
GET /api/v1/differentiation/groups
```

This endpoint returns all differentiation groups in the system.

#### Get Differentiation Group Details

```
GET /api/v1/differentiation/groups/{group_id}
```

This endpoint retrieves detailed information about a differentiation group:
- Includes all molecules in the group and their differences
- Provides access to the full molecule data for each member

#### Get Molecule Differentiation Information

```
GET /api/v1/molecules/{molecule_id}/differentiation
```

This endpoint returns differentiation information for a specific molecule:
- Includes the differentiation group it belongs to
- Lists other molecules in the same group
- Provides the differentiation description explaining how it differs from others

## Using the API with Consolidated Molecules

### Regular vs. Consolidated-Aware Endpoints

All existing API endpoints continue to work as before, but may not provide consolidated molecule information. The new consolidated-aware endpoints provide additional information about molecule relationships.

### Migration Guide

When working with molecules, consider these best practices:

1. **Check for Consolidation**: Always check if a molecule is consolidated using the `is_consolidated` field
2. **Use Primary Molecules**: When storing references to molecules, use the primary molecule ID
3. **Handle Redirections**: Be aware that requests for consolidated molecules will redirect to their primaries
4. **Include Consolidated Information**: When displaying molecule information, include consolidation details to help users understand relationships

### Example: Fetching a Molecule

**Basic Approach**:
```javascript
// Fetch molecule with ID
fetch('/api/v1/molecules/123e4567-e89b-12d3-a456-426614174000')
  .then(response => response.json())
  .then(data => {
    // Process molecule data
    console.log(data);
  });
```

**Consolidated-Aware Approach**:
```javascript
// Fetch molecule with consolidated handling
fetch('/api/v1/consolidated/molecules/123e4567-e89b-12d3-a456-426614174000')
  .then(response => response.json())
  .then(data => {
    const molecule = data.data;
    
    // Check if this is a consolidated molecule
    if (molecule.is_consolidated) {
      console.log(`This molecule is consolidated into: ${molecule.consolidated_to}`);
      // You may want to notify the user about the redirection
    }
    
    // Check if this is a primary molecule with consolidated molecules
    if (molecule.consolidated_molecules && molecule.consolidated_molecules.length > 0) {
      console.log(`This molecule has ${molecule.consolidated_molecules.length} consolidated molecules`);
    }
    
    // Check if this is part of a differentiation group
    if (molecule.differentiation_group) {
      console.log(`This molecule is part of differentiation group: ${molecule.differentiation_group}`);
      console.log(`Differentiation: ${molecule.differentiation_description}`);
    }
    
    // Process molecule data
    console.log(molecule);
  });
```

### Example: Batch Operations

**Basic Approach**:
```javascript
// Fetch multiple molecules
fetch('/api/v1/batch', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({ molecule_ids: ['id1', 'id2', 'id3'] })
})
  .then(response => response.json())
  .then(data => {
    // Process molecule data
    console.log(data);
  });
```

**Consolidated-Aware Approach**:
```javascript
// Fetch multiple molecules with consolidated handling
fetch('/api/v1/consolidated/batch', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({ molecule_ids: ['id1', 'id2', 'id3'] })
})
  .then(response => response.json())
  .then(data => {
    const molecules = data.data.molecules;
    
    // Check for consolidated redirections
    if (data.data.meta && data.data.meta.consolidated_redirections) {
      const redirections = data.data.meta.consolidated_redirections;
      console.log('The following molecules were redirected:');
      for (const [original, primary] of Object.entries(redirections)) {
        console.log(`${original} → ${primary}`);
      }
    }
    
    // Process molecule data
    console.log(molecules);
  });
```

## Response Format

All API responses follow the standardized format:

```json
{
  "status": "success",
  "timestamp": "2023-05-13T12:34:56.789Z",
  "code": 200,
  "message": "Request succeeded",
  "data": {
    // Response data
  }
}
```

### Consolidated Molecule Fields

Molecule responses include these additional fields when using consolidated-aware endpoints:

```json
{
  "id": "123e4567-e89b-12d3-a456-426614174000",
  "name": "Dimethyl sulfoxide",
  "smiles": "CS(=O)C",
  "is_consolidated": false,
  "consolidated_to": null,
  "consolidated_molecules": ["123e4567-e89b-12d3-a456-426614174001"],
  "differentiation_group": null,
  "differentiation_description": null
}
```

### Differentiation Fields

Molecules in differentiation groups include these additional fields:

```json
{
  "id": "123e4567-e89b-12d3-a456-426614174000",
  "name": "Glycerol",
  "smiles": "C(C(CO)O)O",
  "is_consolidated": false,
  "consolidated_to": null,
  "consolidated_molecules": [],
  "differentiation_group": "glycerol-group",
  "differentiation_description": "Anhydrous form"
}
```

## Best Practices

1. **Always Use Primary Molecules**: Store and reference primary molecule IDs, not consolidated (secondary) molecules.

2. **Handle Redirections Gracefully**: When displaying molecule data, handle redirections from consolidated molecules to primaries in a user-friendly way.

3. **Include Differentiation Information**: When working with molecules that are part of differentiation groups, include the differentiation description to help users understand the differences.

4. **Check Consolidation Status**: Before performing operations that modify molecule data, check if the molecule is consolidated and use the primary molecule instead.

5. **Use Batch Operations Efficiently**: When working with multiple molecules, use the batch operations endpoints to handle consolidation efficiently.

## Conclusion

The consolidated molecule API enhancement provides a more consistent and accurate representation of molecules in the CryoProtect database. By properly handling consolidated molecules and differentiation groups, you can provide a better user experience and ensure data consistency in your applications.