# API Standardization Implementation

## Overview

This implementation adds standardized API support for the data quality improvements made in Phase 2, specifically handling consolidated molecules and differentiation groups. The implementation ensures that all API endpoints follow standardized patterns for response formats, error handling, and HTTP status codes.

## Components Implemented

### 1. Consolidated Molecule Utilities

The `consolidated_utils.py` module provides the foundation for working with consolidated molecules:

- `is_consolidated()`: Checks if a molecule is consolidated (secondary)
- `get_primary_molecule()`: Gets the primary molecule for a consolidated molecule
- `get_consolidated_molecules()`: Gets all consolidated molecules for a primary molecule
- `get_differentiation_group()`: Gets the differentiation group for a molecule
- `get_differentiation_group_members()`: Gets all molecules in a differentiation group
- `enrich_molecule_data()`: Adds consolidated and differentiation information to molecule data

### 2. API Decorators

The `consolidated_decorators.py` module provides decorators for handling consolidated molecules in API endpoints:

- `handle_consolidated_molecules()`: Decorator for single molecule endpoints
- `handle_batch_consolidated_molecules()`: Decorator for batch operations with molecules

These decorators automatically:
- Check if requested molecules are consolidated
- Redirect to primary molecules when needed
- Add consolidation information to responses
- Maintain backward compatibility with existing clients

### 3. API Resources

New API resources have been implemented for working with consolidated molecules:

- `ConsolidatedMoleculeResource`: Handles consolidated molecule retrieval with automatic redirection
- `ConsolidatedMoleculeBatchResource`: Handles batch operations with consolidated molecules
- `PrimaryMoleculeResource`: Gets the primary molecule for a molecule
- `ConsolidatedMoleculesListResource`: Lists all consolidated molecule relationships

And for differentiation groups:

- `DifferentiationGroupListResource`: Lists all differentiation groups
- `DifferentiationGroupResource`: Gets information about a specific differentiation group
- `MoleculeDifferentiationResource`: Gets differentiation information for a molecule

### 4. Integration with API Framework

The implementation integrates with the existing API framework:

- Registering new resources with the API
- Adding documentation for new endpoints
- Ensuring consistent response formats
- Following standardized error handling patterns

### 5. Updated Models

The molecule model has been updated to include consolidated molecule information:

- `is_consolidated`: Boolean indicating if a molecule is consolidated
- `consolidated_to`: Reference to the primary molecule (if consolidated)
- `consolidated_molecules`: List of consolidated molecules (for primary molecules)
- `differentiation_group`: Differentiation group ID (if differentiated)
- `differentiation_description`: Description of the molecule's differentiation

## New API Endpoints

The implementation adds the following new API endpoints:

- `/api/v1/consolidated/molecules/{molecule_id}`: Get a molecule with consolidated handling
- `/api/v1/consolidated/batch`: Batch operations with consolidated molecules
- `/api/v1/molecules/{molecule_id}/primary`: Get the primary molecule for a molecule
- `/api/v1/consolidated`: List all consolidated molecule relationships
- `/api/v1/differentiation/groups`: List all differentiation groups
- `/api/v1/differentiation/groups/{group_id}`: Get details about a differentiation group
- `/api/v1/molecules/{molecule_id}/differentiation`: Get differentiation information for a molecule

## Testing

A comprehensive test suite has been implemented in `test_consolidated_api.py` to verify that:

- Consolidated molecule endpoints correctly redirect to primary molecules
- Batch operations properly handle consolidated molecules
- Differentiation group endpoints return correct information
- All responses follow the standardized format

## Example Usage

An example script in `examples/use_consolidated_api.py` demonstrates how to use the consolidated molecule API:

- Fetching a molecule with consolidated handling
- Using batch operations with consolidated molecules
- Working with differentiation groups
- Listing consolidated molecule relationships

## Response Format Standardization

All responses from the new endpoints follow the standardized format:

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

Error responses follow a similar format with appropriate status codes and error details.

## Next Steps

To complete the API standardization process:

1. Apply the consolidated molecule handling to existing endpoints
2. Update client code to handle consolidated molecule information
3. Add documentation for clients on working with consolidated molecules
4. Develop integration tests for the entire API
5. Monitor API performance with the new handlers