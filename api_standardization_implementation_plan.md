# API Standardization Implementation Plan

## Overview

This plan outlines the steps needed to implement API standardization for CryoProtect, ensuring that API endpoints properly handle the data quality improvements from Phase 2, including:

1. Consolidated molecules handling
2. Standardized names
3. Differentiated molecules
4. Molecule properties standardization
5. Database triggers integration

## Implementation Steps

### 1. Update Response Models

1. Create molecule response models that include consolidated molecule information:
   - Add `is_consolidated` flag to molecule response
   - Add `consolidated_to` reference to primary molecule when applicable
   - Add `consolidated_molecules` list for primary molecules
   - Add `differentiation_group` information when applicable

2. Update schema definitions:
   ```python
   # In api/models.py
   molecule_fields = {
       'id': fields.String,
       'name': fields.String,
       'molecular_formula': fields.String,
       'smiles': fields.String,
       'inchi_key': fields.String,
       'is_consolidated': fields.Boolean,
       'consolidated_to': fields.String,
       'consolidated_molecules': fields.List(fields.String),
       'differentiation_group': fields.String,
       # Other existing fields
   }
   ```

### 2. Create Consolidated Molecule Utilities

1. Create a utility module for handling consolidated molecules in the API:
   ```python
   # api/consolidated_utils.py
   
   def get_primary_molecule(molecule_id):
       """
       Gets the primary molecule for a given molecule ID.
       
       Args:
           molecule_id: The molecule ID to check
           
       Returns:
           Primary molecule ID if consolidated, otherwise the original ID
       """
       # Implementation
   
   def get_consolidated_molecules(molecule_id):
       """
       Gets the list of consolidated molecules for a primary molecule.
       
       Args:
           molecule_id: The primary molecule ID
           
       Returns:
           List of consolidated molecule IDs
       """
       # Implementation
   ```

### 3. Update API Endpoint Handlers

1. Add decorators to handle consolidated molecules automatically:
   ```python
   # api/api_decorators.py
   
   def handle_consolidated_molecules(f):
       """
       Decorator to handle consolidated molecules in API responses.
       """
       @functools.wraps(f)
       def decorated(*args, **kwargs):
           # Get molecule ID from request
           molecule_id = kwargs.get('molecule_id') or request.args.get('molecule_id')
           
           if molecule_id:
               # Check if molecule is consolidated
               primary_id = get_primary_molecule(molecule_id)
               
               if primary_id != molecule_id:
                   # Redirect to primary molecule
                   kwargs['molecule_id'] = primary_id
                   kwargs['original_molecule_id'] = molecule_id
           
           # Call the original function
           result = f(*args, **kwargs)
           
           # Add consolidated information to response
           if isinstance(result, dict) and 'id' in result:
               result['is_consolidated'] = is_consolidated(result['id'])
               result['consolidated_to'] = get_primary_molecule(result['id'])
               result['consolidated_molecules'] = get_consolidated_molecules(result['id'])
           
           return result
       
       return decorated
   ```

2. Apply the decorator to relevant endpoints:
   ```python
   # In API resource classes
   
   @handle_consolidated_molecules
   def get(self, molecule_id):
       # Existing implementation
   ```

### 4. Update Response Formatting

1. Modify response creation to include standardized molecule names:
   ```python
   # When returning molecule data
   def get_molecule_data(molecule_id):
       # Get molecule from database
       molecule = query_molecule(molecule_id)
       
       # Use standardized name
       if molecule and molecule.name:
           molecule_data = {
               'id': molecule.id,
               'name': molecule.name,  # This is now standardized
               # Other fields
           }
           
           # Add standardization information
           if hasattr(molecule, 'original_name') and molecule.original_name != molecule.name:
               molecule_data['original_name'] = molecule.original_name
               
           return molecule_data
   ```

### 5. Add Differentiation Group Support

1. Create utilities for handling differentiated molecules:
   ```python
   # api/differentiation_utils.py
   
   def get_differentiation_group(molecule_id):
       """
       Gets the differentiation group for a molecule.
       
       Args:
           molecule_id: The molecule ID to check
           
       Returns:
           Differentiation group ID if differentiated, otherwise None
       """
       # Implementation
   
   def get_differentiation_group_members(group_id):
       """
       Gets all molecules in a differentiation group.
       
       Args:
           group_id: The differentiation group ID
           
       Returns:
           List of molecule IDs in the group
       """
       # Implementation
   ```

2. Create new API endpoints for differentiation groups:
   ```python
   # api/resources.py
   
   class DifferentiationGroupResource(Resource):
       """Resource for differentiation groups."""
       
       @doc(description='Get information about a differentiation group',
            tags=['Molecules', 'Differentiation'])
       def get(self, group_id):
           """Get information about a differentiation group."""
           # Implementation
   
   class DifferentiationGroupListResource(Resource):
       """Resource for listing differentiation groups."""
       
       @doc(description='List all differentiation groups',
            tags=['Molecules', 'Differentiation'])
       def get(self):
           """List all differentiation groups."""
           # Implementation
   ```

### 6. Update Batch Resources

1. Update batch processing endpoints to handle consolidated molecules:
   ```python
   # api/batch_resources.py
   
   class MoleculeBatchResource(Resource):
       """Resource for batch operations on molecules."""
       
       @doc(description='Get multiple molecules by ID',
            tags=['Molecules', 'Batch'])
       def post(self):
           """Get multiple molecules by ID."""
           # Get molecule IDs from request
           molecule_ids = request.json.get('molecule_ids', [])
           
           # Handle consolidated molecules
           primary_ids = [get_primary_molecule(id) for id in molecule_ids]
           
           # Query primary molecules
           molecules = query_molecules(primary_ids)
           
           # Add consolidated information to response
           for molecule in molecules:
               molecule['is_consolidated'] = is_consolidated(molecule['id'])
               molecule['consolidated_to'] = get_primary_molecule(molecule['id'])
               molecule['consolidated_molecules'] = get_consolidated_molecules(molecule['id'])
           
           return molecules
   ```

### 7. Update Database Query Utilities

1. Create a utility to ensure all database queries handle consolidated molecules:
   ```python
   # api/db_utils.py
   
   def get_molecule_query_builder():
       """
       Creates a query builder for molecules that handles consolidated molecules.
       
       Returns:
           QueryBuilder instance
       """
       # Import query builder
       from database.query_builder import QueryBuilder
       
       # Create query builder
       query = QueryBuilder('molecules')
       
       # Add hook to handle consolidated molecules
       query.add_hook('after_select', handle_consolidated_molecule_rows)
       
       return query
   ```

### 8. Update API Tests

1. Create tests for consolidated molecule handling:
   ```python
   # tests/test_api_consolidated_handling.py
   
   def test_consolidated_molecule_redirect():
       """Test that consolidated molecules redirect to primary."""
       # Create a consolidated molecule
       # Test API response
   
   def test_consolidated_molecule_info():
       """Test that molecule responses include consolidated information."""
       # Create a primary molecule with consolidated molecules
       # Test API response
   ```

### 9. Update API Documentation

1. Update OpenAPI documentation to reflect consolidated molecule handling:
   ```python
   # api/openapi.py
   
   # Define consolidated molecule schemas
   consolidated_molecule_schema = {
       'type': 'object',
       'properties': {
           'is_consolidated': {
               'type': 'boolean',
               'description': 'Whether this molecule is consolidated'
           },
           'consolidated_to': {
               'type': 'string',
               'description': 'ID of the primary molecule if this is consolidated'
           },
           'consolidated_molecules': {
               'type': 'array',
               'items': {
                   'type': 'string'
               },
               'description': 'List of consolidated molecule IDs for a primary molecule'
           }
       }
   }
   
   # Update molecule schema
   molecule_schema.update(consolidated_molecule_schema)
   ```

## Deployment Strategy

1. **Development Environment**
   - Implement and test API standardization in development
   - Run automated tests for all endpoints
   - Verify consolidated molecule handling

2. **Staging Environment**
   - Deploy to staging environment
   - Test with real consolidated data
   - Verify performance and database query optimization

3. **Production Environment**
   - Deploy to production during low-traffic period
   - Monitor API response times and errors
   - Ensure backward compatibility for existing clients

## Backward Compatibility

To maintain backward compatibility:

1. Continue to accept molecule IDs for consolidated molecules
2. Automatically redirect to primary molecules
3. Include both primary and legacy fields in responses
4. Add deprecation headers for endpoints that will change in the future

## Next Steps

After implementing API standardization:

1. Update client applications to handle consolidated molecules
2. Create documentation for API consumers explaining the consolidation model
3. Develop additional API endpoints for advanced molecule queries
4. Implement performance monitoring for API endpoints