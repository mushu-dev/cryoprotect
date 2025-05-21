# Phase 3 Implementation Plan

## Overview

This document outlines the implementation plan for Phase 3 of the Molecule Data Quality Enhancement project. Phase 3 focuses on integrating the data quality improvements from Phases 1 and 2 into the application layers, specifically the API and underlying database optimizations.

## Goals

1. Ensure all API endpoints properly handle consolidated and differentiated molecules
2. Optimize database queries for working with consolidated molecules
3. Implement audit trail for tracking molecule consolidation changes
4. Complete naming of molecules with "None" names
5. Maintain backward compatibility for existing integrations

## Implementation Phases

The implementation will be broken down into the following phases, each with specific deliverables:

### Phase 3.1: Core API Middleware and Database Optimizations

#### Tasks

1. **Create API Middleware for Molecule ID Resolution**
   - Implement `api/middleware.py` with `resolve_molecule_id` function
   - Add helper functions for batch molecule ID resolution
   - Create tests for middleware functions

2. **Add Database Indexes for Optimized Queries**
   - Create migration for adding index on `primary_molecule_id`
   - Create migration for adding index on differentiation group fields
   - Implement query optimization functions

3. **Update Database Query Patterns**
   - Create common query patterns for consolidated molecules
   - Implement database helper functions for common operations
   - Create materialized views for frequent complex queries

#### Deliverables

- `api/middleware.py` implementation
- Database migration scripts for indexes
- Unit tests for middleware functions
- Documentation for query patterns

### Phase 3.2: API Endpoint Integration

#### Tasks

1. **Update Existing API Endpoints**
   - Modify molecule endpoints to use middleware
   - Update molecular property endpoints to handle consolidated molecules
   - Update mixture component endpoints to handle consolidated molecules
   - Ensure all endpoints follow standardized response format

2. **Implement New Consolidated Molecule API Endpoints**
   - Create endpoint for molecule consolidation status
   - Create endpoint for secondary molecules
   - Create endpoint for differentiation groups
   - Implement appropriate request validation and error handling

3. **Update API Documentation**
   - Update OpenAPI schema to include new endpoints
   - Add documentation section for working with consolidated molecules
   - Create examples for common consolidated molecule operations

#### Deliverables

- Updated API endpoint implementations
- New API endpoints for consolidated molecules
- Updated OpenAPI documentation
- Integration tests for API endpoints

### Phase 3.3: Audit Trail and Data Completeness

#### Tasks

1. **Implement Consolidated Molecule Audit Trail**
   - Create audit table for tracking consolidation changes
   - Implement database trigger for automatic logging
   - Create API endpoint for retrieving audit trail
   - Add admin interface for viewing audit trail

2. **Implement Batch Processor for 'None' Names**
   - Create batch processor for molecules with 'None' names
   - Implement PubChem name retrieval with error handling
   - Create fallback naming approach using molecular formula
   - Develop command-line interface for manual naming

#### Deliverables

- Database migrations for audit table
- Trigger implementation for audit logging
- Batch processor for 'None' names
- Documentation for audit trail and batch processing

## Implementation Details

### API Middleware Implementation

The core of Phase 3 is the API middleware for molecule ID resolution:

```python
# api/middleware.py

from flask import g, current_app
import logging
from typing import Tuple, List, Dict, Optional, Any

from database.adapter import get_connection

logger = logging.getLogger(__name__)

def resolve_molecule_id(molecule_id: str) -> Tuple[str, bool]:
    """
    Resolve a molecule ID to its primary molecule ID if it's consolidated.
    
    Args:
        molecule_id: The molecule ID to resolve
        
    Returns:
        Tuple of (resolved_id, is_consolidated)
        - resolved_id: The primary molecule ID if consolidated, otherwise the original ID
        - is_consolidated: True if the molecule is consolidated, False otherwise
    """
    with get_connection() as conn:
        cursor = conn.cursor()
        
        try:
            query = """
                SELECT 
                    id, 
                    consolidated_to,
                    consolidated_to IS NOT NULL as is_consolidated
                FROM molecules
                WHERE id = %s
            """
            
            cursor.execute(query, (molecule_id,))
            result = cursor.fetchone()
            
            if not result:
                return molecule_id, False
                
            if result[2]:  # Is consolidated
                return result[1], True
            
            return result[0], False
        except Exception as e:
            logger.error(f"Error resolving molecule ID {molecule_id}: {str(e)}")
            return molecule_id, False
        finally:
            cursor.close()

def resolve_molecule_ids(molecule_ids: List[str]) -> Dict[str, Tuple[str, bool]]:
    """
    Resolve multiple molecule IDs to their primary molecule IDs.
    
    Args:
        molecule_ids: List of molecule IDs to resolve
        
    Returns:
        Dictionary mapping original IDs to tuples of (resolved_id, is_consolidated)
    """
    results = {}
    
    # Default to original IDs if resolution fails
    for molecule_id in molecule_ids:
        results[molecule_id] = (molecule_id, False)
    
    if not molecule_ids:
        return results
    
    with get_connection() as conn:
        cursor = conn.cursor()
        
        try:
            placeholders = ','.join(['%s'] * len(molecule_ids))
            query = f"""
                SELECT 
                    id, 
                    consolidated_to,
                    consolidated_to IS NOT NULL as is_consolidated
                FROM molecules
                WHERE id IN ({placeholders})
            """
            
            cursor.execute(query, molecule_ids)
            for row in cursor.fetchall():
                molecule_id = row[0]
                consolidated_to = row[1]
                is_consolidated = row[2]
                
                if is_consolidated:
                    results[molecule_id] = (consolidated_to, True)
                else:
                    results[molecule_id] = (molecule_id, False)
            
            return results
        except Exception as e:
            logger.error(f"Error resolving molecule IDs {molecule_ids}: {str(e)}")
            return results
        finally:
            cursor.close()

def get_molecule_with_consolidated_info(molecule_id: str) -> Optional[Dict[str, Any]]:
    """
    Get a molecule with consolidated molecule information.
    
    Args:
        molecule_id: The molecule ID to get
        
    Returns:
        Dictionary containing molecule data with consolidated information,
        or None if the molecule doesn't exist
    """
    resolved_id, is_consolidated = resolve_molecule_id(molecule_id)
    
    with get_connection() as conn:
        cursor = conn.cursor()
        
        try:
            query = """
                SELECT 
                    id, 
                    name,
                    molecular_formula,
                    smiles,
                    inchi_key,
                    consolidated_to,
                    (SELECT jsonb_agg(id) 
                     FROM molecules 
                     WHERE consolidated_to = m.id) as consolidated_molecules,
                    (SELECT property_value 
                     FROM molecular_properties mp
                     JOIN property_types pt ON mp.property_type_id = pt.id
                     WHERE mp.molecule_id = m.id AND pt.name = 'differentiationGroup'
                     LIMIT 1) as differentiation_group,
                    (SELECT property_value 
                     FROM molecular_properties mp
                     JOIN property_types pt ON mp.property_type_id = pt.id
                     WHERE mp.molecule_id = m.id AND pt.name = 'differentiationDescription'
                     LIMIT 1) as differentiation_description,
                    created_at,
                    updated_at
                FROM molecules m
                WHERE id = %s
            """
            
            cursor.execute(query, (resolved_id,))
            result = cursor.fetchone()
            
            if not result:
                return None
                
            # Convert to dictionary
            columns = [desc[0] for desc in cursor.description]
            molecule_data = dict(zip(columns, result))
            
            # Add consolidated info
            molecule_data['is_consolidated'] = is_consolidated
            
            # If this is a consolidated molecule (secondary), include the original ID
            if is_consolidated:
                molecule_data['original_molecule_id'] = molecule_id
            
            return molecule_data
        except Exception as e:
            logger.error(f"Error getting molecule with consolidated info {molecule_id}: {str(e)}")
            return None
        finally:
            cursor.close()
```

### Database Indexes

To optimize queries involving consolidated molecules, we'll add these indexes:

```sql
-- Add index for consolidated_to lookups
CREATE INDEX IF NOT EXISTS idx_molecules_consolidated_to 
ON molecules(consolidated_to)
WHERE consolidated_to IS NOT NULL;

-- Add index for differentiation group lookups
CREATE INDEX IF NOT EXISTS idx_molecular_properties_differentiation_group
ON molecular_properties(property_value, molecule_id)
WHERE property_type_id = 'differentiationGroup';
```

### Consolidated Molecule Audit Trail

To track changes to molecule consolidation, we'll implement:

```sql
-- Audit table for tracking consolidation changes
CREATE TABLE IF NOT EXISTS molecule_consolidation_audit (
    id SERIAL PRIMARY KEY,
    operation_type TEXT NOT NULL,
    primary_molecule_id UUID NOT NULL,
    secondary_molecule_id UUID NOT NULL,
    performed_by TEXT,
    performed_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    operation_details JSONB
);

-- Indexes for efficient querying
CREATE INDEX idx_molecule_consolidation_audit_primary
ON molecule_consolidation_audit(primary_molecule_id);

CREATE INDEX idx_molecule_consolidation_audit_secondary
ON molecule_consolidation_audit(secondary_molecule_id);

-- Trigger for automatic logging
CREATE OR REPLACE FUNCTION log_molecule_consolidation()
RETURNS TRIGGER AS $$
BEGIN
    IF (TG_OP = 'UPDATE' AND 
        NEW.consolidated_to IS DISTINCT FROM OLD.consolidated_to) THEN
        -- Log consolidation change
        INSERT INTO molecule_consolidation_audit (
            operation_type, 
            primary_molecule_id, 
            secondary_molecule_id,
            performed_by,
            operation_details
        ) VALUES (
            CASE
                WHEN NEW.consolidated_to IS NULL THEN 'DECONSOLIDATE'
                WHEN OLD.consolidated_to IS NULL THEN 'CONSOLIDATE'
                ELSE 'CHANGE_PRIMARY'
            END,
            COALESCE(NEW.consolidated_to, OLD.consolidated_to),
            NEW.id,
            current_user,
            jsonb_build_object(
                'previous_state', 
                CASE WHEN OLD.consolidated_to IS NULL 
                     THEN 'INDEPENDENT' 
                     ELSE 'CONSOLIDATED TO ' || OLD.consolidated_to::text 
                END,
                'new_state', 
                CASE WHEN NEW.consolidated_to IS NULL
                     THEN 'INDEPENDENT'
                     ELSE 'CONSOLIDATED TO ' || NEW.consolidated_to::text
                END
            )
        );
    END IF;
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER molecule_consolidation_audit_trigger
AFTER UPDATE ON molecules
FOR EACH ROW
WHEN (NEW.consolidated_to IS NOT NULL OR OLD.consolidated_to IS NOT NULL)
EXECUTE FUNCTION log_molecule_consolidation();
```

## Timeline

| Phase | Task | Priority | Estimated Effort | Dependencies |
|------|------|----------|------------------|--------------|
| 3.1 | Create API Middleware | High | 2 days | None |
| 3.1 | Add Database Indexes | High | 1 day | None |
| 3.1 | Update Database Query Patterns | Medium | 2 days | Database Indexes |
| 3.2 | Update Existing API Endpoints | High | 3 days | API Middleware |
| 3.2 | Implement New API Endpoints | Medium | 2 days | API Middleware |
| 3.2 | Update API Documentation | Medium | 1 day | API Implementation |
| 3.3 | Implement Audit Trail | Low | 2 days | None |
| 3.3 | Implement Batch Processor | Low | 3 days | None |

## Testing Strategy

1. **Unit Tests**: Write tests for middleware functions, API endpoints, and database helpers
2. **Integration Tests**: Test end-to-end workflows with consolidated molecules
3. **Performance Tests**: Verify query optimization with large datasets

## Success Criteria

Phase 3 will be considered complete when:

1. All API endpoints properly handle consolidated and differentiated molecules
2. Database queries are optimized for consolidated molecules
3. Audit trail is in place for tracking consolidation changes
4. Batch processor is available for handling 'None' names
5. All tests pass and documentation is updated

## Next Steps

After Phase 3 is complete, focus will shift to UI enhancements to better visualize consolidated and differentiated molecules in the user interface.