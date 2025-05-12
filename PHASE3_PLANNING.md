# Molecule Data Quality Enhancement Plan - Phase 3 Planning

## Overview

Phase 3 of the Molecule Data Quality Enhancement Plan focuses on integrating the data quality improvements from Phases 1 and 2 into the application layers, specifically the API and user interface. This document outlines the key tasks, approaches, and considerations for Phase 3 implementation.

## 1. API Integration

### 1.1 Goals

- Ensure all API endpoints properly handle consolidated and differentiated molecules
- Provide consistent behavior when working with molecule IDs
- Maintain backward compatibility for existing integrations
- Add new endpoints to support consolidated molecule functionality

### 1.2 Implementation Tasks

#### 1.2.1 Core API Middleware

Create a middleware function in `api/middleware.py` that:

```python
def resolve_molecule_id(molecule_id):
    """
    Resolve a molecule ID to its primary, if it's consolidated.
    
    Args:
        molecule_id: The molecule ID to resolve
        
    Returns:
        Tuple of (resolved_id, is_consolidated)
    """
    with get_db_connection() as conn:
        with conn.cursor() as cursor:
            cursor.execute("""
                SELECT 
                    id, 
                    primary_molecule_id,
                    CASE WHEN primary_molecule_id IS NOT NULL THEN true ELSE false END as is_consolidated
                FROM molecules
                WHERE id = %s
            """, (molecule_id,))
            result = cursor.fetchone()
            
            if not result:
                return None, False
                
            if result[1]:  # Has a primary_molecule_id
                return result[1], True
            
            return result[0], False
```

#### 1.2.2 API Endpoint Updates

Update the following API endpoints to use the middleware:

1. **GET /api/molecules/{id}**
   - Add consolidation status information to the response
   - Add primary molecule reference if consolidated
   - Add differentiation group information if differentiated

2. **GET /api/molecular-properties/{molecule_id}**
   - Redirect to primary molecule properties if consolidated
   - Add indicator that properties are from a consolidated molecule

3. **GET /api/mixture-components/{molecule_id}**
   - Redirect to primary molecule components if consolidated
   - Add indicator that results are from a consolidated molecule

4. **POST /api/molecular-properties**
   - Automatically redirect property creation to primary molecule

5. **POST /api/mixture-components**
   - Automatically redirect component creation to primary molecule

#### 1.2.3 New API Endpoints

Add new endpoints to support consolidated molecule functionality:

1. **GET /api/molecules/{id}/consolidation-status**
   - Returns detailed information about consolidation status
   - Includes references to primary and secondary molecules

2. **GET /api/molecules/{id}/secondary-molecules**
   - For primary molecules, returns list of secondary molecules

3. **GET /api/differentiation-groups/{group_id}**
   - Returns all molecules in a differentiation group
   - Includes comparisons of key properties

### 1.3 API Documentation

Update API documentation to reflect changes:

1. Update OpenAPI schema to include new endpoints and response fields
2. Add documentation section for working with consolidated molecules
3. Create examples for common consolidated molecule operations

## 2. UI Enhancement

### 2.1 Goals

- Clearly indicate consolidated and differentiated molecules in the UI
- Provide intuitive navigation between primary and secondary molecules
- Show differentiation information for similar molecules

### 2.2 Implementation Tasks

#### 2.2.1 Molecule Detail View

Update the molecule detail view to:

1. Show a banner for secondary molecules indicating they are consolidated
2. Provide a link to the primary molecule on secondary molecule pages
3. Show a list of secondary molecules on primary molecule pages
4. Display differentiation properties clearly for differentiated molecules

#### 2.2.2 Search Results

Enhance search results to:

1. Prioritize primary molecules in search results
2. Indicate consolidated status with appropriate icons
3. Optionally include or exclude secondary molecules in results
4. Group differentiated molecules and highlight differences

#### 2.2.3 Molecule Creation and Editing

Update molecule creation and editing to:

1. Check for potential duplicates before creation
2. Suggest consolidation for similar molecules
3. Prevent editing of secondary molecules
4. Guide users through the consolidation process for duplicate molecules

## 3. Data Completeness Enhancement

### 3.1 Remaining "None" Names

Continue addressing molecules with "None" names:

1. Create a prioritized list of remaining molecules with "None" names
2. Implement batch processing for PubChem name retrieval with error handling
3. Create a fallback approach using molecular formula and SMILES for naming
4. Develop a UI component for manual naming of molecules

### 3.2 Implementation Approach

Create `complete_molecule_names_batch.py`:

```python
def process_none_name_molecules(batch_size=50, max_batches=None):
    """
    Process molecules with 'None' names in batches.
    
    Args:
        batch_size: Number of molecules to process in each batch
        max_batches: Maximum number of batches to process (None for all)
        
    Returns:
        Dict with statistics about the processing
    """
    # Implementation:
    # 1. Get molecules with None names prioritized by usage
    # 2. Process in batches with PubChem API
    # 3. Use molecular formula and SMILES as fallback
    # 4. Update database and log changes
```

## 4. Performance Optimization

### 4.1 Database Indexes

Add the following indexes to optimize queries for consolidated molecules:

```sql
-- Add index for primary_molecule_id lookups
CREATE INDEX IF NOT EXISTS idx_molecules_primary_molecule_id 
ON molecules(primary_molecule_id);

-- Add index for differentiation group lookups
CREATE INDEX IF NOT EXISTS idx_molecules_differentiation_group
ON molecules((properties->'differentiation'->>'differentiation_group'))
WHERE properties->'differentiation' IS NOT NULL;
```

### 4.2 Query Optimization

Create optimized query patterns:

1. Create common query patterns for consolidated molecules
2. Implement server-side functions for common operations
3. Create materialized views for frequent complex queries

## 5. Consolidated Molecule Audit Trail

### 5.1 Audit Table

Create an audit table to track consolidation operations:

```sql
CREATE TABLE IF NOT EXISTS molecule_consolidation_audit (
    id SERIAL PRIMARY KEY,
    operation_type TEXT NOT NULL,
    primary_molecule_id UUID NOT NULL,
    secondary_molecule_id UUID NOT NULL,
    performed_by TEXT,
    performed_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    operation_details JSONB
);

CREATE INDEX idx_molecule_consolidation_audit_primary
ON molecule_consolidation_audit(primary_molecule_id);

CREATE INDEX idx_molecule_consolidation_audit_secondary
ON molecule_consolidation_audit(secondary_molecule_id);
```

### 5.2 Trigger for Audit Trail

Create a trigger to automatically log consolidation changes:

```sql
CREATE OR REPLACE FUNCTION log_molecule_consolidation()
RETURNS TRIGGER AS $$
BEGIN
    IF (TG_OP = 'UPDATE' AND 
        NEW.primary_molecule_id IS DISTINCT FROM OLD.primary_molecule_id) THEN
        -- Log consolidation change
        INSERT INTO molecule_consolidation_audit (
            operation_type, 
            primary_molecule_id, 
            secondary_molecule_id,
            performed_by,
            operation_details
        ) VALUES (
            'CONSOLIDATE',
            NEW.primary_molecule_id,
            NEW.id,
            current_user,
            jsonb_build_object(
                'previous_state', 
                CASE WHEN OLD.primary_molecule_id IS NULL 
                     THEN 'INDEPENDENT' 
                     ELSE 'CONSOLIDATED TO ' || OLD.primary_molecule_id::text 
                END,
                'new_state', 'CONSOLIDATED TO ' || NEW.primary_molecule_id::text
            )
        );
    END IF;
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER molecule_consolidation_audit_trigger
AFTER UPDATE ON molecules
FOR EACH ROW
WHEN (NEW.primary_molecule_id IS NOT NULL OR OLD.primary_molecule_id IS NOT NULL)
EXECUTE FUNCTION log_molecule_consolidation();
```

### 5.3 Audit Report

Create an API endpoint to retrieve the consolidation audit trail:

```python
@app.route('/api/admin/consolidation-audit', methods=['GET'])
@admin_required
def get_consolidation_audit():
    """Get the consolidation audit trail."""
    page = request.args.get('page', 1, type=int)
    per_page = request.args.get('per_page', 50, type=int)
    
    # Implementation:
    # 1. Query the audit table with pagination
    # 2. Join with molecules to get additional context
    # 3. Return formatted audit records
```

## 6. Implementation Timeline

| Task | Priority | Estimated Effort | Dependencies |
|------|----------|------------------|--------------|
| API Middleware | High | 2 days | None |
| Core API Endpoint Updates | High | 3 days | API Middleware |
| New API Endpoints | Medium | 2 days | API Middleware |
| API Documentation | Medium | 1 day | All API implementation |
| UI Molecule Detail Updates | High | 2 days | API Endpoint Updates |
| UI Search Results Updates | Medium | 2 days | API Endpoint Updates |
| UI Creation/Editing Updates | Medium | 3 days | API Endpoint Updates |
| Database Indexes | High | 1 day | None |
| Query Optimization | Medium | 2 days | Database Indexes |
| Audit Trail Implementation | Low | 2 days | None |
| Remaining "None" Names | Low | 3 days | None |

## 7. Testing Strategy

### 7.1 Unit Tests

Create specific unit tests for:

1. Middleware functions for molecule ID resolution
2. Updated API endpoint behaviors with consolidated molecules
3. UI components that deal with consolidated and differentiated molecules

### 7.2 Integration Tests

Create integration tests for:

1. End-to-end workflows involving consolidated molecules
2. API response validation for consolidated molecule endpoints
3. Database trigger behavior and audit trail creation

### 7.3 Performance Tests

Implement performance tests for:

1. Query response times with consolidated molecules
2. API endpoint performance with large molecule sets
3. UI responsiveness with consolidated molecule data

## 8. Conclusion

Phase 3 builds on the solid foundation established in Phases 1 and 2 by integrating the data quality improvements into the application layers. This phase will significantly enhance the user experience by providing clear visibility into consolidated and differentiated molecules while maintaining data integrity through the database trigger system already implemented.