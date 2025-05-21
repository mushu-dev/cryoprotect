# Consolidated Molecule Audit Trail Guide

This guide explains the consolidated molecule audit trail functionality implemented in Phase 3 of the CryoProtect project. The audit trail provides a comprehensive history of all molecule consolidation operations, making it possible to track when molecules are consolidated, deconsolidated, or have their primary molecule changed.

## Overview

The consolidated molecule audit trail consists of:

1. A database table (`molecule_consolidation_audit`) to store audit records
2. Database triggers to automatically create audit records when consolidation changes occur
3. An API endpoint to retrieve and filter audit records

## Database Components

### Audit Table Schema

The `molecule_consolidation_audit` table has the following structure:

```sql
CREATE TABLE molecule_consolidation_audit (
    id SERIAL PRIMARY KEY,
    operation_type TEXT NOT NULL,          -- Type of operation (CONSOLIDATE, DECONSOLIDATE, CHANGE_PRIMARY)
    primary_molecule_id UUID NOT NULL,     -- ID of the primary molecule
    secondary_molecule_id UUID NOT NULL,   -- ID of the secondary molecule
    performed_by TEXT,                     -- User who performed the operation
    performed_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(), -- When the operation was performed
    operation_details JSONB                -- Additional details about the operation
);
```

### Database Trigger

A trigger on the `molecules` table automatically creates audit records when the `consolidated_to` field changes:

```sql
CREATE TRIGGER molecule_consolidation_audit_trigger
AFTER UPDATE ON molecules
FOR EACH ROW
WHEN (NEW.consolidated_to IS DISTINCT FROM OLD.consolidated_to)
EXECUTE FUNCTION log_molecule_consolidation();
```

The trigger function `log_molecule_consolidation()` determines the type of operation based on the changes:

- **CONSOLIDATE**: When a molecule is consolidated to another molecule (NULL → UUID)
- **DECONSOLIDATE**: When a molecule is deconsolidated (UUID → NULL)
- **CHANGE_PRIMARY**: When a molecule's primary molecule is changed (UUID → different UUID)

## API Components

### Audit Endpoint

The audit records can be accessed through the following API endpoint:

```
GET /api/v1/admin/consolidation-audit
```

This endpoint requires administrative privileges to access.

### Query Parameters

The following query parameters can be used to filter audit records:

- `molecule_id`: Filter by primary or secondary molecule ID
- `operation_type`: Filter by operation type (CONSOLIDATE, DECONSOLIDATE, CHANGE_PRIMARY)
- `performed_by`: Filter by user who performed the operation
- `start_date`: Filter by operations performed on or after this date
- `end_date`: Filter by operations performed on or before this date
- `page`: Page number for pagination (default: 1)
- `per_page`: Number of records per page (default: 20, max: 100)

### Response Format

The endpoint returns a standard response with the following structure:

```json
{
  "status": "success",
  "message": "Consolidation audit records retrieved successfully",
  "data": {
    "audit_records": [
      {
        "id": 1,
        "operation_type": "CONSOLIDATE",
        "primary_molecule_id": "uuid-of-primary-molecule",
        "secondary_molecule_id": "uuid-of-secondary-molecule",
        "performed_by": "username",
        "performed_at": "2025-05-13T12:34:56.789Z",
        "operation_details": {
          "previous_state": "INDEPENDENT",
          "new_state": "CONSOLIDATED TO uuid-of-primary-molecule",
          "timestamp": "2025-05-13T12:34:56.789Z",
          "molecule_name": "Example Molecule"
        },
        "primary_molecule_name": "Primary Molecule Name",
        "secondary_molecule_name": "Secondary Molecule Name"
      }
      // ... more audit records
    ],
    "count": 1
  },
  "pagination": {
    "page": 1,
    "per_page": 20,
    "total_items": 1,
    "total_pages": 1,
    "has_next": false,
    "has_prev": false
  }
}
```

## Use Cases

### Tracking Consolidation History

The audit trail is particularly useful for:

1. **Accountability**: Tracking who made consolidation changes and when
2. **Troubleshooting**: Investigating issues related to molecule consolidation
3. **Regulatory Compliance**: Providing a complete history of changes for audit purposes
4. **Data Quality Assurance**: Verifying that consolidation operations are correctly applied

### Example Queries

To view all consolidation operations for a specific molecule:

```
GET /api/v1/admin/consolidation-audit?molecule_id=uuid-of-molecule
```

To view all consolidations performed by a specific user:

```
GET /api/v1/admin/consolidation-audit?operation_type=CONSOLIDATE&performed_by=username
```

To view all deconsolidations in the last week:

```
GET /api/v1/admin/consolidation-audit?operation_type=DECONSOLIDATE&start_date=2025-05-06T00:00:00Z
```

## Implementation Details

The audit trail implementation is divided into three main components:

1. **Database Migration** (`migrations/025_consolidated_molecule_audit.sql`): Creates the audit table, indexes, and trigger
2. **API Resource** (`api/audit_resources.py`): Provides the endpoint for retrieving audit records
3. **Tests** (`tests/test_consolidation_audit.py`): Verifies the functionality of the audit trail

## Testing

The audit trail implementation includes comprehensive tests that verify:

1. The trigger correctly creates audit records for consolidation operations
2. The API endpoint returns the correct audit records
3. Admin permissions are properly enforced
4. Filtering and pagination work correctly

To run the tests:

```bash
python -m unittest tests/test_consolidation_audit.py
```

## Verification

A verification script is provided to check that all components of the audit trail have been correctly implemented:

```bash
python verify_audit_implementation.py
```

This script checks for the presence of all required files, functions, and configuration entries.

## Conclusion

The consolidated molecule audit trail provides complete transparency into molecule consolidation operations, enhancing data governance and accountability in the CryoProtect system. The implementation leverages database triggers for automatic record creation and provides a secure API for retrieving the audit history.