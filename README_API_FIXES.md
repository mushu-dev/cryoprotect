# CryoProtect v2 API Integration Fixes

This document provides information about the API integration fixes for the CryoProtect v2 project to work with the new standardized database schema.

## Overview

The database tables have been renamed from singular to plural form (e.g., "molecule" to "molecules"), and these fixes ensure the API endpoints reference the correct tables.

### Key Changes

1. **Updated API resource classes** to reference the correct plural table names
2. **Fixed endpoint duplication issues** (especially in ComparisonResource)
3. **Implemented consistent error handling**
4. **Added retry logic for API calls**

## Table Name Mapping

The following table name mapping was applied:

```python
TABLE_NAME_MAPPING = {
    "molecule": "molecules",
    "mixture": "mixtures",
    "prediction": "predictions",
    "experiment": "experiments",
    "molecular_property": "molecular_properties",
    "mixture_component": "mixture_components",
    "calculation_method": "calculation_methods",
    "property_type": "property_types"
}
```

## Files Updated

The following files were updated:

1. `api/resources.py` - Updated resource classes to use plural table names
2. `api/utils.py` - Enhanced error handling and added retry logic
3. `api/__init__.py` - Fixed endpoint registration issues
4. `app.py` - Added proper error handling for API requests

## How to Apply the Fixes

### Automatic Method

Run the provided script to automatically apply all the fixes:

#### Windows

```
run_api_fixes.bat
```

#### Linux/macOS

```
chmod +x run_api_fixes.sh
./run_api_fixes.sh
```

### Manual Method

If you prefer to apply the changes manually, follow these steps:

1. Update `api/resources.py`:
   - Replace all occurrences of singular table names with their plural forms
   - Update view names (e.g., "molecule_with_properties" to "molecules_with_properties")

2. Update `api/utils.py`:
   - Add retry logic for API calls
   - Enhance error handling

3. Update `api/__init__.py`:
   - Fix endpoint duplication in ComparisonResource

4. Update `app.py`:
   - Add proper error handling for API requests

## Verification

After applying the fixes, you can verify that the API is working correctly by:

1. Running the application:
   ```
   python app.py
   ```

2. Accessing the health check endpoint:
   ```
   curl http://localhost:5000/health
   ```

3. Testing API endpoints:
   ```
   curl http://localhost:5000/api/v1/molecules
   curl http://localhost:5000/api/v1/mixtures
   ```

## Rollback

If you need to rollback the changes, backup files were created with the `.bak` extension. You can restore them by:

1. Renaming the backup files to their original names:
   ```
   mv api/resources.py.bak api/resources.py
   mv api/utils.py.bak api/utils.py
   mv api/__init__.py.bak api/__init__.py
   mv app.py.bak app.py
   ```

## Troubleshooting

If you encounter any issues after applying the fixes:

1. Check the log files in the `logs` directory
2. Ensure all database tables have been properly renamed
3. Verify that the database views have been updated to reflect the new table names
4. Check for any custom code that might be referencing the old table names

## Additional Information

The retry logic added to the API calls helps handle transient errors by:

1. Automatically retrying failed API calls up to 3 times
2. Using exponential backoff with jitter to prevent thundering herd problems
3. Logging detailed information about retries and failures

The enhanced error handling provides:

1. More detailed error messages
2. Proper HTTP status codes
3. Consistent error response format
4. Improved logging of errors

## Contact

If you have any questions or need assistance, please contact the CryoProtect development team.
## [2025-04-20] Dashboard API Standardization

- Standardized all methods in `api/dashboard_resources.py` according to recommendations in `API_Standardization_Report.md` and `API_Standardization_Summary.md`.
- Applied consistent error handling using `handle_error` utility.
- Enforced authentication with `@token_required` decorator on all appropriate endpoints.
- Standardized response formatting using `@marshal_with` and removed direct use of `_handle_json_serialization`.
- Updated field schemas for marshalling responses.
- All changes verified by passing `tests/test_resources_standardized.py`.
- This update ensures maintainable, consistent, and robust API behavior for all dashboard endpoints.

## [2025-04-21] Core API Resources Standardization

- Standardized all methods in `api/resources.py` according to the API standardization plan.
- Applied consistent error handling using `handle_error` utility instead of direct `abort` calls.
- Enforced authentication with `@token_required` decorator on all protected endpoints.
- Standardized response formatting using `@marshal_with` and removed direct use of `_handle_json_serialization`.
- Added comprehensive docstrings for all resource classes and methods.
- Improved request validation using class-level `reqparse.RequestParser` definitions.
- Added proper user parameter to all authenticated methods.
- Standardized error responses with consistent context information and status codes.
- All changes verified by passing API endpoint tests.
- This update ensures consistent error handling, authentication, and response formatting across all core API resources.