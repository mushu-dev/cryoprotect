# Roo PM Directive: Complete Phase 1 - Essential Core Functionality

## Current Status and Goals

We're currently in **Phase 1: Essential Core Functionality** of the CryoProtect v2 project. Our analysis shows that we need to complete the following high-priority items:

1. Complete the lab verification workflow by registering API endpoints
2. Remove login requirements from core pages (molecules, mixtures)
3. Finalize remaining API implementation

These tasks are critical for completing Phase 1 and moving on to Phase 2 (Minimal Viable Deployment).

## Task 1: Complete Lab Verification API Registration

### Requirements
- Register the lab verification API endpoints in the API system

### Implementation Details
Most of the lab verification workflow is already implemented:
- ✅ Database schema (`migrations/013_lab_verification_schema.sql`)
- ✅ Data model (`LabVerification` in `api/models.py`)
- ✅ API resources (`api/lab_verification_resources.py`) 
- ✅ Schema definitions (`lab_verification_fields` in `api/schemas.py`)
- ✅ Frontend implementation (`static/js/verification.js`)
- ✅ Tests (`tests/test_lab_verification.py`)
- ❌ **API Route Registration** (missing)

### Implementation Steps
1. Edit `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/__init__.py` to add these routes after line 164:

```python
# Register Lab Verification resources
api.add_resource(LabVerificationResource, 
                '/experiments/<string:experiment_id>/verification', 
                endpoint='experiment_verification')
api.add_resource(LabVerificationResource, 
                '/verifications/<string:verification_id>', 
                endpoint='verification_update')
api.add_resource(VerificationStatsResource, 
                '/verification/stats', 
                endpoint='verification_stats')
```

2. Verify that all lab verification endpoints are working correctly.

## Task 2: Remove Login Requirements for Core Pages

### Requirements
- Allow users to access molecules and mixtures pages without logging in
- Maintain authentication for secured operations (create, update, delete)
- Keep the login page but make it optional

### Implementation Steps

1. Modify `/mnt/c/Users/1edwa/Documents/CryoProtect v2/app.py` to remove login requirements from core route handlers.

2. Update templates to handle unauthenticated users:
   - Modify `/mnt/c/Users/1edwa/Documents/CryoProtect v2/templates/molecules.html`
   - Modify `/mnt/c/Users/1edwa/Documents/CryoProtect v2/templates/mixtures.html`

3. Update JavaScript to check for authentication only for secured operations:
   - Modify `/mnt/c/Users/1edwa/Documents/CryoProtect v2/static/js/app.js`
   - Modify `/mnt/c/Users/1edwa/Documents/CryoProtect v2/static/js/auth.js`

4. Update navigation to show login/logout in header but not force redirect:
   - Modify `/mnt/c/Users/1edwa/Documents/CryoProtect v2/templates/base.html`

5. Ensure read-only operations work without authentication:
   - Update any API endpoint decorators to allow GET requests without authentication
   - Update any frontend code that assumes authentication

## Task 3: Finalize API Implementation

### Requirements
- Ensure all API endpoints are properly registered
- Implement any missing API validation
- Ensure consistent error handling across all endpoints

### Implementation Steps

1. Register any remaining API endpoints:
   - Ensure lab verification endpoints are registered
   - Check for any other unregistered endpoints

2. Standardize API error handling:
   - Verify all endpoints use the standardized error response format
   - Review API calls in `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/resources.py`

3. Improve API validation:
   - Ensure all endpoints validate input parameters
   - Complete any TODOs in validation logic

4. Update API documentation:
   - Ensure all endpoints are documented in OpenAPI format
   - Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/api_docs.py`

## Supporting Files

### Lab Verification Implementation
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/lab_verification_resources.py` (existing)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/models.py` (lines 2604-2740)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/schemas.py` (lab_verification_fields section)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/static/js/verification.js` (existing)

### Authentication-Related Files
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/app.py` (route handlers)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/utils.py` (authentication utils)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/jwt_auth.py` (authentication)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/static/js/auth.js` (frontend authentication)

### API-Related Files
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/__init__.py` (API routes)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/resources.py` (API resources)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/api_standards.py` (API standards)

## Success Criteria

1. Lab verification workflow is fully functional:
   - API endpoints are properly registered
   - Verification data can be created, retrieved, and updated
   - Verification statistics are accessible

2. Molecules and mixtures pages are accessible without login:
   - Users can browse and view molecules/mixtures without authentication
   - The application still prompts for login when attempting to modify data
   - The header shows login/logout options

3. API implementation is complete and consistent:
   - All endpoints are properly registered
   - All endpoints use standardized response format
   - All endpoints validate input properly
   - All endpoints are documented

## Testing Instructions

1. Test lab verification workflow:
   ```bash
   python tests/run_tests.py -t test_lab_verification.py
   ```

2. Test authentication changes:
   - Manual testing with browser in incognito mode
   - Verify read-only access without authentication
   - Verify write operations require authentication

3. Test API endpoints:
   ```bash
   python tests/run_supabase_api_tests.py
   ```

## Next Steps After Completion

After completing Phase 1, we'll move to **Phase 2: Minimal Viable Deployment** which focuses on:

1. Essential Maintenance Utilities
2. Predictive Models Completion
3. Basic CI/CD Pipeline
4. Minimal Monitoring

Please provide a detailed completion report when Phase 1 tasks are finished.