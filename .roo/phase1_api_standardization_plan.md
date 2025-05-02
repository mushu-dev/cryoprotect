# Phase 1: API Standardization Plan

## Objective
Standardize all API resource files in the CryoProtect v2 project to ensure consistent error handling, authentication, response formatting, and documentation.

## Current Status
- Project is ~35-40% complete overall
- API Standardization is ~10% complete (1 of ~10 resource files)
- `api/dashboard_resources.py` has been successfully standardized as an example

## Implementation Pattern
Each API resource file should be updated to include:

1. **Consistent Error Handling**
   - Use the `handle_error` utility for all error responses
   - Example from dashboard_resources.py:
   ```python
   try:
       # Operation code
   except Exception as e:
       return handle_error(str(e), 500)
   ```

2. **Authentication Enforcement**
   - Use `@token_required` decorator on all protected endpoints
   - Example from dashboard_resources.py:
   ```python
   @token_required
   def get(self, user=None):
       # Endpoint code
   ```

3. **Response Formatting**
   - Use `@marshal_with` with defined field schemas
   - Define field schemas at the top of each resource file
   - Example from dashboard_resources.py:
   ```python
   dashboard_fields = {
       'id': fields.String,
       'name': fields.String,
       'data': fields.Raw,
       'created_at': fields.DateTime
   }

   class DashboardResource(Resource):
       @marshal_with(dashboard_fields)
       def get(self):
           # Endpoint code
   ```

4. **Request Validation**
   - Validate incoming request data before processing
   - Check for required fields and validate data types
   - Example from dashboard_resources.py:
   ```python
   parser = reqparse.RequestParser()
   parser.add_argument('name', type=str, required=True)
   parser.add_argument('data', type=dict, required=True)
   
   # In the post method:
   args = parser.parse_args()
   ```

5. **Comprehensive Docstrings**
   - Add detailed docstrings for each resource class and method
   - Include parameter descriptions and return values
   - Example from dashboard_resources.py:
   ```python
   class DashboardResource(Resource):
       """Resource for managing dashboards.
       
       Provides endpoints for creating, retrieving, updating and deleting dashboards.
       """
       
       def get(self):
           """Retrieve a dashboard by ID.
           
           Returns:
               Dashboard object if found, otherwise 404 error
           """
   ```

## Implementation Order
Standardize API resource files in this sequence:

1. `api/resources.py` (Core API functionality)
2. `api/mixture_analysis_resources.py` (Scientific functionality)
3. `api/rdkit_resources.py` (Chemical processing)
4. `api/predictive_models_resources.py` (ML interfaces)
5. `api/scoring_resources.py`
6. `api/export_resources.py`
7. `api/team_resources.py`
8. `api/protocol_designer_resources.py`
9. Remaining resource files

## Implementation Approach for Each File

1. **Analyze the File**
   - Examine each endpoint in the resource file
   - Identify existing error handling, authentication, and response formatting
   - Note any custom patterns that need to be preserved

2. **Update Error Handling**
   - Replace custom error responses with `handle_error`
   - Ensure appropriate HTTP status codes
   - Maintain any specific error information needed

3. **Implement Authentication**
   - Add `@token_required` decorator to protected endpoints
   - Ensure the `user` parameter is properly utilized
   - Verify authorization checks within endpoints

4. **Standardize Response Formatting**
   - Define field schemas for all response types
   - Apply `@marshal_with` to appropriate methods
   - Ensure consistent response structure

5. **Add Request Validation**
   - Implement `RequestParser` for each endpoint with input
   - Add appropriate type checking and required flags
   - Validate request data before processing

6. **Update Docstrings**
   - Add class and method level docstrings
   - Document parameters, return values, and exceptions
   - Include usage examples where helpful

7. **Update Tests**
   - Modify tests to expect standardized responses
   - Ensure all tests pass with the updated code
   - Add tests for edge cases if needed

8. **Document Changes**
   - Update `README_API_FIXES.md` with changes made
   - Note any significant refactoring or pattern changes
   - Document any breaking changes to API contracts

## Verification for Each File

After standardizing each file:

1. Run relevant API tests:
   ```
   python -m unittest tests/test_api_endpoints.py
   ```
   
2. Verify endpoints still function correctly:
   ```
   python -m app.py
   # Test with API client or curl commands
   ```

3. Check for consistent patterns across resources

## Key Files Reference

### 1. `api/resources.py`
Core API functionality with base resources:
- Contains base resource classes used throughout the API
- Includes endpoint registration and routing logic
- Has fundamental CRUD operations for primary entities

### 2. `api/mixture_analysis_resources.py`
Scientific functionality for mixture analysis:
- Endpoints for analyzing cryoprotectant mixtures
- Includes calculation endpoints for scientific properties
- Contains visualization data generation functions

### 3. `api/rdkit_resources.py`
Chemical processing with RDKit integration:
- Molecular structure handling endpoints
- Chemical property calculation functions
- Structure conversion and visualization endpoints

### 4. `api/predictive_models_resources.py`
Machine learning interfaces:
- Prediction generation endpoints
- Model training and management functions
- Prediction evaluation and comparison endpoints

## Utility Functions Reference

### Error Handling
```python
# From api/utils.py
def handle_error(message, status_code=400, details=None):
    """Handle API errors with consistent formatting."""
    response = {
        'error': message,
        'status_code': status_code
    }
    if details:
        response['details'] = details
    return response, status_code
```

### Authentication
```python
# From api/utils.py
def token_required(f):
    """Decorator for endpoints that require authentication."""
    @wraps(f)
    def decorated(*args, **kwargs):
        token = None
        # Extract token from Authorization header
        auth_header = request.headers.get('Authorization')
        if auth_header:
            parts = auth_header.split()
            if len(parts) == 2 and parts[0].lower() == 'bearer':
                token = parts[1]
                
        if not token:
            return {'message': 'Token is missing'}, 401
            
        try:
            # Verify token and get user
            supabase = get_supabase_client()
            user = supabase.auth.get_user(token)
            
            return f(*args, **kwargs, user=user)
        except Exception as e:
            return {'message': 'Invalid token'}, 401
            
    return decorated
```

## Example: StandardizedResource

Here's an example of a fully standardized resource:

```python
from flask_restful import Resource, reqparse, fields, marshal_with
from api.utils import token_required, handle_error

# Define response fields
item_fields = {
    'id': fields.String,
    'name': fields.String,
    'description': fields.String,
    'created_at': fields.DateTime,
    'updated_at': fields.DateTime
}

class ItemResource(Resource):
    """Resource for managing items.
    
    Provides endpoints for creating, retrieving, updating and deleting items.
    """
    
    # Define request parser
    parser = reqparse.RequestParser()
    parser.add_argument('name', type=str, required=True, help='Name cannot be blank')
    parser.add_argument('description', type=str, required=False)
    
    @marshal_with(item_fields)
    @token_required
    def get(self, item_id, user=None):
        """Retrieve an item by ID.
        
        Args:
            item_id: ID of the item to retrieve
            user: User object from token_required decorator
            
        Returns:
            Item object if found
            
        Raises:
            404: If item not found
        """
        try:
            supabase = get_supabase_client()
            result = supabase.table('items').select('*').eq('id', item_id).execute()
            
            if not result.data:
                return handle_error("Item not found", 404)
                
            return result.data[0]
        except Exception as e:
            return handle_error(str(e), 500)
    
    @marshal_with(item_fields)
    @token_required
    def post(self, user=None):
        """Create a new item.
        
        Args:
            user: User object from token_required decorator
            
        Returns:
            Created item object
            
        Raises:
            400: If request data is invalid
            500: If database error occurs
        """
        try:
            args = self.parser.parse_args()
            
            supabase = get_supabase_client()
            result = supabase.table('items').insert({
                'name': args['name'],
                'description': args.get('description', ''),
                'user_id': user.id
            }).execute()
            
            return result.data[0]
        except Exception as e:
            return handle_error(str(e), 500)
```

## Documentation

After each standardization, update `README_API_FIXES.md` with:

1. File standardized
2. Endpoints updated
3. Any breaking changes
4. Patterns implemented
5. Date of implementation