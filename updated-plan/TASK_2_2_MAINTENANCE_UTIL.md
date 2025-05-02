# Task 2.2: Implement Authentication Fix Module

## Objective
Implement the Authentication Fix module in the maintenance utility to enable automated fixing of authentication issues.

## Context
Authentication is a critical component of the CryoProtect application. We have various authentication-related scripts (`fix_auth_service_role.py`, `fix_auth_simple.py`, etc.) that need to be consolidated into a cohesive module within the maintenance utility. The current service role implementation appears to be a workaround and needs a more robust solution.

## Acceptance Criteria
- The `fix_auth_service_role` function in the maintenance utility is fully implemented
- Proper service role authentication is configured for appropriate scenarios
- Authentication-related configuration files are properly updated
- Security best practices are followed
- The fix can be run via the command line and interactive menu

## Implementation Steps

1. Examine the existing authentication fix scripts:
   ```bash
   cat fix_auth_service_role.py
   cat fix_auth_simple.py
   ```

2. Implement the `fix_auth_service_role` function in `maintenance_utils.py`:
   ```python
   def fix_auth_service_role(args):
       """
       Fix authentication issues related to service role.
       
       This function:
       1. Updates authentication configuration
       2. Implements proper service role authentication
       3. Ensures RLS policies work correctly with service role
       
       Args:
           args: Command line arguments (may include dry_run, verify, etc.)
       
       Returns:
           bool: True if successful, False otherwise
       """
       logger.info("Running Auth Service Role Fix...")
       
       try:
           # Backup current authentication configuration
           auth_config_path = "auth_config.py"
           if os.path.exists(auth_config_path):
               backup_file(auth_config_path)
           
           # Update or create auth_config.py
           with open(auth_config_path, 'w') as f:
               f.write("""\"\"\"
               Authentication Configuration for CryoProtect v2
               
               This module provides configuration for authentication mechanisms,
               particularly for service role authentication.
               \"\"\"
               
               import os
               
               # Service role configuration
               USER_ID = os.environ.get('SERVICE_ROLE_USER_ID', '748b5eb7-15dd-4019-b128-ae9d80d9d446')
               USE_SERVICE_ROLE = os.environ.get('USE_SERVICE_ROLE', 'true').lower() in ('true', 'yes', '1')
               
               # Authentication modes
               AUTH_MODE = os.environ.get('AUTH_MODE', 'service_role')  # Options: 'service_role', 'user', 'anon'
               
               # JWT configuration
               JWT_SECRET = os.environ.get('JWT_SECRET', 'your-secret-key-here')
               JWT_EXPIRATION = int(os.environ.get('JWT_EXPIRATION', '3600'))  # 1 hour default
               
               # Service role credentials - secure handling
               SERVICE_ROLE_KEY = os.environ.get('SUPABASE_SERVICE_ROLE_KEY', '')
               
               def get_auth_config():
                   \"\"\"Get the current authentication configuration.\"\"\"
                   return {
                       'user_id': USER_ID,
                       'use_service_role': USE_SERVICE_ROLE,
                       'auth_mode': AUTH_MODE
                   }
               """)
           
           # Create or update service_role_helper.py
           service_helper_path = "service_role_helper.py"
           if os.path.exists(service_helper_path):
               backup_file(service_helper_path)
           
           with open(service_helper_path, 'w') as f:
               f.write("""\"\"\"
               Service Role Helper for CryoProtect v2
               
               This module provides utilities for service role authentication.
               \"\"\"
               
               import os
               from auth_config import USER_ID, USE_SERVICE_ROLE, SERVICE_ROLE_KEY
               
               def get_service_role_auth():
                   \"\"\"
                   Get a service role authentication object.
                   
                   Returns:
                       object: A user-like object for service role authentication
                   \"\"\"
                   if USE_SERVICE_ROLE:
                       # Create a service role authentication object
                       class ServiceRoleAuth:
                           def __init__(self):
                               self.id = USER_ID
                               self.email = os.environ.get("SUPABASE_USER", "service@cryoprotect.com")
                               self.user_metadata = {"display_name": "CryoProtect Service"}
                               self.role = "service_role"
                               
                       return ServiceRoleAuth()
                   else:
                       return None
                       
               def get_service_role_headers():
                   \"\"\"
                   Get headers for service role authentication.
                   
                   Returns:
                       dict: Headers including service role key if enabled
                   \"\"\"
                   headers = {
                       "Content-Type": "application/json"
                   }
                   
                   if USE_SERVICE_ROLE and SERVICE_ROLE_KEY:
                       headers["Authorization"] = f"Bearer {SERVICE_ROLE_KEY}"
                       
                   return headers
               """)
           
           # Update app.py imports if needed
           app_path = "app.py"
           if os.path.exists(app_path):
               with open(app_path, 'r') as f:
                   app_content = f.read()
               
               # Check if auth import exists and update if needed
               if "from auth_config import" not in app_content:
                   app_content = app_content.replace(
                       "# Service role authentication - Added by patch_app_with_service_role.py",
                       "# Service role authentication - Updated by maintenance utility\nfrom auth_config import USER_ID, USE_SERVICE_ROLE"
                   )
                   
                   # Remove old try/except block if it exists
                   import_block = app_content.find("try:")
                   import_end = app_content.find("except ImportError:")
                   if import_block != -1 and import_end != -1:
                       end_block = app_content.find("\n", app_content.find("\n", import_end))
                       if end_block != -1:
                           app_content = app_content[:import_block] + app_content[end_block+1:]
                   
                   with open(app_path, 'w') as f:
                       f.write(app_content)
           
           # Update .env.template with authentication variables
           env_template_path = ".env.template"
           if os.path.exists(env_template_path):
               with open(env_template_path, 'r') as f:
                   env_content = f.read()
               
               # Add authentication variables if not present
               auth_vars = """
               # Authentication Configuration
               SERVICE_ROLE_USER_ID=748b5eb7-15dd-4019-b128-ae9d80d9d446
               USE_SERVICE_ROLE=true
               AUTH_MODE=service_role
               JWT_SECRET=your-secret-key-here
               JWT_EXPIRATION=3600
               """
               
               if "Authentication Configuration" not in env_content:
                   with open(env_template_path, 'a') as f:
                       f.write(auth_vars)
           
           logger.info("Auth Service Role Fix completed successfully!")
           return True
       except Exception as e:
           logger.error(f"Auth Service Role Fix failed: {str(e)}", exc_info=True)
           return False
   ```

3. Update the function mapping in the maintenance utility:
   ```python
   # Ensure the function is listed in the FIX_FUNCTIONS dictionary
   FIX_FUNCTIONS = {
       # Existing mappings...
       "auth_service_role": fix_auth_service_role,
       # Other mappings...
   }
   ```

4. Add verification function:
   ```python
   def verify_auth_service_role():
       """Verify that service role authentication is properly configured."""
       try:
           # Check auth_config.py exists
           if not os.path.exists("auth_config.py"):
               logger.error("auth_config.py does not exist")
               return False
           
           # Check service_role_helper.py exists
           if not os.path.exists("service_role_helper.py"):
               logger.error("service_role_helper.py does not exist")
               return False
           
           # Import and test authentication functions
           try:
               from auth_config import get_auth_config
               config = get_auth_config()
               if not isinstance(config, dict) or 'use_service_role' not in config:
                   logger.error("Authentication configuration is invalid")
                   return False
           except ImportError:
               logger.error("Could not import authentication configuration")
               return False
           
           # Test service role helper
           try:
               from service_role_helper import get_service_role_auth
               auth = get_service_role_auth()
               if USE_SERVICE_ROLE and (auth is None or not hasattr(auth, 'id')):
                   logger.error("Service role authentication object is invalid")
                   return False
           except ImportError:
               logger.error("Could not import service role helper")
               return False
           
           logger.info("Service role authentication verification passed!")
           return True
       except Exception as e:
           logger.error(f"Service role authentication verification failed: {str(e)}", exc_info=True)
           return False
   ```

5. Add verification to the main verification function:
   ```python
   def verify_all_fixes():
       """Run all verification functions."""
       # Existing verifications...
       success = success and verify_auth_service_role()
       # More verifications...
       return success
   ```

## Files to Modify
- `maintenance_utils.py` - Add the authentication fix implementation
- `auth_config.py` - Will be created or modified by the fix
- `service_role_helper.py` - Will be created or modified by the fix
- Potentially `app.py` - May need updates to imports

## Verification
1. Run the utility with the authentication fix:
   ```bash
   python maintenance_utils.py --fix auth_service_role
   ```

2. Check that the fix succeeds without errors

3. Verify the authentication configuration:
   ```bash
   python maintenance_utils.py --verify
   ```

4. Test an API endpoint that requires authentication to ensure it works with the service role:
   ```bash
   curl -H "Authorization: Bearer $SUPABASE_SERVICE_ROLE_KEY" http://localhost:5000/api/v1/molecules
   ```

## Notes for Roo Code Agent
- This implementation creates a more robust and configurable authentication system
- The service role key should be handled securely, preferably through environment variables
- Make sure to test thoroughly as authentication is a critical security component
- Consider implementing a simple test script to verify the authentication works without exposing sensitive keys
- Add appropriate comments explaining the security considerations