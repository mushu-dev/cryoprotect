#!/usr/bin/env python3
"""
CryoProtect v2 - Patch App for Service Role Authentication

This script patches the app.py file to add service role authentication.
"""

import os
import sys

# Check if app.py exists
if not os.path.exists("app.py"):
    print("Error: app.py file not found")
    sys.exit(1)

# Check if auth_config.py exists
if not os.path.exists("auth_config.py"):
    print("Error: auth_config.py not found. Run fix_auth_simple.py first.")
    sys.exit(1)

# Create a backup of the app.py file
with open("app.py", "r") as f:
    original_content = f.read()

with open("app.py.bak", "w") as f:
    f.write(original_content)
    print("Created backup of app.py as app.py.bak")

# Define the patch to add service role authentication
patch = """
# Service role authentication - Added by patch_app_with_service_role.py
try:
    from auth_config import USER_ID, USE_SERVICE_ROLE
except ImportError:
    USER_ID = "748b5eb7-15dd-4019-b128-ae9d80d9d446"  # ID from user creation
    USE_SERVICE_ROLE = True

# Service role authentication helper
def get_service_role_auth():
    \"\"\"Get a fake user for service role authentication.\"\"\"
    if USE_SERVICE_ROLE:
        # Create a fake user object for service role authentication
        class FakeUser:
            def __init__(self):
                self.id = USER_ID
                self.email = os.environ.get("SUPABASE_USER", "user@example.com")
                self.user_metadata = {"display_name": "CryoProtect User"}
                
        return FakeUser()
    else:
        return None
"""

# Define the patch for authenticate_user function
authenticate_user_patch = """
def authenticate_user():
    \"\"\"
    Authenticate user with Supabase.
    
    Returns:
        User object if authenticated, None otherwise.
    \"\"\"
    # Service role authentication bypass
    if USE_SERVICE_ROLE:
        return get_service_role_auth()
        
    # Original authentication code
    supabase = get_supabase_client()
    try:
        # Check if we have a valid session
        session = supabase.auth.get_session()
        if session:
            return session.user
    except Exception as e:
        app.logger.error(f"Error authenticating user: {str(e)}")
    
    # Try to authenticate with credentials from request
    try:
        if request.method == 'POST' and request.content_type == 'application/json':
            data = request.get_json()
            email = data.get('email')
            password = data.get('password')
            
            if email and password:
                response = supabase.auth.sign_in_with_password({
                    'email': email,
                    'password': password
                })
                
                if response and hasattr(response, 'user'):
                    return response.user
    except Exception as e:
        app.logger.error(f"Error authenticating with credentials: {str(e)}")
    
    return None
"""

# Find and append our imports and helpers
modified_content = original_content

# Add imports after other imports
import_section_end = modified_content.find("logger = logging.getLogger(__name__)")
if import_section_end != -1:
    modified_content = modified_content[:import_section_end] + patch + "\n\n" + modified_content[import_section_end:]
else:
    # Insert at the top if we can't find a suitable position
    modified_content = patch + "\n\n" + modified_content

# Replace authenticate_user function if it's imported from api.utils
if "from api.utils import get_supabase_client, authenticate_user" in modified_content:
    print("Found authenticate_user import. Need to modify api/utils.py instead.")
    
    # Check if the utils.py file exists
    if os.path.exists("api/utils.py"):
        with open("api/utils.py", "r") as f:
            utils_content = f.read()
        
        with open("api/utils.py.bak", "w") as f:
            f.write(utils_content)
            print("Created backup of api/utils.py as api/utils.py.bak")
        
        # Find the authenticate_user function in utils.py
        authenticate_user_start = utils_content.find("def authenticate_user()")
        if authenticate_user_start != -1:
            # Find the end of the function
            function_body_start = utils_content.find(":", authenticate_user_start)
            indentation = utils_content[utils_content.rfind("\n", 0, authenticate_user_start) + 1:authenticate_user_start].count(" ")
            
            # Get the function body
            function_body_end = authenticate_user_start
            depth = 0
            i = function_body_start + 1
            while i < len(utils_content):
                if utils_content[i:i+1] == "\n":
                    next_line_start = i + 1
                    next_line_end = utils_content.find("\n", next_line_start)
                    if next_line_end == -1:
                        next_line_end = len(utils_content)
                    
                    next_line = utils_content[next_line_start:next_line_end]
                    
                    if next_line.strip() and next_line.count(" ") <= indentation:
                        function_body_end = i
                        break
                
                i += 1
            
            if function_body_end == authenticate_user_start:
                function_body_end = len(utils_content)
            
            # Replace the function
            service_role_patch = """
def authenticate_user():
    \"\"\"
    Authenticate user with Supabase.
    
    Returns:
        User object if authenticated, None otherwise.
    \"\"\"
    # Service role authentication bypass
    try:
        from auth_config import USER_ID, USE_SERVICE_ROLE
        if USE_SERVICE_ROLE:
            # Create a fake user object for service role authentication
            class FakeUser:
                def __init__(self):
                    self.id = USER_ID
                    self.email = os.environ.get("SUPABASE_USER", "user@example.com")
                    self.user_metadata = {"display_name": "CryoProtect User"}
                    
            return FakeUser()
    except ImportError:
        pass
        
    # Original authentication code
"""
            
            # Add our patch
            indentation_spaces = " " * indentation
            service_role_patch_with_indent = service_role_patch.replace("\n", "\n" + indentation_spaces)
            
            # Get the original function implementation but skip the first function definition line
            original_function_body = utils_content[function_body_start+1:function_body_end]
            
            # Replace the function with our patched version
            new_utils_content = utils_content[:authenticate_user_start] + service_role_patch_with_indent.strip() + original_function_body + utils_content[function_body_end:]
            
            # Write the modified file
            with open("api/utils.py", "w") as f:
                f.write(new_utils_content)
            
            print("Updated authenticate_user function in api/utils.py")
        else:
            print("Could not find authenticate_user function in api/utils.py")
    else:
        print("Could not find api/utils.py file")

# Handle the check at the end of app.py
if "__name__ == '__main__'" in modified_content:
    # Find the authentication check
    auth_check_start = modified_content.find("# Try to authenticate with Supabase")
    if auth_check_start != -1:
        # Find the end of the authentication check
        auth_check_end = modified_content.find("# Run the application", auth_check_start)
        if auth_check_end != -1:
            # Replace the authentication check
            new_auth_check = """    # Try to authenticate with Supabase
    with app.app_context():
        if USE_SERVICE_ROLE:
            user = get_service_role_auth()
            app.logger.info(f"Using service role authentication with user ID: {USER_ID}")
        else:
            user = authenticate_user()
            if user:
                app.logger.info(f"Authenticated as {user.email}")
            else:
                app.logger.warning("No authentication. Some operations may fail due to Row Level Security (RLS) policies.")
    
"""
            modified_content = modified_content[:auth_check_start] + new_auth_check + modified_content[auth_check_end:]
            print("Updated authentication check at app startup")

# Write the modified content back to app.py
with open("app.py", "w") as f:
    f.write(modified_content)

print("\nUpdated app.py with service role authentication")
print("\nNext steps:")
print("1. Run the test script to verify service role authentication:")
print("   python test_service_role_auth.py")
print("2. If the test is successful, run the application:")
print("   python run_app_with_fix.bat")
