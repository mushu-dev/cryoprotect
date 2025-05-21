#!/usr/bin/env python3
"""
CryoProtect v2 - Update App Authentication

This script modifies the app.py file to use the service role authentication approach.
"""

import os
import re
import sys
from dotenv import load_dotenv

print("CryoProtect v2 - Update App Authentication")
print("------------------------------------------")

# Check if app.py exists
if not os.path.exists("app.py"):
    print("Error: app.py file not found")
    sys.exit(1)

# Load environment variables
load_dotenv()

# Create a backup of the app.py file
with open("app.py", "r") as f:
    app_content = f.read()

with open("app.py.backup", "w") as f:
    f.write(app_content)
    print("Created backup of app.py as app.py.backup")

# Add the service role import at the top of the file
if "from auth_config import USER_ID, USE_SERVICE_ROLE" not in app_content:
    # Add after other imports
    import_section_end = app_content.find("\n\n", app_content.find("import"))
    if import_section_end == -1:
        import_section_end = app_content.find("\n", app_content.find("import"))
    
    # Insert our import after the import section
    modified_app_content = app_content[:import_section_end] + "\n\n# Service role authentication\ntry:\n    from auth_config import USER_ID, USE_SERVICE_ROLE\nexcept ImportError:\n    USER_ID = \"748b5eb7-15dd-4019-b128-ae9d80d9d446\"  # ID from user creation\n    USE_SERVICE_ROLE = True\n" + app_content[import_section_end:]
    
    print("Added service role import")
else:
    modified_app_content = app_content
    print("Service role import already exists")

# Modify the authentication functions if they exist
# Look for patterns that might indicate authentication functionality

# Pattern 1: Initialize Supabase client
supabase_init_pattern = r"supabase\s*=\s*create_client\(.*\)"
if re.search(supabase_init_pattern, modified_app_content):
    print("Found Supabase client initialization")
    # No need to modify this, as we still use the same initialization code

# Pattern 2: Authentication routes or functions
auth_patterns = [
    r"@app\.route\s*\(\s*['\"]\/login['\"].*\)",
    r"@app\.route\s*\(\s*['\"]\/auth.*\)",
    r"def\s+login\s*\(",
    r"def\s+authenticate\s*\(",
    r"supabase\.auth\.sign_in",
]

has_auth_code = False
for pattern in auth_patterns:
    if re.search(pattern, modified_app_content):
        has_auth_code = True
        print(f"Found authentication code matching pattern: {pattern}")

# Add our service role authentication helper
if has_auth_code:
    service_role_helper = """
# Service role authentication helper
def get_user_id_service_role():
    \"\"\"Get user ID using service role approach.\"\"\"
    return USER_ID if USE_SERVICE_ROLE else None

def authenticate_service_role(supabase_client):
    \"\"\"Authenticate using service role approach.\"\"\"
    if not USE_SERVICE_ROLE:
        return None
    
    # Return a fake session object
    class FakeUser:
        def __init__(self):
            self.id = USER_ID
    
    class FakeSession:
        def __init__(self):
            self.user = FakeUser()
    
    return FakeSession()

"""
    
    # Find a good place to insert our helper
    # Try to find after imports but before routes
    app_definition = modified_app_content.find("app = Flask(__name__)")
    if app_definition != -1:
        insert_point = modified_app_content.find("\n\n", app_definition)
        if insert_point != -1:
            modified_app_content = modified_app_content[:insert_point] + "\n" + service_role_helper + modified_app_content[insert_point:]
            print("Added service role authentication helper")
        else:
            # If we can't find a good insertion point, append it to the end
            modified_app_content += "\n" + service_role_helper
            print("Added service role authentication helper at the end")
    else:
        # If we can't find the app definition, append it to the end
        modified_app_content += "\n" + service_role_helper
        print("Added service role authentication helper at the end")

    # Now find authentication code and modify it to use our helper
    auth_modified = False
    
    # Pattern for user ID retrieval
    user_id_pattern = r"(user_id|user\.id|session\[['"]user_id['\"]\])"
    
    # Check each line for user ID references and add bypass code
    lines = modified_app_content.split("\n")
    for i, line in enumerate(lines):
        if re.search(user_id_pattern, line) and "get_user_id_service_role()" not in line:
            # Add a comment to indicate we're modifying this line
            indentation = len(line) - len(line.lstrip())
            spaces = " " * indentation
            
            # Check if this is an assignment
            if "=" in line:
                var_name = line.split("=")[0].strip()
                new_line = f"{spaces}# Original: {line.strip()}\n{spaces}if USE_SERVICE_ROLE:\n{spaces}    {var_name} = get_user_id_service_role()\n{spaces}else:\n{spaces}    {line.strip()}"
                lines[i] = new_line
                auth_modified = True
            else:
                # It's a reference but not an assignment, more complex to modify
                # Just add a comment for now
                lines[i] = f"{spaces}# TODO: Check if service role bypass needed: {line.strip()}"
    
    # Authentication function pattern
    auth_func_pattern = r"def\s+(authenticate|login|sign_in).*:"
    for i, line in enumerate(lines):
        if re.search(auth_func_pattern, line):
            # Found an authentication function, add our bypass code
            # Find the function body indentation
            j = i + 1
            while j < len(lines) and (not lines[j].strip() or lines[j].startswith(" ") or lines[j].startswith("\t")):
                j += 1
            
            # Get the indentation of the function body
            func_body_line = None
            for k in range(i+1, j):
                if lines[k].strip():
                    func_body_line = lines[k]
                    break
            
            if func_body_line:
                indentation = len(func_body_line) - len(func_body_line.lstrip())
                spaces = " " * indentation
                
                # Add our bypass code at the beginning of the function
                bypass_code = f"{spaces}# Service role bypass\n{spaces}if USE_SERVICE_ROLE:\n{spaces}    return authenticate_service_role(supabase)"
                
                # Insert after the function definition
                lines.insert(i+1, bypass_code)
                auth_modified = True
    
    if auth_modified:
        modified_app_content = "\n".join(lines)
        print("Modified authentication code to use service role bypass")
    else:
        print("Could not find specific authentication code to modify")

# Write the modified content back to app.py
with open("app.py", "w") as f:
    f.write(modified_app_content)

print("\nUpdated app.py with service role authentication")
print("\nNotes:")
print("1. The app.py file has been backed up to app.py.backup")
print("2. Service role authentication has been added to bypass email confirmation")
print("3. You may need to manually adjust some code if the automatic modifications are not sufficient")
print("\nNext steps:")
print("1. Run the test script to verify service role authentication:")
print("   python test_service_role_auth_compat.py")
print("2. If the test is successful, try running the application:")
print("   python app.py")
