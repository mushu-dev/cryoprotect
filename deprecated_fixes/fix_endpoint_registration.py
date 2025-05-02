#!/usr/bin/env python3
"""
CryoProtect v2 - Fix Endpoint Registration

This script fixes the duplicate endpoint registration issue in api/__init__.py
by modifying the way endpoint aliases are registered.

Usage:
    python fix_endpoint_registration.py
"""

import os
import sys
import re
import logging
import shutil
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("fix_endpoint_registration.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def backup_file(file_path):
    """Create a backup of the file."""
    backup_path = f"{file_path}.bak.{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    try:
        shutil.copy2(file_path, backup_path)
        logger.info(f"Created backup: {backup_path}")
        return True
    except Exception as e:
        logger.error(f"Error creating backup: {e}")
        return False

def fix_endpoint_registration(file_path):
    """Fix the endpoint registration in api/__init__.py."""
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return False
    
    # Create a backup
    if not backup_file(file_path):
        logger.error("Failed to create backup. Aborting.")
        return False
    
    try:
        # Read the file
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Find the problematic lines
        pattern = r"api\.add_resource\(ComparisonResource, '/mixtures/<string:mixture_id>/comparisons'\)\s*api\.add_resource\(ComparisonResource, '/mixtures/<string:mixture_id>/compare'\)"
        
        # Replace with the fixed version using add_url_rule
        replacement = """api.add_resource(ComparisonResource, '/mixtures/<string:mixture_id>/comparisons', endpoint='comparison_resource')
# Add an alias for the comparison endpoint using add_url_rule instead of add_resource
api_bp.add_url_rule('/api/v1/mixtures/<string:mixture_id>/compare', view_func=ComparisonResource.as_view('comparison_alias'))"""
        
        # Check if the pattern exists
        if re.search(pattern, content, re.DOTALL):
            # Replace the pattern
            new_content = re.sub(pattern, replacement, content, flags=re.DOTALL)
            
            # Write the modified content back to the file
            with open(file_path, 'w') as f:
                f.write(new_content)
            
            logger.info(f"Successfully fixed endpoint registration in {file_path}")
            return True
        else:
            # Try a more flexible approach
            pattern1 = r"api\.add_resource\(ComparisonResource, '/mixtures/<string:mixture_id>/comparisons'\)"
            pattern2 = r"api\.add_resource\(ComparisonResource, '/mixtures/<string:mixture_id>/compare'\)"
            
            if re.search(pattern1, content) and re.search(pattern2, content):
                # Remove the second line
                new_content = re.sub(pattern2, "", content)
                
                # Add the new line after the first pattern
                new_content = re.sub(
                    pattern1, 
                    "api.add_resource(ComparisonResource, '/mixtures/<string:mixture_id>/comparisons', endpoint='comparison_resource')\n# Add an alias for the comparison endpoint using add_url_rule instead of add_resource\napi_bp.add_url_rule('/api/v1/mixtures/<string:mixture_id>/compare', view_func=ComparisonResource.as_view('comparison_alias'))", 
                    new_content
                )
                
                # Write the modified content back to the file
                with open(file_path, 'w') as f:
                    f.write(new_content)
                
                logger.info(f"Successfully fixed endpoint registration in {file_path} using flexible approach")
                return True
            else:
                logger.error(f"Could not find the expected patterns in {file_path}")
                return False
    except Exception as e:
        logger.error(f"Error fixing endpoint registration: {e}")
        return False

def main():
    """Main function to fix endpoint registration."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Fix Endpoint Registration")
    print("=" * 80)
    
    file_path = "api/__init__.py"
    
    print(f"\nFixing endpoint registration in {file_path}...")
    if fix_endpoint_registration(file_path):
        print(f"\nSuccessfully fixed endpoint registration in {file_path}")
        print("\nThe fix addresses the following issue:")
        print("  - Duplicate endpoint registration for '/mixtures/<string:mixture_id>/compare'")
        print("  - Error: 'View function mapping is overwriting an existing endpoint function'")
        print("\nThe fix uses Flask's add_url_rule instead of add_resource for the alias.")
        print("This allows multiple URL patterns to point to the same resource without conflicts.")
        
        print("\nNext steps:")
        print("1. Run the API verification script to confirm the fix")
        print("2. Test the API endpoints to ensure they work correctly")
        print("3. Update documentation to reflect the changes")
        
        return 0
    else:
        print(f"\nFailed to fix endpoint registration in {file_path}")
        print("Please check the logs for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())