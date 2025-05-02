#!/usr/bin/env python3
"""
CryoProtect v2 - Update API Resources

This script updates the API resource classes to work with both singular and plural table names
during the database migration. It modifies the existing API resources to use the compatibility
layer, ensuring backward compatibility while encouraging migration to the new schema.

Tables affected:
1. molecule → molecules
2. mixture → mixtures
3. experiment → experiments
4. prediction → predictions
5. project → projects

Author: CryoProtect Team
Date: 2025-04-18
"""

import os
import sys
import re
import logging
import argparse
import importlib
import inspect
from datetime import datetime
from typing import Dict, List, Any, Tuple, Optional
import shutil

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"update_api_resources_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Table mappings
TABLE_MAPPINGS = {
    "molecule": "molecules",
    "mixture": "mixtures",
    "experiment": "experiments",
    "prediction": "predictions",
    "project": "projects"
}

# Files to update
API_FILES = [
    "api/resources.py",
    "api/models.py",
    "api/__init__.py",
    "api/batch_resources.py",
    "api/mixture_analysis_resources.py",
    "api/predictive_models_resources.py",
    "api/protocol_designer_resources.py",
    "api/rdkit_resources.py",
    "api/scoring_resources.py",
    "api/team_resources.py",
    "api/user_profile_resources.py",
    "api/export_api_resources.py"
]

def backup_file(file_path: str) -> str:
    """
    Create a backup of a file.
    
    Args:
        file_path (str): Path to the file
        
    Returns:
        str: Path to the backup file
    """
    backup_path = f"{file_path}.bak.{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    try:
        shutil.copy2(file_path, backup_path)
        logger.info(f"Created backup: {backup_path}")
        return backup_path
    except Exception as e:
        logger.error(f"Error creating backup of {file_path}: {str(e)}")
        return ""

def update_imports(file_content: str) -> str:
    """
    Update import statements to include the compatibility layer.
    
    Args:
        file_content (str): Original file content
        
    Returns:
        str: Updated file content
    """
    # Check if compatibility layer is already imported
    if "from api_compatibility_layer import" in file_content:
        return file_content
    
    # Add import for compatibility layer
    import_pattern = r"(from api\.utils import.*?(?:\n|$))"
    compatibility_import = "from api_compatibility_layer import compatibility_layer, log_deprecated_endpoint_usage, inject_deprecation_notice, get_table_with_compatibility, track_deprecated_usage\n"
    
    if re.search(import_pattern, file_content):
        return re.sub(import_pattern, r"\1" + compatibility_import, file_content)
    else:
        # If the specific import pattern isn't found, add after the last import
        import_section_end = 0
        for match in re.finditer(r"^import .*?$|^from .*? import .*?$", file_content, re.MULTILINE):
            import_section_end = max(import_section_end, match.end())
        
        if import_section_end > 0:
            return file_content[:import_section_end] + "\n" + compatibility_import + file_content[import_section_end:]
        else:
            # If no imports found, add at the beginning after any docstrings
            docstring_end = file_content.find('"""', file_content.find('"""') + 3) + 3 if '"""' in file_content else 0
            return file_content[:docstring_end] + "\n\n" + compatibility_import + file_content[docstring_end:]

def update_resource_classes(file_content: str) -> str:
    """
    Update API resource classes to use the compatibility layer.
    
    Args:
        file_content (str): Original file content
        
    Returns:
        str: Updated file content
    """
    # Pattern to find class definitions for resources
    class_pattern = r"(class\s+(\w+)Resource\(Resource\):.*?(?=class|\Z))"
    
    # Find all resource classes
    for match in re.finditer(class_pattern, file_content, re.DOTALL):
        class_content = match.group(1)
        class_name = match.group(2)
        
        # Check if this class is related to one of our target tables
        is_target_class = False
        for singular in TABLE_MAPPINGS.keys():
            if singular.lower() in class_name.lower():
                is_target_class = True
                break
        
        if not is_target_class:
            continue
        
        # Check if compatibility_layer decorator is already applied
        if "@compatibility_layer" in class_content:
            continue
        
        # Add compatibility_layer decorator to methods
        updated_class_content = class_content
        
        # Pattern to find methods in the class
        method_pattern = r"(\s+@.*?\n)?\s+def\s+(get|post|put|delete|patch)\s*\("
        
        # Add compatibility_layer decorator to methods
        for method_match in re.finditer(method_pattern, class_content):
            decorator_pos = method_match.start()
            existing_decorator = method_match.group(1)
            
            # If there's already a decorator, add our decorator before it
            if existing_decorator:
                updated_class_content = (
                    updated_class_content[:decorator_pos] +
                    "    @compatibility_layer\n" +
                    updated_class_content[decorator_pos:]
                )
            else:
                # If no decorator, add ours before the method definition
                method_pos = method_match.start()
                updated_class_content = (
                    updated_class_content[:method_pos] +
                    "    @compatibility_layer\n" +
                    updated_class_content[method_pos:]
                )
        
        # Replace the original class content with the updated one
        file_content = file_content.replace(class_content, updated_class_content)
    
    return file_content

def update_table_references(file_content: str) -> str:
    """
    Update direct table name references to use the compatibility function.
    
    Args:
        file_content (str): Original file content
        
    Returns:
        str: Updated file content
    """
    updated_content = file_content
    
    # Pattern to find Supabase table queries
    table_query_pattern = r"supabase\.table\(['\"](\w+)['\"]\)"
    
    # Replace with get_table_with_compatibility function
    for match in re.finditer(table_query_pattern, file_content):
        table_name = match.group(1)
        if table_name in TABLE_MAPPINGS or table_name in TABLE_MAPPINGS.values():
            original = match.group(0)
            replacement = f"supabase.table(get_table_with_compatibility('{table_name}'))"
            updated_content = updated_content.replace(original, replacement)
    
    # Pattern to find Supabase from_ queries
    from_query_pattern = r"supabase\.from_\(['\"](\w+)['\"]\)"
    
    # Replace with get_table_with_compatibility function
    for match in re.finditer(from_query_pattern, file_content):
        table_name = match.group(1)
        if table_name in TABLE_MAPPINGS or table_name in TABLE_MAPPINGS.values():
            original = match.group(0)
            replacement = f"supabase.from_(get_table_with_compatibility('{table_name}'))"
            updated_content = updated_content.replace(original, replacement)
    
    return updated_content

def add_deprecation_tracking(file_content: str) -> str:
    """
    Add deprecation tracking to API endpoints.
    
    Args:
        file_content (str): Original file content
        
    Returns:
        str: Updated file content
    """
    updated_content = file_content
    
    # Pattern to find methods in resource classes
    method_pattern = r"(\s+def\s+(get|post|put|delete|patch)\s*\(self(?:,\s*[\w_]+)*\):.*?)(return.*?(?:\n\s+\w|\n\s*$|\n\s*except))"
    
    # Add tracking code before return statements
    for match in re.finditer(method_pattern, file_content, re.DOTALL):
        method_body = match.group(1)
        return_stmt = match.group(3)
        
        # Check if this is a singular endpoint by looking at the class and method
        context = file_content[max(0, match.start() - 200):match.start()]
        class_match = re.search(r"class\s+(\w+)Resource", context)
        
        if class_match:
            class_name = class_match.group(1)
            
            # Check if this class is related to one of our target tables
            is_target_class = False
            for singular in TABLE_MAPPINGS.keys():
                if singular.lower() in class_name.lower():
                    is_target_class = True
                    break
            
            if is_target_class:
                # Add tracking code before the return statement
                tracking_code = "\n        # Track usage of deprecated endpoints\n        if is_singular_endpoint(request.path):\n            track_deprecated_usage(request.path, get_user_id() if hasattr(g, 'user_id') else None)\n        "
                
                # Insert tracking code before the return statement
                updated_return = tracking_code + return_stmt
                updated_method = method_body + updated_return
                
                # Replace in the updated content
                updated_content = updated_content.replace(method_body + return_stmt, updated_method)
    
    return updated_content

def update_model_classes(file_content: str) -> str:
    """
    Update model classes to use the compatibility layer.
    
    Args:
        file_content (str): Original file content
        
    Returns:
        str: Updated file content
    """
    updated_content = file_content
    
    # Pattern to find table_name class variables
    table_name_pattern = r"(class\s+(\w+)\(BaseModel\):.*?table_name\s*=\s*['\"](\w+)['\"])"
    
    # Update table_name to use get_table_with_compatibility
    for match in re.finditer(table_name_pattern, file_content, re.DOTALL):
        class_def = match.group(1)
        class_name = match.group(2)
        table_name = match.group(3)
        
        if table_name in TABLE_MAPPINGS:
            # Replace the static table_name with a property
            updated_class_def = class_def.replace(
                f"table_name = '{table_name}'",
                f"_table_name = '{table_name}'\n    \n    @property\n    def table_name(self):\n        return get_table_with_compatibility(self._table_name)"
            )
            updated_content = updated_content.replace(class_def, updated_class_def)
    
    return updated_content

def update_api_file(file_path: str, dry_run: bool = False) -> bool:
    """
    Update an API file to use the compatibility layer.
    
    Args:
        file_path (str): Path to the file
        dry_run (bool): If True, don't write changes to file
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        logger.info(f"Updating file: {file_path}")
        
        # Read the file
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Create a backup
        if not dry_run:
            backup_file(file_path)
        
        # Apply updates
        updated_content = content
        updated_content = update_imports(updated_content)
        updated_content = update_resource_classes(updated_content)
        updated_content = update_table_references(updated_content)
        updated_content = add_deprecation_tracking(updated_content)
        
        # For models.py, update model classes
        if file_path.endswith('models.py'):
            updated_content = update_model_classes(updated_content)
        
        # Check if any changes were made
        if content == updated_content:
            logger.info(f"No changes needed for {file_path}")
            return True
        
        # Write the updated content
        if not dry_run:
            with open(file_path, 'w') as f:
                f.write(updated_content)
            logger.info(f"Updated {file_path}")
        else:
            logger.info(f"Would update {file_path} (dry run)")
        
        return True
    except Exception as e:
        logger.error(f"Error updating {file_path}: {str(e)}")
        return False

def update_api_init(file_path: str, dry_run: bool = False) -> bool:
    """
    Update the API initialization file to register both singular and plural endpoints.
    
    Args:
        file_path (str): Path to the __init__.py file
        dry_run (bool): If True, don't write changes to file
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        logger.info(f"Updating API initialization: {file_path}")
        
        # Read the file
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Create a backup
        if not dry_run:
            backup_file(file_path)
        
        # Add imports
        updated_content = update_imports(content)
        
        # Pattern to find resource registrations
        registration_pattern = r"api\.add_resource\((\w+Resource),\s*['\"]\/(\w+)(?:\/.*?)?['\"](,\s*endpoint=['\"]\w+['\"])?\)"
        
        # Find all resource registrations
        for match in re.finditer(registration_pattern, updated_content):
            resource_class = match.group(1)
            endpoint_name = match.group(2)
            endpoint_param = match.group(3) or ""
            
            # Check if this is a plural endpoint that needs a singular alias
            if endpoint_name in TABLE_MAPPINGS.values():
                singular_name = next((k for k, v in TABLE_MAPPINGS.items() if v == endpoint_name), None)
                
                if singular_name:
                    # Create a new registration for the singular endpoint
                    original = match.group(0)
                    singular_endpoint = f"/{singular_name}"
                    
                    if "/<" in original:
                        # Handle endpoints with parameters
                        param_part = original.split("/<", 1)[1].split("'", 1)[0]
                        singular_endpoint = f"/{singular_name}/<{param_part}"
                    
                    # Add the new registration after the original one
                    new_registration = f"\n# Register compatibility endpoint for singular table name\napi.add_resource({resource_class}, '{singular_endpoint}', endpoint='{singular_name}_compat')"
                    
                    # Add the new registration after the original one
                    updated_content = updated_content.replace(original, original + new_registration)
        
        # Check if any changes were made
        if content == updated_content:
            logger.info(f"No changes needed for {file_path}")
            return True
        
        # Write the updated content
        if not dry_run:
            with open(file_path, 'w') as f:
                f.write(updated_content)
            logger.info(f"Updated {file_path}")
        else:
            logger.info(f"Would update {file_path} (dry run)")
        
        return True
    except Exception as e:
        logger.error(f"Error updating {file_path}: {str(e)}")
        return False

def update_all_api_files(dry_run: bool = False) -> Dict[str, Any]:
    """
    Update all API files to use the compatibility layer.
    
    Args:
        dry_run (bool): If True, don't write changes to file
        
    Returns:
        Dict[str, Any]: Report of updated files
    """
    logger.info(f"Starting API resource update (dry_run={dry_run})")
    
    updated_files = []
    failed_files = []
    
    # Update API initialization file first
    init_file = "api/__init__.py"
    if os.path.exists(init_file):
        success = update_api_init(init_file, dry_run)
        if success:
            updated_files.append(init_file)
        else:
            failed_files.append(init_file)
    
    # Update other API files
    for file_path in API_FILES:
        if file_path == init_file:
            continue  # Already processed
        
        if os.path.exists(file_path):
            success = update_api_file(file_path, dry_run)
            if success:
                updated_files.append(file_path)
            else:
                failed_files.append(file_path)
        else:
            logger.warning(f"File not found: {file_path}")
    
    report = {
        "success": len(failed_files) == 0,
        "timestamp": datetime.now().isoformat(),
        "updated_files": updated_files,
        "failed_files": failed_files,
        "dry_run": dry_run
    }
    
    # Save report to file
    if not dry_run:
        report_file = f"update_api_resources_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        import json
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {report_file}")
    
    return report

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Update API resources for CryoProtect v2")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without writing changes")
    args = parser.parse_args()
    
    report = update_all_api_files(dry_run=args.dry_run)
    
    if report["success"]:
        logger.info("Successfully updated API resources")
        if args.dry_run:
            logger.info("Dry run completed successfully, no changes were written")
    else:
        logger.error("Failed to update some API resources")
        for file_path in report["failed_files"]:
            logger.error(f"Failed to update: {file_path}")
    
    return 0 if report["success"] else 1

if __name__ == "__main__":
    sys.exit(main())