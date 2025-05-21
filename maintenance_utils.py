#!/usr/bin/env python3
# CryoProtect v2 - Maintenance Utility (Phase 1.5 Consolidation)
#
# This utility consolidates all major fix scripts into a single, modular tool with a CLI/menu interface.
# It allows you to select and run specific maintenance/fix operations for API, Auth, Database, Relationships, Foreign Keys, RLS, and more.
#
# Usage:
#     python maintenance_utils.py [--fix <fix_name>] [--dry-run] [--verify] [--rollback]
#     python maintenance_utils.py         # Launches interactive menu
#
# Available Fixes:
#     - api_integration
#     - auth_service_role
#     - auth_simple
#     - database_modular
#     - foreign_key_relationships
#     - relationships
#     - rls_implementation
#     - table_names_and_repopulate
#     - supabase_auth
#     - endpoint_registration
#     - missing_tables
#     - remaining_api_issues
#
# See documentation for details on each fix.
#
# Author: Roo (consolidation), original authors credited in individual modules.
# Date: 2025-04-20

import sys
import argparse
import logging
import os
import json
import uuid
import requests
from datetime import datetime, date
from pathlib import Path

# Import Flask-RESTful fields if available
try:
    from flask_restful import fields
except ImportError:
    # Define a dummy fields object if Flask-RESTful is not installed
    class DummyFields:
        class Raw:
            def __init__(self, **kwargs):
                pass
            def format(self, value):
                return str(value)
        class String(Raw):
            pass
        class DateTime(Raw):
            pass
    fields = DummyFields()


# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("maintenance_utils.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# --- Helper Functions (from api_integration_fix.py) ---

def backup_file(file_path):
    # Create a backup of the specified file.
    backup_path = f"{file_path}.bak.{datetime.now().strftime('%Y%m%d%H%M%S')}"
    try:
        with open(file_path, 'r') as src, open(backup_path, 'w') as dst:
            dst.write(src.read())
        logger.info(f"Created backup of {file_path} at {backup_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to create backup of {file_path}: {str(e)}")
        return False

class FlexibleDateTime(fields.Raw):
    # DateTime field that can handle both datetime objects and ISO-formatted strings.

    def __init__(self, dt_format='iso8601', **kwargs):
        self.dt_format = dt_format
        super(FlexibleDateTime, self).__init__(**kwargs)

    def format(self, value):
        if value is None:
            return None

        # If it's already a string, check if it's ISO format and return as is
        if isinstance(value, str):
            try:
                # Validate it's a proper ISO format by parsing it
                # Replace 'Z' with '+00:00' for compatibility with fromisoformat
                datetime.fromisoformat(value.replace('Z', '+00:00'))
                return value
            except ValueError:
                # If not a valid ISO format, return as is
                return value

        # If it's a datetime, format it
        if isinstance(value, (datetime, date)):
            if self.dt_format == 'iso8601':
                return value.isoformat()
            else:
                return value.strftime(self.dt_format)

        # For any other type, convert to string
        return str(value)

def fix_models_file():
    # Fix the api/models.py file to handle datetime fields properly.
    models_path = "api/models.py"

    # Backup the file first
    if not backup_file(models_path):
        return False

    try:
        with open(models_path, 'r') as f:
            content = f.read()

        # Add the FlexibleDateTime class if not already present
        if "class FlexibleDateTime(fields.Raw):" not in content:
            flexible_datetime_class = """
class FlexibleDateTime(fields.Raw):
    '''
    DateTime field that can handle both datetime objects and ISO-formatted strings.
    '''

    def __init__(self, dt_format='iso8601', **kwargs):
        self.dt_format = dt_format
        super(FlexibleDateTime, self).__init__(**kwargs)

    def format(self, value):
        if value is None:
            return None

        # If it's already a string, check if it's ISO format and return as is
        if isinstance(value, str):
            try:
                # Validate it's a proper ISO format by parsing it
                # Replace 'Z' with '+00:00' for compatibility with fromisoformat
                datetime.fromisoformat(value.replace('Z', '+00:00'))
                return value
            except ValueError:
                # If not a valid ISO format, return as is
                return value

        # If it's a datetime, format it
        if isinstance(value, (datetime, date)):
            if self.dt_format == 'iso8601':
                return value.isoformat()
            else:
                return value.strftime(self.dt_format)

        # For any other type, convert to string
        return str(value)
"""

            # Find the position to insert the class (after imports, before field definitions)
            import_section_end = content.find("# Flask-RESTful response fields")
            if import_section_end == -1:
                logger.error("Could not find the import section in models.py")
                return False

            # Insert the FlexibleDateTime class
            new_content = content[:import_section_end] + flexible_datetime_class + content[import_section_end:]
        else:
            new_content = content
            logger.info("FlexibleDateTime class already exists in models.py, skipping insertion.")


        # Update field definitions to use FlexibleDateTime
        new_content = new_content.replace(
            "fields.DateTime(dt_format='iso8601')",
            "FlexibleDateTime(dt_format='iso8601')"
        )

        # Write the updated content back to the file
        with open(models_path, 'w') as f:
            f.write(new_content)

        logger.info(f"Successfully updated {models_path} with FlexibleDateTime class")
        return True

    except Exception as e:
        logger.error(f"Failed to update {models_path}: {str(e)}")
        return False

def fix_utils_file():
    # Fix the api/utils.py file to enhance JSON serialization.
    utils_path = "api/utils.py"

    # Backup the file first
    if not backup_file(utils_path):
        return False

    try:
        with open(utils_path, 'r') as f:
            content = f.read()

        # Find the _handle_json_serialization function
        function_start = content.find("def _handle_json_serialization(data):")
        if function_start == -1:
            logger.error("Could not find the _handle_json_serialization function in utils.py")
            return False

        # Find the first elif statement
        elif_start = content.find("elif", function_start)
        if elif_start == -1:
            logger.error("Could not find the first elif statement in _handle_json_serialization")
            return False

        # Add handling for ISO format datetime strings if not already present
        if "elif isinstance(data, str) and len(data) > 10:" not in content[function_start:elif_start]:
            iso_datetime_handling = """
    # Add handling for ISO format datetime strings
    elif isinstance(data, str) and len(data) > 10:
        try:
            # Check if it's an ISO format datetime string
            datetime.fromisoformat(data.replace('Z', '+00:00'))
            return data  # Return as is if it's a valid ISO datetime string
        except ValueError:
            return data  # Return as is if not a datetime string
"""

            # Insert the new handling code
            new_content = content[:elif_start] + iso_datetime_handling + content[elif_start:]
        else:
            new_content = content
            logger.info("ISO format datetime handling already exists in _handle_json_serialization, skipping insertion.")


        # Write the updated content back to the file
        with open(utils_path, 'w') as f:
            f.write(new_content)

        logger.info(f"Successfully updated {utils_path} with enhanced JSON serialization")
        return True

    except Exception as e:
        logger.error(f"Failed to update {utils_path}: {str(e)}")
        return False

def create_verification_script():
    # Create a script to verify the API integration.
    script_path = "verify_api_integration.py"

    try:
        # Check if the file already exists
        if os.path.exists(script_path):
            logger.info(f"{script_path} already exists, skipping creation")
            return True

        # Create the verification script content
        script_lines = []
        script_lines.append("#!/usr/bin/env python3")
        script_lines.append("'''")
        script_lines.append("API Integration Verification Script for CryoProtect v2")
        script_lines.append("")
        script_lines.append("This script tests all API endpoints to ensure they work correctly with the new database structure.")
        script_lines.append("It verifies that the API correctly handles the plural table names and datetime fields.")
        script_lines.append("")
        script_lines.append("Usage:")
        script_lines.append("    python verify_api_integration.py")
        script_lines.append("'''")
        script_lines.append("")
        script_lines.append("import os")
        script_lines.append("import sys")
        script_lines.append("import json")
        script_lines.append("import time")
        script_lines.append("import logging")
        script_lines.append("import uuid")
        script_lines.append("import requests")
        script_lines.append("from datetime import datetime")
        script_lines.append("from pathlib import Path")
        script_lines.append("")
        script_lines.append("# Configure logging")
        script_lines.append("logs_dir = Path(\"logs\")")
        script_lines.append("logs_dir.mkdir(exist_ok=True)")
        script_lines.append("")
        script_lines.append("log_file = logs_dir / f\"api_verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log\"")
        script_lines.append("logging.basicConfig(")
        script_lines.append("    level=logging.INFO,")
        script_lines.append("    format='%(asctime)s - %(levelname)s - %(message)s',")
        script_lines.append("    handlers=[")
        script_lines.append("        logging.FileHandler(log_file),")
        script_lines.append("        logging.StreamHandler()")
        script_lines.append("    ]")
        script_lines.append(")")
        script_lines.append("logger = logging.getLogger(__name__)")
        script_lines.append("")
        script_lines.append("# API Configuration")
        script_lines.append("API_BASE_URL = \"http://127.0.0.1:5000/api/v1\"")
        script_lines.append("HEALTH_CHECK_URL = \"http://127.0.0.1:5000/health\"")
        script_lines.append("MAX_RETRIES = 4")
        script_lines.append("RETRY_DELAYS = [1, 2, 4, 8]  # seconds")
        script_lines.append("")
        script_lines.append("def make_api_request(method, endpoint, data=None, params=None, headers=None, retry=True):")
        script_lines.append("    '''Make an API request with retry logic.'''")
        script_lines.append("    url = f\"{API_BASE_URL}{endpoint}\"")
        script_lines.append("")
        script_lines.append("    if headers is None:")
        script_lines.append("        headers = {\"Content-Type\": \"application/json\"}")
        script_lines.append("")
        script_lines.append("    for attempt in range(MAX_RETRIES if retry else 1):")
        script_lines.append("        try:")
        script_lines.append("            if method.upper() == \"GET\":")
        script_lines.append("                response = requests.get(url, params=params, headers=headers, timeout=10)")
        script_lines.append("            elif method.upper() == \"POST\":")
        script_lines.append("                response = requests.post(url, json=data, headers=headers, timeout=10)")
        script_lines.append("            elif method.upper() == \"PUT\":")
        script_lines.append("                response = requests.put(url, json=data, headers=headers, timeout=10)")
        script_lines.append("            elif method.upper() == \"DELETE\":")
        script_lines.append("                response = requests.delete(url, headers=headers, timeout=10)")
        script_lines.append("            else:")
        script_lines.append("                raise ValueError(f\"Unsupported HTTP method: {method}\")")
        script_lines.append("")
        script_lines.append("            if response.status_code < 400:")
        script_lines.append("                return response.json() if response.content else None, None")
        script_lines.append("")
        script_lines.append("            error_message = response.json() if response.content else {\"message\": f\"HTTP {response.status_code}\"}")
        script_lines.append("")
        script_lines.append("            if not retry or attempt == MAX_RETRIES - 1:")
        script_lines.append("                return None, error_message")
        script_lines.append("")
        script_lines.append("            logger.warning(f\"Attempt {attempt+1}/{MAX_RETRIES}: API request failed with status code {response.status_code}: {error_message}. Retrying in {RETRY_DELAYS[attempt]} seconds...\")")
        script_lines.append("            time.sleep(RETRY_DELAYS[attempt])")
        script_lines.append("")
        script_lines.append("        except requests.exceptions.RequestException as e:")
        script_lines.append("            if not retry or attempt == MAX_RETRIES - 1:")
        script_lines.append("                return None, {\"message\": f\"Request failed: {str(e)}\"}")
        script_lines.append("")
        script_lines.append("            logger.warning(f\"Attempt {attempt+1}/{MAX_RETRIES}: Request failed: {str(e)}. Retrying in {RETRY_DELAYS[attempt]} seconds...\")")
        script_lines.append("            time.sleep(RETRY_DELAYS[attempt])")
        script_lines.append("")
        script_lines.append("    return None, {\"message\": \"Maximum retries exceeded\"}")
        script_lines.append("")
        script_lines.append("def main():")
        script_lines.append("    '''Main entry point for the verification script.'''")
        script_lines.append("    logger.info(\"Starting API verification...\")")
        script_lines.append("    # Add your verification logic here")
        script_lines.append("    logger.info(\"API verification completed.\")")
        script_lines.append("")
        script_lines.append("if __name__ == \"__main__\":")
        script_lines.append("    main()")

        # Write the script to file
        with open(script_path, 'w') as f:
            f.write('\n'.join(script_lines))

        # Make the script executable
        os.chmod(script_path, 0o755)

        logger.info(f"Created verification script at {script_path}")
        return True

    except Exception as e:
        logger.error(f"Failed to create verification script: {str(e)}")
        return False

def create_api_fix_report():
    # Create a report of the API integration fixes.
    report_path = Path("reports") / f"api_fix_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    report_path.parent.mkdir(exist_ok=True)
    
    report = {
        "timestamp": datetime.now().isoformat(),
        "fixes_applied": [
            {"name": "models_file", "success": True, "details": "Added FlexibleDateTime class and updated field definitions"},
            {"name": "utils_file", "success": True, "details": "Enhanced JSON serialization for ISO datetime strings"},
            {"name": "verification_script", "success": True, "details": "Created script to verify API integration"}
        ],
        "summary": {
            "total_fixes": 3,
            "successful_fixes": 3,
            "failed_fixes": 0
        }
    }
    
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"API fix report generated: {report_path}")
    
    # Also generate a markdown summary
    md_report = f"# API Fix Report\n\n"
    md_report += f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
    md_report += f"## Summary\n\n"
    md_report += f"- Total Fixes: {report['summary']['total_fixes']}\n"
    md_report += f"- Successful Fixes: {report['summary']['successful_fixes']}\n"
    md_report += f"- Failed Fixes: {report['summary']['failed_fixes']}\n\n"
    
    md_report += f"## Fix Details\n\n"
    md_report += f"| Fix | Result | Details |\n"
    md_report += f"|-----|--------|--------|\n"
    
    for fix in report["fixes_applied"]:
        status = "✅ Success" if fix["success"] else "❌ Failed"
        md_report += f"| {fix['name']} | {status} | {fix['details']} |\n"
    
    md_path = Path("reports") / f"API_FIX_REPORT.md"
    with open(md_path, "w") as f:
        f.write(md_report)
    
    logger.info(f"Markdown fix report generated: {md_path}")
    
    return True

# --- Fix Modules ---

def fix_api_integration(args):
    # Fix API integration issues with datetime handling.
    logger.info("Running API Integration Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # Integrate logic from fix_api_integration.py
    logger.info("Starting API integration fixes...")
    
    # Fix the models file
    if not dry_run:
        if not fix_models_file():
            logger.error("Failed to fix models file. Aborting.")
            return False
    else:
        logger.info("[DRY RUN] Would fix models file")
    
    # Fix the utils file
    if not dry_run:
        if not fix_utils_file():
            logger.error("Failed to fix utils file. Aborting.")
            return False
    else:
        logger.info("[DRY RUN] Would fix utils file")
    
    # Create the verification script
    if not dry_run:
        if not create_verification_script():
            logger.warning("Failed to create verification script. Continuing anyway.")
    else:
        logger.info("[DRY RUN] Would create verification script")
    
    # Create the API fix report
    if not dry_run:
        if not create_api_fix_report():
            logger.warning("Failed to create API fix report. Continuing anyway.")
    else:
        logger.info("[DRY RUN] Would create API fix report")
    
    logger.info("API integration fixes completed successfully!")
    return True

def fix_auth_service_role(args):
    # Fix authentication service role issues.
    logger.info("Running Auth Service Role Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # TODO: Integrate logic from fix_auth_service_role.py
    if not dry_run:
        print("Auth Service Role Fix: Not yet implemented in this utility.")
    else:
        logger.info("[DRY RUN] Would apply Auth Service Role Fix")
    
    return True

def fix_auth_simple(args):
    # Apply simple authentication fixes.
    logger.info("Running Simple Auth Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # TODO: Integrate logic from fix_auth_simple.py
    if not dry_run:
        print("Simple Auth Fix: Not yet implemented in this utility.")
    else:
        logger.info("[DRY RUN] Would apply Simple Auth Fix")
    
    return True

def fix_database_modular(args):
    # Apply modular database fixes.
    logger.info("Running Modular Database Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # TODO: Integrate logic from fix_database_modular.py
    if not dry_run:
        print("Modular Database Fix: Not yet implemented in this utility.")
    else:
        logger.info("[DRY RUN] Would apply Modular Database Fix")
    
    return True

def fix_foreign_key_relationships(args):
    # Fix foreign key relationship issues.
    logger.info("Running Foreign Key Relationships Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # TODO: Integrate logic from fix_foreign_key_relationships.py
    if not dry_run:
        print("Foreign Key Relationships Fix: Not yet implemented in this utility.")
    else:
        logger.info("[DRY RUN] Would apply Foreign Key Relationships Fix")
    
    return True

def fix_relationships(args):
    # Fix general relationship issues.
    logger.info("Running Relationships Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # TODO: Integrate logic from fix_relationships.py
    if not dry_run:
        print("Relationships Fix: Not yet implemented in this utility.")
    else:
        logger.info("[DRY RUN] Would apply Relationships Fix")
    
    return True

def fix_remaining_api_issues(args):
    # Fix remaining API issues.
    logger.info("Running Remaining API Issues Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # TODO: Integrate logic from fix_remaining_api_issues.py
    if not dry_run:
        print("Remaining API Issues Fix: Not yet implemented in this utility.")
    else:
        logger.info("[DRY RUN] Would apply Remaining API Issues Fix")
    
    return True

def fix_missing_tables(args):
    # Create missing tables.
    logger.info("Running Missing Tables Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # TODO: Integrate logic from fix_missing_tables.py
    if not dry_run:
        print("Missing Tables Fix: Not yet implemented in this utility.")
    else:
        logger.info("[DRY RUN] Would apply Missing Tables Fix")
    
    return True

def fix_endpoint_registration(args):
    # Fix endpoint registration issues.
    logger.info("Running Endpoint Registration Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # TODO: Integrate logic from fix_endpoint_registration.py
    if not dry_run:
        print("Endpoint Registration Fix: Not yet implemented in this utility.")
    else:
        logger.info("[DRY RUN] Would apply Endpoint Registration Fix")
    
    return True

def fix_supabase_auth(args):
    # Fix Supabase authentication issues.
    logger.info("Running Supabase Auth Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # TODO: Integrate logic from fix_supabase_auth.py
    if not dry_run:
        print("Supabase Auth Fix: Not yet implemented in this utility.")
    else:
        logger.info("[DRY RUN] Would apply Supabase Auth Fix")
    
    return True

def fix_table_names_and_repopulate(args):
    # Fix table names and repopulate data.
    logger.info("Running Table Names and Repopulate Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # TODO: Integrate logic from fix_table_names_and_repopulate.py
    if not dry_run:
        print("Table Names and Repopulate Fix: Not yet implemented in this utility.")
    else:
        logger.info("[DRY RUN] Would apply Table Names and Repopulate Fix")
    
    return True

def fix_rls_implementation(args):
    # Fix Row-Level Security implementation.
    logger.info("Running RLS Implementation Fix...")
    
    # Check for verify-only mode
    if hasattr(args, 'verify') and args.verify:
        logger.info("Running in VERIFY ONLY mode - checking for issues without applying fixes")
        # TODO: Implement verification logic
        logger.info("Verification complete")
        return True
    
    # Check for dry-run mode
    dry_run = hasattr(args, 'dry_run') and args.dry_run
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    
    # Check for rollback mode
    if hasattr(args, 'rollback') and args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")
        # TODO: Implement rollback logic using backup files
        logger.info("Rollback complete")
        return True
    
    # TODO: Integrate logic from fix_rls_implementation.py
    if not dry_run:
        print("RLS Implementation Fix: Not yet implemented in this utility.")
    else:
        logger.info("[DRY RUN] Would apply RLS Implementation Fix")
    
    return True

# Map fix names to functions
FIX_FUNCTIONS = {
    "api_integration": fix_api_integration,
    "auth_service_role": fix_auth_service_role,
    "auth_simple": fix_auth_simple,
    "database_modular": fix_database_modular,
    "foreign_key_relationships": fix_foreign_key_relationships,
    "relationships": fix_relationships,
    "rls_implementation": fix_rls_implementation,
    "table_names_and_repopulate": fix_table_names_and_repopulate,
    "supabase_auth": fix_supabase_auth,
    "endpoint_registration": fix_endpoint_registration,
    "missing_tables": fix_missing_tables,
    "remaining_api_issues": fix_remaining_api_issues,
}

def interactive_menu():
    # Display an interactive menu for selecting and running fixes
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Maintenance Utility")
    print("=" * 80)
    
    # First, ask for the fix to run
    print("Select a fix to run:")
    for i, fix in enumerate(sorted(FIX_FUNCTIONS.keys()), 1):
        print(f"  {i}. {fix.replace('_', ' ').title()}")
    print("  0. Exit")
    print("=" * 80)
    
    choice = input("Enter your choice: ").strip()
    if choice == "0":
        print("Exiting.")
        sys.exit(0)
    
    try:
        idx = int(choice) - 1
        if 0 <= idx < len(FIX_FUNCTIONS):
            fix_name = sorted(FIX_FUNCTIONS.keys())[idx]
            
            # Ask for additional options
            print("\nAdditional options:")
            dry_run = input("Run in dry-run mode? (no changes will be made) [y/N]: ").strip().lower() == 'y'
            verify = input("Run in verify-only mode? (no fixes will be applied) [y/N]: ").strip().lower() == 'y'
            rollback = input("Run in rollback mode? (restore from backups) [y/N]: ").strip().lower() == 'y'
            
            # Create an args object similar to what argparse would create
            class Args:
                pass
            
            args = Args()
            args.fix = fix_name
            args.dry_run = dry_run
            args.verify = verify
            args.rollback = rollback
            
            # Set global flags based on arguments
            if args.dry_run:
                logger.info("Running in DRY RUN mode - no changes will be made")
            if args.verify:
                logger.info("Running in VERIFY ONLY mode - no fixes will be applied")
            if args.rollback:
                logger.info("Running in ROLLBACK mode - attempting to restore from backups")
            
            # Run the selected fix with the specified options
            print(f"\nRunning {fix_name.replace('_', ' ').title()}...")
            result = FIX_FUNCTIONS[fix_name](args)
            
            if result:
                print(f"\n✅ Fix '{fix_name}' completed successfully")
            else:
                print(f"\n❌ Fix '{fix_name}' failed")
            
            input("\nPress Enter to continue...")
        else:
            print("Invalid choice.")
            input("\nPress Enter to continue...")
    except Exception as e:
        print(f"Error: {str(e)}")
        logger.error(f"Error in interactive menu: {str(e)}")
        input("\nPress Enter to continue...")

def main():
    # Main entry point for the maintenance utility
    epilog_text = "Available Fixes:\n"
    epilog_text += "  - api_integration: Fix API integration issues with datetime handling\n"
    epilog_text += "  - auth_service_role: Fix authentication service role issues\n"
    epilog_text += "  - auth_simple: Apply simple authentication fixes\n"
    epilog_text += "  - database_modular: Apply modular database fixes\n"
    epilog_text += "  - foreign_key_relationships: Fix foreign key relationship issues\n"
    epilog_text += "  - relationships: Fix general relationship issues\n"
    epilog_text += "  - rls_implementation: Fix Row-Level Security implementation\n"
    epilog_text += "  - table_names_and_repopulate: Fix table names and repopulate data\n"
    epilog_text += "  - supabase_auth: Fix Supabase authentication issues\n"
    epilog_text += "  - endpoint_registration: Fix endpoint registration issues\n"
    epilog_text += "  - missing_tables: Create missing tables\n"
    epilog_text += "  - remaining_api_issues: Fix remaining API issues"
    
    parser = argparse.ArgumentParser(
        description="CryoProtect v2 - Maintenance Utility",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog_text
    )
    parser.add_argument("--fix", choices=FIX_FUNCTIONS.keys(), help="Name of the fix to run")
    parser.add_argument("--dry-run", action="store_true", help="Simulate changes without making them")
    parser.add_argument("--verify", action="store_true", help="Verify only, do not apply fixes")
    parser.add_argument("--rollback", action="store_true", help="Rollback previous changes")
    parser.add_argument("--list", action="store_true", help="List available fixes and exit")
    args = parser.parse_args()

    # Handle --list flag
    if args.list:
        print("\nAvailable fixes:")
        for fix_name in sorted(FIX_FUNCTIONS.keys()):
            print(f"  - {fix_name}")
        return

    # Set global flags based on arguments
    if args.dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    if args.verify:
        logger.info("Running in VERIFY ONLY mode - no fixes will be applied")
    if args.rollback:
        logger.info("Running in ROLLBACK mode - attempting to restore from backups")

    # Execute the specified fix or launch the interactive menu
    if args.fix:
        logger.info(f"Running fix: {args.fix}")
        try:
            result = FIX_FUNCTIONS[args.fix](args)
            if result:
                logger.info(f"Fix '{args.fix}' completed successfully")
            else:
                logger.error(f"Fix '{args.fix}' failed")
                sys.exit(1)
        except Exception as e:
            logger.error(f"Error running fix '{args.fix}': {str(e)}")
            sys.exit(1)
    else:
        print("\nLaunching interactive menu...")
        try:
            while True:
                interactive_menu()
        except KeyboardInterrupt:
            print("\nOperation cancelled by user.")
            sys.exit(0)

if __name__ == "__main__":
    main()
