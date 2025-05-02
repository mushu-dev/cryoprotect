#!/usr/bin/env python3
"""
Environment Variable Validation Script for CryoProtect v2

This script validates that all required environment variables are set by:
1. Parsing the .env.template file to extract variable names
2. Checking if each variable is set in the current environment
3. Reporting any missing variables
4. Exiting with appropriate status code

Usage:
  python validate_env.py [options]

Options:
  --template PATH     Path to the .env.template file (default: .env.template)
  --section SECTION   Only check variables in the specified section
  --all               Check all variables, not just required ones
  --verbose           Show more detailed output
  --help              Show this help message and exit
"""

import os
import sys
import re
import argparse
from collections import defaultdict

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate environment variables for CryoProtect v2",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--template",
        default=".env.template",
        help="Path to the .env.template file (default: .env.template)"
    )
    parser.add_argument(
        "--section",
        help="Only check variables in the specified section"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Check all variables, not just required ones"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show more detailed output"
    )
    
    return parser.parse_args()

def parse_env_template(template_path='.env.template'):
    """
    Parse the .env.template file to extract environment variable names.
    
    Args:
        template_path (str): Path to the .env.template file
        
    Returns:
        dict: Dictionary of variable names with their required status
    """
    if not os.path.exists(template_path):
        print(f"Error: Template file '{template_path}' not found.")
        sys.exit(1)
    
    variables = {}
    current_section = "GENERAL"
    
    # Read the entire file content
    with open(template_path, 'r') as file:
        lines = file.readlines()
    
    # Process lines with context awareness
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Skip empty lines
        if not line:
            i += 1
            continue
            
        # Check if this is a section header
        if line.startswith('# ====='):
            current_section = line.strip('# =').strip()
            i += 1
            continue
            
        # Check if this is a variable definition
        if '=' in line and not line.startswith('#'):
            var_name = line.split('=', 1)[0].strip()
            
            # Look at preceding comments to determine if required
            is_required = False
            description = ""
            
            # Look back up to 5 lines for comments about this variable
            for j in range(max(0, i-5), i):
                prev_line = lines[j].strip()
                if prev_line.startswith('#'):
                    comment_text = prev_line[1:].strip()
                    description += comment_text + " "
                    if 'required' in comment_text.lower() and not 'optional' in comment_text.lower():
                        is_required = True
            
            # If no explicit requirement found, use heuristics
            if not is_required and not 'optional' in description.lower():
                # Consider variables in certain sections as required by default
                if current_section in ["DATABASE CONFIGURATION"]:
                    # In database section, consider URL and key as required
                    if var_name in ["SUPABASE_URL", "SUPABASE_KEY"]:
                        is_required = True
                
                # Consider variables with "password" or "key" in critical sections as required
                if current_section in ["DATABASE CONFIGURATION", "APPLICATION CONFIGURATION"]:
                    if "key" in var_name.lower() and "secret" in description.lower():
                        is_required = True
            
            variables[var_name] = {
                'required': is_required,
                'section': current_section,
                'description': description.strip()
            }
        
        i += 1
    
    return variables

def validate_environment_variables(variables, check_all=False, section_filter=None):
    """
    Check if required environment variables are set.
    
    Args:
        variables (dict): Dictionary of variable names with their required status
        check_all (bool): Whether to check all variables or just required ones
        section_filter (str): Only check variables in this section if provided
        
    Returns:
        tuple: Lists of missing required variables and missing optional variables
    """
    missing_required = []
    missing_optional = []
    
    for var_name, info in variables.items():
        # Skip if we're filtering by section and this isn't in the target section
        if section_filter and info['section'] != section_filter:
            continue
            
        # Check if the variable is set in the environment
        if not os.environ.get(var_name):
            if info['required']:
                missing_required.append((var_name, info['section'], info.get('description', '')))
            elif check_all:
                missing_optional.append((var_name, info['section'], info.get('description', '')))
    
    return missing_required, missing_optional

def main():
    """Main function to validate environment variables."""
    args = parse_arguments()
    
    print(f"Validating environment variables for CryoProtect v2 using template: {args.template}")
    if args.section:
        print(f"Checking only variables in section: {args.section}")
    if args.all:
        print("Checking all variables (required and optional)")
    
    # Parse .env.template file
    variables = parse_env_template(args.template)
    
    if args.verbose:
        print(f"\nFound {len(variables)} variables in template file")
        required_count = sum(1 for info in variables.values() if info['required'])
        print(f"  - {required_count} required variables")
        print(f"  - {len(variables) - required_count} optional variables")
    
    # Validate environment variables
    missing_required, missing_optional = validate_environment_variables(
        variables,
        check_all=args.all,
        section_filter=args.section
    )
    
    # Report results
    exit_code = 0
    
    if missing_required:
        print("\nERROR: The following required environment variables are missing:")
        exit_code = 1
        
        # Group by section
        by_section = defaultdict(list)
        for var_name, section, description in missing_required:
            by_section[section].append((var_name, description))
        
        # Print missing variables by section
        for section, vars in by_section.items():
            print(f"\n{section}:")
            for var, desc in vars:
                if desc and args.verbose:
                    print(f"  - {var} ({desc})")
                else:
                    print(f"  - {var}")
    
    if missing_optional and args.all:
        print("\nWARNING: The following optional environment variables are not set:")
        
        # Group by section
        by_section = defaultdict(list)
        for var_name, section, description in missing_optional:
            by_section[section].append((var_name, description))
        
        # Print missing variables by section
        for section, vars in by_section.items():
            print(f"\n{section}:")
            for var, desc in vars:
                if desc and args.verbose:
                    print(f"  - {var} ({desc})")
                else:
                    print(f"  - {var}")
    
    if exit_code == 1:
        print("\nPlease set these required variables in your environment or .env file before running the application.")
    else:
        print("All required environment variables are set.")
    
    sys.exit(exit_code)

if __name__ == "__main__":
    main()