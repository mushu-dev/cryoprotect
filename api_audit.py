#!/usr/bin/env python3
"""
CryoProtect Analyzer API - Endpoint Audit Tool

This script audits all API endpoints for consistency and generates a report
of inconsistencies found. It checks for:
- Consistent response formats
- Proper HTTP status code usage
- Proper error handling
- Documentation completeness

Usage:
    python api_audit.py [--output-file FILENAME] [--verbose]

Options:
    --output-file FILENAME    Output file for the audit report (default: api_audit_report.md)
    --verbose                 Enable verbose output
"""

import os
import sys
import json
import inspect
import importlib
import argparse
import logging
from typing import Dict, List, Any, Tuple, Set
from datetime import datetime
from http import HTTPStatus

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Import Flask and related modules
try:
    from flask import Flask, Blueprint
    from flask_restful import Api, Resource
except ImportError:
    logger.error("Flask and Flask-RESTful are required. Install with: pip install flask flask-restful")
    sys.exit(1)

# Define audit categories
AUDIT_CATEGORIES = {
    'response_format': 'Response Format',
    'status_codes': 'HTTP Status Codes',
    'error_handling': 'Error Handling',
    'documentation': 'Documentation'
}

def find_api_modules() -> List[str]:
    """
    Find all API modules in the project.
    
    Returns:
        List of module names
    """
    api_modules = []
    api_dir = os.path.join(os.path.dirname(__file__), 'api')
    
    if not os.path.exists(api_dir):
        logger.warning(f"API directory not found: {api_dir}")
        return api_modules
    
    for filename in os.listdir(api_dir):
        if filename.endswith('.py') and not filename.startswith('__'):
            module_name = f"api.{filename[:-3]}"
            api_modules.append(module_name)
    
    return api_modules

def find_resource_classes(module_name: str) -> List[Tuple[str, Any]]:
    """
    Find all Resource classes in a module.
    
    Args:
        module_name: Name of the module to search
        
    Returns:
        List of (class_name, class_object) tuples
    """
    resource_classes = []
    
    try:
        module = importlib.import_module(module_name)
        
        for name, obj in inspect.getmembers(module):
            if inspect.isclass(obj) and issubclass(obj, Resource) and obj != Resource:
                resource_classes.append((name, obj))
    except ImportError as e:
        logger.warning(f"Failed to import module {module_name}: {e}")
    except Exception as e:
        logger.warning(f"Error processing module {module_name}: {e}")
    
    return resource_classes

def check_response_format(resource_class: Any) -> List[Dict[str, Any]]:
    """
    Check if a resource class uses consistent response formats.
    
    Args:
        resource_class: Resource class to check
        
    Returns:
        List of inconsistencies found
    """
    inconsistencies = []
    
    # Check HTTP methods
    http_methods = ['get', 'post', 'put', 'delete', 'patch']
    
    for method_name in http_methods:
        method = getattr(resource_class, method_name, None)
        
        if method and method.__qualname__.startswith(resource_class.__name__):
            # Check if method uses standardized response format
            source = inspect.getsource(method)
            
            # Check for direct return of data without standardization
            if 'return ' in source and not any(x in source for x in [
                'create_standard_response',
                'create_success_response',
                'create_error_response',
                'jsonify_standard_response',
                '@standardize_response'
            ]):
                inconsistencies.append({
                    'resource': resource_class.__name__,
                    'method': method_name.upper(),
                    'issue': 'Does not use standardized response format',
                    'location': f"{resource_class.__module__}.{resource_class.__name__}.{method_name}"
                })
    
    return inconsistencies

def check_status_codes(resource_class: Any) -> List[Dict[str, Any]]:
    """
    Check if a resource class uses proper HTTP status codes.
    
    Args:
        resource_class: Resource class to check
        
    Returns:
        List of inconsistencies found
    """
    inconsistencies = []
    
    # Check HTTP methods
    http_methods = ['get', 'post', 'put', 'delete', 'patch']
    
    for method_name in http_methods:
        method = getattr(resource_class, method_name, None)
        
        if method and method.__qualname__.startswith(resource_class.__name__):
            # Check if method uses proper HTTP status codes
            source = inspect.getsource(method)
            
            # Check for hardcoded status codes
            if 'return ' in source and any(x in source for x in [
                ', 200',
                ', 201',
                ', 204',
                ', 400',
                ', 401',
                ', 403',
                ', 404',
                ', 409',
                ', 429',
                ', 500'
            ]) and not any(x in source for x in [
                'HTTPStatus.',
                'HTTP_STATUS_CODES'
            ]):
                inconsistencies.append({
                    'resource': resource_class.__name__,
                    'method': method_name.upper(),
                    'issue': 'Uses hardcoded HTTP status codes',
                    'location': f"{resource_class.__module__}.{resource_class.__name__}.{method_name}"
                })
    
    return inconsistencies

def check_error_handling(resource_class: Any) -> List[Dict[str, Any]]:
    """
    Check if a resource class uses proper error handling.
    
    Args:
        resource_class: Resource class to check
        
    Returns:
        List of inconsistencies found
    """
    inconsistencies = []
    
    # Check HTTP methods
    http_methods = ['get', 'post', 'put', 'delete', 'patch']
    
    for method_name in http_methods:
        method = getattr(resource_class, method_name, None)
        
        if method and method.__qualname__.startswith(resource_class.__name__):
            # Check if method uses proper error handling
            source = inspect.getsource(method)
            
            # Check for try-except blocks
            if 'try:' not in source:
                inconsistencies.append({
                    'resource': resource_class.__name__,
                    'method': method_name.upper(),
                    'issue': 'Missing error handling (no try-except block)',
                    'location': f"{resource_class.__module__}.{resource_class.__name__}.{method_name}"
                })
            
            # Check for standardized error handling
            if 'except' in source and not any(x in source for x in [
                'create_error_response',
                'handle_error',
                'handle_supabase_error'
            ]):
                inconsistencies.append({
                    'resource': resource_class.__name__,
                    'method': method_name.upper(),
                    'issue': 'Does not use standardized error handling',
                    'location': f"{resource_class.__module__}.{resource_class.__name__}.{method_name}"
                })
    
    return inconsistencies

def check_documentation(resource_class: Any) -> List[Dict[str, Any]]:
    """
    Check if a resource class has proper documentation.
    
    Args:
        resource_class: Resource class to check
        
    Returns:
        List of inconsistencies found
    """
    inconsistencies = []
    
    # Check class docstring
    if not resource_class.__doc__:
        inconsistencies.append({
            'resource': resource_class.__name__,
            'method': 'CLASS',
            'issue': 'Missing class docstring',
            'location': f"{resource_class.__module__}.{resource_class.__name__}"
        })
    
    # Check HTTP methods
    http_methods = ['get', 'post', 'put', 'delete', 'patch']
    
    for method_name in http_methods:
        method = getattr(resource_class, method_name, None)
        
        if method and method.__qualname__.startswith(resource_class.__name__):
            # Check method docstring
            if not method.__doc__:
                inconsistencies.append({
                    'resource': resource_class.__name__,
                    'method': method_name.upper(),
                    'issue': 'Missing method docstring',
                    'location': f"{resource_class.__module__}.{resource_class.__name__}.{method_name}"
                })
            
            # Check for OpenAPI documentation
            source = inspect.getsource(method)
            if not any(x in source for x in [
                '@doc(',
                '@document_endpoint(',
                '@apispec_marshal_with(',
                '@use_kwargs('
            ]):
                inconsistencies.append({
                    'resource': resource_class.__name__,
                    'method': method_name.upper(),
                    'issue': 'Missing OpenAPI documentation',
                    'location': f"{resource_class.__module__}.{resource_class.__name__}.{method_name}"
                })
    
    return inconsistencies

def audit_api_endpoints() -> Dict[str, List[Dict[str, Any]]]:
    """
    Audit all API endpoints for consistency.
    
    Returns:
        Dictionary of inconsistencies by category
    """
    inconsistencies = {
        'response_format': [],
        'status_codes': [],
        'error_handling': [],
        'documentation': []
    }
    
    # Find all API modules
    api_modules = find_api_modules()
    logger.info(f"Found {len(api_modules)} API modules")
    
    # Find all resource classes
    resource_classes = []
    for module_name in api_modules:
        module_resources = find_resource_classes(module_name)
        resource_classes.extend(module_resources)
    
    logger.info(f"Found {len(resource_classes)} resource classes")
    
    # Check each resource class
    for class_name, resource_class in resource_classes:
        logger.debug(f"Checking resource class: {class_name}")
        
        # Check response format
        response_format_issues = check_response_format(resource_class)
        inconsistencies['response_format'].extend(response_format_issues)
        
        # Check status codes
        status_code_issues = check_status_codes(resource_class)
        inconsistencies['status_codes'].extend(status_code_issues)
        
        # Check error handling
        error_handling_issues = check_error_handling(resource_class)
        inconsistencies['error_handling'].extend(error_handling_issues)
        
        # Check documentation
        documentation_issues = check_documentation(resource_class)
        inconsistencies['documentation'].extend(documentation_issues)
    
    return inconsistencies

def generate_report(inconsistencies: Dict[str, List[Dict[str, Any]]], output_file: str) -> None:
    """
    Generate a report of inconsistencies found.
    
    Args:
        inconsistencies: Dictionary of inconsistencies by category
        output_file: Output file for the report
    """
    total_issues = sum(len(issues) for issues in inconsistencies.values())
    
    with open(output_file, 'w') as f:
        f.write("# CryoProtect API Endpoint Audit Report\n\n")
        f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Summary\n\n")
        f.write(f"Total issues found: {total_issues}\n\n")
        
        for category, issues in inconsistencies.items():
            f.write(f"- {AUDIT_CATEGORIES[category]}: {len(issues)} issues\n")
        
        f.write("\n## Issues by Category\n\n")
        
        for category, issues in inconsistencies.items():
            if not issues:
                continue
                
            f.write(f"### {AUDIT_CATEGORIES[category]}\n\n")
            
            f.write("| Resource | Method | Issue | Location |\n")
            f.write("|----------|--------|-------|----------|\n")
            
            for issue in issues:
                f.write(f"| {issue['resource']} | {issue['method']} | {issue['issue']} | {issue['location']} |\n")
            
            f.write("\n")
        
        f.write("## Recommendations\n\n")
        
        f.write("1. **Standardize Response Formats**\n")
        f.write("   - Use the `create_standard_response` function for all responses\n")
        f.write("   - Apply the `@standardize_response` decorator to all endpoint methods\n\n")
        
        f.write("2. **Use Proper HTTP Status Codes**\n")
        f.write("   - Use the `HTTPStatus` enum for all status codes\n")
        f.write("   - Follow REST conventions for status code usage\n\n")
        
        f.write("3. **Implement Consistent Error Handling**\n")
        f.write("   - Use try-except blocks in all endpoint methods\n")
        f.write("   - Use the `create_error_response` function for all error responses\n\n")
        
        f.write("4. **Complete API Documentation**\n")
        f.write("   - Add docstrings to all classes and methods\n")
        f.write("   - Use the `@document_endpoint` decorator for OpenAPI documentation\n")
    
    logger.info(f"Report generated: {output_file}")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Audit API endpoints for consistency')
    parser.add_argument('--output-file', default='api_audit_report.md', help='Output file for the audit report')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info("Starting API endpoint audit")
    
    # Audit API endpoints
    inconsistencies = audit_api_endpoints()
    
    # Generate report
    generate_report(inconsistencies, args.output_file)
    
    logger.info("API endpoint audit completed")

if __name__ == '__main__':
    main()