#!/usr/bin/env python3
"""
CryoProtect Analyzer API - Endpoint Standardization Tool

This script applies standardization to all API endpoints based on the audit results.
It updates the code to use standardized response formats, proper HTTP status codes,
consistent error handling, and complete documentation.

Usage:
    python apply_api_standardization.py [--audit-file FILENAME] [--dry-run] [--verbose]

Options:
    --audit-file FILENAME    Audit report file (default: api_audit_report.md)
    --dry-run                Show changes without applying them
    --verbose                Enable verbose output
"""

import os
import sys
import re
import json
import inspect
import importlib
import argparse
import logging
import ast
from typing import Dict, List, Any, Tuple, Set
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_audit_report(audit_file: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Parse the audit report to extract inconsistencies.
    
    Args:
        audit_file: Path to the audit report file
        
    Returns:
        Dictionary of inconsistencies by category
    """
    inconsistencies = {
        'response_format': [],
        'status_codes': [],
        'error_handling': [],
        'documentation': []
    }
    
    if not os.path.exists(audit_file):
        logger.warning(f"Audit file not found: {audit_file}")
        return inconsistencies
    
    try:
        with open(audit_file, 'r') as f:
            content = f.read()
        
        # Extract issues by category
        for category in inconsistencies.keys():
            pattern = rf"### .*?{category.replace('_', ' ').title()}.*?\n\n\| Resource \| Method \| Issue \| Location \|\n\|.*?\|\n(.*?)(?:\n\n|\Z)"
            match = re.search(pattern, content, re.DOTALL | re.IGNORECASE)
            
            if match:
                issues_table = match.group(1)
                
                # Parse table rows
                for row in issues_table.strip().split('\n'):
                    if not row.startswith('|'):
                        continue
                        
                    parts = [part.strip() for part in row.split('|')]
                    if len(parts) < 6:
                        continue
                        
                    resource = parts[1]
                    method = parts[2]
                    issue = parts[3]
                    location = parts[4]
                    
                    inconsistencies[category].append({
                        'resource': resource,
                        'method': method,
                        'issue': issue,
                        'location': location
                    })
    except Exception as e:
        logger.error(f"Error parsing audit report: {e}")
    
    return inconsistencies

def find_module_path(module_name: str) -> str:
    """
    Find the file path for a module.
    
    Args:
        module_name: Name of the module
        
    Returns:
        File path for the module
    """
    try:
        module = importlib.import_module(module_name)
        file_path = inspect.getfile(module)
        return file_path
    except (ImportError, TypeError) as e:
        logger.warning(f"Failed to find module path for {module_name}: {e}")
        
        # Try to construct the path manually
        parts = module_name.split('.')
        base_path = os.path.dirname(__file__)
        
        for i in range(len(parts)):
            path = os.path.join(base_path, *parts[i:]) + '.py'
            if os.path.exists(path):
                return path
                
        return None

def extract_method_source(file_path: str, class_name: str, method_name: str) -> Tuple[str, int, int]:
    """
    Extract the source code for a method from a file.
    
    Args:
        file_path: Path to the file
        class_name: Name of the class
        method_name: Name of the method
        
    Returns:
        Tuple of (source_code, start_line, end_line)
    """
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Parse the file
    tree = ast.parse(content)
    
    # Find the class
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            # Find the method
            for child in node.body:
                if isinstance(child, ast.FunctionDef) and child.name == method_name:
                    start_line = child.lineno
                    end_line = child.end_lineno if hasattr(child, 'end_lineno') else 0
                    
                    # If end_line is not available, estimate it
                    if end_line == 0:
                        lines = content.split('\n')
                        indent = 0
                        
                        # Find the indentation of the method
                        for i in range(start_line - 1, len(lines)):
                            line = lines[i]
                            if i == start_line - 1:
                                # First line of the method
                                indent = len(line) - len(line.lstrip())
                            elif line.strip() and len(line) - len(line.lstrip()) <= indent:
                                # Found a line with less indentation
                                end_line = i
                                break
                        
                        if end_line == 0:
                            # If we couldn't find the end, use the end of the file
                            end_line = len(lines)
                    
                    # Extract the source code
                    lines = content.split('\n')
                    source_code = '\n'.join(lines[start_line - 1:end_line])
                    
                    return source_code, start_line, end_line
    
    return None, 0, 0

def standardize_response_format(source_code: str) -> str:
    """
    Standardize the response format in a method.
    
    Args:
        source_code: Source code of the method
        
    Returns:
        Updated source code
    """
    # Add imports if needed
    imports_added = False
    if 'from api.api_standards import' not in source_code:
        # Add imports at the beginning of the method
        lines = source_code.split('\n')
        indent = len(lines[0]) - len(lines[0].lstrip())
        indent_str = ' ' * indent
        
        # Find the first non-docstring line
        start_index = 0
        in_docstring = False
        for i, line in enumerate(lines):
            if i == 0:
                continue
                
            if '"""' in line or "'''" in line:
                in_docstring = not in_docstring
                continue
                
            if not in_docstring and line.strip():
                start_index = i
                break
        
        # Add imports
        import_line = f"{indent_str}from api.api_standards import create_success_response, jsonify_standard_response"
        lines.insert(start_index, import_line)
        source_code = '\n'.join(lines)
        imports_added = True
    
    # Replace direct returns with standardized responses
    pattern = r'return\s+([^,]+),\s*(\d+)'
    replacement = r'return jsonify_standard_response(*create_success_response(data=\1, status_code=\2))'
    source_code = re.sub(pattern, replacement, source_code)
    
    # Add decorator if needed
    if '@standardize_response' not in source_code:
        lines = source_code.split('\n')
        indent = len(lines[0]) - len(lines[0].lstrip())
        indent_str = ' ' * indent
        
        # Add decorator before the method definition
        decorator_line = f"{indent_str}@standardize_response"
        lines.insert(0, decorator_line)
        
        # Add import if needed
        if not imports_added:
            import_line = f"{indent_str}from api.api_decorators import standardize_response"
            lines.insert(0, import_line)
            
        source_code = '\n'.join(lines)
    
    return source_code

def standardize_status_codes(source_code: str) -> str:
    """
    Standardize HTTP status codes in a method.
    
    Args:
        source_code: Source code of the method
        
    Returns:
        Updated source code
    """
    # Add imports if needed
    if 'from http import HTTPStatus' not in source_code:
        # Add imports at the beginning of the method
        lines = source_code.split('\n')
        indent = len(lines[0]) - len(lines[0].lstrip())
        indent_str = ' ' * indent
        
        # Find the first non-docstring line
        start_index = 0
        in_docstring = False
        for i, line in enumerate(lines):
            if i == 0:
                continue
                
            if '"""' in line or "'''" in line:
                in_docstring = not in_docstring
                continue
                
            if not in_docstring and line.strip():
                start_index = i
                break
        
        # Add imports
        import_line = f"{indent_str}from http import HTTPStatus"
        lines.insert(start_index, import_line)
        source_code = '\n'.join(lines)
    
    # Replace hardcoded status codes with HTTPStatus enum
    status_code_map = {
        '200': 'HTTPStatus.OK',
        '201': 'HTTPStatus.CREATED',
        '204': 'HTTPStatus.NO_CONTENT',
        '400': 'HTTPStatus.BAD_REQUEST',
        '401': 'HTTPStatus.UNAUTHORIZED',
        '403': 'HTTPStatus.FORBIDDEN',
        '404': 'HTTPStatus.NOT_FOUND',
        '409': 'HTTPStatus.CONFLICT',
        '429': 'HTTPStatus.TOO_MANY_REQUESTS',
        '500': 'HTTPStatus.INTERNAL_SERVER_ERROR',
        '503': 'HTTPStatus.SERVICE_UNAVAILABLE'
    }
    
    for code, enum in status_code_map.items():
        pattern = r'(\W)' + code + r'(\W)'
        replacement = r'\1' + enum + r'\2'
        source_code = re.sub(pattern, replacement, source_code)
    
    return source_code

def standardize_error_handling(source_code: str) -> str:
    """
    Standardize error handling in a method.
    
    Args:
        source_code: Source code of the method
        
    Returns:
        Updated source code
    """
    # Add imports if needed
    imports_added = False
    if 'from api.api_standards import' not in source_code:
        # Add imports at the beginning of the method
        lines = source_code.split('\n')
        indent = len(lines[0]) - len(lines[0].lstrip())
        indent_str = ' ' * indent
        
        # Find the first non-docstring line
        start_index = 0
        in_docstring = False
        for i, line in enumerate(lines):
            if i == 0:
                continue
                
            if '"""' in line or "'''" in line:
                in_docstring = not in_docstring
                continue
                
            if not in_docstring and line.strip():
                start_index = i
                break
        
        # Add imports
        import_line = f"{indent_str}from api.api_standards import create_error_response, jsonify_standard_response"
        lines.insert(start_index, import_line)
        source_code = '\n'.join(lines)
        imports_added = True
    
    # Add try-except block if needed
    if 'try:' not in source_code:
        lines = source_code.split('\n')
        indent = len(lines[0]) - len(lines[0].lstrip())
        indent_str = ' ' * indent
        
        # Find the first non-docstring line
        start_index = 0
        end_index = len(lines) - 1
        in_docstring = False
        for i, line in enumerate(lines):
            if '"""' in line or "'''" in line:
                in_docstring = not in_docstring
                continue
                
            if not in_docstring and line.strip() and 'def ' not in line:
                start_index = i
                break
        
        # Add try-except block
        try_line = f"{indent_str}try:"
        except_line = f"{indent_str}except Exception as e:"
        error_line = f"{indent_str}    return jsonify_standard_response(*create_error_response(error=e, context=\"Error in {lines[0].strip()}\"))"
        
        # Indent the method body
        for i in range(start_index, end_index + 1):
            lines[i] = indent_str + "    " + lines[i][len(indent_str):]
        
        # Add try-except block
        lines.insert(start_index, try_line)
        lines.append(except_line)
        lines.append(error_line)
        
        source_code = '\n'.join(lines)
    
    # Replace error handling with standardized error handling
    pattern = r'except\s+([^\:]+):\s*\n\s+return\s+([^,]+),\s*(\d+)'
    replacement = r'except \1:\n        return jsonify_standard_response(*create_error_response(error=e, status_code=\3, context="Error in method"))'
    source_code = re.sub(pattern, replacement, source_code)
    
    return source_code

def standardize_documentation(source_code: str, class_name: str, method_name: str) -> str:
    """
    Standardize documentation in a method.
    
    Args:
        source_code: Source code of the method
        class_name: Name of the class
        method_name: Name of the method
        
    Returns:
        Updated source code
    """
    # Add docstring if needed
    if '"""' not in source_code and "'''" not in source_code:
        lines = source_code.split('\n')
        indent = len(lines[0]) - len(lines[0].lstrip())
        indent_str = ' ' * indent
        
        # Add docstring after the method definition
        docstring_start = f'{indent_str}"""'
        docstring_line1 = f'{indent_str}Handle {method_name.upper()} request for {class_name}.'
        docstring_line2 = f'{indent_str}'
        docstring_line3 = f'{indent_str}Returns:'
        docstring_line4 = f'{indent_str}    Standardized API response'
        docstring_end = f'{indent_str}"""'
        
        # Find the method definition line
        for i, line in enumerate(lines):
            if 'def ' in line:
                # Add docstring after the method definition
                lines.insert(i + 1, docstring_start)
                lines.insert(i + 2, docstring_line1)
                lines.insert(i + 3, docstring_line2)
                lines.insert(i + 4, docstring_line3)
                lines.insert(i + 5, docstring_line4)
                lines.insert(i + 6, docstring_end)
                break
        
        source_code = '\n'.join(lines)
    
    # Add OpenAPI documentation if needed
    if not any(x in source_code for x in ['@doc(', '@document_endpoint(']):
        lines = source_code.split('\n')
        indent = len(lines[0]) - len(lines[0].lstrip())
        indent_str = ' ' * indent
        
        # Add OpenAPI documentation before the method definition
        doc_line1 = f'{indent_str}@document_endpoint('
        doc_line2 = f'{indent_str}    summary="{method_name.capitalize()} {class_name.replace("Resource", "").lower()}",'
        doc_line3 = f'{indent_str}    description="Handle {method_name.upper()} request for {class_name}.",'
        doc_line4 = f'{indent_str}    tags=["{class_name.replace("Resource", "")}"]'
        doc_line5 = f'{indent_str})'
        
        # Find the method definition line
        for i, line in enumerate(lines):
            if 'def ' in line:
                # Add OpenAPI documentation before the method definition
                lines.insert(i, doc_line1)
                lines.insert(i + 1, doc_line2)
                lines.insert(i + 2, doc_line3)
                lines.insert(i + 3, doc_line4)
                lines.insert(i + 4, doc_line5)
                break
        
        # Add import if needed
        if 'from api.api_docs import' not in source_code:
            # Find the first non-docstring line
            start_index = 0
            in_docstring = False
            for i, line in enumerate(lines):
                if i == 0:
                    continue
                    
                if '"""' in line or "'''" in line:
                    in_docstring = not in_docstring
                    continue
                    
                if not in_docstring and line.strip():
                    start_index = i
                    break
            
            # Add imports
            import_line = f"{indent_str}from api.api_docs import document_endpoint"
            lines.insert(start_index, import_line)
        
        source_code = '\n'.join(lines)
    
    return source_code

def apply_standardization(inconsistencies: Dict[str, List[Dict[str, Any]]], dry_run: bool = False) -> Dict[str, int]:
    """
    Apply standardization to API endpoints based on inconsistencies.
    
    Args:
        inconsistencies: Dictionary of inconsistencies by category
        dry_run: Whether to show changes without applying them
        
    Returns:
        Dictionary of counts of changes applied by category
    """
    changes_applied = {
        'response_format': 0,
        'status_codes': 0,
        'error_handling': 0,
        'documentation': 0
    }
    
    # Process each category
    for category, issues in inconsistencies.items():
        for issue in issues:
            try:
                # Extract location information
                location = issue['location']
                parts = location.split('.')
                
                if len(parts) < 3:
                    logger.warning(f"Invalid location format: {location}")
                    continue
                
                module_name = '.'.join(parts[:-2])
                class_name = parts[-2]
                method_name = parts[-1]
                
                # Find the file path
                file_path = find_module_path(module_name)
                if not file_path:
                    logger.warning(f"Could not find file path for module: {module_name}")
                    continue
                
                # Extract the method source
                source_code, start_line, end_line = extract_method_source(file_path, class_name, method_name)
                if not source_code:
                    logger.warning(f"Could not extract source code for method: {location}")
                    continue
                
                # Apply standardization based on category
                updated_source = source_code
                if category == 'response_format':
                    updated_source = standardize_response_format(source_code)
                elif category == 'status_codes':
                    updated_source = standardize_status_codes(source_code)
                elif category == 'error_handling':
                    updated_source = standardize_error_handling(source_code)
                elif category == 'documentation':
                    updated_source = standardize_documentation(source_code, class_name, method_name)
                
                # Check if changes were made
                if updated_source != source_code:
                    changes_applied[category] += 1
                    
                    if dry_run:
                        logger.info(f"Would update {location} ({category})")
                        logger.debug(f"Original:\n{source_code}\n\nUpdated:\n{updated_source}")
                    else:
                        # Update the file
                        with open(file_path, 'r') as f:
                            content = f.read()
                        
                        # Replace the method source
                        lines = content.split('\n')
                        updated_lines = updated_source.split('\n')
                        
                        # Replace the lines
                        lines[start_line - 1:end_line] = updated_lines
                        
                        # Write the updated content
                        with open(file_path, 'w') as f:
                            f.write('\n'.join(lines))
                        
                        logger.info(f"Updated {location} ({category})")
            except Exception as e:
                logger.error(f"Error applying standardization to {issue['location']}: {e}")
    
    return changes_applied

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Apply API endpoint standardization')
    parser.add_argument('--audit-file', default='api_audit_report.md', help='Audit report file')
    parser.add_argument('--dry-run', action='store_true', help='Show changes without applying them')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info("Starting API endpoint standardization")
    
    # Parse audit report
    inconsistencies = parse_audit_report(args.audit_file)
    
    # Apply standardization
    changes_applied = apply_standardization(inconsistencies, args.dry_run)
    
    # Print summary
    total_changes = sum(changes_applied.values())
    logger.info(f"Standardization completed: {total_changes} changes applied")
    
    for category, count in changes_applied.items():
        logger.info(f"- {category.replace('_', ' ').title()}: {count} changes")
    
    if args.dry_run:
        logger.info("Dry run completed. No changes were applied.")

if __name__ == '__main__':
    main()