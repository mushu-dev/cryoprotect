#!/usr/bin/env python
"""
Backend Configuration Consolidation Script

This script merges multiple backend configuration files into a single unified config.py file.
It analyzes the configuration values across different files to ensure consistency.
"""

import os
import sys
import re
import ast
import argparse
from pathlib import Path


def extract_config_values(file_path):
    """Extract configuration values and their types from a Python file."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Parse the file
    tree = ast.parse(content)
    
    # Extract different types of variables
    assigned_vars = {}
    env_vars = {}
    class_vars = {}
    
    # Find top-level assignments (e.g., DB_HOST = os.getenv(...))
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    var_name = target.id
                    
                    # Check if it's an os.getenv call
                    if isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Attribute):
                        func_name = f"{node.value.func.value.id}.{node.value.func.attr}" if hasattr(node.value.func.value, 'id') else ''
                        if func_name == 'os.getenv':
                            env_vars[var_name] = {
                                'source': 'env',
                                'env_var': ast.literal_eval(node.value.args[0]) if node.value.args else var_name,
                                'default': ast.literal_eval(node.value.args[1]) if len(node.value.args) > 1 else None,
                                'line': node.lineno
                            }
                        else:
                            try:
                                assigned_vars[var_name] = {
                                    'source': 'assignment',
                                    'value': ast.literal_eval(node.value) if isinstance(node.value, (ast.Str, ast.Num, ast.NameConstant, ast.List, ast.Dict)) else str(node.value),
                                    'line': node.lineno
                                }
                            except (ValueError, TypeError, SyntaxError):
                                assigned_vars[var_name] = {
                                    'source': 'assignment',
                                    'value': 'complex_value',
                                    'line': node.lineno
                                }
    
    # Find class variables in config classes
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and (node.name.endswith('Config') or 'Config' in node.name):
            class_name = node.name
            class_vars[class_name] = {}
            
            for item in node.body:
                if isinstance(item, ast.AnnAssign):  # Type-annotated assignment (e.g., DEBUG: bool = True)
                    var_name = item.target.id if isinstance(item.target, ast.Name) else str(item.target)
                    type_hint = ast.unparse(item.annotation) if hasattr(ast, 'unparse') else str(item.annotation)
                    
                    try:
                        if item.value:
                            value = ast.literal_eval(item.value) if isinstance(item.value, (ast.Str, ast.Num, ast.NameConstant, ast.List, ast.Dict)) else str(item.value)
                        else:
                            value = None
                    except (ValueError, TypeError, SyntaxError):
                        value = 'complex_value'
                    
                    class_vars[class_name][var_name] = {
                        'type': type_hint,
                        'value': value,
                        'line': item.lineno
                    }
                elif isinstance(item, ast.Assign):  # Regular assignment (e.g., DEBUG = True)
                    for target in item.targets:
                        if isinstance(target, ast.Name):
                            var_name = target.id
                            
                            try:
                                value = ast.literal_eval(item.value) if isinstance(item.value, (ast.Str, ast.Num, ast.NameConstant, ast.List, ast.Dict)) else str(item.value)
                            except (ValueError, TypeError, SyntaxError):
                                value = 'complex_value'
                            
                            class_vars[class_name][var_name] = {
                                'type': None,
                                'value': value,
                                'line': item.lineno
                            }
    
    return {
        'assigned_vars': assigned_vars,
        'env_vars': env_vars,
        'class_vars': class_vars
    }


def consolidate_configs(main_config_path, secondary_config_paths):
    """Consolidate multiple configuration files into the main config."""
    # Extract values from the main config
    main_config = extract_config_values(main_config_path)
    
    # Extract values from secondary configs
    secondary_configs = []
    for path in secondary_config_paths:
        secondary_configs.append({
            'path': path,
            'values': extract_config_values(path)
        })
    
    # Print a summary of what was found
    print(f"\nConfiguration Analysis Summary")
    print(f"==============================")
    print(f"Main Config: {main_config_path}")
    print(f"  - {len(main_config['env_vars'])} environment variables")
    print(f"  - {len(main_config['assigned_vars'])} assigned variables")
    print(f"  - {len(main_config['class_vars'])} configuration classes")
    
    for secondary in secondary_configs:
        print(f"\nSecondary Config: {secondary['path']}")
        print(f"  - {len(secondary['values']['env_vars'])} environment variables")
        print(f"  - {len(secondary['values']['assigned_vars'])} assigned variables")
        print(f"  - {len(secondary['values']['class_vars'])} configuration classes")
    
    # Identify overlaps and differences
    all_env_vars = set(main_config['env_vars'].keys())
    all_assigned_vars = set(main_config['assigned_vars'].keys())
    all_class_vars = set(main_config['class_vars'].keys())
    
    for secondary in secondary_configs:
        # Env vars
        sec_env_vars = set(secondary['values']['env_vars'].keys())
        overlap_env = all_env_vars.intersection(sec_env_vars)
        unique_sec_env = sec_env_vars - all_env_vars
        
        # Assigned vars
        sec_assigned_vars = set(secondary['values']['assigned_vars'].keys())
        overlap_assigned = all_assigned_vars.intersection(sec_assigned_vars)
        unique_sec_assigned = sec_assigned_vars - all_assigned_vars
        
        # Class vars
        sec_class_vars = set(secondary['values']['class_vars'].keys())
        overlap_class = all_class_vars.intersection(sec_class_vars)
        unique_sec_class = sec_class_vars - all_class_vars
        
        print(f"\nOverlap Analysis for {secondary['path']}")
        print(f"  - Environment Variables: {len(overlap_env)} overlapping, {len(unique_sec_env)} unique")
        print(f"  - Assigned Variables: {len(overlap_assigned)} overlapping, {len(unique_sec_assigned)} unique")
        print(f"  - Config Classes: {len(overlap_class)} overlapping, {len(unique_sec_class)} unique")
        
        if overlap_env:
            print(f"    Overlapping env vars: {', '.join(sorted(overlap_env))}")
        if unique_sec_env:
            print(f"    Unique env vars: {', '.join(sorted(unique_sec_env))}")
        
        # Update the master sets
        all_env_vars.update(sec_env_vars)
        all_assigned_vars.update(sec_assigned_vars)
        all_class_vars.update(sec_class_vars)
    
    # Generate a consolidated configuration file
    output_path = Path(main_config_path).with_suffix('.consolidated.py')
    
    with open(main_config_path, 'r') as f:
        main_content = f.read()
    
    # Find import section end
    import_end = 0
    in_imports = False
    import_section = []
    
    for i, line in enumerate(main_content.splitlines()):
        if line.startswith('import ') or line.startswith('from '):
            in_imports = True
            import_section.append(line)
            import_end = i + 1
        elif in_imports and line.strip() and not line.startswith('#') and not line.startswith('import ') and not line.startswith('from '):
            # End of import section
            in_imports = False
    
    # Collect unique imports from secondary files
    secondary_imports = set()
    for secondary in secondary_configs:
        with open(secondary['path'], 'r') as f:
            secondary_content = f.readlines()
            
        for line in secondary_content:
            line = line.strip()
            if line.startswith('import ') or line.startswith('from '):
                secondary_imports.add(line)
    
    # Add unique imports to the import section
    for imp in sorted(list(secondary_imports)):
        if imp not in import_section:
            import_section.append(imp)
    
    # Insert a divider for added sections
    import_section.append("\n# Additional imports from consolidated files")
    
    # Construct the new file content
    lines = main_content.splitlines()
    
    # Replace import section
    new_content = import_section + lines[import_end:]
    
    # Write the consolidated file
    with open(output_path, 'w') as f:
        f.write('\n'.join(new_content))
    
    print(f"\nConsolidated configuration written to: {output_path}")
    print("Please review and edit as needed.")
    return output_path


def main():
    parser = argparse.ArgumentParser(description="Consolidate backend configuration files")
    parser.add_argument("--main", default="../config.py", help="Path to main config.py")
    parser.add_argument("--secondary", nargs="+", default=["../frontend/db_config.py"], help="Paths to secondary config files")
    
    args = parser.parse_args()
    
    # Resolve paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    main_config_path = os.path.abspath(os.path.join(script_dir, args.main))
    secondary_config_paths = [os.path.abspath(os.path.join(script_dir, path)) for path in args.secondary]
    
    # Ensure main config file exists
    if not os.path.exists(main_config_path):
        print(f"Error: Main config file not found at {main_config_path}")
        return 1
    
    # Filter out non-existent secondary configs
    valid_secondary_paths = []
    for path in secondary_config_paths:
        if os.path.exists(path):
            valid_secondary_paths.append(path)
        else:
            print(f"Warning: Secondary config file not found at {path}")
    
    if not valid_secondary_paths:
        print("No valid secondary config files found")
        return 1
    
    # Consolidate configs
    output_path = consolidate_configs(main_config_path, valid_secondary_paths)
    
    print(f"\nNext Steps:")
    print(f"1. Review {output_path} for any inconsistencies or errors")
    print(f"2. Update references to secondary config files to use the consolidated version")
    print(f"3. Test the application with the new configuration")
    print(f"4. Once satisfied, rename {output_path} to {main_config_path}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())