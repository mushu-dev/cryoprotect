#!/usr/bin/env python3
"""
Update existing code to use the rdkit_wrapper module instead of direct RDKit imports.
This script scans Python files for direct RDKit imports and suggests changes.
It can either print suggestions or apply them directly to the files.
"""

import os
import re
import sys
import argparse
from typing import List, Dict, Set, Tuple, Optional

# Patterns to search for
DIRECT_IMPORT_PATTERNS = [
    r'from\s+rdkit\s+import\s+(\w+)',
    r'from\s+rdkit\.Chem\s+import\s+(\w+)',
    r'from\s+rdkit\.Chem\s+import\s+([^#\n]+)',
    r'import\s+rdkit',
    r'import\s+rdkit\.Chem'
]

RDKIT_USAGE_PATTERNS = [
    r'Chem\.MolFromSmiles',
    r'Chem\.MolToSmiles',
    r'Descriptors\.MolWt',
    r'Descriptors\.MolLogP',
    r'MolSurf\.TPSA',
    r'Lipinski\.NumHDonors',
    r'Lipinski\.NumHAcceptors',
    r'Descriptors\.NumRotatableBonds',
    r'Chem\.rdMolDescriptors\.CalcNumRings',
    r'MolDraw2DSVG',
    r'AllChem\.GetMorganFingerprintAsBitVect',
    r'DataStructs\.TanimotoSimilarity'
]

# Mapping from RDKit functions to wrapper functions
FUNCTION_MAPPING = {
    'Chem.MolFromSmiles': 'rdkit_wrapper.create_molecule_from_smiles',
    'Chem.MolToSmiles': 'mol_to_smiles(mol)',  # Will need special handling
    'Descriptors.MolWt': 'properties["molecular_weight"]',  # Use with calculate_properties
    'Descriptors.MolLogP': 'properties["logp"]',  # Use with calculate_properties
    'MolSurf.TPSA': 'properties["tpsa"]',  # Use with calculate_properties
    'Lipinski.NumHDonors': 'properties["h_donors"]',  # Use with calculate_properties
    'Lipinski.NumHAcceptors': 'properties["h_acceptors"]',  # Use with calculate_properties
    'Descriptors.NumRotatableBonds': 'properties["rotatable_bonds"]',  # Use with calculate_properties
    'Chem.rdMolDescriptors.CalcNumRings': 'properties["ring_count"]',  # Use with calculate_properties
    'MolDraw2DSVG': 'rdkit_wrapper.generate_molecule_svg',
    'AllChem.GetMorganFingerprintAsBitVect': 'rdkit_wrapper.generate_fingerprint',
    'DataStructs.TanimotoSimilarity': 'rdkit_wrapper.calculate_similarity'
}

def scan_file(file_path: str) -> Dict[str, List[int]]:
    """
    Scan a file for RDKit import and usage patterns.
    Returns a dictionary mapping patterns to line numbers.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    results = {}
    
    # Check for import patterns
    for pattern in DIRECT_IMPORT_PATTERNS:
        matches = []
        for i, line in enumerate(lines):
            if re.search(pattern, line):
                matches.append(i + 1)  # Line numbers are 1-based
        if matches:
            results[pattern] = matches
    
    # Check for usage patterns
    for pattern in RDKIT_USAGE_PATTERNS:
        matches = []
        for i, line in enumerate(lines):
            if re.search(pattern, line):
                matches.append(i + 1)  # Line numbers are 1-based
        if matches:
            results[pattern] = matches
    
    return results

def generate_suggestions(file_path: str, results: Dict[str, List[int]]) -> List[str]:
    """
    Generate suggestions based on the scan results.
    """
    suggestions = []
    
    if any(pattern in results for pattern in DIRECT_IMPORT_PATTERNS):
        suggestions.append(f"Replace direct RDKit imports with: from rdkit_wrapper import *")
    
    for pattern in RDKIT_USAGE_PATTERNS:
        if pattern in results:
            function_name = pattern.replace('\\', '')
            if function_name in FUNCTION_MAPPING:
                replacement = FUNCTION_MAPPING[function_name]
                if 'properties' in replacement:
                    suggestions.append(
                        f"Replace {function_name} with: properties = rdkit_wrapper.calculate_properties(mol_or_smiles) and then use {replacement}"
                    )
                else:
                    suggestions.append(f"Replace {function_name} with: {replacement}")
    
    return suggestions

def apply_changes(file_path: str, backup: bool = True) -> Tuple[bool, List[str]]:
    """
    Apply suggested changes directly to the file.
    Returns a tuple of (success, changes_made).
    """
    with open(file_path, 'r') as f:
        content = f.read()
    
    if backup:
        backup_path = f"{file_path}.bak"
        with open(backup_path, 'w') as f:
            f.write(content)
    
    changes_made = []
    
    # Replace imports
    for pattern in DIRECT_IMPORT_PATTERNS:
        if re.search(pattern, content):
            new_content = re.sub(pattern, '# \\g<0> # Replaced by rdkit_wrapper', content)
            if new_content != content:
                changes_made.append(f"Commented out direct RDKit imports")
                content = new_content
    
    # Add wrapper import
    if changes_made:
        import_line = 'from rdkit_wrapper import *  # Unified RDKit wrapper\n'
        if 'import ' in content:
            # Find the last import line
            import_lines = re.findall(r'^.*import.*$', content, re.MULTILINE)
            if import_lines:
                last_import = import_lines[-1]
                content = content.replace(last_import, last_import + '\n' + import_line)
                changes_made.append("Added rdkit_wrapper import")
        else:
            # Add import at the top
            content = import_line + content
            changes_made.append("Added rdkit_wrapper import at the top")
    
    # Write modified content
    with open(file_path, 'w') as f:
        f.write(content)
    
    return (True, changes_made)

def main():
    parser = argparse.ArgumentParser(description='Update code to use rdkit_wrapper')
    parser.add_argument('path', help='File or directory to scan')
    parser.add_argument('--apply', action='store_true', help='Apply suggested changes')
    parser.add_argument('--no-backup', action='store_true', help='Do not create backup files')
    args = parser.parse_args()
    
    path = args.path
    apply_changes_flag = args.apply
    create_backup = not args.no_backup
    
    if os.path.isfile(path):
        files = [path]
    elif os.path.isdir(path):
        files = []
        for root, _, filenames in os.walk(path):
            for filename in filenames:
                if filename.endswith('.py'):
                    files.append(os.path.join(root, filename))
    else:
        print(f"Error: {path} is not a valid file or directory")
        sys.exit(1)
    
    total_files = len(files)
    files_with_rdkit = 0
    files_changed = 0
    
    for file_path in files:
        print(f"Scanning {file_path}...")
        results = scan_file(file_path)
        
        if results:
            files_with_rdkit += 1
            print(f"  Found RDKit usage in {file_path}")
            
            suggestions = generate_suggestions(file_path, results)
            for suggestion in suggestions:
                print(f"  - {suggestion}")
            
            if apply_changes_flag:
                success, changes = apply_changes(file_path, backup=create_backup)
                if success:
                    files_changed += 1
                    print(f"  Applied changes:")
                    for change in changes:
                        print(f"    - {change}")
                    if create_backup:
                        print(f"  Backup created at {file_path}.bak")
                else:
                    print(f"  Failed to apply changes to {file_path}")
    
    print(f"\nSummary:")
    print(f"  Total files scanned: {total_files}")
    print(f"  Files with RDKit usage: {files_with_rdkit}")
    
    if apply_changes_flag:
        print(f"  Files modified: {files_changed}")
    else:
        print("\nTo apply changes, run with --apply flag")

if __name__ == "__main__":
    main()