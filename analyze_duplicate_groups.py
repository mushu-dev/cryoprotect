#!/usr/bin/env python3
"""
Analyze duplicate groups in detail to determine consolidation strategy.

This script:
1. Retrieves molecules that have been marked as duplicates
2. Analyzes each duplicate group in detail
3. Checks properties, structure, and relationships
4. Classifies duplicates as true duplicates or variants
5. Generates a detailed report for consolidation planning
"""

import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
from collections import defaultdict

# Load environment variables
load_dotenv()

# Constants
DUPLICATES_REPORT_FILE = "duplicate_molecules_report.json"
DUPLICATE_ANALYSIS_FILE = "duplicate_groups_analysis.json"
CONSOLIDATION_PLAN_FILE = "duplicate_consolidation_plan.json"

# Database connection
def connect_to_db():
    """Connect to the database."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }
    return psycopg2.connect(**db_params)

def get_duplicate_molecules_by_property(conn):
    """
    Get molecules that have been marked as duplicates based on the
    properties field containing a duplicate_group value.
    """
    duplicate_molecules = []
    duplicate_groups = {}
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, molecular_formula, pubchem_cid, smiles,
                  data_source, properties, created_at, updated_at
            FROM molecules
            WHERE properties ? 'duplicate_group'
            ORDER BY (properties->>'duplicate_group'), name
        """)
        
        for molecule in cursor.fetchall():
            duplicate_molecules.append(molecule)
            
            # Organize by group
            group_id = molecule['properties']['duplicate_group']
            if group_id not in duplicate_groups:
                duplicate_groups[group_id] = {
                    'group_id': group_id,
                    'duplicate_type': molecule['properties'].get('duplicate_type'),
                    'duplicate_value': molecule['properties'].get('duplicate_value'),
                    'molecules': []
                }
            
            duplicate_groups[group_id]['molecules'].append(molecule)
    
    return duplicate_molecules, duplicate_groups

def get_molecule_relationships(conn, molecule_id):
    """Get all relationships for a molecule."""
    relationships = {
        'properties': [],
        'mixture_components': [],
        'experiments': [],
        'predictions': []
    }
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Check molecular properties
        cursor.execute("""
            SELECT id, property_type, property_value, unit, source
            FROM molecular_properties
            WHERE molecule_id = %s
        """, (molecule_id,))
        relationships['properties'] = cursor.fetchall()
        
        # Check mixture components
        cursor.execute("""
            SELECT id, mixture_id, concentration, units
            FROM mixture_components
            WHERE molecule_id = %s
        """, (molecule_id,))
        relationships['mixture_components'] = cursor.fetchall()

        # Check experiment results linked to this molecule
        cursor.execute("""
            SELECT id FROM experiment_results
            WHERE molecule_id = %s
        """, (molecule_id,))

        experiment_ids = [row['id'] for row in cursor.fetchall()]
        experiment_properties = []

        # If there are related experiment results, get their details
        if experiment_ids:
            # Use IN clause with tuple format
            cursor.execute("""
                SELECT id, protocol_id, molecule_id, mixture_id, tissue_type_id,
                        concentration, concentration_unit, viability_percentage,
                        recovery_rate, functionality_score
                FROM experiment_results
                WHERE id IN %s
            """, (tuple(experiment_ids),))
            experiment_properties = cursor.fetchall()

        relationships['experiments'] = experiment_properties

        # Check predictions
        cursor.execute("""
            SELECT id, property_type_id, calculation_method_id, numeric_value, text_value, boolean_value, confidence
            FROM predictions
            WHERE molecule_id = %s
        """, (molecule_id,))
        relationships['predictions'] = cursor.fetchall()
    
    return relationships

def analyze_duplicate_group(conn, group):
    """
    Analyze a duplicate group to determine if they are true duplicates
    or variants, and identify the best strategy for consolidation.
    """
    molecules = group['molecules']
    analysis = {
        'group_id': group['group_id'],
        'duplicate_type': group['duplicate_type'],
        'duplicate_value': group['duplicate_value'],
        'molecule_count': len(molecules),
        'pubchem_cids': set(),
        'has_formula_differences': False,
        'has_smiles_differences': False,
        'has_property_differences': False,
        'has_relationship_conflicts': False,
        'molecules_with_relationships': 0,
        'relationship_counts': {
            'properties': 0,
            'mixture_components': 0,
            'experiments': 0,
            'predictions': 0
        },
        'primary_candidate': None,
        'consolidation_type': None,
        'molecules_detail': []
    }
    
    # Check for differences in chemical identity
    formulas = set()
    smiles_strings = set()
    
    for molecule in molecules:
        # Track PubChem CIDs
        if molecule['pubchem_cid']:
            analysis['pubchem_cids'].add(molecule['pubchem_cid'])
        
        # Track formulas and SMILES
        if molecule['molecular_formula']:
            formulas.add(molecule['molecular_formula'])
        
        if molecule['smiles']:
            smiles_strings.add(molecule['smiles'])
        
        # Get relationships
        relationships = get_molecule_relationships(conn, molecule['id'])
        
        # Track relationship counts
        total_relationships = sum(len(rel) for rel in relationships.values())
        has_relationships = total_relationships > 0
        
        if has_relationships:
            analysis['molecules_with_relationships'] += 1
            for rel_type, rel_items in relationships.items():
                analysis['relationship_counts'][rel_type] += len(rel_items)
        
        # Add detailed molecule information
        molecule_detail = {
            'id': molecule['id'],
            'name': molecule['name'],
            'molecular_formula': molecule['molecular_formula'],
            'pubchem_cid': molecule['pubchem_cid'],
            'has_relationships': has_relationships,
            'relationship_count': total_relationships,
            'relationships': relationships,
            'data_source': molecule['data_source'],
            'created_at': molecule['created_at'],
            'updated_at': molecule['updated_at']
        }
        
        analysis['molecules_detail'].append(molecule_detail)
    
    # Check for differences
    analysis['has_formula_differences'] = len(formulas) > 1
    analysis['has_smiles_differences'] = len(smiles_strings) > 1
    
    # Property differences are already captured in the has_property_differences flag
    
    # Determine consolidation type
    if analysis['has_formula_differences'] or analysis['has_smiles_differences']:
        # Different chemical structures - these are variants, not duplicates
        analysis['consolidation_type'] = 'DIFFERENTIATE'
    elif len(analysis['pubchem_cids']) > 1:
        # Multiple PubChem CIDs - may need special handling due to uniqueness constraint
        analysis['consolidation_type'] = 'SELECTIVE_MERGE'
    elif analysis['molecules_with_relationships'] == 0:
        # No relationships - can safely merge or delete duplicates
        analysis['consolidation_type'] = 'SAFE_MERGE'
    elif analysis['molecules_with_relationships'] == 1:
        # Only one molecule has relationships - can designate it as primary
        analysis['consolidation_type'] = 'PRIMARY_SELECTION'
    else:
        # Multiple molecules have relationships - need to analyze conflicts
        analysis['consolidation_type'] = 'COMPLEX_MERGE'
    
    # Identify primary candidate based on completeness and relationships
    primary_candidates = sorted(
        analysis['molecules_detail'],
        key=lambda m: (
            m['pubchem_cid'] is not None,  # Prefer molecules with PubChem CIDs
            m['molecular_formula'] is not None,  # Prefer molecules with formulas
            m['relationship_count'],  # Prefer molecules with more relationships
            m['created_at']  # Prefer older molecules as a tiebreaker
        ),
        reverse=True
    )
    
    if primary_candidates:
        analysis['primary_candidate'] = primary_candidates[0]['id']
    
    return analysis

def generate_consolidation_plan(analyses):
    """
    Generate a plan for consolidating duplicate groups based on the analyses.
    """
    plan = {
        'safe_merge': [],
        'selective_merge': [],
        'primary_selection': [],
        'complex_merge': [],
        'differentiate': [],
        'summary': {
            'total_groups': len(analyses),
            'safe_merge_count': 0,
            'selective_merge_count': 0,
            'primary_selection_count': 0,
            'complex_merge_count': 0,
            'differentiate_count': 0
        }
    }
    
    for group_id, analysis in analyses.items():
        consolidation_type = analysis['consolidation_type'].lower()
        
        # Convert to snake_case for dictionary keys
        if consolidation_type == 'safe_merge':
            plan['safe_merge'].append(analysis)
            plan['summary']['safe_merge_count'] += 1
        elif consolidation_type == 'selective_merge':
            plan['selective_merge'].append(analysis)
            plan['summary']['selective_merge_count'] += 1
        elif consolidation_type == 'primary_selection':
            plan['primary_selection'].append(analysis)
            plan['summary']['primary_selection_count'] += 1
        elif consolidation_type == 'complex_merge':
            plan['complex_merge'].append(analysis)
            plan['summary']['complex_merge_count'] += 1
        elif consolidation_type == 'differentiate':
            plan['differentiate'].append(analysis)
            plan['summary']['differentiate_count'] += 1
    
    return plan

def main():
    """
    Analyze duplicate groups and generate a consolidation plan.
    """
    import argparse
    parser = argparse.ArgumentParser(description="Analyze duplicate groups for consolidation planning")
    parser.add_argument("--output-dir", default=".", help="Directory to save output files")
    args = parser.parse_args()
    
    output_dir = args.output_dir
    
    # Connect to database
    print("Connecting to database...")
    conn = connect_to_db()
    
    try:
        # Get duplicate molecules and groups
        print("Retrieving marked duplicate molecules...")
        duplicate_molecules, duplicate_groups = get_duplicate_molecules_by_property(conn)
        
        print(f"Found {len(duplicate_molecules)} molecules in {len(duplicate_groups)} duplicate groups")
        
        # Analyze each duplicate group
        print("Analyzing duplicate groups...")
        group_analyses = {}
        
        for group_id, group in duplicate_groups.items():
            print(f"  Analyzing group {group_id} ({len(group['molecules'])} molecules)")
            analysis = analyze_duplicate_group(conn, group)
            group_analyses[group_id] = analysis
        
        # Generate consolidation plan
        print("Generating consolidation plan...")
        consolidation_plan = generate_consolidation_plan(group_analyses)
        
        # Save results
        analyses_path = os.path.join(output_dir, DUPLICATE_ANALYSIS_FILE)
        plan_path = os.path.join(output_dir, CONSOLIDATION_PLAN_FILE)
        
        with open(analyses_path, 'w') as f:
            json.dump(group_analyses, f, indent=2, default=str)
        
        with open(plan_path, 'w') as f:
            json.dump(consolidation_plan, f, indent=2, default=str)
        
        print(f"Analysis saved to {analyses_path}")
        print(f"Consolidation plan saved to {plan_path}")
        
        # Print summary
        print("\nConsolidation Plan Summary:")
        print(f"Total duplicate groups: {consolidation_plan['summary']['total_groups']}")
        print(f"  Safe merge groups: {consolidation_plan['summary']['safe_merge_count']}")
        print(f"  Selective merge groups: {consolidation_plan['summary']['selective_merge_count']}")
        print(f"  Primary selection groups: {consolidation_plan['summary']['primary_selection_count']}")
        print(f"  Complex merge groups: {consolidation_plan['summary']['complex_merge_count']}")
        print(f"  Differentiate groups: {consolidation_plan['summary']['differentiate_count']}")
    
    except Exception as e:
        print(f"Error analyzing duplicate groups: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()