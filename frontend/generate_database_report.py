#!/usr/bin/env python3
"""
Script to generate a comprehensive report on the CryoProtect database.
This helps identify issues and track improvements.
"""

import os
import sys
import json
import argparse
from datetime import datetime
import logging
from psycopg2.extras import RealDictCursor
import db_utils

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def get_table_counts():
    """Get row counts for all tables in the database."""
    query = """
    SELECT 
        tablename as table_name,
        pg_total_relation_size(quote_ident(tablename)) as table_size_bytes,
        (SELECT COUNT(*) FROM information_schema.columns WHERE table_name=tables.tablename) as column_count,
        (SELECT COUNT(*) FROM pg_indexes WHERE tablename=tables.tablename) as index_count,
        (SELECT reltuples::bigint FROM pg_class WHERE oid = (quote_ident(tablename))::regclass) as estimated_row_count
    FROM 
        pg_tables as tables
    WHERE 
        schemaname = 'public'
    ORDER BY 
        table_name;
    """
    try:
        return db_utils.execute_query(query, cursor_factory=RealDictCursor)
    except Exception as e:
        logger.error(f"Error getting table counts: {e}")
        return []

def get_table_schemas():
    """Get schema information for all tables."""
    query = """
    SELECT 
        table_name,
        column_name,
        data_type,
        character_maximum_length,
        column_default,
        is_nullable,
        (
            SELECT 
                pg_description.description 
            FROM 
                pg_description 
            JOIN 
                pg_class ON pg_description.objoid = pg_class.oid 
            JOIN 
                pg_attribute ON pg_attribute.attrelid = pg_class.oid AND pg_attribute.attnum = pg_description.objsubid 
            WHERE 
                pg_class.relname = columns.table_name AND pg_attribute.attname = columns.column_name
        ) as column_description
    FROM 
        information_schema.columns
    WHERE 
        table_schema = 'public'
    ORDER BY 
        table_name, 
        ordinal_position;
    """
    try:
        columns = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        # Group by table_name
        tables = {}
        for col in columns:
            table_name = col['table_name']
            if table_name not in tables:
                tables[table_name] = []
            tables[table_name].append(col)
        
        return tables
    except Exception as e:
        logger.error(f"Error getting table schemas: {e}")
        return {}

def get_foreign_keys():
    """Get foreign key relationships between tables."""
    query = """
    SELECT
        tc.table_schema, 
        tc.constraint_name, 
        tc.table_name, 
        kcu.column_name, 
        ccu.table_schema AS foreign_table_schema,
        ccu.table_name AS foreign_table_name,
        ccu.column_name AS foreign_column_name 
    FROM 
        information_schema.table_constraints AS tc 
    JOIN 
        information_schema.key_column_usage AS kcu ON tc.constraint_name = kcu.constraint_name
        AND tc.table_schema = kcu.table_schema
    JOIN 
        information_schema.constraint_column_usage AS ccu ON ccu.constraint_name = tc.constraint_name
        AND ccu.table_schema = tc.table_schema
    WHERE 
        tc.constraint_type = 'FOREIGN KEY' 
        AND tc.table_schema = 'public'
    ORDER BY 
        tc.table_name, 
        kcu.column_name;
    """
    try:
        return db_utils.execute_query(query, cursor_factory=RealDictCursor)
    except Exception as e:
        logger.error(f"Error getting foreign keys: {e}")
        return []

def get_indexes():
    """Get index information for all tables."""
    query = """
    SELECT
        t.relname AS table_name,
        i.relname AS index_name,
        a.attname AS column_name,
        ix.indisunique AS is_unique,
        ix.indisprimary AS is_primary,
        pg_get_indexdef(ix.indexrelid) AS index_definition
    FROM
        pg_index ix
    JOIN
        pg_class i ON i.oid = ix.indexrelid
    JOIN
        pg_class t ON t.oid = ix.indrelid
    JOIN
        pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(ix.indkey)
    JOIN
        pg_namespace n ON n.oid = t.relnamespace
    WHERE
        n.nspname = 'public'
    ORDER BY
        t.relname,
        i.relname,
        a.attnum;
    """
    try:
        indexes = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        # Group by table_name and index_name
        tables = {}
        for idx in indexes:
            table_name = idx['table_name']
            index_name = idx['index_name']
            
            if table_name not in tables:
                tables[table_name] = {}
            
            if index_name not in tables[table_name]:
                tables[table_name][index_name] = {
                    'index_name': index_name,
                    'is_unique': idx['is_unique'],
                    'is_primary': idx['is_primary'],
                    'index_definition': idx['index_definition'],
                    'columns': []
                }
            
            tables[table_name][index_name]['columns'].append(idx['column_name'])
        
        return tables
    except Exception as e:
        logger.error(f"Error getting indexes: {e}")
        return {}

def get_missing_indexes():
    """Identify potential missing indexes based on foreign keys and common query patterns."""
    potential_missing_indexes = []
    
    # Get foreign keys - these should generally have indexes
    fks = get_foreign_keys()
    indexes = get_indexes()
    
    for fk in fks:
        table_name = fk['table_name']
        column_name = fk['column_name']
        
        # Check if there's an index on this FK column
        has_index = False
        if table_name in indexes:
            for idx_name in indexes[table_name]:
                idx = indexes[table_name][idx_name]
                if column_name in idx['columns'] and idx['columns'][0] == column_name:
                    has_index = True
                    break
        
        if not has_index:
            potential_missing_indexes.append({
                'table_name': table_name,
                'column_name': column_name,
                'reason': 'Foreign key without index',
                'suggested_sql': f'CREATE INDEX idx_{table_name}_{column_name} ON {table_name}({column_name});'
            })
    
    # Check for specific columns that are likely to be queried frequently
    common_columns = [
        {'table': 'molecules', 'column': 'name'},
        {'table': 'molecules', 'column': 'pubchem_cid'},
        {'table': 'molecules', 'column': 'smiles'},
        {'table': 'consolidated_molecules', 'column': 'name'},
        {'table': 'consolidated_molecules', 'column': 'pubchem_cid'},
        {'table': 'consolidated_molecules', 'column': 'smiles'},
        {'table': 'molecular_properties', 'column': 'property_name'},
        {'table': 'molecular_properties', 'column': 'property_type'},
        {'table': 'mixtures', 'column': 'name'},
    ]
    
    for col in common_columns:
        table_name = col['table']
        column_name = col['column']
        
        # Skip if table doesn't exist in indexes dict
        if table_name not in indexes:
            continue
            
        # Check if there's an index on this column
        has_index = False
        for idx_name in indexes[table_name]:
            idx = indexes[table_name][idx_name]
            if column_name in idx['columns'] and idx['columns'][0] == column_name:
                has_index = True
                break
        
        if not has_index:
            potential_missing_indexes.append({
                'table_name': table_name,
                'column_name': column_name,
                'reason': 'Frequently queried column without index',
                'suggested_sql': f'CREATE INDEX idx_{table_name}_{column_name} ON {table_name}({column_name});'
            })
    
    return potential_missing_indexes

def get_duplicate_stats():
    """Get statistics on duplicate rows in key tables."""
    duplicate_stats = {}
    
    # Check for duplicate SMILES in consolidated_molecules
    query = """
    WITH duplicate_check AS (
      SELECT 
        smiles,
        COUNT(*) as count
      FROM 
        consolidated_molecules
      WHERE 
        smiles IS NOT NULL
      GROUP BY 
        smiles
      HAVING 
        COUNT(*) > 1
    )
    SELECT 
      COUNT(*) as total_duplicates,
      SUM(count) as total_duplicate_records,
      MAX(count) as max_duplicates_for_one_molecule
    FROM 
      duplicate_check;
    """
    try:
        duplicate_stats['consolidated_molecules_smiles'] = db_utils.execute_query(
            query, cursor_factory=RealDictCursor)[0]
    except Exception as e:
        logger.error(f"Error getting duplicate SMILES stats: {e}")
    
    # Check for duplicate pubchem_cid in consolidated_molecules
    query = """
    WITH duplicate_check AS (
      SELECT 
        pubchem_cid,
        COUNT(*) as count
      FROM 
        consolidated_molecules
      WHERE 
        pubchem_cid IS NOT NULL
      GROUP BY 
        pubchem_cid
      HAVING 
        COUNT(*) > 1
    )
    SELECT 
      COUNT(*) as total_duplicates,
      SUM(count) as total_duplicate_records,
      MAX(count) as max_duplicates_for_one_molecule
    FROM 
      duplicate_check;
    """
    try:
        duplicate_stats['consolidated_molecules_pubchem_cid'] = db_utils.execute_query(
            query, cursor_factory=RealDictCursor)[0]
    except Exception as e:
        logger.error(f"Error getting duplicate pubchem_cid stats: {e}")
    
    return duplicate_stats

def get_data_completeness():
    """Get statistics on data completeness in key tables."""
    completeness_stats = {}
    
    # Check completeness of consolidated_molecules
    query = """
    SELECT 
      COUNT(*) as total_molecules,
      SUM(CASE WHEN smiles IS NULL THEN 1 ELSE 0 END) as missing_smiles,
      SUM(CASE WHEN formula IS NULL THEN 1 ELSE 0 END) as missing_formula,
      SUM(CASE WHEN molecular_weight IS NULL THEN 1 ELSE 0 END) as missing_weight,
      SUM(CASE WHEN pubchem_cid IS NULL THEN 1 ELSE 0 END) as missing_pubchem_cid,
      SUM(CASE WHEN primary_molecule_id IS NULL THEN 1 ELSE 0 END) as missing_primary_id
    FROM 
      consolidated_molecules;
    """
    try:
        completeness_stats['consolidated_molecules'] = db_utils.execute_query(
            query, cursor_factory=RealDictCursor)[0]
    except Exception as e:
        logger.error(f"Error getting molecule completeness stats: {e}")
    
    # Check completeness of molecular_properties
    query = """
    SELECT 
      COUNT(*) as total_properties,
      SUM(CASE WHEN property_type_id IS NULL THEN 1 ELSE 0 END) as missing_property_type,
      SUM(CASE WHEN property_name IS NULL THEN 1 ELSE 0 END) as missing_property_name,
      SUM(CASE WHEN numeric_value IS NULL AND text_value IS NULL AND boolean_value IS NULL THEN 1 ELSE 0 END) 
        as missing_all_values,
      COUNT(DISTINCT molecule_id) as distinct_molecules,
      AVG(COUNT(*)) OVER (PARTITION BY molecule_id) as avg_properties_per_molecule
    FROM 
      molecular_properties
    GROUP BY 
      molecule_id;
    """
    try:
        property_stats = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        if property_stats:
            # Take the first row for most stats, but calculate avg_properties_per_molecule correctly
            completeness_stats['molecular_properties'] = property_stats[0]
            
            # Calculate the average properties per molecule
            if property_stats and len(property_stats) > 0:
                avg_props = sum(row['avg_properties_per_molecule'] for row in property_stats) / len(property_stats)
                completeness_stats['molecular_properties']['avg_properties_per_molecule'] = avg_props
    except Exception as e:
        logger.error(f"Error getting property completeness stats: {e}")
    
    return completeness_stats

def get_rls_policy_stats():
    """Get statistics on Row Level Security policies."""
    query = """
    SELECT 
        schemaname, 
        tablename, 
        policyname, 
        roles, 
        cmd, 
        qual
    FROM 
        pg_policies
    WHERE 
        schemaname = 'public'
    ORDER BY 
        tablename, 
        policyname;
    """
    try:
        policies = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        # Analyze policies for potential issues
        analysis = {
            'total_policies': len(policies),
            'tables_with_policies': len(set(p['tablename'] for p in policies)),
            'potential_duplicate_policies': [],
            'tables_without_service_role': [],
            'tables_with_complex_policies': []
        }
        
        # Group policies by table
        policies_by_table = {}
        for policy in policies:
            table = policy['tablename']
            if table not in policies_by_table:
                policies_by_table[table] = []
            policies_by_table[table].append(policy)
        
        # Find potential duplicate policies (same table, same cmd, similar qual)
        for table, table_policies in policies_by_table.items():
            for i, policy1 in enumerate(table_policies):
                for j, policy2 in enumerate(table_policies):
                    if i >= j:
                        continue
                    
                    if policy1['cmd'] == policy2['cmd'] and ('service_role' in policy1['qual'] and 'service_role' in policy2['qual']):
                        analysis['potential_duplicate_policies'].append({
                            'table': table,
                            'policy1': policy1['policyname'],
                            'policy2': policy2['policyname'],
                            'cmd': policy1['cmd']
                        })
            
            # Check if this table has service_role policy
            has_service_role = any('service_role' in p['qual'] for p in table_policies)
            if not has_service_role:
                analysis['tables_without_service_role'].append(table)
            
            # Check for complex policies
            complex_policies = [p for p in table_policies if p['qual'] and len(p['qual']) > 200]
            if complex_policies:
                analysis['tables_with_complex_policies'].append({
                    'table': table,
                    'complex_policies': [p['policyname'] for p in complex_policies]
                })
        
        return {
            'policies': policies,
            'analysis': analysis
        }
    except Exception as e:
        logger.error(f"Error getting RLS policy stats: {e}")
        return {
            'policies': [],
            'analysis': {
                'total_policies': 0,
                'tables_with_policies': 0,
                'potential_duplicate_policies': [],
                'tables_without_service_role': [],
                'tables_with_complex_policies': []
            }
        }

def save_report(report, filename_prefix="database_analysis_report"):
    """Save the report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.json"
    
    try:
        with open(filename, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {filename}")
        return filename
    except Exception as e:
        logger.error(f"Error saving report: {e}")
        return None

def generate_recommendations(report):
    """Generate recommendations based on the database report."""
    recommendations = []
    
    # Duplicate molecules recommendations
    if report.get('duplicate_stats', {}).get('consolidated_molecules_smiles', {}).get('total_duplicates', 0) > 0:
        recommendations.append({
            'priority': 'high',
            'category': 'data_quality',
            'issue': 'Duplicate molecules exist in consolidated_molecules table',
            'recommendation': 'Run the molecule consolidation script to merge duplicates',
            'suggested_action': 'python3 apply_molecule_consolidation.py'
        })
    
    # Missing data recommendations
    data_completeness = report.get('data_completeness', {}).get('consolidated_molecules', {})
    if data_completeness.get('missing_formula', 0) > 0 or data_completeness.get('missing_weight', 0) > 0:
        recommendations.append({
            'priority': 'medium',
            'category': 'data_quality',
            'issue': f"Missing molecular data: {data_completeness.get('missing_formula', 0)} missing formulas, {data_completeness.get('missing_weight', 0)} missing molecular weights",
            'recommendation': 'Run a script to calculate missing properties from SMILES using RDKit',
            'suggested_action': 'Create and run a property calculation script'
        })
    
    # Foreign key recommendations
    if len(report.get('missing_indexes', [])) > 0:
        fk_indexes = [idx for idx in report.get('missing_indexes', []) if idx['reason'] == 'Foreign key without index']
        if fk_indexes:
            recommendations.append({
                'priority': 'high',
                'category': 'performance',
                'issue': f"Missing indexes on {len(fk_indexes)} foreign key columns",
                'recommendation': 'Create indexes on foreign key columns to improve join performance',
                'suggested_action': '\n'.join(idx['suggested_sql'] for idx in fk_indexes)
            })
    
    # Duplicate RLS policies
    rls_analysis = report.get('rls_policy_stats', {}).get('analysis', {})
    if len(rls_analysis.get('potential_duplicate_policies', [])) > 0:
        recommendations.append({
            'priority': 'medium',
            'category': 'security',
            'issue': f"Found {len(rls_analysis.get('potential_duplicate_policies', []))} potentially duplicate RLS policies",
            'recommendation': 'Consolidate duplicate RLS policies to improve maintainability and performance',
            'suggested_action': 'Create a migration to consolidate RLS policies'
        })
    
    # Complex RLS policies
    if len(rls_analysis.get('tables_with_complex_policies', [])) > 0:
        recommendations.append({
            'priority': 'low',
            'category': 'security',
            'issue': f"Found complex RLS policies on {len(rls_analysis.get('tables_with_complex_policies', []))} tables",
            'recommendation': 'Simplify complex RLS policies by creating reusable functions',
            'suggested_action': 'Review and refactor complex RLS policies'
        })
    
    return recommendations

def main():
    """Main function to generate a comprehensive database report."""
    parser = argparse.ArgumentParser(description="Generate a comprehensive database report.")
    parser.add_argument("--no-save", action="store_true", help="Don't save the report to a file")
    parser.add_argument("--summary", action="store_true", help="Print a summary instead of the full report")
    args = parser.parse_args()
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    # Generate report
    logger.info("Generating database report...")
    report = {
        'timestamp': datetime.now().isoformat(),
        'table_counts': get_table_counts(),
        'table_schemas': get_table_schemas(),
        'foreign_keys': get_foreign_keys(),
        'indexes': get_indexes(),
        'missing_indexes': get_missing_indexes(),
        'duplicate_stats': get_duplicate_stats(),
        'data_completeness': get_data_completeness(),
        'rls_policy_stats': get_rls_policy_stats()
    }
    
    # Generate recommendations
    report['recommendations'] = generate_recommendations(report)
    
    # Save report
    if not args.no_save:
        filename = save_report(report)
        if filename:
            print(f"Report saved to {filename}")
    
    # Print summary or full report
    if args.summary:
        print("\n=== Database Report Summary ===")
        print(f"Tables: {len(report['table_counts'])}")
        print(f"Foreign Keys: {len(report['foreign_keys'])}")
        print(f"Potential Missing Indexes: {len(report['missing_indexes'])}")
        
        duplicate_smiles = report['duplicate_stats'].get('consolidated_molecules_smiles', {})
        print(f"Duplicate Molecules: {duplicate_smiles.get('total_duplicates', 0)} molecules with duplicates, {duplicate_smiles.get('total_duplicate_records', 0)} total duplicate records")
        
        data_completeness = report['data_completeness'].get('consolidated_molecules', {})
        print(f"Data Completeness: {data_completeness.get('total_molecules', 0)} molecules, {data_completeness.get('missing_formula', 0)} missing formulas, {data_completeness.get('missing_weight', 0)} missing weights")
        
        rls_analysis = report['rls_policy_stats']['analysis']
        print(f"RLS Policies: {rls_analysis['total_policies']} policies on {rls_analysis['tables_with_policies']} tables, {len(rls_analysis['potential_duplicate_policies'])} potential duplicates")
        
        print("\n=== Recommendations ===")
        for rec in report['recommendations']:
            print(f"[{rec['priority'].upper()}] {rec['issue']}")
            print(f"  => {rec['recommendation']}")
    else:
        print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()