#!/usr/bin/env python3
"""
Script to analyze database query performance and identify opportunities for indexing.
This tool helps identify slow queries and suggests indexes to improve performance.
"""

import os
import sys
import logging
import json
from datetime import datetime
import re
from dotenv import load_dotenv
import db_utils

# Load environment variables
load_dotenv()

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def get_table_column_info():
    """Get information about tables and columns in the database."""
    try:
        query = """
        SELECT 
            t.table_schema,
            t.table_name,
            c.column_name,
            c.data_type,
            tc.constraint_type
        FROM 
            information_schema.tables t
        LEFT JOIN 
            information_schema.columns c 
            ON t.table_schema = c.table_schema 
            AND t.table_name = c.table_name
        LEFT JOIN 
            information_schema.key_column_usage kcu
            ON c.table_schema = kcu.table_schema
            AND c.table_name = kcu.table_name
            AND c.column_name = kcu.column_name
        LEFT JOIN 
            information_schema.table_constraints tc
            ON kcu.constraint_schema = tc.constraint_schema
            AND kcu.constraint_name = tc.constraint_name
        WHERE 
            t.table_schema = 'public'
            AND t.table_type = 'BASE TABLE'
        ORDER BY 
            t.table_schema,
            t.table_name,
            c.ordinal_position;
        """
        
        results = db_utils.execute_query(query, cursor_factory=db_utils.RealDictCursor)
        
        # Organize by table
        tables = {}
        for row in results:
            table_name = row['table_name']
            if table_name not in tables:
                tables[table_name] = {
                    'name': table_name,
                    'schema': row['table_schema'],
                    'columns': []
                }
            
            # Add column info
            constraint_type = row['constraint_type']
            is_primary_key = constraint_type == 'PRIMARY KEY'
            is_foreign_key = constraint_type == 'FOREIGN KEY'
            
            tables[table_name]['columns'].append({
                'name': row['column_name'],
                'data_type': row['data_type'],
                'is_primary_key': is_primary_key,
                'is_foreign_key': is_foreign_key,
                'has_constraint': constraint_type is not None
            })
        
        return tables
    except Exception as e:
        logger.error(f"Error getting table schema: {e}")
        return None

def get_existing_indexes():
    """Get information about existing indexes in the database."""
    try:
        query = """
        SELECT 
            n.nspname AS schema_name,
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
            pg_namespace n ON n.oid = t.relnamespace
        JOIN 
            pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(ix.indkey)
        WHERE 
            n.nspname = 'public'
        ORDER BY 
            n.nspname,
            t.relname,
            i.relname,
            a.attnum;
        """
        
        results = db_utils.execute_query(query, cursor_factory=db_utils.RealDictCursor)
        
        # Organize by table and index
        indexes = {}
        for row in results:
            table_name = row['table_name']
            index_name = row['index_name']
            
            if table_name not in indexes:
                indexes[table_name] = {}
            
            if index_name not in indexes[table_name]:
                indexes[table_name][index_name] = {
                    'name': index_name,
                    'columns': [],
                    'is_unique': row['is_unique'],
                    'is_primary': row['is_primary'],
                    'definition': row['index_definition']
                }
            
            indexes[table_name][index_name]['columns'].append(row['column_name'])
        
        return indexes
    except Exception as e:
        logger.error(f"Error getting existing indexes: {e}")
        return None

def analyze_query_patterns():
    """Analyze common query patterns using execution plans."""
    try:
        # Define common query patterns to analyze
        test_queries = [
            {
                'name': 'Molecule lookup by ID',
                'query': 'EXPLAIN ANALYZE SELECT * FROM molecules WHERE id = uuid_generate_v4();',
                'tables': ['molecules'],
                'columns': ['id']
            },
            {
                'name': 'Molecule lookup by name',
                'query': 'EXPLAIN ANALYZE SELECT * FROM molecules WHERE name ILIKE \'%glycerol%\';',
                'tables': ['molecules'],
                'columns': ['name']
            },
            {
                'name': 'Property lookup by molecule',
                'query': 'EXPLAIN ANALYZE SELECT * FROM molecular_properties WHERE molecule_id = uuid_generate_v4();',
                'tables': ['molecular_properties'],
                'columns': ['molecule_id']
            },
            {
                'name': 'Property lookup by name',
                'query': 'EXPLAIN ANALYZE SELECT * FROM molecular_properties WHERE property_name = \'molecular_weight\';',
                'tables': ['molecular_properties'],
                'columns': ['property_name']
            },
            {
                'name': 'Molecule search by formula',
                'query': 'EXPLAIN ANALYZE SELECT * FROM molecules WHERE formula = \'C3H8O3\';',
                'tables': ['molecules'],
                'columns': ['formula']
            },
            {
                'name': 'Search molecules by type',
                'query': 'EXPLAIN ANALYZE SELECT * FROM molecules WHERE type = \'cryoprotectant\';',
                'tables': ['molecules'],
                'columns': ['type']
            },
            {
                'name': 'Find molecules by property range',
                'query': 'EXPLAIN ANALYZE SELECT m.* FROM molecules m JOIN molecular_properties mp ON m.id = mp.molecule_id WHERE mp.property_name = \'molecular_weight\' AND mp.numeric_value BETWEEN 50 AND 200;',
                'tables': ['molecules', 'molecular_properties'],
                'columns': ['id', 'molecule_id', 'property_name', 'numeric_value']
            },
            {
                'name': 'Lookup consolidated molecule',
                'query': 'EXPLAIN ANALYZE SELECT * FROM consolidated_molecules WHERE primary_id = uuid_generate_v4();',
                'tables': ['consolidated_molecules'],
                'columns': ['primary_id']
            },
            {
                'name': 'Search by pubchem_cid',
                'query': 'EXPLAIN ANALYZE SELECT * FROM molecules WHERE pubchem_cid = 1;',
                'tables': ['molecules'],
                'columns': ['pubchem_cid']
            },
            {
                'name': 'Join molecules and properties with filtering',
                'query': 'EXPLAIN ANALYZE SELECT m.*, mp.numeric_value FROM molecules m JOIN molecular_properties mp ON m.id = mp.molecule_id WHERE mp.property_name = \'melting_point\' AND m.is_public = TRUE ORDER BY mp.numeric_value DESC LIMIT 10;',
                'tables': ['molecules', 'molecular_properties'],
                'columns': ['id', 'molecule_id', 'property_name', 'is_public', 'numeric_value']
            }
        ]
        
        # Execute each test query to get the execution plan
        query_results = []
        for test in test_queries:
            try:
                result = db_utils.execute_query(test['query'], cursor_factory=db_utils.RealDictCursor)
                
                # Extract execution time from EXPLAIN ANALYZE output
                execution_time = None
                for row in result:
                    for key, value in row.items():
                        if isinstance(value, str) and "execution time" in value.lower():
                            match = re.search(r"execution time: (\d+\.\d+)", value.lower())
                            if match:
                                execution_time = float(match.group(1))
                
                query_results.append({
                    'name': test['name'],
                    'tables': test['tables'],
                    'columns': test['columns'],
                    'execution_time_ms': execution_time,
                    'plan': result
                })
            except Exception as e:
                logger.warning(f"Error executing test query '{test['name']}': {e}")
                query_results.append({
                    'name': test['name'],
                    'tables': test['tables'],
                    'columns': test['columns'],
                    'error': str(e)
                })
        
        return query_results
    except Exception as e:
        logger.error(f"Error analyzing query patterns: {e}")
        return None

def identify_missing_indexes(tables, existing_indexes, query_patterns):
    """Identify columns that would benefit from indexes based on query patterns."""
    # Track potential index candidates
    index_candidates = []
    
    # Tables that we've analyzed
    analyzed_tables = set()
    
    # Check each query pattern for potential index opportunities
    for pattern in query_patterns:
        if 'error' in pattern:
            continue
        
        for table_name in pattern.get('tables', []):
            if table_name not in tables:
                continue
            
            analyzed_tables.add(table_name)
            table = tables[table_name]
            
            # Get columns mentioned in the query
            query_columns = [col for col in pattern.get('columns', []) 
                             if any(c['name'] == col for c in table['columns'])]
            
            # For each column, check if it's already indexed
            table_indexes = existing_indexes.get(table_name, {})
            for column in query_columns:
                is_indexed = any(column in idx['columns'] for idx in table_indexes.values())
                
                # If not indexed, consider adding an index
                if not is_indexed:
                    # Check if column exists in this table
                    column_info = next((c for c in table['columns'] if c['name'] == column), None)
                    if column_info:
                        # Don't suggest indexes on primary keys (they're already indexed)
                        if column_info.get('is_primary_key'):
                            continue
                        
                        # Add to candidates
                        candidate = {
                            'table': table_name,
                            'column': column,
                            'data_type': column_info.get('data_type'),
                            'reason': f"Used in {pattern['name']} query",
                            'is_foreign_key': column_info.get('is_foreign_key', False),
                            'query_name': pattern['name']
                        }
                        
                        # Check if we already have this candidate
                        if not any(c['table'] == table_name and c['column'] == column 
                                  for c in index_candidates):
                            index_candidates.append(candidate)
    
    # Analyze foreign keys that may not have been caught by query patterns
    for table_name, table in tables.items():
        if table_name not in analyzed_tables:
            continue
        
        table_indexes = existing_indexes.get(table_name, {})
        
        for column_info in table['columns']:
            if column_info.get('is_foreign_key') and not column_info.get('is_primary_key'):
                column = column_info['name']
                is_indexed = any(column in idx['columns'] for idx in table_indexes.values())
                
                if not is_indexed:
                    candidate = {
                        'table': table_name,
                        'column': column,
                        'data_type': column_info.get('data_type'),
                        'reason': "Foreign key relationship",
                        'is_foreign_key': True,
                        'query_name': None
                    }
                    
                    if not any(c['table'] == table_name and c['column'] == column 
                              for c in index_candidates):
                        index_candidates.append(candidate)
    
    return index_candidates

def generate_index_creation_sql(index_candidates):
    """Generate SQL statements to create the suggested indexes."""
    sql_statements = []
    
    for candidate in index_candidates:
        table = candidate['table']
        column = candidate['column']
        
        # Generate a suitable index name
        index_name = f"idx_{table}_{column}"
        
        # Different index types for different column types
        data_type = candidate.get('data_type', '').lower()
        
        # For text columns, use a more specific index type
        if 'text' in data_type or 'char' in data_type:
            # For case-insensitive searches
            sql = f"CREATE INDEX {index_name} ON {table} USING btree (lower({column}));"
            sql_statements.append({
                'sql': sql,
                'table': table,
                'column': column,
                'name': index_name,
                'reason': candidate['reason'],
                'type': 'btree-lower'
            })
            
            # If this might be used for pattern matching, add a trigram index
            sql = f"CREATE INDEX {index_name}_trgm ON {table} USING gin ({column} gin_trgm_ops);"
            sql_statements.append({
                'sql': sql,
                'table': table,
                'column': column,
                'name': f"{index_name}_trgm",
                'reason': candidate['reason'] + " (for pattern matching)",
                'type': 'gin-trigram'
            })
        else:
            # Regular B-tree index for most columns
            sql = f"CREATE INDEX {index_name} ON {table} ({column});"
            sql_statements.append({
                'sql': sql,
                'table': table,
                'column': column,
                'name': index_name,
                'reason': candidate['reason'],
                'type': 'btree'
            })
    
    return sql_statements

def save_analysis_report(tables, existing_indexes, query_patterns, index_candidates, sql_statements):
    """Save the analysis results to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"index_analysis_report_{timestamp}.json"
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'database_info': {
            'tables_analyzed': len(tables) if tables else 0,
            'existing_indexes': sum(len(indexes) for indexes in existing_indexes.values()) if existing_indexes else 0
        },
        'query_patterns': query_patterns,
        'index_candidates': index_candidates,
        'sql_statements': sql_statements,
        'summary': {
            'total_candidates': len(index_candidates),
            'total_sql_statements': len(sql_statements)
        }
    }
    
    try:
        with open(filename, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {filename}")
        return filename
    except Exception as e:
        logger.error(f"Error saving report: {e}")
        return None

def create_index_migration_file(sql_statements, filename=None):
    """Create a migration file with the index creation statements."""
    if not filename:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"migrations/033_add_performance_indexes.sql"
    
    sql_content = """-- Migration: 033_add_performance_indexes.sql
-- Purpose: Add performance-optimizing indexes to frequently queried columns
-- This migration adds indexes to improve query performance

"""
    
    # Add SQL statements with comments
    for statement in sql_statements:
        sql_content += f"-- {statement['reason']}\n"
        sql_content += f"{statement['sql']}\n\n"
    
    try:
        # Ensure migrations directory exists
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        
        with open(filename, 'w') as f:
            f.write(sql_content)
        logger.info(f"Migration file created at {filename}")
        return filename
    except Exception as e:
        logger.error(f"Error creating migration file: {e}")
        return None

def main():
    """Main function to run the analysis."""
    parser = argparse.ArgumentParser(description="Analyze database query performance and suggest indexes.")
    parser.add_argument("--dry-run", action="store_true", help="Analyze only, don't create migration files")
    args = parser.parse_args()
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    # Get table and column information
    logger.info("Getting table and column information...")
    tables = get_table_column_info()
    
    if not tables:
        logger.error("Failed to get table information. Exiting.")
        sys.exit(1)
    
    # Get existing indexes
    logger.info("Getting existing index information...")
    existing_indexes = get_existing_indexes()
    
    if not existing_indexes:
        logger.warning("Failed to get existing index information.")
        existing_indexes = {}
    
    # Analyze query patterns
    logger.info("Analyzing query patterns...")
    query_patterns = analyze_query_patterns()
    
    if not query_patterns:
        logger.error("Failed to analyze query patterns. Exiting.")
        sys.exit(1)
    
    # Identify missing indexes
    logger.info("Identifying missing indexes...")
    index_candidates = identify_missing_indexes(tables, existing_indexes, query_patterns)
    
    if not index_candidates:
        logger.info("No missing indexes identified.")
        sys.exit(0)
    
    # Generate SQL statements
    logger.info("Generating SQL statements...")
    sql_statements = generate_index_creation_sql(index_candidates)
    
    # Save analysis report
    report_file = save_analysis_report(tables, existing_indexes, query_patterns, 
                                       index_candidates, sql_statements)
    
    if not report_file:
        logger.error("Failed to save analysis report.")
        sys.exit(1)
    
    # Create migration file if not dry run
    if not args.dry_run:
        migration_file = create_index_migration_file(sql_statements)
        
        if not migration_file:
            logger.error("Failed to create migration file.")
            sys.exit(1)
    
    # Print summary
    print("\nIndex Analysis Summary:")
    print("======================")
    print(f"Tables analyzed: {len(tables)}")
    print(f"Existing indexes: {sum(len(indexes) for indexes in existing_indexes.values())}")
    print(f"Query patterns analyzed: {len(query_patterns)}")
    print(f"Index candidates identified: {len(index_candidates)}")
    print(f"SQL statements generated: {len(sql_statements)}")
    print(f"\nAnalysis report saved to: {report_file}")
    
    if not args.dry_run:
        print(f"Migration file created at: {migration_file}")
    else:
        print("\nDry run mode - no migration file was created.")
    
    if index_candidates:
        print("\nTop Index Recommendations:")
        for i, candidate in enumerate(index_candidates[:5], 1):
            print(f"{i}. Table: {candidate['table']}, Column: {candidate['column']}")
            print(f"   Reason: {candidate['reason']}")
    
    sys.exit(0)

if __name__ == "__main__":
    import argparse
    main()