#!/usr/bin/env python3
"""
Analyze query patterns and recommend database indexes.

This script connects to the database, monitors queries, and suggests indexes
that would improve performance for common query patterns. It applies the recommended
indexes if specified in the command-line arguments.
"""

import sys
import argparse
import json
import time
import re
from pathlib import Path
from datetime import datetime
from psycopg2.extras import RealDictCursor

# Import our database modules
from database import db, db_service_role

# Configure argument parsing
parser = argparse.ArgumentParser(description="Analyze query patterns and recommend indexes")
parser.add_argument('--monitor-time', type=int, default=60, 
                   help='Time in seconds to monitor queries (default: 60)')
parser.add_argument('--min-occurrences', type=int, default=5,
                   help='Minimum query occurrences to consider for indexing (default: 5)')
parser.add_argument('--apply', action='store_true',
                   help='Apply the recommended indexes')
parser.add_argument('--report-file', type=str, default='index_recommendations.json',
                   help='File to store the index recommendations (default: index_recommendations.json)')
args = parser.parse_args()

def load_config():
    """Load database configuration."""
    config_path = Path('config/config.json')
    if not config_path.exists():
        print(f"Error: Configuration file not found at {config_path}")
        return None
    
    with open(config_path, 'r') as f:
        config_data = json.load(f)
    
    if 'database' in config_data and 'connection' in config_data['database']:
        connection_config = config_data['database']['connection']
        if 'supabase' in connection_config:
            config = connection_config['supabase'].copy()
            # Add pooling settings if available
            if 'pooling' in config_data['database']:
                pooling = config_data['database']['pooling']
                if 'min_connections' in pooling:
                    config['min_connections'] = pooling['min_connections']
                if 'max_connections' in pooling:
                    config['max_connections'] = pooling['max_connections']
            return config
    
    print("Error: No valid database configuration found")
    return None

def init_database():
    """Initialize database connection."""
    config = load_config()
    if not config:
        print("Failed to load configuration")
        return False
    
    # Initialize database modules
    print("Initializing database modules...")
    db.init_connection_pool(config=config)
    
    # Create service role config
    service_role_config = config.copy()
    service_role_config['options'] = "-c role=service_role"
    db_service_role.init_connection_pool(config=service_role_config)
    
    return True

def get_existing_indexes():
    """Get existing indexes in the database."""
    print("Getting existing indexes...")
    indexes = {}
    
    try:
        conn = db_service_role.get_connection()
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Query for all indexes in the database
            cursor.execute("""
                SELECT
                    t.relname AS table_name,
                    i.relname AS index_name,
                    a.attname AS column_name,
                    ix.indisunique AS is_unique,
                    ix.indisprimary AS is_primary
                FROM
                    pg_class t,
                    pg_class i,
                    pg_index ix,
                    pg_attribute a
                WHERE
                    t.oid = ix.indrelid
                    AND i.oid = ix.indexrelid
                    AND a.attrelid = t.oid
                    AND a.attnum = ANY(ix.indkey)
                    AND t.relkind = 'r'
                    AND t.relnamespace = (SELECT oid FROM pg_namespace WHERE nspname = 'public')
                ORDER BY
                    t.relname,
                    i.relname,
                    a.attnum
            """)
            
            for row in cursor.fetchall():
                table_name = row['table_name']
                index_name = row['index_name']
                column_name = row['column_name']
                is_unique = row['is_unique']
                is_primary = row['is_primary']
                
                if table_name not in indexes:
                    indexes[table_name] = []
                
                # Check if this index is already in our list
                found = False
                for idx in indexes[table_name]:
                    if idx['name'] == index_name:
                        idx['columns'].append(column_name)
                        found = True
                        break
                
                if not found:
                    indexes[table_name].append({
                        'name': index_name,
                        'columns': [column_name],
                        'unique': is_unique,
                        'primary': is_primary
                    })
        
        db_service_role.release_connection(conn)
    except Exception as e:
        print(f"Error getting existing indexes: {e}")
        return {}
    
    return indexes

def analyze_schema():
    """Analyze the database schema and recommend indexes for foreign keys and common query patterns."""
    print("Analyzing database schema...")
    tables_info = {}
    foreign_keys = []
    
    try:
        conn = db_service_role.get_connection()
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Get tables and column information
            cursor.execute("""
                SELECT
                    t.table_name,
                    c.column_name,
                    c.data_type,
                    c.is_nullable,
                    tc.constraint_type
                FROM
                    information_schema.tables t
                JOIN
                    information_schema.columns c ON t.table_name = c.table_name
                LEFT JOIN
                    information_schema.constraint_column_usage ccu ON c.table_name = ccu.table_name AND c.column_name = ccu.column_name
                LEFT JOIN
                    information_schema.table_constraints tc ON ccu.constraint_name = tc.constraint_name
                WHERE
                    t.table_schema = 'public'
                    AND t.table_type = 'BASE TABLE'
                ORDER BY
                    t.table_name,
                    c.ordinal_position
            """)
            
            for row in cursor.fetchall():
                table_name = row['table_name']
                column_name = row['column_name']
                data_type = row['data_type']
                is_nullable = row['is_nullable']
                constraint_type = row['constraint_type']
                
                if table_name not in tables_info:
                    tables_info[table_name] = {
                        'columns': [],
                        'primary_key': None,
                        'row_count': 0
                    }
                
                tables_info[table_name]['columns'].append({
                    'name': column_name,
                    'data_type': data_type,
                    'is_nullable': is_nullable,
                    'constraint_type': constraint_type
                })
                
                if constraint_type == 'PRIMARY KEY':
                    tables_info[table_name]['primary_key'] = column_name
            
            # Get foreign key relationships
            cursor.execute("""
                SELECT
                    kcu.table_name,
                    kcu.column_name,
                    ccu.table_name AS foreign_table_name,
                    ccu.column_name AS foreign_column_name
                FROM
                    information_schema.table_constraints AS tc
                JOIN
                    information_schema.key_column_usage AS kcu ON tc.constraint_name = kcu.constraint_name
                JOIN
                    information_schema.constraint_column_usage AS ccu ON ccu.constraint_name = tc.constraint_name
                WHERE
                    tc.constraint_type = 'FOREIGN KEY'
                    AND tc.table_schema = 'public'
            """)
            
            for row in cursor.fetchall():
                foreign_keys.append({
                    'table': row['table_name'],
                    'column': row['column_name'],
                    'foreign_table': row['foreign_table_name'],
                    'foreign_column': row['foreign_column_name']
                })
            
            # Get table row counts
            for table_name in tables_info:
                try:
                    cursor.execute(f"SELECT COUNT(*) AS count FROM {table_name}")
                    result = cursor.fetchone()
                    if result:
                        tables_info[table_name]['row_count'] = result['count']
                except Exception as e:
                    print(f"  Error getting row count for {table_name}: {e}")
        
        db_service_role.release_connection(conn)
    except Exception as e:
        print(f"Error analyzing schema: {e}")
        return {}, []
    
    return tables_info, foreign_keys

def generate_common_queries(tables_info):
    """Generate common query patterns based on schema analysis."""
    print("Generating common query patterns...")
    queries = {}
    
    # Typical query patterns for each table
    for table_name, info in tables_info.items():
        if info['row_count'] < 100:
            # Skip tables with few rows
            continue
        
        primary_key = info['primary_key']
        if primary_key:
            # Lookup by primary key
            query = f"SELECT * FROM {table_name} WHERE {primary_key} = $1"
            queries[query] = {
                'calls': 100,  # Estimated frequency
                'tables': [(table_name, table_name)],
                'columns': [primary_key],
                'pattern': 'primary_key_lookup'
            }
        
        # Find text columns for LIKE queries
        text_columns = [col['name'] for col in info['columns'] 
                      if col['data_type'] in ('character varying', 'text', 'char') 
                      and col['name'] not in ('id', 'uuid', primary_key)]
        
        for column in text_columns:
            # Text search queries
            query = f"SELECT * FROM {table_name} WHERE {column} LIKE $1"
            queries[query] = {
                'calls': 20,  # Estimated frequency
                'tables': [(table_name, table_name)],
                'columns': [column],
                'pattern': 'text_search'
            }
        
        # Find numeric columns for range queries
        numeric_columns = [col['name'] for col in info['columns'] 
                         if col['data_type'] in ('integer', 'numeric', 'decimal', 'real', 'double precision') 
                         and col['name'] not in ('id', 'uuid', primary_key)]
        
        for column in numeric_columns:
            # Range queries
            query = f"SELECT * FROM {table_name} WHERE {column} BETWEEN $1 AND $2"
            queries[query] = {
                'calls': 10,  # Estimated frequency
                'tables': [(table_name, table_name)],
                'columns': [column],
                'pattern': 'numeric_range'
            }
        
        # ORDER BY queries for sortable columns
        sortable_columns = text_columns + numeric_columns
        if sortable_columns:
            for column in sortable_columns:
                query = f"SELECT * FROM {table_name} ORDER BY {column} LIMIT 20"
                queries[query] = {
                    'calls': 5,  # Estimated frequency
                    'tables': [(table_name, table_name)],
                    'columns': [column],
                    'pattern': 'sorting'
                }
    
    return queries

def analyze_foreign_keys(foreign_keys, tables_info):
    """Generate query patterns based on foreign key relationships."""
    print("Analyzing foreign key relationships...")
    queries = {}
    
    for fk in foreign_keys:
        table = fk['table']
        column = fk['column']
        foreign_table = fk['foreign_table']
        foreign_column = fk['foreign_column']
        
        # Skip tables with few rows
        if table in tables_info and tables_info[table]['row_count'] < 100:
            continue
        
        # JOIN queries
        query = f"SELECT * FROM {table} JOIN {foreign_table} ON {table}.{column} = {foreign_table}.{foreign_column}"
        queries[query] = {
            'calls': 15,  # Estimated frequency
            'tables': [(table, table), (foreign_table, foreign_table)],
            'columns': [column],
            'pattern': 'foreign_key_join'
        }
        
        # Lookup by foreign key
        query = f"SELECT * FROM {table} WHERE {column} = $1"
        queries[query] = {
            'calls': 30,  # Estimated frequency
            'tables': [(table, table)],
            'columns': [column],
            'pattern': 'foreign_key_lookup'
        }
    
    return queries

def parse_query(query):
    """Parse a SQL query to extract tables and columns."""
    tables = []
    columns = []
    
    # Extract tables from FROM clause
    from_match = re.search(r'FROM\s+([^WHERE;]+)', query, re.IGNORECASE)
    if from_match:
        from_clause = from_match.group(1).strip()
        table_list = re.split(r',\s*', from_clause)
        for table_item in table_list:
            # Extract table name and alias
            table_match = re.search(r'([^\s(]+)(?:\s+(?:AS\s+)?([^\s]+))?', table_item.strip(), re.IGNORECASE)
            if table_match:
                table_name = table_match.group(1).strip()
                alias = table_match.group(2).strip() if table_match.group(2) else table_name
                tables.append((table_name, alias))
    
    # Extract columns from WHERE clause
    where_match = re.search(r'WHERE\s+([^;]+)', query, re.IGNORECASE)
    if where_match:
        where_clause = where_match.group(1).strip()
        # Extract column comparisons
        column_matches = re.findall(r'([^\s=<>!]+)\s*(?:[=<>!]+|IS|LIKE|IN)', where_clause, re.IGNORECASE)
        for column in column_matches:
            # Remove table alias if present
            if '.' in column:
                parts = column.split('.')
                columns.append(parts[1])
            else:
                columns.append(column)
    
    # Extract columns from ORDER BY clause
    order_match = re.search(r'ORDER\s+BY\s+([^;]+)', query, re.IGNORECASE)
    if order_match:
        order_clause = order_match.group(1).strip()
        # Extract columns
        order_columns = re.split(r',\s*', order_clause)
        for column in order_columns:
            # Remove DESC or ASC if present
            column = re.sub(r'\s+(?:ASC|DESC)$', '', column.strip(), flags=re.IGNORECASE)
            # Remove table alias if present
            if '.' in column:
                parts = column.split('.')
                columns.append(parts[1])
            else:
                columns.append(column)
    
    return tables, columns

def recommend_indexes(queries, existing_indexes, min_occurrences=5):
    """Recommend indexes based on query patterns."""
    print("Analyzing query patterns and recommending indexes...")
    
    # Count occurrences of columns in WHERE and ORDER BY clauses
    column_occurrences = {}
    for query_info in queries.values():
        calls = query_info.get('calls', 0)
        tables = query_info.get('tables', [])
        columns = query_info.get('columns', [])
        pattern = query_info.get('pattern', 'unknown')
        
        # Skip queries with low occurrences
        if calls < min_occurrences:
            continue
        
        # Weight by query pattern
        weight = 1.0
        if pattern == 'primary_key_lookup':
            weight = 0.2  # Lower weight as these are already indexed
        elif pattern == 'foreign_key_lookup':
            weight = 2.0  # Higher weight as these are important for performance
        elif pattern == 'foreign_key_join':
            weight = 1.5  # Moderate weight for joins
        elif pattern == 'text_search':
            weight = 1.2  # Moderate weight for text search
        elif pattern == 'numeric_range':
            weight = 1.3  # Moderate weight for range queries
        elif pattern == 'sorting':
            weight = 1.0  # Normal weight for sorting
        
        # Count column occurrences
        for table_name, _ in tables:
            if table_name not in column_occurrences:
                column_occurrences[table_name] = {}
            
            for column in columns:
                if column not in column_occurrences[table_name]:
                    column_occurrences[table_name][column] = 0
                column_occurrences[table_name][column] += calls * weight
    
    # Generate index recommendations
    recommendations = []
    for table_name, columns in column_occurrences.items():
        # Skip tables that don't exist in the database
        if table_name not in existing_indexes:
            continue
        
        # Get existing indexes for this table
        table_indexes = existing_indexes[table_name]
        
        # Sort columns by occurrence
        sorted_columns = sorted(columns.items(), key=lambda x: x[1], reverse=True)
        
        # Generate recommendations for frequently used columns
        for column, occurrences in sorted_columns:
            # Skip columns that are already indexed
            already_indexed = False
            for index in table_indexes:
                if column in index['columns']:
                    already_indexed = True
                    break
            
            if already_indexed:
                continue
            
            # Recommend an index for this column
            index_name = f"idx_{table_name}_{column}"
            recommendations.append({
                'table': table_name,
                'column': column,
                'occurrences': occurrences,
                'index_name': index_name,
                'sql': f"CREATE INDEX {index_name} ON {table_name} ({column})"
            })
    
    return recommendations

def apply_indexes(recommendations):
    """Apply the recommended indexes to the database."""
    print("Applying recommended indexes...")
    results = []
    
    try:
        conn = db_service_role.get_connection()
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            for recommendation in recommendations:
                try:
                    print(f"Creating index {recommendation['index_name']}...")
                    start_time = time.time()
                    cursor.execute(recommendation['sql'])
                    conn.commit()
                    end_time = time.time()
                    
                    results.append({
                        'status': 'success',
                        'index_name': recommendation['index_name'],
                        'table': recommendation['table'],
                        'column': recommendation['column'],
                        'execution_time': end_time - start_time
                    })
                    
                    print(f"  ✓ Index created successfully ({end_time - start_time:.2f} seconds)")
                except Exception as e:
                    print(f"  ✗ Error creating index: {e}")
                    results.append({
                        'status': 'error',
                        'index_name': recommendation['index_name'],
                        'table': recommendation['table'],
                        'column': recommendation['column'],
                        'error': str(e)
                    })
        
        db_service_role.release_connection(conn)
    except Exception as e:
        print(f"Error applying indexes: {e}")
    
    return results

def create_report(recommendations, results=None):
    """Create a report of index recommendations and results."""
    report = {
        'timestamp': datetime.now().isoformat(),
        'recommendations': recommendations
    }
    
    if results:
        report['results'] = results
    
    return report

def main():
    """Main entry point."""
    print("=" * 60)
    print(" DATABASE INDEX ANALYZER ")
    print("=" * 60)
    
    # Initialize database connection
    if not init_database():
        return 1
    
    # Get existing indexes
    existing_indexes = get_existing_indexes()
    print(f"Found {sum(len(indexes) for indexes in existing_indexes.values())} existing indexes.")
    
    # Analyze schema
    tables_info, foreign_keys = analyze_schema()
    print(f"Analyzed {len(tables_info)} tables and {len(foreign_keys)} foreign key relationships.")
    
    # Generate common query patterns
    schema_queries = generate_common_queries(tables_info)
    print(f"Generated {len(schema_queries)} common query patterns.")
    
    # Generate query patterns from foreign keys
    fk_queries = analyze_foreign_keys(foreign_keys, tables_info)
    print(f"Generated {len(fk_queries)} query patterns from foreign keys.")
    
    # Combine all query patterns
    queries = {**schema_queries, **fk_queries}
    print(f"Total query patterns: {len(queries)}")
    
    # Recommend indexes
    recommendations = recommend_indexes(queries, existing_indexes, args.min_occurrences)
    print(f"Generated {len(recommendations)} index recommendations.")
    
    # Apply indexes if requested
    results = None
    if args.apply and recommendations:
        results = apply_indexes(recommendations)
        print(f"Applied {len([r for r in results if r['status'] == 'success'])} indexes.")
    
    # Create report
    report = create_report(recommendations, results)
    report['tables_analyzed'] = len(tables_info)
    report['foreign_keys_analyzed'] = len(foreign_keys)
    report['query_patterns_analyzed'] = len(queries)
    
    # Save report
    with open(args.report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"Report saved to {args.report_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())