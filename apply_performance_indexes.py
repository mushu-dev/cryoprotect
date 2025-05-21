#!/usr/bin/env python3
"""
Apply performance indexes to the database.

This script applies the recommended indexes from the index_recommendations.json file.
It also creates a migration SQL file for version control.
"""

import sys
import json
import time
import argparse
from pathlib import Path
from datetime import datetime
from psycopg2.extras import RealDictCursor

# Import database modules
from database import db, db_service_role

# Configure argument parsing
parser = argparse.ArgumentParser(description="Apply performance indexes to the database")
parser.add_argument('--recommendations-file', type=str, default='index_recommendations.json',
                    help='File containing index recommendations (default: index_recommendations.json)')
parser.add_argument('--migration-file', type=str, default=None,
                    help='File to write the migration SQL to (default: migrations/NNN_performance_indexes.sql)')
parser.add_argument('--dry-run', action='store_true',
                    help='Only print what would be done, do not apply indexes')
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

def load_recommendations(file_path):
    """Load index recommendations from a file."""
    print(f"Loading recommendations from {file_path}...")
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        if 'recommendations' not in data:
            print("Error: Recommendations file has no 'recommendations' field")
            return []
        
        recommendations = data['recommendations']
        print(f"Loaded {len(recommendations)} recommendations")
        return recommendations
    except FileNotFoundError:
        print(f"Error: Recommendations file not found at {file_path}")
        return []
    except json.JSONDecodeError:
        print(f"Error: Recommendations file is not valid JSON")
        return []
    except Exception as e:
        print(f"Error loading recommendations: {e}")
        return []

def generate_migration_file(recommendations, file_path=None):
    """Generate a migration SQL file with all index creation statements."""
    if not file_path:
        # Find the next migration number
        migration_dir = Path('migrations')
        if not migration_dir.exists():
            migration_dir.mkdir(parents=True)
        
        existing_migrations = [p.name for p in migration_dir.glob('*.sql') if p.name.startswith('0')]
        if existing_migrations:
            # Extract the highest migration number
            highest_num = max(int(m.split('_')[0]) for m in existing_migrations)
            next_num = highest_num + 1
        else:
            next_num = 1
        
        # Format with leading zeros
        file_path = migration_dir / f"{next_num:03d}_performance_indexes.sql"
    
    print(f"Generating migration file at {file_path}...")
    
    # Create the migration SQL
    sql = f"""-- Migration: Performance Indexes
-- Created: {datetime.now().isoformat()}
-- Description: Adds performance indexes for common query patterns

BEGIN;

-- Performance indexes for frequently queried columns
"""
    
    # Add each index creation statement
    for rec in recommendations:
        sql += f"{rec['sql']};\n"
    
    sql += """
-- End of migration
COMMIT;
"""
    
    # Write the migration file
    with open(file_path, 'w') as f:
        f.write(sql)
    
    print(f"Migration file generated at {file_path}")
    return file_path

def apply_indexes(recommendations, dry_run=False):
    """Apply the recommended indexes to the database."""
    print("Applying indexes to the database...")
    
    if dry_run:
        print("DRY RUN - Not actually applying indexes")
        for rec in recommendations:
            print(f"Would execute: {rec['sql']}")
        return []
    
    results = []
    
    try:
        conn = db_service_role.get_connection()
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            for rec in recommendations:
                try:
                    print(f"Creating index {rec['index_name']}...")
                    start_time = time.time()
                    cursor.execute(rec['sql'])
                    conn.commit()
                    end_time = time.time()
                    
                    results.append({
                        'status': 'success',
                        'index_name': rec['index_name'],
                        'table': rec['table'],
                        'column': rec['column'],
                        'execution_time': end_time - start_time
                    })
                    
                    print(f"  ✓ Index created successfully ({end_time - start_time:.2f} seconds)")
                except Exception as e:
                    conn.rollback()
                    print(f"  ✗ Error creating index: {e}")
                    results.append({
                        'status': 'error',
                        'index_name': rec['index_name'],
                        'table': rec['table'],
                        'column': rec['column'],
                        'error': str(e)
                    })
        
        db_service_role.release_connection(conn)
    except Exception as e:
        print(f"Error applying indexes: {e}")
    
    return results

def main():
    """Main entry point."""
    print("=" * 60)
    print(" DATABASE PERFORMANCE INDEX APPLIER ")
    print("=" * 60)
    
    # Initialize database connection
    if not init_database():
        return 1
    
    # Load recommendations
    recommendations = load_recommendations(args.recommendations_file)
    if not recommendations:
        print("No recommendations to apply")
        return 1
    
    # Generate migration file
    migration_file = generate_migration_file(recommendations, args.migration_file)
    
    # Apply indexes if not dry run
    if not args.dry_run:
        results = apply_indexes(recommendations)
        success_count = len([r for r in results if r['status'] == 'success'])
        error_count = len([r for r in results if r['status'] == 'error'])
        
        print(f"\nIndex application summary:")
        print(f"  Total indexes: {len(recommendations)}")
        print(f"  Successfully applied: {success_count}")
        print(f"  Failed: {error_count}")
        
        # Save results
        results_file = f"index_application_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(results_file, 'w') as f:
            json.dump({
                'timestamp': datetime.now().isoformat(),
                'migration_file': str(migration_file),
                'results': results
            }, f, indent=2)
        
        print(f"Results saved to {results_file}")
    else:
        print("\nDRY RUN - No indexes were applied")
        print(f"Migration file generated at {migration_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())