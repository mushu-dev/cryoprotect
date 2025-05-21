#!/usr/bin/env python3
"""
Benchmark database performance after adding indexes.

This script runs a series of common queries and measures their performance.
"""

import time
import json
from pathlib import Path
from datetime import datetime
from psycopg2.extras import RealDictCursor

# Import database modules
from database import db

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
    
    # Initialize database module
    print("Initializing database module...")
    db.init_connection_pool(config=config)
    return True

def run_benchmark_queries():
    """Run benchmark queries and measure performance."""
    print(f"Starting benchmark at {datetime.now().isoformat()}")
    results = []
    
    # Get a database connection
    conn = db.get_connection()
    if not conn:
        print("Failed to get database connection")
        return None
    
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Define benchmark queries
            benchmark_queries = [
                {
                    'name': 'Get all molecule count',
                    'query': "SELECT COUNT(*) AS count FROM molecules",
                    'params': None
                },
                {
                    'name': 'Get molecule by name pattern',
                    'query': "SELECT * FROM molecules WHERE name LIKE %s LIMIT 10",
                    'params': ('%glyc%',)
                },
                {
                    'name': 'Get molecule by formula',
                    'query': "SELECT * FROM molecules WHERE formula = %s LIMIT 10",
                    'params': ('C6H12O6',)
                },
                {
                    'name': 'Join molecules with properties',
                    'query': """
                        SELECT 
                            m.id, 
                            m.name, 
                            pt.name as property_name, 
                            mp.numeric_value, 
                            mp.text_value
                        FROM 
                            molecules m
                        JOIN 
                            molecular_properties mp ON m.id = mp.molecule_id
                        JOIN 
                            property_types pt ON mp.property_type_id = pt.id
                        LIMIT 20
                    """,
                    'params': None
                },
                {
                    'name': 'Order by molecular weight',
                    'query': "SELECT * FROM molecules ORDER BY molecular_weight DESC LIMIT 10",
                    'params': None
                },
                {
                    'name': 'Get properties by source',
                    'query': "SELECT * FROM molecular_properties WHERE source = %s LIMIT 10",
                    'params': ('PubChem',)
                },
                {
                    'name': 'Get property calculation queue by status',
                    'query': "SELECT * FROM property_calculation_queue WHERE status = %s LIMIT 10",
                    'params': ('pending',)
                },
                {
                    'name': 'Scientific data audit query',
                    'query': "SELECT * FROM scientific_data_audit WHERE operation = %s LIMIT 10",
                    'params': ('UPDATE',)
                }
            ]
            
            # Run each query 5 times and measure performance
            for query_info in benchmark_queries:
                query_name = query_info['name']
                query_sql = query_info['query']
                query_params = query_info['params']
                
                print(f"Running query: {query_name}")
                times = []
                
                for i in range(5):
                    # Clear the cache first with a dummy query
                    cursor.execute("SELECT 1")
                    
                    # Run the actual query with timing
                    start_time = time.time()
                    cursor.execute(query_sql, query_params)
                    cursor.fetchall()
                    end_time = time.time()
                    
                    execution_time = (end_time - start_time) * 1000  # Convert to ms
                    times.append(execution_time)
                    print(f"  Run {i+1}: {execution_time:.2f} ms")
                
                # Calculate average time
                avg_time = sum(times) / len(times)
                print(f"  Average time: {avg_time:.2f} ms")
                
                # Store result
                results.append({
                    'query_name': query_name,
                    'execution_times': times,
                    'average_time': avg_time
                })
    except Exception as e:
        print(f"Error running benchmark: {e}")
    finally:
        db.release_connection(conn)
    
    print(f"Benchmark completed at {datetime.now().isoformat()}")
    return results

def main():
    """Main entry point."""
    print("=" * 60)
    print(" DATABASE PERFORMANCE BENCHMARK ")
    print("=" * 60)
    
    # Initialize database
    if not init_database():
        return 1
    
    # Run benchmark queries
    results = run_benchmark_queries()
    if not results:
        print("No benchmark results")
        return 1
    
    # Save benchmark results
    results_file = f"db_benchmark_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(results_file, 'w') as f:
        json.dump({
            'timestamp': datetime.now().isoformat(),
            'results': results
        }, f, indent=2)
    
    print(f"Benchmark results saved to {results_file}")
    
    # Display summary
    print("\nBenchmark Summary:")
    for result in results:
        print(f"  {result['query_name']}: {result['average_time']:.2f} ms")
    
    return 0

if __name__ == "__main__":
    main()