#!/usr/bin/env python3
"""
Generate database load for index analysis.

This script generates database queries to help with index analysis.
"""

import time
import random
from datetime import datetime
from database import db

def run_benchmark_queries():
    """Run a series of benchmark queries to simulate application load."""
    print(f"Starting benchmark queries at {datetime.now().isoformat()}")
    
    # Connect to the database
    conn = db.get_connection()
    if not conn:
        print("Failed to get database connection")
        return
    
    try:
        cursor = conn.cursor()
        
        # Run queries on molecules table
        print("Running molecule queries...")
        for _ in range(10):
            # Random ID lookup
            cursor.execute("SELECT * FROM molecules ORDER BY RANDOM() LIMIT 1")
            molecule = cursor.fetchone()
            if molecule:
                molecule_id = molecule[0]  # Assuming ID is the first column
                
                # Query by ID
                cursor.execute("SELECT * FROM molecules WHERE id = %s", (molecule_id,))
                cursor.fetchone()
                
                # Query by name pattern
                cursor.execute("SELECT * FROM molecules WHERE name LIKE %s", (f"%{chr(random.randint(65, 90))}%",))
                cursor.fetchall()
                
                # Query with sorting
                cursor.execute("SELECT * FROM molecules ORDER BY name LIMIT 10")
                cursor.fetchall()
                
                # Query with multiple conditions
                cursor.execute(
                    "SELECT * FROM molecules WHERE molecular_weight > %s AND formula IS NOT NULL LIMIT 10",
                    (random.randint(100, 500),)
                )
                cursor.fetchall()
        
        # Run queries on molecular_properties table
        print("Running property queries...")
        for _ in range(10):
            # Query properties by molecule ID
            cursor.execute(
                "SELECT * FROM molecular_properties WHERE molecule_id = %s",
                (molecule_id,)
            )
            cursor.fetchall()
            
            # Query properties by type
            cursor.execute(
                "SELECT mp.*, m.name FROM molecular_properties mp JOIN molecules m ON mp.molecule_id = m.id JOIN property_types pt ON mp.property_type_id = pt.id WHERE pt.name = %s LIMIT 10",
                (random.choice(['LogP', 'Molecular Weight', 'Hydrogen Bond Donors', 'Hydrogen Bond Acceptors']),)
            )
            cursor.fetchall()
            
            # Query with numeric range
            cursor.execute(
                "SELECT mp.*, m.name FROM molecular_properties mp JOIN molecules m ON mp.molecule_id = m.id WHERE mp.numeric_value BETWEEN %s AND %s LIMIT 10",
                (random.randint(0, 100), random.randint(101, 200))
            )
            cursor.fetchall()
        
        # Run queries on property_types table
        print("Running property type queries...")
        for _ in range(5):
            cursor.execute("SELECT * FROM property_types WHERE data_type = %s",
                         (random.choice(['numeric', 'text', 'boolean']),))
            cursor.fetchall()
        
        # Run complex join queries
        print("Running complex join queries...")
        for _ in range(5):
            cursor.execute("""
                SELECT 
                    m.id, 
                    m.name, 
                    m.formula, 
                    pt.name as property_name, 
                    mp.numeric_value,
                    mp.text_value
                FROM 
                    molecules m
                JOIN 
                    molecular_properties mp ON m.id = mp.molecule_id
                JOIN 
                    property_types pt ON mp.property_type_id = pt.id
                WHERE 
                    pt.data_type = %s
                ORDER BY 
                    m.name
                LIMIT 20
            """, (random.choice(['numeric', 'text']),))
            cursor.fetchall()
            
            # Query with subquery
            cursor.execute("""
                SELECT 
                    m.id, 
                    m.name, 
                    m.formula,
                    (SELECT COUNT(*) FROM molecular_properties mp WHERE mp.molecule_id = m.id) as property_count
                FROM 
                    molecules m
                WHERE 
                    m.name LIKE %s
                ORDER BY 
                    property_count DESC
                LIMIT 10
            """, (f"%{chr(random.randint(65, 90))}%",))
            cursor.fetchall()
    
    except Exception as e:
        print(f"Error running benchmark queries: {e}")
    finally:
        db.release_connection(conn)
    
    print(f"Completed benchmark queries at {datetime.now().isoformat()}")

if __name__ == "__main__":
    # Initialize database connection
    from database import db
    
    # Load configuration
    import json
    from pathlib import Path
    
    config_path = Path('config/config.json')
    if config_path.exists():
        with open(config_path, 'r') as f:
            config_data = json.load(f)
        
        if 'database' in config_data and 'connection' in config_data['database']:
            connection_config = config_data['database']['connection']
            if 'supabase' in connection_config:
                config = connection_config['supabase']
                # Add pooling settings if available
                if 'pooling' in config_data['database']:
                    pooling = config_data['database']['pooling']
                    if 'min_connections' in pooling:
                        config['min_connections'] = pooling['min_connections']
                    if 'max_connections' in pooling:
                        config['max_connections'] = pooling['max_connections']
                
                # Initialize database
                db.init_connection_pool(config=config)
                
                # Run benchmark queries
                run_benchmark_queries()
            else:
                print("No supabase configuration found")
        else:
            print("No database connection configuration found")
    else:
        print(f"Configuration file not found at {config_path}")