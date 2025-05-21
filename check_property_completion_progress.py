#!/usr/bin/env python3
"""
Monitor the progress of the property completion background job.
Provides statistics on how many molecules have properties and what's remaining.
"""

import os
import sys
import psycopg2
from psycopg2.extras import RealDictCursor
from datetime import datetime
import logging
from dotenv import load_dotenv

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

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
    
    logger.info(f"Connecting to database")
    return psycopg2.connect(**db_params)

def check_background_job():
    """Check if the property completion background job is running."""
    import subprocess
    result = subprocess.run(
        ["ps", "aux"], 
        capture_output=True, 
        text=True
    )
    
    processes = result.stdout.split('\n')
    running_jobs = [p for p in processes if 'complete_missing_properties' in p and 'grep' not in p]
    
    if running_jobs:
        for job in running_jobs:
            parts = job.split()
            pid = parts[1] if len(parts) > 1 else "Unknown"
            start_time = parts[9] if len(parts) > 9 else "Unknown"
            logger.info(f"Property completion job running: PID={pid}, started at {start_time}")
        return True
    else:
        logger.info("No property completion background job found running")
        return False

def get_property_statistics(conn):
    """Get statistics on molecular properties in the database."""
    queries = {
        "total_molecules": """
            SELECT COUNT(*) as count FROM molecules
        """,
        "molecules_with_smiles": """
            SELECT COUNT(*) as count FROM molecules WHERE smiles IS NOT NULL
        """,
        "molecules_with_formula": """
            SELECT COUNT(*) as count 
            FROM molecules m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE pt.name = 'molecular_formula'
        """,
        "molecules_with_molecular_weight": """
            SELECT COUNT(*) as count 
            FROM molecules m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE pt.name = 'molecular_weight'
        """,
        "molecules_with_logp": """
            SELECT COUNT(*) as count 
            FROM molecules m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE pt.name = 'logP'
        """,
        "molecules_with_tpsa": """
            SELECT COUNT(*) as count 
            FROM molecules m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE pt.name = 'TPSA'
        """,
        "molecules_with_all_properties": """
            SELECT COUNT(DISTINCT m.id) as count
            FROM molecules m
            WHERE m.id IN (
                SELECT molecule_id 
                FROM molecular_properties mp
                JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE pt.name IN ('molecular_formula', 'molecular_weight', 'logP', 'TPSA')
                GROUP BY molecule_id
                HAVING COUNT(DISTINCT pt.name) >= 4
            )
        """,
        "property_counts_by_source": """
            SELECT mp.source, COUNT(*) as count
            FROM molecular_properties mp
            GROUP BY mp.source
            ORDER BY count DESC
        """,
        "recently_updated_properties": """
            SELECT COUNT(*) as count
            FROM molecular_properties
            WHERE updated_at > NOW() - INTERVAL '1 hour'
        """
    }
    
    results = {}
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        for key, query in queries.items():
            try:
                cursor.execute(query)
                result = cursor.fetchone()
                results[key] = result['count'] if result else 0
            except Exception as e:
                logger.error(f"Error executing query {key}: {e}")
                results[key] = "Error"
    
    return results

def calculate_progress(stats):
    """Calculate completion percentages based on statistics."""
    if stats["total_molecules"] > 0:
        progress = {
            "molecules_with_smiles": (stats["molecules_with_smiles"] / stats["total_molecules"]) * 100,
            "molecules_with_formula": (stats["molecules_with_formula"] / stats["total_molecules"]) * 100,
            "molecules_with_molecular_weight": (stats["molecules_with_molecular_weight"] / stats["total_molecules"]) * 100,
            "molecules_with_logp": (stats["molecules_with_logp"] / stats["total_molecules"]) * 100,
            "molecules_with_tpsa": (stats["molecules_with_tpsa"] / stats["total_molecules"]) * 100,
            "molecules_with_all_properties": (stats["molecules_with_all_properties"] / stats["total_molecules"]) * 100
        }
        return progress
    return {}

def estimate_completion_time(stats, progress):
    """Estimate time to completion based on recent update rate."""
    if stats["recently_updated_properties"] > 0 and "molecules_with_all_properties" in progress:
        # Properties updated in the last hour
        hourly_rate = stats["recently_updated_properties"]
        
        # Remaining properties to add (rough estimate)
        remaining_molecules = stats["total_molecules"] - stats["molecules_with_all_properties"]
        remaining_properties = remaining_molecules * 4  # Assuming 4 main properties per molecule
        
        # Estimated hours to completion
        hours_remaining = remaining_properties / hourly_rate if hourly_rate > 0 else float('inf')
        
        if hours_remaining == float('inf'):
            return "Cannot calculate - no recent updates"
        elif hours_remaining > 48:
            return f"Approximately {hours_remaining/24:.1f} days"
        else:
            return f"Approximately {hours_remaining:.1f} hours"
    
    return "Cannot calculate - insufficient data"

def main():
    """Main function to check property completion progress."""
    # Check if background job is running
    job_running = check_background_job()
    
    # Connect to the database
    conn = connect_to_db()
    
    try:
        # Get property statistics
        stats = get_property_statistics(conn)
        
        # Calculate progress percentages
        progress = calculate_progress(stats)
        
        # Estimate completion time
        estimated_time = estimate_completion_time(stats, progress)
        
        # Print report
        print("\nProperty Completion Progress Report")
        print("=================================")
        print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Background Job Running: {'Yes' if job_running else 'No'}")
        
        print("\nDatabase Statistics:")
        print(f"Total Molecules: {stats['total_molecules']}")
        print(f"Molecules with SMILES: {stats['molecules_with_smiles']} ({progress.get('molecules_with_smiles', 0):.2f}%)")
        
        print("\nProperty Coverage:")
        print(f"Molecules with Formula: {stats['molecules_with_formula']} ({progress.get('molecules_with_formula', 0):.2f}%)")
        print(f"Molecules with Molecular Weight: {stats['molecules_with_molecular_weight']} ({progress.get('molecules_with_molecular_weight', 0):.2f}%)")
        print(f"Molecules with LogP: {stats['molecules_with_logp']} ({progress.get('molecules_with_logp', 0):.2f}%)")
        print(f"Molecules with TPSA: {stats['molecules_with_tpsa']} ({progress.get('molecules_with_tpsa', 0):.2f}%)")
        print(f"Molecules with All Properties: {stats['molecules_with_all_properties']} ({progress.get('molecules_with_all_properties', 0):.2f}%)")
        
        print("\nRecent Activity:")
        print(f"Properties Updated in Last Hour: {stats['recently_updated_properties']}")
        print(f"Estimated Time to Completion: {estimated_time}")
        
        print("\nProperty Sources:")
        if "property_counts_by_source" in stats and isinstance(stats["property_counts_by_source"], list):
            for source_stat in stats["property_counts_by_source"]:
                print(f"  {source_stat['source']}: {source_stat['count']} properties")
        
        # Return code based on completion status
        if progress.get('molecules_with_all_properties', 0) >= 99.5:
            print("\nStatus: COMPLETE - Property completion is essentially finished!")
            return 0
        elif job_running:
            print("\nStatus: RUNNING - Property completion job is still running")
            return 0
        else:
            print("\nStatus: STOPPED - Property completion job is not running")
            return 1
    
    except Exception as e:
        logger.error(f"Error checking property completion: {e}", exc_info=True)
        return 1
    
    finally:
        conn.close()

if __name__ == "__main__":
    sys.exit(main())