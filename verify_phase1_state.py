#!/usr/bin/env python3
"""
Verify Phase 1 State

This script checks the current state of the database after Phase 1 implementation attempts,
identifying what's working and what needs to be fixed.
"""

import os
import sys
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
import json
from datetime import datetime
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"phase1_verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# Load environment variables
load_dotenv()

def connect_to_db():
    """Connect to the database."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD')
    }
    
    try:
        conn = psycopg2.connect(**db_params)
        logger.info("Connected to database")
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        return None

def check_molecule_counts(conn):
    """Check basic molecule counts in the database."""
    with conn.cursor() as cursor:
        cursor.execute("SELECT COUNT(*) FROM molecules")
        total_molecules = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM consolidated_molecules")
        consolidated_count = cursor.fetchone()[0]
        
        # Check for duplicates based on InChIKey
        cursor.execute("""
            SELECT COUNT(*) FROM (
                SELECT inchikey, COUNT(*) 
                FROM molecules 
                WHERE inchikey IS NOT NULL 
                GROUP BY inchikey 
                HAVING COUNT(*) > 1
            ) AS duplicates
        """)
        duplicate_inchikeys = cursor.fetchone()[0]
        
        # Count molecules with the same InChIKey
        cursor.execute("""
            SELECT SUM(molecule_count) - COUNT(*) AS duplicate_molecules FROM (
                SELECT inchikey, COUNT(*) AS molecule_count
                FROM molecules
                WHERE inchikey IS NOT NULL
                GROUP BY inchikey
                HAVING COUNT(*) > 1
            ) AS duplicate_groups
        """)
        duplicate_molecules = cursor.fetchone()[0] or 0
        
        return {
            "total_molecules": total_molecules,
            "consolidated_molecules": consolidated_count,
            "duplicate_inchikeys": duplicate_inchikeys,
            "duplicate_molecules": duplicate_molecules,
            "consolidation_ratio": consolidated_count / total_molecules if total_molecules > 0 else 0
        }

def check_property_standardization(conn):
    """Check property standardization status."""
    with conn.cursor() as cursor:
        # Get all property types
        cursor.execute("SELECT name FROM property_types ORDER BY name")
        property_types = [row[0] for row in cursor.fetchall()]
        
        # Look for potential duplicates (case-insensitive)
        cursor.execute("""
            SELECT LOWER(name), COUNT(*) 
            FROM property_types 
            GROUP BY LOWER(name) 
            HAVING COUNT(*) > 1
            ORDER BY COUNT(*) DESC
        """)
        duplicate_properties = {row[0]: row[1] for row in cursor.fetchall()}
        
        # Check for standardized vs non-standardized names
        standard_patterns = [
            ("LogP", ["logp", "log_p"]),
            ("Molecular Weight", ["molecular_weight", "mol_weight", "mw"]),
            ("TPSA", ["tpsa", "topological_polar_surface_area"]),
            ("Hydrogen Bond Donor Count", ["h_donors", "hbd", "donor_count"]),
            ("Hydrogen Bond Acceptor Count", ["h_acceptors", "hba", "acceptor_count"])
        ]
        
        standardization_issues = {}
        for standard, variants in standard_patterns:
            found_variants = []
            for variant in variants:
                cursor.execute("SELECT COUNT(*) FROM property_types WHERE name ILIKE %s", (f"%{variant}%",))
                if cursor.fetchone()[0] > 0:
                    found_variants.append(variant)
            
            if len(found_variants) > 0:
                standardization_issues[standard] = found_variants
        
        return {
            "total_property_types": len(property_types),
            "duplicate_properties": duplicate_properties,
            "standardization_issues": standardization_issues,
            "property_list": property_types[:30]  # Show first 30 properties
        }

def check_rdkit_properties(conn):
    """Check RDKit property calculation status."""
    with conn.cursor() as cursor:
        # Get calculation methods
        cursor.execute("SELECT id, name, version FROM calculation_methods WHERE name LIKE '%RDKit%'")
        rdkit_methods = cursor.fetchall()
        
        # Check key properties
        key_properties = ["molecular_weight", "logp", "tpsa", "h_donors", "h_acceptors"]
        property_coverage = {}
        
        for prop in key_properties:
            # Find property types matching this property
            cursor.execute("""
                SELECT id FROM property_types 
                WHERE name ILIKE %s OR name ILIKE %s
            """, (prop, f"%{prop}%"))
            property_ids = [row[0] for row in cursor.fetchall()]
            
            if property_ids:
                property_id_list = ','.join([f"'{pid}'" for pid in property_ids])
                cursor.execute(f"""
                    SELECT COUNT(DISTINCT molecule_id) 
                    FROM molecular_properties
                    WHERE property_type_id IN ({property_id_list})
                """)
                molecules_with_property = cursor.fetchone()[0]
                
                cursor.execute("SELECT COUNT(*) FROM molecules")
                total_molecules = cursor.fetchone()[0]
                
                property_coverage[prop] = {
                    "molecules_with_property": molecules_with_property,
                    "total_molecules": total_molecules,
                    "coverage_percentage": (molecules_with_property / total_molecules) * 100 if total_molecules > 0 else 0
                }
        
        return {
            "rdkit_methods": [{
                "id": method[0], 
                "name": method[1],
                "version": method[2]
            } for method in rdkit_methods],
            "property_coverage": property_coverage
        }

def check_cryoprotection_scores(conn):
    """Check cryoprotection scores status."""
    with conn.cursor() as cursor:
        # Find cryoprotection score property type
        cursor.execute("""
            SELECT id FROM property_types 
            WHERE name ILIKE '%cryoprotect%' OR name ILIKE '%cryo%score%'
        """)
        property_ids = cursor.fetchall()
        
        if not property_ids:
            return {
                "score_property_found": False,
                "molecules_with_scores": 0,
                "total_molecules": 0,
                "coverage_percentage": 0,
                "top_molecules": []
            }
        
        score_property_id = property_ids[0][0]
        
        # Count molecules with scores
        cursor.execute("""
            SELECT COUNT(DISTINCT molecule_id) 
            FROM molecular_properties
            WHERE property_type_id = %s
        """, (score_property_id,))
        molecules_with_scores = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM molecules")
        total_molecules = cursor.fetchone()[0]
        
        # Get top molecules by score
        cursor.execute("""
            SELECT m.name, mp.numeric_value
            FROM molecular_properties mp
            JOIN molecules m ON mp.molecule_id = m.id
            WHERE mp.property_type_id = %s
            ORDER BY mp.numeric_value DESC
            LIMIT 10
        """, (score_property_id,))
        top_molecules = [(name, float(score) if score is not None else 0.0) for name, score in cursor.fetchall()]
        
        return {
            "score_property_found": True,
            "score_property_id": score_property_id,
            "molecules_with_scores": molecules_with_scores,
            "total_molecules": total_molecules,
            "coverage_percentage": (molecules_with_scores / total_molecules) * 100 if total_molecules > 0 else 0,
            "top_molecules": top_molecules
        }

def main():
    """Main function."""
    conn = connect_to_db()
    if not conn:
        logger.error("Failed to connect to database")
        return 1
    
    try:
        results = {
            "timestamp": datetime.now().isoformat(),
            "molecule_counts": check_molecule_counts(conn),
            "property_standardization": check_property_standardization(conn),
            "rdkit_properties": check_rdkit_properties(conn),
            "cryoprotection_scores": check_cryoprotection_scores(conn)
        }
        
        # Output results
        logger.info("=== PHASE 1 VERIFICATION RESULTS ===")
        
        # Molecule counts
        logger.info("\n== MOLECULE COUNTS ==")
        logger.info(f"Total molecules: {results['molecule_counts']['total_molecules']}")
        logger.info(f"Consolidated molecules: {results['molecule_counts']['consolidated_molecules']}")
        logger.info(f"Duplicate InChIKeys: {results['molecule_counts']['duplicate_inchikeys']}")
        logger.info(f"Duplicate molecules: {results['molecule_counts']['duplicate_molecules']}")
        logger.info(f"Consolidation ratio: {results['molecule_counts']['consolidation_ratio']:.2%}")
        
        # Property standardization
        logger.info("\n== PROPERTY STANDARDIZATION ==")
        logger.info(f"Total property types: {results['property_standardization']['total_property_types']}")
        if results['property_standardization']['duplicate_properties']:
            logger.info("Duplicate properties found:")
            for name, count in results['property_standardization']['duplicate_properties'].items():
                logger.info(f"  {name}: {count} variants")
        else:
            logger.info("No duplicate properties found")
        
        if results['property_standardization']['standardization_issues']:
            logger.info("Standardization issues found:")
            for standard, variants in results['property_standardization']['standardization_issues'].items():
                logger.info(f"  {standard}: {', '.join(variants)}")
        else:
            logger.info("No standardization issues found")
        
        # RDKit properties
        logger.info("\n== RDKIT PROPERTIES ==")
        if results['rdkit_properties']['rdkit_methods']:
            for method in results['rdkit_properties']['rdkit_methods']:
                logger.info(f"RDKit method: {method['name']} (Version: {method['version']})")
        else:
            logger.info("No RDKit calculation methods found")
        
        logger.info("Property coverage:")
        for prop, coverage in results['rdkit_properties']['property_coverage'].items():
            logger.info(f"  {prop}: {coverage['coverage_percentage']:.2f}% ({coverage['molecules_with_property']}/{coverage['total_molecules']})")
        
        # Cryoprotection scores
        logger.info("\n== CRYOPROTECTION SCORES ==")
        if results['cryoprotection_scores']['score_property_found']:
            logger.info(f"Cryoprotection score property found")
            logger.info(f"Coverage: {results['cryoprotection_scores']['coverage_percentage']:.2f}% ({results['cryoprotection_scores']['molecules_with_scores']}/{results['cryoprotection_scores']['total_molecules']})")
            
            if results['cryoprotection_scores']['top_molecules']:
                logger.info("Top cryoprotectants:")
                for i, (name, score) in enumerate(results['cryoprotection_scores']['top_molecules'], 1):
                    logger.info(f"  {i}. {name}: {score:.2f}" if score is not None else f"  {i}. {name}: N/A")
        else:
            logger.info("Cryoprotection score property not found")
        
        # Convert any Decimal values to float for JSON serialization
        def convert_decimal(obj):
            import decimal
            if isinstance(obj, decimal.Decimal):
                return float(obj)
            raise TypeError
                
        # Save detailed results to file
        output_file = f"phase1_verification_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, default=convert_decimal)
        
        logger.info(f"\nDetailed results saved to {output_file}")
        return 0
        
    except Exception as e:
        logger.error(f"Error during verification: {e}")
        return 1
    finally:
        conn.close()

if __name__ == "__main__":
    sys.exit(main())