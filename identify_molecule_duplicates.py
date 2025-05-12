#!/usr/bin/env python3
"""
Identify molecules that are likely duplicates and create a report.

This script:
1. Identifies molecules that are likely duplicates (same name, formula, or structure)
2. Creates groups of duplicate molecules
3. Generates a report for manual review
"""

import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

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
    return psycopg2.connect(**db_params)

def find_duplicate_names(conn):
    """Find molecules with duplicate names."""
    duplicates = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Find names that appear multiple times
        cursor.execute("""
            SELECT lower(name) as name, COUNT(*) as count
            FROM molecules
            WHERE name IS NOT NULL
            GROUP BY lower(name)
            HAVING COUNT(*) > 1
            ORDER BY count DESC
        """)
        
        duplicate_names = cursor.fetchall()
        
        # For each duplicate name, get the molecules
        for dup in duplicate_names:
            name = dup["name"]
            cursor.execute("""
                SELECT id, name, molecular_formula, pubchem_cid, smiles, 
                       data_source, created_at, is_public
                FROM molecules
                WHERE lower(name) = %s
                ORDER BY pubchem_cid NULLS LAST, created_at
            """, (name,))
            
            molecules = cursor.fetchall()
            
            duplicates.append({
                "name": name,
                "count": dup["count"],
                "molecules": molecules
            })
    
    return duplicates

def find_duplicate_formulas(conn):
    """Find molecules with duplicate molecular formulas."""
    duplicates = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Find formulas that appear multiple times
        cursor.execute("""
            SELECT molecular_formula, COUNT(*) as count
            FROM molecules
            WHERE molecular_formula IS NOT NULL
            GROUP BY molecular_formula
            HAVING COUNT(*) > 1
            ORDER BY count DESC
        """)
        
        duplicate_formulas = cursor.fetchall()
        
        # For each duplicate formula, get the molecules
        for dup in duplicate_formulas:
            formula = dup["molecular_formula"]
            cursor.execute("""
                SELECT id, name, molecular_formula, pubchem_cid, smiles, 
                       data_source, created_at, is_public
                FROM molecules
                WHERE molecular_formula = %s
                ORDER BY pubchem_cid NULLS LAST, created_at
            """, (formula,))
            
            molecules = cursor.fetchall()
            
            duplicates.append({
                "formula": formula,
                "count": dup["count"],
                "molecules": molecules
            })
    
    return duplicates

def find_molecules_missing_cids(conn):
    """Find molecules with missing PubChem CIDs."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, molecular_formula, pubchem_cid, smiles,
                   data_source, created_at, is_public
            FROM molecules
            WHERE pubchem_cid IS NULL
            ORDER BY name
        """)
        
        return cursor.fetchall()

def analyze_test_molecules(conn):
    """Analyze test molecules."""
    test_molecules = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, molecular_formula, pubchem_cid, smiles,
                   data_source, created_at, is_public
            FROM molecules
            WHERE name ILIKE '%test%' OR name ILIKE '%dummy%' OR name ILIKE '%example%'
            ORDER BY name
        """)
        
        test_molecules = cursor.fetchall()
    
    return test_molecules

def generate_report(name_duplicates, formula_duplicates, missing_cids, test_molecules):
    """Generate a comprehensive report on molecule duplicates and issues."""
    total_name_dupes = sum(d["count"] for d in name_duplicates)
    total_formula_dupes = sum(d["count"] for d in formula_duplicates)
    
    report = {
        "summary": {
            "total_molecules_analyzed": sum(d["count"] for d in name_duplicates) + len(missing_cids),
            "duplicate_names": {
                "groups": len(name_duplicates),
                "total_molecules": total_name_dupes
            },
            "duplicate_formulas": {
                "groups": len(formula_duplicates),
                "total_molecules": total_formula_dupes
            },
            "missing_pubchem_cids": len(missing_cids),
            "test_molecules": len(test_molecules)
        },
        "duplicate_names": name_duplicates,
        "duplicate_formulas": formula_duplicates,
        "missing_pubchem_cids": missing_cids,
        "test_molecules": test_molecules,
        "recommendations": {
            "general": "Review duplicate groups and identify molecules that should be merged or renamed",
            "missing_cids": "Prioritize adding PubChem CIDs to non-test molecules with complete data",
            "test_molecules": "Consider adding a 'TEST_' prefix to test molecule names for clarity"
        }
    }
    
    return report

def main():
    """Identify and report on molecule duplicates."""
    print("Identifying molecule duplicates and issues...")
    
    # Connect to database
    conn = connect_to_db()
    
    try:
        # Find duplicate names
        print("Finding molecules with duplicate names...")
        name_duplicates = find_duplicate_names(conn)
        print(f"Found {len(name_duplicates)} groups of molecules with duplicate names")
        
        # Find duplicate formulas
        print("Finding molecules with duplicate molecular formulas...")
        formula_duplicates = find_duplicate_formulas(conn)
        print(f"Found {len(formula_duplicates)} groups of molecules with duplicate formulas")
        
        # Find molecules missing CIDs
        print("Finding molecules with missing PubChem CIDs...")
        missing_cids = find_molecules_missing_cids(conn)
        print(f"Found {len(missing_cids)} molecules with missing PubChem CIDs")
        
        # Analyze test molecules
        print("Analyzing test molecules...")
        test_molecules = analyze_test_molecules(conn)
        print(f"Found {len(test_molecules)} test molecules")
        
        # Generate report
        print("\nGenerating comprehensive report...")
        report = generate_report(name_duplicates, formula_duplicates, missing_cids, test_molecules)
        
        # Save report to file
        with open("molecule_duplicates_report.json", "w") as f:
            json.dump(report, f, indent=2, default=str)
        
        print("\nSummary:")
        print(f"Total molecules analyzed: {report['summary']['total_molecules_analyzed']}")
        print(f"Duplicate name groups: {report['summary']['duplicate_names']['groups']}")
        print(f"Duplicate formula groups: {report['summary']['duplicate_formulas']['groups']}")
        print(f"Molecules missing PubChem CIDs: {report['summary']['missing_pubchem_cids']}")
        print(f"Test molecules: {report['summary']['test_molecules']}")
        
        print("\nDetailed report saved to 'molecule_duplicates_report.json'")
        print("\nRecommendations:")
        print("1. Review the report to identify molecules that should be merged or renamed")
        print("2. Add PubChem CIDs to non-test molecules with complete data")
        print("3. Consider adding a 'TEST_' prefix to test molecule names for clarity")
        print("4. Create a data cleanup plan based on the findings")
    
    except Exception as e:
        print(f"Error identifying duplicates: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()