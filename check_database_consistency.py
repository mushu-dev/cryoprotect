#!/usr/bin/env python3
"""
CryoProtect v2 - Automated Database Consistency Checks

This script performs a series of automated consistency checks on the CryoProtect v2 Supabase/Postgres database
to ensure data integrity after population. It outputs a human-readable report and a machine-readable (JSON) report.

Checks performed:
- Referential integrity (foreign keys: molecular_properties, mixture_components, etc.)
- No orphaned records (e.g., properties without a parent molecule)
- No duplicate molecules (by inchikey)
- All required fields are non-null and valid
- Logical consistency (e.g., mixture components reference valid molecules, mixture totals make sense)

Usage:
    python check_database_consistency.py

Environment variables required (can be set in .env):
    SUPABASE_URL, SUPABASE_KEY

Outputs:
    - Human-readable report to stdout
    - JSON report to 'consistency_report.json'

Author: CryoProtect v2 Team
"""

import os
import sys
import json
import logging
from dotenv import load_dotenv
from supabase import create_client, Client
from service_role_helper import get_supabase_client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

def get_supabase_connection():
    """Get a Supabase client connection."""
    try:
        # Try to use the service role helper first
        supabase = get_supabase_client()
        logger.info("Connected to Supabase using service role")
        return supabase
    except Exception as e:
        logger.warning(f"Error using service role helper: {str(e)}")
        # Fall back to direct connection
        supabase_url = os.getenv("SUPABASE_URL")
        supabase_key = os.getenv("SUPABASE_KEY")
        if not supabase_url or not supabase_key:
            logger.error("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
            sys.exit(1)
        
        supabase = create_client(supabase_url, supabase_key)
        logger.info("Connected to Supabase directly")
        return supabase

def check_referential_integrity(supabase):
    """
    Checks that all foreign keys reference valid records using Supabase client.
    """
    results = []

    # 1. molecular_properties.molecule_id -> molecules.id
    response = supabase.rpc(
        'exec_sql',
        {
            'query': """
                SELECT mp.id
                FROM molecular_properties mp
                LEFT JOIN molecules m ON mp.molecule_id = m.id
                WHERE m.id IS NULL
            """
        }
    ).execute()
    
    if hasattr(response, 'data') and response.data:
        orphans = response.data
        results.append({
            "check": "molecular_properties.molecule_id references molecules.id",
            "status": "PASS" if not orphans else "FAIL",
            "details": f"{len(orphans)} orphaned molecular_properties" if orphans else "All OK",
            "orphans": [row["id"] for row in orphans]
        })
    else:
        logger.error("Failed to check molecular_properties.molecule_id references")
        results.append({
            "check": "molecular_properties.molecule_id references molecules.id",
            "status": "ERROR",
            "details": "Failed to execute query",
            "orphans": []
        })

    # 2. molecular_properties.property_type_id -> property_types.id
    response = supabase.rpc(
        'exec_sql',
        {
            'query': """
                SELECT mp.id
                FROM molecular_properties mp
                LEFT JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE pt.id IS NULL
            """
        }
    ).execute()
    
    if hasattr(response, 'data') and response.data:
        orphans = response.data
        results.append({
            "check": "molecular_properties.property_type_id references property_types.id",
            "status": "PASS" if not orphans else "FAIL",
            "details": f"{len(orphans)} orphaned molecular_properties" if orphans else "All OK",
            "orphans": [row["id"] for row in orphans]
        })
    else:
        logger.error("Failed to check molecular_properties.property_type_id references")
        results.append({
            "check": "molecular_properties.property_type_id references property_types.id",
            "status": "ERROR",
            "details": "Failed to execute query",
            "orphans": []
        })

    # 3. mixture_components.mixture_id -> mixtures.id
    response = supabase.rpc(
        'exec_sql',
        {
            'query': """
                SELECT mc.id
                FROM mixture_components mc
                LEFT JOIN mixtures m ON mc.mixture_id = m.id
                WHERE m.id IS NULL
            """
        }
    ).execute()
    
    if hasattr(response, 'data') and response.data:
        orphans = response.data
        results.append({
            "check": "mixture_components.mixture_id references mixtures.id",
            "status": "PASS" if not orphans else "FAIL",
            "details": f"{len(orphans)} orphaned mixture_components" if orphans else "All OK",
            "orphans": [row["id"] for row in orphans]
        })
    else:
        logger.error("Failed to check mixture_components.mixture_id references")
        results.append({
            "check": "mixture_components.mixture_id references mixtures.id",
            "status": "ERROR",
            "details": "Failed to execute query",
            "orphans": []
        })

    # 4. mixture_components.molecule_id -> molecules.id
    response = supabase.rpc(
        'exec_sql',
        {
            'query': """
                SELECT mc.id
                FROM mixture_components mc
                LEFT JOIN molecules m ON mc.molecule_id = m.id
                WHERE m.id IS NULL
            """
        }
    ).execute()
    
    if hasattr(response, 'data') and response.data:
        orphans = response.data
        results.append({
            "check": "mixture_components.molecule_id references molecules.id",
            "status": "PASS" if not orphans else "FAIL",
            "details": f"{len(orphans)} orphaned mixture_components" if orphans else "All OK",
            "orphans": [row["id"] for row in orphans]
        })
    else:
        logger.error("Failed to check mixture_components.molecule_id references")
        results.append({
            "check": "mixture_components.molecule_id references molecules.id",
            "status": "ERROR",
            "details": "Failed to execute query",
            "orphans": []
        })

    return results

def check_no_duplicate_molecules(supabase):
    """
    Checks for duplicate molecules by inchikey using Supabase client.
    """
    response = supabase.rpc(
        'exec_sql',
        {
            'query': """
                SELECT inchikey, COUNT(*) as count
                FROM molecules
                GROUP BY inchikey
                HAVING COUNT(*) > 1
            """
        }
    ).execute()
    
    if hasattr(response, 'data') and response.data:
        dups = response.data
        return {
            "check": "No duplicate molecules by inchikey",
            "status": "PASS" if not dups else "FAIL",
            "details": f"{len(dups)} duplicate inchikeys" if dups else "All OK",
            "duplicates": [{"inchikey": row["inchikey"], "count": row["count"]} for row in dups]
        }
    else:
        logger.error("Failed to check duplicate molecules")
        return {
            "check": "No duplicate molecules by inchikey",
            "status": "ERROR",
            "details": "Failed to execute query",
            "duplicates": []
        }

def check_required_fields(supabase):
    """
    Checks that all required fields are non-null and valid using Supabase client.
    """
    results = []

    # 1. molecules: inchikey, name, cid
    response = supabase.rpc(
        'exec_sql',
        {
            'query': """
                SELECT id FROM molecules
                WHERE inchikey IS NULL OR name IS NULL OR cid IS NULL
            """
        }
    ).execute()
    
    if hasattr(response, 'data') and response.data:
        missing = response.data
        results.append({
            "check": "Molecules required fields (inchikey, name, cid) non-null",
            "status": "PASS" if not missing else "FAIL",
            "details": f"{len(missing)} molecules missing required fields" if missing else "All OK",
            "molecule_ids": [row["id"] for row in missing]
        })
    else:
        logger.error("Failed to check molecules required fields")
        results.append({
            "check": "Molecules required fields (inchikey, name, cid) non-null",
            "status": "ERROR",
            "details": "Failed to execute query",
            "molecule_ids": []
        })

    # 2. molecular_properties: molecule_id, property_type_id, at least one value field
    response = supabase.rpc(
        'exec_sql',
        {
            'query': """
                SELECT id FROM molecular_properties
                WHERE molecule_id IS NULL OR property_type_id IS NULL
                    OR (numeric_value IS NULL AND text_value IS NULL AND boolean_value IS NULL)
            """
        }
    ).execute()
    
    if hasattr(response, 'data') and response.data:
        missing = response.data
        results.append({
            "check": "MolecularProperty required fields (molecule_id, property_type_id, value) non-null",
            "status": "PASS" if not missing else "FAIL",
            "details": f"{len(missing)} molecular_properties missing required fields" if missing else "All OK",
            "property_ids": [row["id"] for row in missing]
        })
    else:
        logger.error("Failed to check molecular_properties required fields")
        results.append({
            "check": "MolecularProperty required fields (molecule_id, property_type_id, value) non-null",
            "status": "ERROR",
            "details": "Failed to execute query",
            "property_ids": []
        })

    # 3. mixture_components: mixture_id, molecule_id, concentration, concentration_unit
    response = supabase.rpc(
        'exec_sql',
        {
            'query': """
                SELECT id FROM mixture_components
                WHERE mixture_id IS NULL OR molecule_id IS NULL OR concentration IS NULL OR concentration_unit IS NULL
            """
        }
    ).execute()
    
    if hasattr(response, 'data') and response.data:
        missing = response.data
        results.append({
            "check": "MixtureComponent required fields (mixture_id, molecule_id, concentration, concentration_unit) non-null",
            "status": "PASS" if not missing else "FAIL",
            "details": f"{len(missing)} mixture_components missing required fields" if missing else "All OK",
            "component_ids": [row["id"] for row in missing]
        })
    else:
        logger.error("Failed to check mixture_components required fields")
        results.append({
            "check": "MixtureComponent required fields (mixture_id, molecule_id, concentration, concentration_unit) non-null",
            "status": "ERROR",
            "details": "Failed to execute query",
            "component_ids": []
        })

    return results

def check_logical_consistency(supabase):
    """
    Checks logical consistency, e.g., mixture totals using Supabase client.
    """
    results = []

    # 1. Mixture components: sum of concentrations per mixture (if units are percent, should be ~100)
    response = supabase.rpc(
        'exec_sql',
        {
            'query': """
                SELECT mixture_id, SUM(concentration) as total, MIN(concentration_unit) as unit
                FROM mixture_components
                GROUP BY mixture_id
            """
        }
    ).execute()
    
    if hasattr(response, 'data') and response.data:
        rows = response.data
        bad_totals = []
        for row in rows:
            if row["unit"] and "percent" in row["unit"].lower():
                if row["total"] is None or abs(row["total"] - 100) > 2:  # Allow small error
                    bad_totals.append({"mixture_id": row["mixture_id"], "total": row["total"], "unit": row["unit"]})
            elif row["total"] is None or row["total"] <= 0:
                bad_totals.append({"mixture_id": row["mixture_id"], "total": row["total"], "unit": row["unit"]})

        results.append({
            "check": "Mixture component totals (percent mixtures sum to ~100%)",
            "status": "PASS" if not bad_totals else "FAIL",
            "details": f"{len(bad_totals)} mixtures with inconsistent totals" if bad_totals else "All OK",
            "mixtures": bad_totals
        })
    else:
        logger.error("Failed to check mixture component totals")
        results.append({
            "check": "Mixture component totals (percent mixtures sum to ~100%)",
            "status": "ERROR",
            "details": "Failed to execute query",
            "mixtures": []
        })

    # 2. Mixture components reference valid molecules (already checked in referential integrity)
    # 3. No negative concentrations
    response = supabase.rpc(
        'exec_sql',
        {
            'query': """
                SELECT id, mixture_id, molecule_id, concentration
                FROM mixture_components
                WHERE concentration < 0
            """
        }
    ).execute()
    
    if hasattr(response, 'data') and response.data:
        negs = response.data
        results.append({
            "check": "No negative concentrations in mixture_components",
            "status": "PASS" if not negs else "FAIL",
            "details": f"{len(negs)} negative concentrations" if negs else "All OK",
            "components": [{"id": row["id"], "mixture_id": row["mixture_id"], "molecule_id": row["molecule_id"], "concentration": row["concentration"]} for row in negs]
        })
    else:
        logger.error("Failed to check negative concentrations")
        results.append({
            "check": "No negative concentrations in mixture_components",
            "status": "ERROR",
            "details": "Failed to execute query",
            "components": []
        })

    return results

def run_all_checks(supabase):
    report = {
        "referential_integrity": check_referential_integrity(supabase),
        "no_duplicate_molecules": check_no_duplicate_molecules(supabase),
        "required_fields": check_required_fields(supabase),
        "logical_consistency": check_logical_consistency(supabase),
    }
    return report

def print_human_readable(report):
    print("="*60)
    print("CryoProtect v2 Database Consistency Check Report")
    print("="*60)
    for section, results in report.items():
        print(f"\n--- {section.replace('_', ' ').title()} ---")
        if isinstance(results, list):
            for res in results:
                print(f"[{res['status']}] {res['check']}: {res['details']}")
        else:
            print(f"[{results['status']}] {results['check']}: {results['details']}")
    print("\nSee 'consistency_report.json' for full details.")

def main():
    supabase = get_supabase_connection()
    report = run_all_checks(supabase)
    print_human_readable(report)
    with open("consistency_report.json", "w") as f:
        json.dump(report, f, indent=2)
    logger.info("Consistency check complete. Report saved to consistency_report.json")

if __name__ == "__main__":
    main()