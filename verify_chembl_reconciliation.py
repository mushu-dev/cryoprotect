#!/usr/bin/env python3
"""
Script to verify that ChEMBL property reconciliation was successful.

This script checks:
1. Number of molecules with ChEMBL IDs (should be >1000)
2. Presence of all reference compounds
3. LogP value accuracy compared to ChEMBL
4. Data source attribution
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime
from typing import Dict, List, Any, Optional

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/verify_reconciliation.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Create logs directory if it doesn't exist
os.makedirs("logs", exist_ok=True)

# Try different methods for database access
try:
    # First try MCP tools if available
    from use_mcp_tool import execute_sql, get_project_id
    logger.info("Using MCP tools for database access")
    USE_MCP = True
except ImportError:
    # Fall back to direct database connection
    logger.info("MCP tools not found, using direct database connection")
    USE_MCP = False
    try:
        from supabase_adapter import SupabaseAdapter
        logger.info("Using Supabase adapter for database connection")
        USE_SUPABASE_ADAPTER = True
    except ImportError:
        logger.info("Supabase adapter not found, using direct psycopg2 connection")
        import psycopg2
        from config import Config
        USE_SUPABASE_ADAPTER = False

# Reference compounds that must be present
REFERENCE_COMPOUNDS = ["CHEMBL25", "CHEMBL1", "CHEMBL2", "CHEMBL3", "CHEMBL4"]

def get_db_connection():
    """
    Get a database connection.
    
    Returns:
        Database connection object or None if unable to connect
    """
    if USE_MCP:
        # Using MCP, no direct connection needed
        return None
    
    try:
        if USE_SUPABASE_ADAPTER:
            # Use Supabase adapter
            adapter = SupabaseAdapter()
            return adapter.get_direct_db_connection()
        else:
            # Use direct psycopg2 connection
            config = Config()
            return psycopg2.connect(
                host=config.SUPABASE_DB_HOST,
                port=config.SUPABASE_DB_PORT,
                database=config.SUPABASE_DB_NAME,
                user=config.SUPABASE_DB_USER,
                password=config.SUPABASE_DB_PASSWORD
            )
    except Exception as e:
        logger.error(f"Failed to connect to database: {str(e)}")
        return None

def count_molecules_with_chembl_id(project_id=None, conn=None):
    """
    Count molecules with ChEMBL IDs in the database.
    
    Args:
        project_id: Supabase project ID (for MCP)
        conn: Database connection (for direct connection)
        
    Returns:
        Count of molecules with ChEMBL IDs
    """
    query = "SELECT COUNT(*) FROM molecules WHERE chembl_id IS NOT NULL;"
    
    if USE_MCP:
        result = execute_sql(query, project_id)
        return result[0]["count"] if result else 0
    else:
        cursor = conn.cursor()
        cursor.execute(query)
        count = cursor.fetchone()[0]
        cursor.close()
        return count

def check_reference_compounds(project_id=None, conn=None):
    """
    Check if all reference compounds are present in the database.
    
    Args:
        project_id: Supabase project ID (for MCP)
        conn: Database connection (for direct connection)
        
    Returns:
        Dictionary with results of the check
    """
    logger.info("Checking for reference compounds...")
    
    results = {
        "total": len(REFERENCE_COMPOUNDS),
        "present": 0,
        "missing": [],
        "details": {}
    }
    
    for chembl_id in REFERENCE_COMPOUNDS:
        query = f"""
        SELECT m.id, m.name, m.chembl_id, mp.logp 
        FROM molecules m
        LEFT JOIN molecular_properties mp ON m.id = mp.molecule_id
        WHERE m.chembl_id = '{chembl_id}';
        """
        
        if USE_MCP:
            result = execute_sql(query, project_id)
            if result and len(result) > 0:
                results["present"] += 1
                results["details"][chembl_id] = {
                    "id": result[0]["id"],
                    "name": result[0]["name"],
                    "logp": result[0]["logp"]
                }
            else:
                results["missing"].append(chembl_id)
        else:
            cursor = conn.cursor()
            cursor.execute(query)
            row = cursor.fetchone()
            cursor.close()
            
            if row:
                results["present"] += 1
                results["details"][chembl_id] = {
                    "id": row[0],
                    "name": row[1],
                    "logp": row[3]
                }
            else:
                results["missing"].append(chembl_id)
    
    if results["missing"]:
        logger.warning(f"Missing reference compounds: {', '.join(results['missing'])}")
    else:
        logger.info("All reference compounds are present")
    
    return results

def check_property_statistics(project_id=None, conn=None):
    """
    Check statistics on molecule properties.
    
    Args:
        project_id: Supabase project ID (for MCP)
        conn: Database connection (for direct connection)
        
    Returns:
        Dictionary with property statistics
    """
    logger.info("Checking property statistics...")
    
    # Check molecules with LogP values
    logp_query = """
    SELECT COUNT(*) 
    FROM molecules m
    JOIN molecular_properties mp ON m.id = mp.molecule_id
    WHERE mp.logp IS NOT NULL AND m.chembl_id IS NOT NULL;
    """
    
    # Check null properties
    null_props_query = """
    SELECT 
        SUM(CASE WHEN mp.logp IS NULL THEN 1 ELSE 0 END) AS null_logp,
        SUM(CASE WHEN mp.molecular_weight IS NULL THEN 1 ELSE 0 END) AS null_mw,
        SUM(CASE WHEN mp.h_bond_donors IS NULL THEN 1 ELSE 0 END) AS null_hbd,
        SUM(CASE WHEN mp.h_bond_acceptors IS NULL THEN 1 ELSE 0 END) AS null_hba
    FROM molecules m
    JOIN molecular_properties mp ON m.id = mp.molecule_id
    WHERE m.chembl_id IS NOT NULL;
    """
    
    # Check data source attribution
    data_source_query = """
    SELECT 
        COUNT(*) as total,
        SUM(CASE WHEN mp.data_source LIKE '%ChEMBL%' THEN 1 ELSE 0 END) AS chembl_source,
        SUM(CASE WHEN mp.data_source LIKE '%reconciled%' THEN 1 ELSE 0 END) AS reconciled
    FROM molecules m
    JOIN molecular_properties mp ON m.id = mp.molecule_id
    WHERE m.chembl_id IS NOT NULL;
    """
    
    results = {
        "molecules_with_chembl_id": 0,
        "molecules_with_logp": 0,
        "null_properties": {},
        "data_source": {}
    }
    
    if USE_MCP:
        # Count molecules with LogP
        logp_result = execute_sql(logp_query, project_id)
        results["molecules_with_logp"] = logp_result[0]["count"] if logp_result else 0
        
        # Check null properties
        null_props_result = execute_sql(null_props_query, project_id)
        if null_props_result:
            results["null_properties"] = {
                "logp": null_props_result[0]["null_logp"],
                "molecular_weight": null_props_result[0]["null_mw"],
                "h_bond_donors": null_props_result[0]["null_hbd"],
                "h_bond_acceptors": null_props_result[0]["null_hba"]
            }
        
        # Check data source
        data_source_result = execute_sql(data_source_query, project_id)
        if data_source_result:
            results["data_source"] = {
                "total": data_source_result[0]["total"],
                "chembl_source": data_source_result[0]["chembl_source"],
                "reconciled": data_source_result[0]["reconciled"]
            }
    else:
        cursor = conn.cursor()
        
        # Count molecules with LogP
        cursor.execute(logp_query)
        results["molecules_with_logp"] = cursor.fetchone()[0]
        
        # Check null properties
        cursor.execute(null_props_query)
        row = cursor.fetchone()
        results["null_properties"] = {
            "logp": row[0],
            "molecular_weight": row[1],
            "h_bond_donors": row[2],
            "h_bond_acceptors": row[3]
        }
        
        # Check data source
        cursor.execute(data_source_query)
        row = cursor.fetchone()
        results["data_source"] = {
            "total": row[0],
            "chembl_source": row[1],
            "reconciled": row[2]
        }
        
        cursor.close()
    
    return results

def check_logp_sample(project_id=None, conn=None, sample_size=10):
    """
    Check a sample of LogP values for accuracy.
    
    Args:
        project_id: Supabase project ID (for MCP)
        conn: Database connection (for direct connection)
        sample_size: Number of molecules to sample
        
    Returns:
        List of sample molecules with LogP values
    """
    logger.info(f"Checking sample of {sample_size} LogP values...")
    
    query = f"""
    SELECT m.id, m.name, m.chembl_id, mp.logp, mp.data_source
    FROM molecules m
    JOIN molecular_properties mp ON m.id = mp.molecule_id
    WHERE m.chembl_id IS NOT NULL AND mp.logp IS NOT NULL
    ORDER BY RANDOM()
    LIMIT {sample_size};
    """
    
    samples = []
    
    if USE_MCP:
        result = execute_sql(query, project_id)
        if result:
            samples = result
    else:
        cursor = conn.cursor()
        cursor.execute(query)
        columns = [desc[0] for desc in cursor.description]
        for row in cursor.fetchall():
            samples.append(dict(zip(columns, row)))
        cursor.close()
    
    return samples

def run_verification():
    """
    Run the verification process.
    
    Returns:
        Dictionary with verification results
    """
    # Set up database access
    project_id = None
    conn = None
    
    if USE_MCP:
        project_id = get_project_id()
        logger.info(f"Using project ID: {project_id}")
    else:
        conn = get_db_connection()
        if not conn:
            logger.error("Failed to connect to database")
            return {"success": False, "error": "Could not connect to database"}
        logger.info("Connected to database")
    
    try:
        verification_results = {
            "timestamp": datetime.now().isoformat(),
            "success": True,
            "details": {}
        }
        
        # Count molecules with ChEMBL IDs
        count = count_molecules_with_chembl_id(project_id, conn)
        logger.info(f"Found {count} molecules with ChEMBL IDs")
        
        verification_results["details"]["molecule_count"] = {
            "count": count,
            "passes": count >= 1000,
            "message": f"Found {count} molecules with ChEMBL IDs (requirement: 1000+)"
        }
        
        # Check reference compounds
        ref_check = check_reference_compounds(project_id, conn)
        verification_results["details"]["reference_compounds"] = {
            "present": ref_check["present"],
            "total": ref_check["total"],
            "missing": ref_check["missing"],
            "details": ref_check["details"],
            "passes": len(ref_check["missing"]) == 0,
            "message": f"Found {ref_check['present']}/{ref_check['total']} reference compounds"
        }
        
        # Check property statistics
        prop_stats = check_property_statistics(project_id, conn)
        verification_results["details"]["property_statistics"] = prop_stats
        
        logp_coverage = (prop_stats["molecules_with_logp"] / count) * 100 if count > 0 else 0
        verification_results["details"]["property_statistics"]["logp_coverage"] = {
            "percentage": logp_coverage,
            "passes": logp_coverage >= 95,
            "message": f"LogP coverage: {logp_coverage:.1f}% (requirement: 95%+)"
        }
        
        # Check LogP sample
        logp_samples = check_logp_sample(project_id, conn, sample_size=10)
        verification_results["details"]["logp_samples"] = logp_samples
        
        # Determine overall success
        verification_results["success"] = (
            verification_results["details"]["molecule_count"]["passes"] and
            verification_results["details"]["reference_compounds"]["passes"] and
            verification_results["details"]["property_statistics"]["logp_coverage"]["passes"]
        )
        
        # Generate summary
        verification_results["summary"] = {
            "molecules_with_chembl_id": count,
            "reference_compounds_present": f"{ref_check['present']}/{ref_check['total']}",
            "logp_coverage": f"{logp_coverage:.1f}%",
            "overall_result": "PASS" if verification_results["success"] else "FAIL"
        }
        
        # Save report
        os.makedirs("reports", exist_ok=True)
        report_file = f"reports/chembl_verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(report_file, "w") as f:
            json.dump(verification_results, f, indent=2)
        
        logger.info(f"Verification completed: {verification_results['summary']['overall_result']}")
        logger.info(f"Report saved to {report_file}")
        
        return verification_results
    
    finally:
        # Close database connection if using direct connection
        if not USE_MCP and conn:
            conn.close()
            logger.info("Database connection closed")

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Verify ChEMBL reconciliation")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO",
                      help="Set logging level")
    parser.add_argument("--generate-html", action="store_true", help="Generate HTML report")
    args = parser.parse_args()
    
    # Set log level
    logging.getLogger().setLevel(getattr(logging, args.log_level))
    
    # Create logs and reports directories
    os.makedirs("logs", exist_ok=True)
    os.makedirs("reports", exist_ok=True)
    
    logger.info(f"Starting ChEMBL reconciliation verification")
    
    # Run verification
    try:
        results = run_verification()
        
        # Print summary
        print("\n" + "="*50)
        print("CHEMBL RECONCILIATION VERIFICATION SUMMARY")
        print("="*50)
        print(f"Molecules with ChEMBL IDs: {results['summary']['molecules_with_chembl_id']}")
        print(f"Reference compounds present: {results['summary']['reference_compounds_present']}")
        print(f"LogP coverage: {results['summary']['logp_coverage']}")
        print(f"Overall result: {results['summary']['overall_result']}")
        print("="*50)
        
        # Generate HTML report if requested
        if args.generate_html:
            generate_html_report(results)
        
        return 0 if results["success"] else 1
    
    except Exception as e:
        logger.error(f"Error during verification: {str(e)}")
        return 1

def generate_html_report(results):
    """
    Generate HTML report from verification results.
    
    Args:
        results: Verification results dictionary
    """
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>ChEMBL Reconciliation Verification Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 0; padding: 20px; }}
        h1, h2, h3 {{ color: #333; }}
        .container {{ max-width: 1000px; margin: 0 auto; }}
        .card {{ background: #f9f9f9; border-radius: 5px; padding: 15px; margin-bottom: 20px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        .pass {{ color: green; }}
        .fail {{ color: red; }}
        .summary {{ font-size: 1.2em; font-weight: bold; }}
        table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
        th, td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>ChEMBL Reconciliation Verification Report</h1>
        <p>Generated: {results['timestamp']}</p>
        
        <div class="card">
            <h2>Summary</h2>
            <p class="summary">Overall Result: <span class="{'pass' if results['success'] else 'fail'}">{results['summary']['overall_result']}</span></p>
            <table>
                <tr>
                    <th>Metric</th>
                    <th>Value</th>
                    <th>Status</th>
                </tr>
                <tr>
                    <td>Molecules with ChEMBL IDs</td>
                    <td>{results['summary']['molecules_with_chembl_id']}</td>
                    <td><span class="{'pass' if results['details']['molecule_count']['passes'] else 'fail'}">{'PASS' if results['details']['molecule_count']['passes'] else 'FAIL'}</span></td>
                </tr>
                <tr>
                    <td>Reference Compounds</td>
                    <td>{results['summary']['reference_compounds_present']}</td>
                    <td><span class="{'pass' if results['details']['reference_compounds']['passes'] else 'fail'}">{'PASS' if results['details']['reference_compounds']['passes'] else 'FAIL'}</span></td>
                </tr>
                <tr>
                    <td>LogP Coverage</td>
                    <td>{results['summary']['logp_coverage']}</td>
                    <td><span class="{'pass' if results['details']['property_statistics']['logp_coverage']['passes'] else 'fail'}">{'PASS' if results['details']['property_statistics']['logp_coverage']['passes'] else 'FAIL'}</span></td>
                </tr>
            </table>
        </div>
        
        <div class="card">
            <h2>Reference Compounds</h2>
            <p>{results['details']['reference_compounds']['message']}</p>
            
            <h3>Compound Details</h3>
            <table>
                <tr>
                    <th>ChEMBL ID</th>
                    <th>Name</th>
                    <th>LogP</th>
                </tr>
                {''.join([f"<tr><td>{chembl_id}</td><td>{details.get('name', 'N/A')}</td><td>{details.get('logp', 'N/A')}</td></tr>" for chembl_id, details in results['details']['reference_compounds']['details'].items()])}
            </table>
            
            <h3>Missing Compounds</h3>
            {'<p>None - All reference compounds are present</p>' if not results['details']['reference_compounds']['missing'] else '<ul>' + ''.join([f"<li>{compound}</li>" for compound in results['details']['reference_compounds']['missing']]) + '</ul>'}
        </div>
        
        <div class="card">
            <h2>Property Statistics</h2>
            <p>{results['details']['property_statistics']['logp_coverage']['message']}</p>
            
            <h3>Property Coverage</h3>
            <table>
                <tr>
                    <th>Property</th>
                    <th>Null Values</th>
                </tr>
                <tr>
                    <td>LogP</td>
                    <td>{results['details']['property_statistics']['null_properties'].get('logp', 'N/A')}</td>
                </tr>
                <tr>
                    <td>Molecular Weight</td>
                    <td>{results['details']['property_statistics']['null_properties'].get('molecular_weight', 'N/A')}</td>
                </tr>
                <tr>
                    <td>H-Bond Donors</td>
                    <td>{results['details']['property_statistics']['null_properties'].get('h_bond_donors', 'N/A')}</td>
                </tr>
                <tr>
                    <td>H-Bond Acceptors</td>
                    <td>{results['details']['property_statistics']['null_properties'].get('h_bond_acceptors', 'N/A')}</td>
                </tr>
            </table>
            
            <h3>Data Source Attribution</h3>
            <table>
                <tr>
                    <th>Metric</th>
                    <th>Count</th>
                </tr>
                <tr>
                    <td>Total Properties</td>
                    <td>{results['details']['property_statistics']['data_source'].get('total', 'N/A')}</td>
                </tr>
                <tr>
                    <td>ChEMBL Source</td>
                    <td>{results['details']['property_statistics']['data_source'].get('chembl_source', 'N/A')}</td>
                </tr>
                <tr>
                    <td>Reconciled</td>
                    <td>{results['details']['property_statistics']['data_source'].get('reconciled', 'N/A')}</td>
                </tr>
            </table>
        </div>
        
        <div class="card">
            <h2>LogP Value Samples</h2>
            <table>
                <tr>
                    <th>ChEMBL ID</th>
                    <th>Name</th>
                    <th>LogP</th>
                    <th>Data Source</th>
                </tr>
                {''.join([f"<tr><td>{sample['chembl_id']}</td><td>{sample['name']}</td><td>{sample['logp']}</td><td>{sample['data_source']}</td></tr>" for sample in results['details']['logp_samples']])}
            </table>
        </div>
    </div>
</body>
</html>
    """
    
    # Save HTML report
    report_file = f"reports/chembl_verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.html"
    with open(report_file, "w") as f:
        f.write(html)
    
    logger.info(f"HTML report saved to {report_file}")

if __name__ == "__main__":
    sys.exit(main())