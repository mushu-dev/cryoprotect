#!/usr/bin/env python3
"""
Verify Database Population

This script verifies that the database has been properly populated with
scientific data, checking for completeness and quality of the data using
the connection factory.
"""

import os
import sys
import logging
import argparse
import json
import time
import traceback
import decimal
from datetime import datetime
from typing import Dict, List, Any, Optional, Union, Tuple

# Import database connection utilities
from database.connection import get_db_connection, close_all_db_connections
from sql_executor import (
    execute_query,
    with_retry,
    process_in_batches,
    get_db
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class VerificationError(Exception):
    """Exception raised for verification failures."""
    pass

def count_molecules() -> Dict[str, int]:
    """
    Count molecules in the database by various criteria.
    
    Returns:
        Dictionary with counts
    """
    # Use a single optimized query to get all counts
    combined_query = """
    SELECT
        (SELECT COUNT(*) FROM molecules) as total_count,
        (SELECT COUNT(*) FROM molecules WHERE pubchem_cid IS NOT NULL) as pubchem_count,
        (SELECT COUNT(*) FROM molecules WHERE chembl_id IS NOT NULL) as chembl_count,
        (SELECT COUNT(*) FROM molecules WHERE pubchem_cid IS NOT NULL AND chembl_id IS NOT NULL) as cross_ref_count,
        (SELECT COUNT(DISTINCT molecule_id) FROM molecular_properties) as with_props_count
    """
    
    # Execute the combined query
    result = execute_query(combined_query)
    
    if not result:
        logger.warning("Failed to get molecule counts")
        return {
            "total_molecules": 0,
            "with_pubchem_cid": 0,
            "with_chembl_id": 0,
            "with_cross_references": 0,
            "with_properties": 0
        }
    
    # Extract counts from the result
    counts = result[0]
    total_count = counts.get('total_count', 0)
    pubchem_count = counts.get('pubchem_count', 0)
    chembl_count = counts.get('chembl_count', 0)
    cross_ref_count = counts.get('cross_ref_count', 0)
    with_props_count = counts.get('with_props_count', 0)
    
    return {
        "total_molecules": total_count,
        "with_pubchem_cid": pubchem_count,
        "with_chembl_id": chembl_count,
        "with_cross_references": cross_ref_count,
        "with_properties": with_props_count
    }

def verify_reference_compounds() -> Dict[str, Any]:
    """
    Verify that all reference compounds are present with complete properties.
    
    Returns:
        Dictionary with verification results
    """
    # Reference compound ChEMBL IDs
    reference_ids = [
        "CHEMBL388978",     # Glycerol
        "CHEMBL1098659",    # DMSO
        "CHEMBL66195",      # beta-Alanine
        "CHEMBL500033",     # tert-Butanol
        "CHEMBL1487",       # Urea
        "CHEMBL6196",       # Ethylene glycol
        "CHEMBL967",        # Propylene glycol
        "CHEMBL262548",     # Trehalose
        "CHEMBL6752"        # Glycine
    ]
    
    # Check each reference compound
    results = {
        "total_reference_compounds": len(reference_ids),
        "found_reference_compounds": 0,
        "missing_reference_compounds": [],
        "incomplete_reference_compounds": [],
        "complete_reference_compounds": [],
        "details": {}
    }
    
    # Get all reference molecules in a single query
    reference_ids_str = "', '".join(reference_ids)
    molecules_query = f"""
    SELECT id, name, chembl_id, smiles, inchi, inchikey, formula, molecular_weight, pubchem_cid
    FROM molecules
    WHERE chembl_id IN ('{reference_ids_str}')
    """
    molecules_result = execute_query(molecules_query)
    
    # Create a map of chembl_id to molecule
    molecule_map = {}
    for molecule in molecules_result:
        molecule_map[molecule["chembl_id"]] = molecule
    
    # Find missing molecules
    for chembl_id in reference_ids:
        if chembl_id not in molecule_map:
            results["missing_reference_compounds"].append(chembl_id)
            results["details"][chembl_id] = {
                "status": "missing",
                "message": "Molecule not found in database"
            }
    
    # If we found any molecules, get their properties in a single query
    if molecule_map:
        # Get molecule IDs
        molecule_ids = [molecule["id"] for molecule in molecules_result]
        molecule_ids_str = "', '".join(molecule_ids)
        
        # Get properties for all molecules in a single query
        properties_query = f"""
        SELECT mp.molecule_id, pt.name, mp.numeric_value, mp.text_value
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE mp.molecule_id IN ('{molecule_ids_str}')
        """
        properties_result = execute_query(properties_query)
        
        # Create a map of molecule_id to properties
        properties_map = {}
        for prop in properties_result:
            molecule_id = prop["molecule_id"]
            if molecule_id not in properties_map:
                properties_map[molecule_id] = {}
            
            prop_name = prop["name"]
            prop_value = prop["numeric_value"] if prop["numeric_value"] is not None else prop["text_value"]
            properties_map[molecule_id][prop_name] = prop_value
        
        # Check each molecule for required properties
        required_properties = ["logP", "h_bond_donors", "h_bond_acceptors"]
        
        for chembl_id, molecule in molecule_map.items():
            molecule_id = molecule["id"]
            property_map = properties_map.get(molecule_id, {})
            
            # Check for required properties
            missing_properties = []
            for required_prop in required_properties:
                if required_prop not in property_map:
                    missing_properties.append(required_prop)
            
            if not missing_properties:
                results["complete_reference_compounds"].append(chembl_id)
                results["details"][chembl_id] = {
                    "status": "complete",
                    "name": molecule["name"],
                    "properties": property_map
                }
            else:
                results["incomplete_reference_compounds"].append(chembl_id)
                results["details"][chembl_id] = {
                    "status": "incomplete",
                    "name": molecule["name"],
                    "missing_properties": missing_properties,
                    "properties": property_map
                }
            
            results["found_reference_compounds"] += 1
    
    return results

def verify_property_completeness() -> Dict[str, Any]:
    """
    Verify the completeness of properties for all molecules.
    
    Returns:
        Dictionary with verification results
    """
    import traceback
    
    # Count molecules with different properties
    results = {
        "property_counts": {},
        "molecules_with_complete_properties": 0,
        "molecules_with_incomplete_properties": 0,
        "total_molecules_with_properties": 0,
        "property_completeness_percentage": 0
    }
    
    try:
        # Use a single optimized query to get all property counts and completeness metrics
        # Modified to handle empty result sets and ensure consistent data types
        optimized_query = """
        WITH
        property_counts AS (
            SELECT pt.name, COUNT(DISTINCT mp.molecule_id) as molecule_count
            FROM molecular_properties mp
            JOIN property_types pt ON mp.property_type_id = pt.id
            GROUP BY pt.name
        ),
        molecule_props AS (
            SELECT
                mp.molecule_id,
                COUNT(DISTINCT pt.name) as property_count,
                BOOL_OR(pt.name = 'logP') as has_logp,
                BOOL_OR(pt.name = 'h_bond_donors') as has_hbd,
                BOOL_OR(pt.name = 'h_bond_acceptors') as has_hba
            FROM molecular_properties mp
            JOIN property_types pt ON mp.property_type_id = pt.id
            GROUP BY mp.molecule_id
        ),
        completeness_metrics AS (
            SELECT
                COUNT(*) as total,
                COUNT(*) FILTER(WHERE has_logp AND has_hbd AND has_hba) as complete
            FROM molecule_props
        ),
        total_molecules AS (
            SELECT COUNT(DISTINCT molecule_id) as count
            FROM molecular_properties
        )
        SELECT
            COALESCE(json_agg(pc), '[]'::json) as property_counts,
            COALESCE((SELECT total FROM completeness_metrics), 0) as total_molecules,
            COALESCE((SELECT complete FROM completeness_metrics), 0) as complete_molecules,
            COALESCE((SELECT count FROM total_molecules), 0) as molecules_with_properties
        FROM property_counts pc
        """
        
        query_result = execute_query(optimized_query)
        
        # Handle different types of query_result
        if query_result is None:
            logger.warning("Query result is None")
            return results
        elif isinstance(query_result, int):
            logger.warning(f"Query result is an integer: {query_result}")
            # Handle integer result (could be a count)
            results["total_molecules_with_properties"] = query_result
            return results
        elif not isinstance(query_result, list):
            logger.warning(f"Query result has unexpected type: {type(query_result)}")
            return results
        elif len(query_result) == 0:
            logger.warning("Query result is an empty list")
            return results
            
        # Extract property counts
        property_counts = query_result[0].get('property_counts', [])
            
        # Handle different types of property_counts that might be returned
        if property_counts is None:
            # Handle None case
            logger.warning("property_counts is None in query result")
        elif isinstance(property_counts, list):
            # Handle list case (expected format)
            for prop in property_counts:
                if isinstance(prop, dict) and "name" in prop and "molecule_count" in prop:
                    results["property_counts"][prop["name"]] = prop["molecule_count"]
        elif isinstance(property_counts, dict):
            # Handle dictionary case
            for name, count in property_counts.items():
                results["property_counts"][name] = count
        elif isinstance(property_counts, int):
            # Handle integer case (unexpected but possible)
            logger.warning(f"property_counts is an integer ({property_counts}) instead of a list or dict")
            # No iteration needed, but we can still record this as a total count
            results["property_counts"]["total"] = property_counts
        else:
            # Handle any other unexpected type
            logger.warning(f"property_counts has unexpected type: {type(property_counts)}")
        
        # Extract completeness metrics
        total_molecules = query_result[0].get('total_molecules', 0)
        complete_molecules = query_result[0].get('complete_molecules', 0)
        molecules_with_properties = query_result[0].get('molecules_with_properties', 0)
            
        # Ensure all values are proper integers
        try:
            total_molecules = int(total_molecules) if total_molecules is not None else 0
        except (TypeError, ValueError):
            logger.warning(f"Could not convert total_molecules to int: {total_molecules}")
            total_molecules = 0
            
        try:
            complete_molecules = int(complete_molecules) if complete_molecules is not None else 0
        except (TypeError, ValueError):
            logger.warning(f"Could not convert complete_molecules to int: {complete_molecules}")
            complete_molecules = 0
            
        try:
            molecules_with_properties = int(molecules_with_properties) if molecules_with_properties is not None else 0
        except (TypeError, ValueError):
            logger.warning(f"Could not convert molecules_with_properties to int: {molecules_with_properties}")
            molecules_with_properties = 0
        
        results["molecules_with_complete_properties"] = complete_molecules
        results["molecules_with_incomplete_properties"] = total_molecules - complete_molecules
        results["total_molecules_with_properties"] = molecules_with_properties
            
        # Calculate property completeness percentage
        if total_molecules is None:
            logger.warning("total_molecules is None in query result")
            results["property_completeness_percentage"] = 0
        elif not isinstance(total_molecules, (int, float)):
            logger.warning(f"total_molecules has unexpected type: {type(total_molecules)}")
            try:
                total_molecules = int(total_molecules)
            except (TypeError, ValueError):
                logger.error(f"Could not convert total_molecules to int: {total_molecules}")
                total_molecules = 0
        
        if complete_molecules is None:
            logger.warning("complete_molecules is None in query result")
            complete_molecules = 0
        elif not isinstance(complete_molecules, (int, float)):
            logger.warning(f"complete_molecules has unexpected type: {type(complete_molecules)}")
            try:
                complete_molecules = int(complete_molecules)
            except (TypeError, ValueError):
                logger.error(f"Could not convert complete_molecules to int: {complete_molecules}")
                complete_molecules = 0
            
        # Now safely calculate the percentage
        if total_molecules > 0:
            results["property_completeness_percentage"] = (complete_molecules / total_molecules * 100)
    except TypeError as e:
        logger.error(f"TypeError in verify_property_completeness: {str(e)}")
        logger.error("Traceback:")
        traceback.print_exc()
    except Exception as e:
        logger.error(f"Error in verify_property_completeness: {str(e)}")
        logger.error("Traceback:")
        traceback.print_exc()
    
    except TypeError as e:
        logger.error(f"TypeError in verify_property_completeness: {str(e)}")
        logger.error("Traceback:")
        traceback.print_exc()
    except Exception as e:
        logger.error(f"Error in verify_property_completeness: {str(e)}")
        logger.error("Traceback:")
        traceback.print_exc()
    
    return results

def verify_query_performance() -> Dict[str, Any]:
    """
    Verify the performance of typical database queries.
    
    Returns:
        Dictionary with verification results
    """
    results = {
        "queries": []
    }
    
    # Define optimized test queries
    test_queries = [
        {
            "name": "Fetch single molecule by ID",
            "query": """
            SELECT id, name, smiles, inchi, inchikey, formula, molecular_weight, pubchem_cid, chembl_id
            FROM molecules
            ORDER BY created_at DESC
            LIMIT 1
            """,
            "params": None
        },
        {
            "name": "Fetch molecule with properties",
            "query": """
            WITH molecule AS (
                SELECT id, name, smiles, inchi, inchikey, formula, molecular_weight, pubchem_cid, chembl_id
                FROM molecules
                ORDER BY created_at DESC
                LIMIT 1
            )
            SELECT
                m.id, m.name, m.smiles, m.inchi, m.inchikey, m.formula, m.molecular_weight, m.pubchem_cid, m.chembl_id,
                array_agg(pt.name) as property_names,
                array_agg(
                    CASE
                        WHEN mp.numeric_value IS NOT NULL THEN mp.numeric_value::text
                        ELSE mp.text_value
                    END
                ) as property_values
            FROM molecule m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            GROUP BY m.id, m.name, m.smiles, m.inchi, m.inchikey, m.formula, m.molecular_weight, m.pubchem_cid, m.chembl_id
            """,
            "params": None
        },
        {
            "name": "Search molecules by name",
            "query": """
            SELECT id, name, smiles, inchi, inchikey, formula, molecular_weight, pubchem_cid, chembl_id
            FROM molecules
            WHERE LOWER(name) LIKE LOWER(%s)
            LIMIT 10
            """,
            "params": ("%glyc%",)
        },
        {
            "name": "Count molecules by property value range",
            "query": """
            SELECT COUNT(*)
            FROM molecules m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE pt.name = 'logP' AND mp.numeric_value BETWEEN -1 AND 3
            """,
            "params": None
        }
    ]
    
    # Warm up the connection and query cache
    logger.info("Warming up database connection and query cache")
    for test in test_queries:
        try:
            execute_query(test["query"], test["params"])
        except Exception as e:
            logger.warning(f"Warm-up query failed: {str(e)}")
    
    # Run each test query multiple times and measure performance
    for test in test_queries:
        query_results = {
            "name": test["name"],
            "query": test["query"],
            "runs": [],
            "average_ms": 0,
            "min_ms": 0,
            "max_ms": 0
        }
        
        # Run the query 5 times
        times = []
        for i in range(5):
            # Clear any previous query cache effects
            execute_query("SELECT 1")
            
            # Execute and time the query
            start_time = time.time()
            execute_query(test["query"], test["params"])
            end_time = time.time()
            
            elapsed_ms = (end_time - start_time) * 1000
            times.append(elapsed_ms)
            query_results["runs"].append(elapsed_ms)
            
            # Small delay between runs to reduce connection contention
            time.sleep(0.1)
        
        # Calculate statistics
        query_results["average_ms"] = sum(times) / len(times)
        query_results["min_ms"] = min(times)
        query_results["max_ms"] = max(times)
        
        results["queries"].append(query_results)
    
    # Calculate overall performance metrics
    overall_times = [q["average_ms"] for q in results["queries"]]
    results["overall_average_ms"] = sum(overall_times) / len(overall_times)
    results["performance_acceptable"] = results["overall_average_ms"] < 50  # Less than 50ms is acceptable
    
    return results

def perform_full_verification(output_report: Optional[str] = None) -> Dict[str, Any]:
    """
    Perform a full verification of the database population.
    
    Args:
        output_report: Optional path to save the verification report
        
    Returns:
        Dictionary with verification results
    """
    logger.info("Starting full database verification")
    
    # Connect to database using connection factory
    db = get_db_connection()
    
    # Test connection
    try:
        success, message = db.test_connection()
        if not success:
            logger.error(f"Failed to connect to database: {message}")
            close_all_db_connections()
            return {
                "success": False,
                "error": f"Failed to connect to database: {message}"
            }
    except Exception as e:
        logger.error(f"Database connection error: {str(e)}")
        close_all_db_connections()
        return {
            "success": False,
            "error": f"Database connection error: {str(e)}"
        }
    
    # Collect verification results
    results = {
        "timestamp": datetime.now().isoformat(),
        "success": True,
        "molecule_counts": None,
        "reference_compounds": None,
        "property_completeness": None,
        "query_performance": None
    }
    
    try:
        # Count molecules
        logger.info("Counting molecules")
        results["molecule_counts"] = count_molecules()
        
        # Verify reference compounds
        logger.info("Verifying reference compounds")
        results["reference_compounds"] = verify_reference_compounds()
        
        # Verify property completeness
        logger.info("Verifying property completeness")
        results["property_completeness"] = verify_property_completeness()
        
        # Verify query performance
        logger.info("Verifying query performance")
        results["query_performance"] = verify_query_performance()
        
        # Determine overall success
        results["success"] = (
            results["molecule_counts"]["total_molecules"] >= 5000 and
            len(results["reference_compounds"]["complete_reference_compounds"]) == 
            results["reference_compounds"]["total_reference_compounds"] and
            results["property_completeness"]["property_completeness_percentage"] >= 90 and
            results["query_performance"]["performance_acceptable"]
        )
        
        # Generate summary
        results["summary"] = {
            "total_molecules": results["molecule_counts"]["total_molecules"],
            "reference_compounds_complete": len(results["reference_compounds"]["complete_reference_compounds"]),
            "reference_compounds_total": results["reference_compounds"]["total_reference_compounds"],
            "property_completeness_percentage": results["property_completeness"]["property_completeness_percentage"],
            "average_query_time_ms": results["query_performance"]["overall_average_ms"]
        }
        
        # Save report if requested
        if output_report:
            os.makedirs(os.path.dirname(output_report), exist_ok=True)
            with open(output_report, "w", encoding="utf-8") as f:
                # Custom JSON encoder to handle Decimal objects
                class DecimalEncoder(json.JSONEncoder):
                    def default(self, obj):
                        if isinstance(obj, decimal.Decimal):
                            return float(obj)
                        return super(DecimalEncoder, self).default(obj)
                
                json.dump(results, f, indent=2, cls=DecimalEncoder)
            logger.info(f"Verification report saved to {output_report}")
        
        # Close all database connections
        close_all_db_connections()
        
        return results
        
    except Exception as e:
        logger.error(f"Error during verification: {str(e)}")
        logger.error("Traceback:")
        traceback.print_exc()
        
        results["success"] = False
        results["error"] = str(e)
        
        # Ensure the summary is populated even with error
        if "summary" not in results:
            results["summary"] = {
                "total_molecules": results.get("molecule_counts", {}).get("total_molecules", 0),
                "reference_compounds_complete": 0,
                "reference_compounds_total": 0,
                "property_completeness_percentage": 0,
                "average_query_time_ms": 0
            }
        
        # Update project state before returning
        try:
            update_project_state(results)
        except Exception as state_error:
            logger.error(f"Error updating project state: {str(state_error)}")
        
        # Close all database connections
        close_all_db_connections()
        
        return results

def update_project_state(results: Dict[str, Any]) -> None:
    """
    Update the project state with verification results.
    
    Args:
        results: Verification results dictionary
    """
    try:
        # Load current project state
        project_state_path = "project_state.json"
        if not os.path.exists(project_state_path):
            logger.error(f"Project state file not found: {project_state_path}")
            return
            
        with open(project_state_path, "r", encoding="utf-8") as f:
            project_state = json.load(f)
            
        # Find the verification phase - check for both "Verification Script Enhancement" and "Direct PostgreSQL Database Population - Verification"
        verification_phase = None
        for phase in project_state.get("highLevelPlan", []):
            if phase.get("phase") in ["Verification Script Enhancement", "Direct PostgreSQL Database Population - Verification"]:
                verification_phase = phase
                break
                
        if verification_phase:
            # Update the phase status based on verification results
            verification_phase["status"] = "Done" if results.get("success", False) else "Failed"
            
            # Add verification results to the phase
            verification_phase["verification_results"] = {
                "timestamp": results.get("timestamp", datetime.now().isoformat()),
                "success": results.get("success", False),
                "summary": results.get("summary", {}),
                "error": results.get("error", None)
            }
            
            # Also check for task-directdb-verify-01 in tasks
            if "tasks" in project_state:
                for task_id, task in project_state["tasks"].items():
                    if task_id == "task-directdb-verify-01":
                        task["status"] = "Done" if results.get("success", False) else "Failed"
                        logger.info(f"Updated task {task_id} status to {task['status']}")
            
            # Save updated project state
            with open(project_state_path, "w", encoding="utf-8") as f:
                json.dump(project_state, f, indent=2)
                
            logger.info(f"Project state updated: verification status = {verification_phase['status']}")
        else:
            # If verification phase not found, still try to update the project state
            logger.warning("Verification phase not found in project state")
            # Try to find any relevant section for database verification
            for section in project_state.get("sections", []):
                if "verification" in section.get("title", "").lower() or "database" in section.get("title", "").lower():
                    section["status"] = "Failed" if not results.get("success", False) else "Done"
                    section["last_updated"] = datetime.now().isoformat()
                    section["error"] = results.get("error", None)
                    
                    # Save updated project state
                    with open(project_state_path, "w", encoding="utf-8") as f:
                        json.dump(project_state, f, indent=2)
                    logger.info(f"Updated alternative section: {section.get('title', 'Unknown')}")
                    break
            
    except Exception as e:
        logger.error(f"Error updating project state: {str(e)}")
        logger.error("Traceback:")
        traceback.print_exc()

def generate_markdown_report(results: Dict[str, Any], output_path: str) -> None:
    """
    Generate a detailed Markdown report from verification results.
    
    Args:
        results: Verification results dictionary
        output_path: Path to save the Markdown report
    """
    try:
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("# Database Population Verification Report\n\n")
            f.write(f"Generated: {results['timestamp']}\n\n")
            
            # Overall status
            status_emoji = "PASS" if results["success"] else "FAIL"
            f.write(f"## Overall Status: {status_emoji}\n\n")
            
            if "error" in results:
                f.write(f"**Error:** {results['error']}\n\n")
                
            # Summary
            if "summary" in results:
                f.write("## Summary\n\n")
                f.write("| Metric | Value | Requirement | Status |\n")
                f.write("|--------|-------|-------------|--------|\n")
                
                summary = results["summary"]
                
                # Total molecules
                molecule_status = "PASS" if summary["total_molecules"] >= 5000 else "FAIL"
                f.write(f"| Total molecules | {summary['total_molecules']:,} | ≥ 5,000 | {molecule_status} |\n")
                
                # Reference compounds
                ref_status = "PASS" if summary["reference_compounds_complete"] == summary["reference_compounds_total"] else "FAIL"
                f.write(f"| Reference compounds | {summary['reference_compounds_complete']}/{summary['reference_compounds_total']} | All complete | {ref_status} |\n")
                
                # Property completeness
                prop_status = "PASS" if summary["property_completeness_percentage"] >= 90 else "FAIL"
                f.write(f"| Property completeness | {summary['property_completeness_percentage']:.1f}% | ≥ 90% | {prop_status} |\n")
                
                # Query performance
                perf_status = "PASS" if summary["average_query_time_ms"] < 50 else "FAIL"
                f.write(f"| Average query time | {summary['average_query_time_ms']:.2f} ms | < 50 ms | {perf_status} |\n\n")
            
            # Molecule counts
            if "molecule_counts" in results:
                f.write("## Molecule Counts\n\n")
                counts = results["molecule_counts"]
                f.write(f"- **Total molecules:** {counts['total_molecules']:,}\n")
                f.write(f"- **With PubChem CID:** {counts['with_pubchem_cid']:,}\n")
                f.write(f"- **With ChEMBL ID:** {counts['with_chembl_id']:,}\n")
                f.write(f"- **With cross-references:** {counts['with_cross_references']:,}\n")
                f.write(f"- **With properties:** {counts['with_properties']:,}\n\n")
            
            # Reference compounds
            if "reference_compounds" in results:
                f.write("## Reference Compounds\n\n")
                ref = results["reference_compounds"]
                f.write(f"- **Total reference compounds:** {ref['total_reference_compounds']}\n")
                f.write(f"- **Found reference compounds:** {ref['found_reference_compounds']}\n")
                f.write(f"- **Complete reference compounds:** {len(ref['complete_reference_compounds'])}\n")
                f.write(f"- **Incomplete reference compounds:** {len(ref['incomplete_reference_compounds'])}\n")
                f.write(f"- **Missing reference compounds:** {len(ref['missing_reference_compounds'])}\n\n")
                
                if ref['missing_reference_compounds']:
                    f.write("### Missing Reference Compounds\n\n")
                    for chembl_id in ref['missing_reference_compounds']:
                        f.write(f"- {chembl_id}\n")
                    f.write("\n")
                    
                if ref['incomplete_reference_compounds']:
                    f.write("### Incomplete Reference Compounds\n\n")
                    f.write("| ChEMBL ID | Name | Missing Properties |\n")
                    f.write("|-----------|------|-------------------|\n")
                    for chembl_id in ref['incomplete_reference_compounds']:
                        details = ref['details'][chembl_id]
                        missing_props = ", ".join(details['missing_properties'])
                        f.write(f"| {chembl_id} | {details['name']} | {missing_props} |\n")
                    f.write("\n")
            
            # Property completeness
            if "property_completeness" in results:
                f.write("## Property Completeness\n\n")
                prop = results["property_completeness"]
                f.write(f"- **Total molecules with properties:** {prop['total_molecules_with_properties']:,}\n")
                f.write(f"- **Molecules with complete properties:** {prop['molecules_with_complete_properties']:,}\n")
                f.write(f"- **Molecules with incomplete properties:** {prop['molecules_with_incomplete_properties']:,}\n")
                f.write(f"- **Property completeness percentage:** {prop['property_completeness_percentage']:.1f}%\n\n")
                
                if prop['property_counts']:
                    f.write("### Property Distribution\n\n")
                    f.write("| Property | Molecule Count |\n")
                    f.write("|----------|---------------|\n")
                    for prop_name, count in prop['property_counts'].items():
                        f.write(f"| {prop_name} | {count:,} |\n")
                    f.write("\n")
            
            # Query performance
            if "query_performance" in results:
                f.write("## Query Performance\n\n")
                perf = results["query_performance"]
                f.write(f"- **Overall average query time:** {perf['overall_average_ms']:.2f} ms\n")
                f.write(f"- **Performance acceptable:** {'Yes' if perf['performance_acceptable'] else 'No'}\n\n")
                
                if perf['queries']:
                    f.write("### Individual Query Performance\n\n")
                    f.write("| Query | Average (ms) | Min (ms) | Max (ms) |\n")
                    f.write("|-------|-------------|----------|----------|\n")
                    for query in perf['queries']:
                        f.write(f"| {query['name']} | {query['average_ms']:.2f} | {query['min_ms']:.2f} | {query['max_ms']:.2f} |\n")
                    f.write("\n")
            
        logger.info(f"Markdown report generated: {output_path}")
    except Exception as e:
        logger.error(f"Error generating Markdown report: {str(e)}")

def main():
    """CLI entry point for verification."""
    parser = argparse.ArgumentParser(description="Verify database population")
    parser.add_argument("--full-verification", action="store_true", help="Perform full verification")
    parser.add_argument("--json-report", 
                      default=f"reports/verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                      help="Output file path for JSON report")
    parser.add_argument("--markdown-report", 
                      default=f"reports/verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md",
                      help="Output file path for Markdown report")
    parser.add_argument("--update-project-state", action="store_true", help="Update project state with results")
    
    args = parser.parse_args()
    
    try:
        if args.full_verification:
            results = perform_full_verification(args.json_report)
            
            # Generate Markdown report
            try:
                generate_markdown_report(results, args.markdown_report)
            except Exception as report_error:
                logger.error(f"Error generating Markdown report: {str(report_error)}")
                logger.error(traceback.format_exc())
            
            # Always update project state, regardless of the --update-project-state flag
            # This ensures the project state is updated even when errors occur
            try:
                update_project_state(results)
            except Exception as state_error:
                logger.error(f"Error updating project state: {str(state_error)}")
                logger.error(traceback.format_exc())
            
            if results.get("success", False):
                logger.info("Verification PASSED!")
                logger.info(f"Summary: {results.get('summary', {})}")
                return 0
            else:
                logger.error("Verification FAILED!")
                if "error" in results:
                    logger.error(f"Error: {results['error']}")
                if "summary" in results:
                    logger.info(f"Summary: {results['summary']}")
                return 1
        else:
            logger.error("No verification type specified. Use --full-verification")
            return 1
    except Exception as e:
        logger.error(f"Verification failed: {str(e)}")
        logger.error("Traceback:")
        traceback.print_exc()
        
        # Create minimal results structure for project state update
        failure_results = {
            "timestamp": datetime.now().isoformat(),
            "success": False,
            "error": str(e),
            "summary": {
                "total_molecules": 0,
                "reference_compounds_complete": 0,
                "reference_compounds_total": 0,
                "property_completeness_percentage": 0,
                "average_query_time_ms": 0
            }
        }
        
        # Always try to update project state on failure
        try:
            update_project_state(failure_results)
        except Exception as state_error:
            logger.error(f"Error updating project state after verification failure: {str(state_error)}")
        
        return 1

if __name__ == "__main__":
    sys.exit(main())