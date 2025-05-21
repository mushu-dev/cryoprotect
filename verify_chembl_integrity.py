#!/usr/bin/env python3
"""
ChEMBL Data Integrity Spot-Check Script

This script implements an optional spot-check for data integrity by programmatically
verifying a sample of ChEMBL records for expected field values. It:

1. Retrieves ChEMBL records from the database
2. Fetches the same records from the ChEMBL API
3. Compares key properties (molecular weight, LogP, etc.) between the two sources
4. Reports any discrepancies and logs the results

Based on specifications in .specs/chembl_data_verification.md
Related Task: task-imp-wv-1-1-integrity-check
"""

import os
import json
import logging
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Set

# Import the official ChEMBL client
try:
    from chembl_webresource_client.new_client import new_client
    chembl_client_available = True
except ImportError:
    chembl_client_available = False

from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/chembl_integrity_check.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Ensure logs directory exists
Path("logs").mkdir(exist_ok=True)

# Ensure reports directory exists
Path("reports").mkdir(exist_ok=True)

def update_project_state(task_id: str, status: str, log_message: str) -> None:
    """
    Update the project state file with the current task status.
    
    Args:
        task_id: The ID of the task to update
        status: The new status of the task
        log_message: The log message to add
    """
    try:
        # Read the current project state as text
        with open("project_state.json", "r") as f:
            content = f.read()
        
        # Find the task in the file content
        task_pattern = f'"{task_id}"\\s*:\\s*\\{{'
        if task_id not in content:
            logger.error(f"Task {task_id} not found in project state")
            return
            
        # Find the status field
        status_pattern = f'"status"\\s*:\\s*"[^"]*"'
        content = content.replace(f'"status": "Running"', f'"status": "{status}"', 1)
        
        # Find the log field and add the new message
        log_message_formatted = f'{datetime.now().strftime("%Y-%m-%d")}: {log_message}'
        if '"log": [' in content:
            # Add to existing log
            content = content.replace('"log": [', f'"log": [\n        "{log_message_formatted}",', 1)
        else:
            # Create new log
            content = content.replace('"outputs": []', f'"outputs": ["verify_chembl_integrity.py"],\n      "log": [\n        "{log_message_formatted}"\n      ]', 1)
        
        # Write the updated content
        with open("project_state.json", "w") as f:
            f.write(content)
        
        logger.info(f"Updated project state: {task_id} -> {status}")
    
    except Exception as e:
        logger.error(f"Error updating project state: {str(e)}")

def initialize_supabase() -> Optional[Client]:
    """
    Initialize the Supabase client using credentials from .env file.
    
    Returns:
        Supabase client or None if initialization fails
    """
    try:
        # Load environment variables from .env file
        load_dotenv()
        
        # Get Supabase credentials from environment
        supabase_url = os.getenv("SUPABASE_URL")
        
        # Look for service role key in different environment variables
        service_role_key = None
        for key_name in ['SUPABASE_SERVICE_ROLE_KEY', 'SUPABASE_SERVICE_KEY', 'SUPABASE_KEY']:
            potential_key = os.getenv(key_name)
            if potential_key:
                # Check if it's likely a service role key (less strict check)
                if potential_key.startswith("eyJ"):
                    service_role_key = potential_key
                    logger.info(f"Using key from {key_name}")
                    break
        
        if not supabase_url or not service_role_key:
            logger.error("Missing Supabase credentials in .env file")
            return None
        
        # Create Supabase client with service role key
        logger.info(f"Initializing Supabase client with URL: {supabase_url}")
        supabase = create_client(supabase_url, service_role_key)
        
        # Test connectivity with a simple query
        try:
            response = supabase.from_("property_types").select("id").limit(1).execute()
            logger.info("Successfully connected to Supabase")
            return supabase
        except Exception as e:
            logger.error(f"Failed to connect to Supabase: {str(e)}")
            return None
    
    except Exception as e:
        logger.error(f"Error initializing Supabase client: {str(e)}")
        return None

def initialize_chembl_client() -> bool:
    """
    Initialize the ChEMBL client.
    
    Returns:
        True if initialization is successful, False otherwise
    """
    if not chembl_client_available:
        logger.error("ChEMBL client not available. Please install chembl_webresource_client package.")
        return False
    
    try:
        # Test connectivity with a simple query
        molecule = new_client.molecule
        test_query = molecule.get('CHEMBL25')
        if test_query:
            logger.info("Successfully connected to ChEMBL API")
            return True
        else:
            logger.error("Failed to connect to ChEMBL API")
            return False
    except Exception as e:
        logger.error(f"Error initializing ChEMBL client: {str(e)}")
        return False

def get_chembl_ids_from_database(supabase: Client, limit: int = 10) -> List[str]:
    """
    Get ChEMBL IDs from the database.
    
    Args:
        supabase: Supabase client
        limit: Maximum number of ChEMBL IDs to retrieve
        
    Returns:
        List of ChEMBL IDs
    """
    try:
        # Query the database for molecules with ChEMBL IDs in data_source column
        response = supabase.from_("molecules").select("data_source").like("data_source", "ChEMBL ID:%").limit(limit).execute()
        
        if not hasattr(response, 'data') or not response.data:
            logger.warning("No ChEMBL IDs found in database")
            return []
        
        # Extract ChEMBL IDs from data_source column
        chembl_ids = []
        for item in response.data:
            data_source = item.get("data_source", "")
            if data_source.startswith("ChEMBL ID:"):
                chembl_id = data_source.replace("ChEMBL ID:", "").strip()
                chembl_ids.append(chembl_id)
        
        logger.info(f"Found {len(chembl_ids)} ChEMBL IDs in database")
        return chembl_ids
    
    except Exception as e:
        logger.error(f"Error getting ChEMBL IDs from database: {str(e)}")
        return []

def get_molecule_from_database(supabase: Client, chembl_id: str) -> Dict[str, Any]:
    """
    Get molecule details from the database.
    
    Args:
        supabase: Supabase client
        chembl_id: ChEMBL ID
        
    Returns:
        Dictionary with molecule details
    """
    try:
        # Query the database for the molecule with the given ChEMBL ID
        response = supabase.from_("molecules").select("*").like("data_source", f"ChEMBL ID: {chembl_id}%").execute()
        
        if not hasattr(response, 'data') or not response.data:
            logger.warning(f"Molecule with ChEMBL ID {chembl_id} not found in database")
            return {}
        
        molecule = response.data[0]
        molecule_id = molecule.get("id")
        
        # Get molecular properties
        properties_response = supabase.from_("molecular_properties").select("*").eq("molecule_id", molecule_id).execute()
        
        properties = {}
        if hasattr(properties_response, 'data') and properties_response.data:
            # Get property types for each property
            property_types = {}
            for prop in properties_response.data:
                property_type_id = prop.get("property_type_id")
                if property_type_id:
                    if property_type_id not in property_types:
                        property_type_response = supabase.from_("property_types").select("name,units").eq("id", property_type_id).execute()
                        if hasattr(property_type_response, 'data') and property_type_response.data:
                            property_types[property_type_id] = property_type_response.data[0]
                    
                    # Add property to properties dictionary
                    if property_type_id in property_types:
                        property_name = property_types[property_type_id].get("name", "").lower()
                        if property_name:
                            # Determine the value based on the property type
                            if prop.get("numeric_value") is not None:
                                properties[property_name] = prop.get("numeric_value")
                            elif prop.get("text_value") is not None:
                                properties[property_name] = prop.get("text_value")
                            elif prop.get("boolean_value") is not None:
                                properties[property_name] = prop.get("boolean_value")
        
        # Add properties to molecule
        molecule["properties"] = properties
        
        return molecule
    
    except Exception as e:
        logger.error(f"Error getting molecule from database: {str(e)}")
        return {}

def get_molecule_from_chembl(chembl_id: str) -> Dict[str, Any]:
    """
    Get molecule details from ChEMBL API.
    
    Args:
        chembl_id: ChEMBL ID
        
    Returns:
        Dictionary with molecule details
    """
    try:
        # Get molecule details from ChEMBL API
        molecule = new_client.molecule
        chembl_molecule = molecule.get(chembl_id)
        
        if not chembl_molecule:
            logger.warning(f"Molecule with ChEMBL ID {chembl_id} not found in ChEMBL API")
            return {}
        
        # Extract properties from molecule_properties
        properties = {}
        if "molecule_properties" in chembl_molecule:
            mol_props = chembl_molecule.get("molecule_properties", {})
            for prop_key, prop_value in mol_props.items():
                if prop_value is not None:
                    properties[prop_key.lower()] = prop_value
        
        # Add properties to molecule
        chembl_molecule["properties"] = properties
        
        return chembl_molecule
    
    except Exception as e:
        logger.error(f"Error getting molecule from ChEMBL API: {str(e)}")
        return {}

def compare_molecules(db_molecule: Dict[str, Any], chembl_molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compare molecule details from database and ChEMBL API.
    
    Args:
        db_molecule: Molecule details from database
        chembl_molecule: Molecule details from ChEMBL API
        
    Returns:
        Dictionary with comparison results
    """
    # Define the properties to compare
    # Map database property names to ChEMBL property names
    property_mapping = {
        "molecular weight": "full_mwt",
        "logp": "alogp",
        "hydrogen bond acceptor count": "hba",
        "hydrogen bond donor count": "hbd",
        "topological polar surface area": "psa",
        "rotatable bond count": "rtb"
    }
    
    # Initialize results
    results = {
        "matching_properties": [],
        "mismatched_properties": [],
        "missing_properties": [],
        "match_percentage": 0.0
    }
    
    # Get properties from database and ChEMBL
    db_properties = db_molecule.get("properties", {})
    chembl_properties = chembl_molecule.get("properties", {})
    
    # Compare properties
    total_properties = 0
    matching_properties = 0
    
    # Check each database property against ChEMBL
    for db_prop_name, db_prop_value in db_properties.items():
        # Skip if property is not in the mapping
        if db_prop_name.lower() not in property_mapping:
            continue
        
        # Get corresponding ChEMBL property name
        chembl_prop_name = property_mapping[db_prop_name.lower()]
        
        # Check if property exists in ChEMBL
        if chembl_prop_name in chembl_properties:
            total_properties += 1
            chembl_prop_value = chembl_properties[chembl_prop_name]
            
            # Compare values with tolerance for floating point values
            if isinstance(db_prop_value, (int, float)) and isinstance(chembl_prop_value, (int, float)):
                # Use relative tolerance for floating point comparison
                if abs(db_prop_value - chembl_prop_value) <= 0.01 * max(abs(db_prop_value), abs(chembl_prop_value)):
                    results["matching_properties"].append({
                        "property": db_prop_name,
                        "db_value": db_prop_value,
                        "chembl_value": chembl_prop_value
                    })
                    matching_properties += 1
                else:
                    results["mismatched_properties"].append({
                        "property": db_prop_name,
                        "db_value": db_prop_value,
                        "chembl_value": chembl_prop_value,
                        "difference": abs(db_prop_value - chembl_prop_value),
                        "percent_difference": abs(db_prop_value - chembl_prop_value) / max(abs(db_prop_value), abs(chembl_prop_value)) * 100
                    })
            elif db_prop_value == chembl_prop_value:
                results["matching_properties"].append({
                    "property": db_prop_name,
                    "db_value": db_prop_value,
                    "chembl_value": chembl_prop_value
                })
                matching_properties += 1
            else:
                results["mismatched_properties"].append({
                    "property": db_prop_name,
                    "db_value": db_prop_value,
                    "chembl_value": chembl_prop_value
                })
        else:
            results["missing_properties"].append({
                "property": db_prop_name,
                "db_value": db_prop_value
            })
    
    # Calculate match percentage
    if total_properties > 0:
        results["match_percentage"] = (matching_properties / total_properties) * 100
    
    return results

def check_for_placeholder_data(molecule: Dict[str, Any]) -> bool:
    """
    Check if a molecule contains placeholder or test data.
    
    Args:
        molecule: Molecule details
        
    Returns:
        True if placeholder data is found, False otherwise
    """
    # Define patterns that indicate placeholder or test data
    placeholder_patterns = [
        "test", "placeholder", "dummy", "example", "sample", "todo", "fixme"
    ]
    
    # Check molecule name
    name = molecule.get("name", "").lower()
    for pattern in placeholder_patterns:
        if pattern in name:
            return True
    
    # Check data_source
    data_source = molecule.get("data_source", "").lower()
    for pattern in placeholder_patterns:
        if pattern in data_source:
            return True
    
    # Check SMILES and InChI for very simple structures that might be placeholders
    smiles = molecule.get("smiles", "")
    if smiles in ["C", "CC", "CCC", "CCCC", "c1ccccc1"]:
        return True
    
    return False

def spot_check_integrity(supabase: Client, limit: int = 10) -> Dict[str, Any]:
    """
    Spot-check data integrity by comparing database records with ChEMBL API.
    
    Args:
        supabase: Supabase client
        limit: Maximum number of records to check
        
    Returns:
        Dictionary with integrity check results
    """
    logger.info(f"Spot-checking data integrity for up to {limit} ChEMBL records")
    
    # Initialize results
    results = {
        "timestamp": datetime.now().isoformat(),
        "records_checked": 0,
        "records_matched": 0,
        "records_mismatched": 0,
        "placeholder_data_found": False,
        "overall_match_percentage": 0.0,
        "details": []
    }
    
    # Check if ChEMBL client is available
    if not chembl_client_available:
        logger.error("ChEMBL client not available. Please install chembl_webresource_client package.")
        results["error"] = "ChEMBL client not available"
        return results
    
    # Initialize ChEMBL client
    if not initialize_chembl_client():
        results["error"] = "Failed to initialize ChEMBL client"
        return results
    
    # Get ChEMBL IDs from database
    chembl_ids = get_chembl_ids_from_database(supabase, limit)
    
    if not chembl_ids:
        results["error"] = "No ChEMBL IDs found in database"
        return results
    
    # Check each ChEMBL ID
    total_match_percentage = 0.0
    
    for chembl_id in chembl_ids:
        try:
            # Get molecule from database
            db_molecule = get_molecule_from_database(supabase, chembl_id)
            
            if not db_molecule:
                logger.warning(f"Molecule with ChEMBL ID {chembl_id} not found in database")
                continue
            
            # Check for placeholder data
            if check_for_placeholder_data(db_molecule):
                logger.warning(f"Placeholder data found for ChEMBL ID {chembl_id}")
                results["placeholder_data_found"] = True
            
            # Get molecule from ChEMBL API
            chembl_molecule = get_molecule_from_chembl(chembl_id)
            
            if not chembl_molecule:
                logger.warning(f"Molecule with ChEMBL ID {chembl_id} not found in ChEMBL API")
                continue
            
            # Compare molecules
            comparison = compare_molecules(db_molecule, chembl_molecule)
            
            # Add to results
            record_result = {
                "chembl_id": chembl_id,
                "name": db_molecule.get("name", ""),
                "smiles": db_molecule.get("smiles", ""),
                "comparison": comparison,
                "placeholder_data": check_for_placeholder_data(db_molecule)
            }
            
            results["details"].append(record_result)
            results["records_checked"] += 1
            
            # Update match statistics
            if comparison["match_percentage"] >= 80.0:
                results["records_matched"] += 1
            else:
                results["records_mismatched"] += 1
            
            total_match_percentage += comparison["match_percentage"]
            
            # Slight delay to be gentle on the API
            time.sleep(0.5)
            
        except Exception as e:
            logger.error(f"Error checking ChEMBL ID {chembl_id}: {str(e)}")
    
    # Calculate overall match percentage
    if results["records_checked"] > 0:
        results["overall_match_percentage"] = total_match_percentage / results["records_checked"]
    
    return results

def save_report(report: Dict[str, Any], filename: str = "chembl_integrity_check_report.json") -> str:
    """
    Save the integrity check report to a file.
    
    Args:
        report: Report dictionary
        filename: Name of the report file
        
    Returns:
        Path to the saved report file
    """
    logger.info(f"Saving report to reports/{filename}")
    
    # Create reports directory if it doesn't exist
    Path("reports").mkdir(exist_ok=True)
    
    # Save report to file
    report_path = f"reports/{filename}"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    
    return report_path

def main():
    """
    Main function to spot-check ChEMBL data integrity.
    """
    logger.info("Starting ChEMBL data integrity spot-check")
    
    # Update project state to Running
    update_project_state("task-imp-wv-1-1-integrity-check", "Running", "Started ChEMBL data integrity spot-check")
    
    # Initialize Supabase client
    supabase = initialize_supabase()
    if not supabase:
        logger.error("Failed to initialize Supabase client")
        update_project_state("task-imp-wv-1-1-integrity-check", "Error", "Failed to initialize Supabase client")
        return
    
    try:
        # Spot-check data integrity
        results = spot_check_integrity(supabase)
        
        # Save report
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_path = save_report(results, f"chembl_integrity_check_{timestamp}.json")
        
        # Update project state based on verification results
        if "error" in results:
            update_project_state(
                "task-imp-wv-1-1-integrity-check",
                "Error",
                f"ChEMBL data integrity spot-check failed: {results['error']}"
            )
        elif results["records_checked"] == 0:
            update_project_state(
                "task-imp-wv-1-1-integrity-check",
                "Done",
                f"ChEMBL data integrity spot-check completed, but no records were checked. Report saved to {report_path}."
            )
        elif results["placeholder_data_found"]:
            update_project_state(
                "task-imp-wv-1-1-integrity-check",
                "Done",
                f"ChEMBL data integrity spot-check completed. Report saved to {report_path}. " +
                f"Placeholder data was found in some records. " +
                f"Overall match percentage: {results['overall_match_percentage']:.2f}%."
            )
        elif results["overall_match_percentage"] >= 80.0:
            update_project_state(
                "task-imp-wv-1-1-integrity-check",
                "Done",
                f"ChEMBL data integrity spot-check completed successfully. Report saved to {report_path}. " +
                f"Checked {results['records_checked']} records with {results['records_matched']} matches. " +
                f"Overall match percentage: {results['overall_match_percentage']:.2f}%."
            )
        else:
            update_project_state(
                "task-imp-wv-1-1-integrity-check",
                "Done",
                f"ChEMBL data integrity spot-check completed with discrepancies. Report saved to {report_path}. " +
                f"Checked {results['records_checked']} records with {results['records_matched']} matches and {results['records_mismatched']} mismatches. " +
                f"Overall match percentage: {results['overall_match_percentage']:.2f}%."
            )
        
        # Print summary
        logger.info("ChEMBL data integrity spot-check completed")
        logger.info(f"Report saved to {report_path}")
        logger.info(f"Records checked: {results['records_checked']}")
        logger.info(f"Records matched: {results['records_matched']}")
        logger.info(f"Records mismatched: {results['records_mismatched']}")
        logger.info(f"Placeholder data found: {results['placeholder_data_found']}")
        logger.info(f"Overall match percentage: {results['overall_match_percentage']:.2f}%")
        
        return results["overall_match_percentage"] >= 80.0 and not results["placeholder_data_found"]
    
    except Exception as e:
        logger.error(f"Error during ChEMBL data integrity spot-check: {str(e)}")
        update_project_state("task-imp-wv-1-1-integrity-check", "Error", f"Error during ChEMBL data integrity spot-check: {str(e)}")
        return False

if __name__ == "__main__":
    main()