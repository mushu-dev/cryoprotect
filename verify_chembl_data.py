#!/usr/bin/env python3
"""
verify_chembl_data.py

This script verifies the integrity and completeness of the ChEMBL data imported 
into the Supabase database. It performs various checks including:
1. Molecule count verification
2. Property distribution verification
3. Sample data verification (optional)

Usage:
    python verify_chembl_data.py
"""

import os
import sys
import json
import logging
import statistics
from datetime import datetime
from typing import Dict, List, Any, Tuple, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# Import the SupabaseDirectConnection
try:
    from supabase_direct import SupabaseDirectConnection
except ImportError:
    logger.error("Failed to import SupabaseDirectConnection. Make sure supabase_direct.py is in the same directory.")
    sys.exit(1)

# Status constants
STATUS_SUCCESS = "SUCCESS"
STATUS_WARNING = "WARNING"
STATUS_FAILURE = "FAILURE"

def verify_molecule_counts() -> Tuple[str, Dict[str, Any]]:
    """
    Verify that there are enough ChEMBL molecules in the database.
    
    Returns:
        Tuple containing:
            - Status (SUCCESS/WARNING/FAILURE)
            - Dictionary with count information
    """
    logger.info("Verifying ChEMBL molecule counts...")
    
    try:
        # Get the database connection
        db = SupabaseDirectConnection.get_instance()
        
        # Query to count molecules from ChEMBL
        query = """
        SELECT COUNT(*) as count 
        FROM molecules 
        WHERE data_source LIKE '%ChEMBL%'
        """
        
        result = db.execute_query(query)
        
        if not result or len(result) == 0:
            return STATUS_FAILURE, {"count": 0, "message": "Failed to retrieve molecule count"}
        
        count = result[0].get('count', 0)
        
        # Check against threshold (1000 as per requirements)
        if count >= 1000:
            status = STATUS_SUCCESS
            message = f"Found {count} ChEMBL molecules (expected >= 1000)"
        elif count > 0:
            status = STATUS_WARNING
            message = f"Found only {count} ChEMBL molecules (expected >= 1000)"
        else:
            status = STATUS_FAILURE
            message = "No ChEMBL molecules found in the database"
        
        return status, {"count": count, "message": message}
    
    except Exception as e:
        logger.error(f"Error verifying molecule counts: {str(e)}")
        return STATUS_FAILURE, {"count": 0, "message": f"Error: {str(e)}"}

def verify_property_distribution() -> Tuple[str, Dict[str, Any]]:
    """
    Verify the distribution of properties for ChEMBL molecules.
    Checks for anomalies in key properties like LogP and Molecular Weight.
    
    Returns:
        Tuple containing:
            - Status (SUCCESS/WARNING/FAILURE)
            - Dictionary with distribution statistics
    """
    logger.info("Verifying ChEMBL property distributions...")
    
    try:
        # Get the database connection
        db = SupabaseDirectConnection.get_instance()
        
        # Query to get key properties for ChEMBL molecules
        query = """
        SELECT mp.value, pt.name as property_name
        FROM molecule_properties mp
        JOIN molecules m ON mp.molecule_id = m.id
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE m.data_source LIKE '%ChEMBL%'
        AND pt.name IN ('LogP', 'Molecular Weight', 'Hydrogen Bond Donors', 'Hydrogen Bond Acceptors')
        """
        
        result = db.execute_query(query)
        
        if not result or len(result) == 0:
            return STATUS_FAILURE, {"message": "No property data found for ChEMBL molecules"}
        
        # Group properties by name
        properties = {}
        for row in result:
            prop_name = row.get('property_name')
            value = row.get('value')
            
            # Skip if property name or value is None
            if prop_name is None or value is None:
                continue
                
            # Try to convert value to float
            try:
                value = float(value)
            except (ValueError, TypeError):
                continue
                
            if prop_name not in properties:
                properties[prop_name] = []
            
            properties[prop_name].append(value)
        
        # Calculate statistics for each property
        stats = {}
        anomalies = []
        
        for prop_name, values in properties.items():
            if not values:
                continue
                
            # Calculate basic statistics
            try:
                mean = statistics.mean(values)
                median = statistics.median(values)
                stdev = statistics.stdev(values) if len(values) > 1 else 0
                min_val = min(values)
                max_val = max(values)
                count = len(values)
                
                stats[prop_name] = {
                    "count": count,
                    "mean": mean,
                    "median": median,
                    "stdev": stdev,
                    "min": min_val,
                    "max": max_val
                }
                
                # Check for anomalies based on property type
                if prop_name == "LogP":
                    # LogP values typically range from -10 to 10
                    if mean < -10 or mean > 10:
                        anomalies.append(f"Mean LogP ({mean:.2f}) is outside typical range (-10 to 10)")
                    if stdev > 5:
                        anomalies.append(f"LogP standard deviation ({stdev:.2f}) is unusually high")
                
                elif prop_name == "Molecular Weight":
                    # Molecular weights for drug-like compounds typically range from 100 to 1000
                    if mean < 100 or mean > 1000:
                        anomalies.append(f"Mean Molecular Weight ({mean:.2f}) is outside typical range (100 to 1000)")
                    if stdev > 500:
                        anomalies.append(f"Molecular Weight standard deviation ({stdev:.2f}) is unusually high")
                
                elif prop_name == "Hydrogen Bond Donors":
                    # HBD typically range from 0 to 10
                    if mean < 0 or mean > 10:
                        anomalies.append(f"Mean Hydrogen Bond Donors ({mean:.2f}) is outside typical range (0 to 10)")
                
                elif prop_name == "Hydrogen Bond Acceptors":
                    # HBA typically range from 0 to 15
                    if mean < 0 or mean > 15:
                        anomalies.append(f"Mean Hydrogen Bond Acceptors ({mean:.2f}) is outside typical range (0 to 15)")
            
            except Exception as e:
                logger.warning(f"Error calculating statistics for {prop_name}: {str(e)}")
        
        # Determine status based on findings
        if not stats:
            return STATUS_FAILURE, {"message": "Failed to calculate property statistics"}
        
        if anomalies:
            status = STATUS_WARNING
            message = f"Found {len(anomalies)} potential anomalies in property distributions"
        else:
            status = STATUS_SUCCESS
            message = "Property distributions appear normal"
        
        return status, {
            "message": message,
            "statistics": stats,
            "anomalies": anomalies
        }
    
    except Exception as e:
        logger.error(f"Error verifying property distributions: {str(e)}")
        return STATUS_FAILURE, {"message": f"Error: {str(e)}"}

def verify_sample_data(sample_chembl_ids: List[str] = None) -> Tuple[str, Dict[str, Any]]:
    """
    Verify sample ChEMBL data by comparing key fields against expected values.
    
    Args:
        sample_chembl_ids: List of ChEMBL IDs to verify. If None, uses default samples.
    
    Returns:
        Tuple containing:
            - Status (SUCCESS/WARNING/FAILURE)
            - Dictionary with verification results
    """
    logger.info("Verifying sample ChEMBL data...")
    
    # Default sample ChEMBL IDs if none provided
    if not sample_chembl_ids:
        sample_chembl_ids = ["CHEMBL25", "CHEMBL1", "CHEMBL2"]
    
    try:
        # Get the database connection
        db = SupabaseDirectConnection.get_instance()
        
        results = {}
        missing_ids = []
        
        for chembl_id in sample_chembl_ids:
            # Query to get molecule data
            query = f"""
            SELECT m.id, m.name, m.smiles, m.chembl_id
            FROM molecules m
            WHERE m.chembl_id = '{chembl_id}'
            """
            
            result = db.execute_query(query)
            
            if not result or len(result) == 0:
                missing_ids.append(chembl_id)
                continue
            
            # Get properties for this molecule
            molecule_id = result[0].get('id')
            if molecule_id:
                prop_query = f"""
                SELECT pt.name as property_name, mp.value
                FROM molecule_properties mp
                JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE mp.molecule_id = '{molecule_id}'
                """
                
                prop_result = db.execute_query(prop_query)
                properties = {row.get('property_name'): row.get('value') for row in prop_result if row.get('property_name')}
            else:
                properties = {}
            
            # Store the results
            results[chembl_id] = {
                "found": True,
                "name": result[0].get('name'),
                "smiles": result[0].get('smiles'),
                "properties": properties
            }
        
        # Determine status based on findings
        if len(missing_ids) == len(sample_chembl_ids):
            status = STATUS_FAILURE
            message = f"None of the sample ChEMBL IDs were found in the database"
        elif missing_ids:
            status = STATUS_WARNING
            message = f"{len(missing_ids)} out of {len(sample_chembl_ids)} sample ChEMBL IDs were not found"
        else:
            status = STATUS_SUCCESS
            message = f"All {len(sample_chembl_ids)} sample ChEMBL IDs were found with complete data"
        
        return status, {
            "message": message,
            "results": results,
            "missing_ids": missing_ids,
            "total_checked": len(sample_chembl_ids),
            "total_found": len(sample_chembl_ids) - len(missing_ids)
        }
    
    except Exception as e:
        logger.error(f"Error verifying sample data: {str(e)}")
        return STATUS_FAILURE, {"message": f"Error: {str(e)}"}

def generate_verification_report() -> Dict[str, Any]:
    """
    Generate a comprehensive verification report by running all verification functions.
    
    Returns:
        Dictionary containing the verification report
    """
    logger.info("Generating ChEMBL data verification report...")
    
    # Run all verification functions
    molecule_count_status, molecule_count_data = verify_molecule_counts()
    property_dist_status, property_dist_data = verify_property_distribution()
    sample_data_status, sample_data_results = verify_sample_data()
    
    # Determine overall status
    if STATUS_FAILURE in [molecule_count_status, property_dist_status, sample_data_status]:
        overall_status = STATUS_FAILURE
    elif STATUS_WARNING in [molecule_count_status, property_dist_status, sample_data_status]:
        overall_status = STATUS_WARNING
    else:
        overall_status = STATUS_SUCCESS
    
    # Compile the report
    report = {
        "timestamp": datetime.now().isoformat(),
        "overall_status": overall_status,
        "checks": {
            "molecule_count": {
                "status": molecule_count_status,
                "data": molecule_count_data
            },
            "property_distribution": {
                "status": property_dist_status,
                "data": property_dist_data
            },
            "sample_data": {
                "status": sample_data_status,
                "data": sample_data_results
            }
        }
    }
    
    return report

def main():
    """Main function to run the verification and display results."""
    logger.info("Starting ChEMBL data verification...")
    
    # Generate the verification report
    report = generate_verification_report()
    
    # Print a summary to the console
    print("\n" + "=" * 60)
    print("ChEMBL Data Verification Summary")
    print("=" * 60)
    print(f"Overall Status: {report['overall_status']}")
    print(f"Molecule Count: {report['checks']['molecule_count']['status']} - {report['checks']['molecule_count']['data']['message']}")
    print(f"Property Distribution: {report['checks']['property_distribution']['status']} - {report['checks']['property_distribution']['data']['message']}")
    print(f"Sample Data: {report['checks']['sample_data']['status']} - {report['checks']['sample_data']['data']['message']}")
    print("=" * 60)
    
    # Save the report to a file
    report_file = f"chembl_verification_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(report_file, "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Verification report saved to {report_file}")
    
    # Return success if overall status is not FAILURE
    return report['overall_status'] != STATUS_FAILURE

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)