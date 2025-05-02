#!/usr/bin/env python3
"""
Simple Cross-Reference Validation Script

This script verifies that molecule identifiers (PubChem CID, ChEMBL ID, etc.)
are correctly cross-referenced and consistent within the database.
"""

import os
import sys
import json
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
)
logger = logging.getLogger(__name__)

# Status constants
STATUS_SUCCESS = "SUCCESS"
STATUS_WARNING = "WARNING"
STATUS_FAILURE = "FAILURE"

def main():
    """Main function"""
    # Create results structure
    results = {
        "timestamp": datetime.now().isoformat(),
        "cross_references": {},
        "assessments": {},
        "overall_assessment": "ERROR"  # Default to ERROR
    }
    
    # Ensure output directory exists
    os.makedirs("reports", exist_ok=True)
    
    try:
        # Use the MCP tool directly
        from use_mcp_tool import supabase
        
        # Execute the query to get molecule counts by identifier source
        project_id = "tsdlmynydfuypiugmkev"
        query = """
        SELECT
          count(*) as total_molecules,
          sum(CASE WHEN pubchem_cid IS NOT NULL AND chembl_id IS NOT NULL THEN 1 ELSE 0 END) as molecules_with_both_sources,
          sum(CASE WHEN pubchem_cid IS NOT NULL AND chembl_id IS NULL THEN 1 ELSE 0 END) as pubchem_only,
          sum(CASE WHEN pubchem_cid IS NULL AND chembl_id IS NOT NULL THEN 1 ELSE 0 END) as chembl_only
        FROM molecules
        """
        
        result = supabase.execute_sql(project_id, query)
        
        if not result:
            logger.error("Failed to execute cross-reference query")
            results["assessments"]["cross_references"] = "ERROR"
        else:
            # Store source overlap statistics
            results["cross_references"]["source_overlap"] = {
                "total_molecules": result[0]["total_molecules"],
                "molecules_with_both_sources": result[0]["molecules_with_both_sources"],
                "pubchem_only": result[0]["pubchem_only"],
                "chembl_only": result[0]["chembl_only"]
            }
            
            # Calculate percentages
            total = result[0]["total_molecules"]
            if total > 0:
                results["cross_references"]["source_overlap"]["both_percent"] = (result[0]["molecules_with_both_sources"] / total) * 100
                results["cross_references"]["source_overlap"]["pubchem_only_percent"] = (result[0]["pubchem_only"] / total) * 100
                results["cross_references"]["source_overlap"]["chembl_only_percent"] = (result[0]["chembl_only"] / total) * 100
            
            # Determine assessment
            both_percent = results["cross_references"]["source_overlap"].get("both_percent", 0)
            
            if both_percent >= 20:
                results["assessments"]["cross_references"] = STATUS_SUCCESS
                logger.info(f"Good cross-reference coverage: {both_percent:.1f}% of molecules have both PubChem and ChEMBL identifiers")
            elif both_percent >= 5:
                results["assessments"]["cross_references"] = STATUS_WARNING
                logger.warning(f"Limited cross-reference coverage: {both_percent:.1f}% of molecules have both PubChem and ChEMBL identifiers")
            else:
                results["assessments"]["cross_references"] = STATUS_FAILURE
                logger.warning(f"Poor cross-reference coverage: {both_percent:.1f}% of molecules have both PubChem and ChEMBL identifiers")
            
            # Set overall assessment
            results["overall_assessment"] = results["assessments"]["cross_references"]
    except Exception as e:
        logger.error(f"Error verifying cross-references: {str(e)}")
        results["assessments"]["cross_references"] = "ERROR"
    
    # Check for identifier conflicts in cryoprotectant_master_list.json
    try:
        master_list_path = os.path.join('data', 'cryoprotectant_master_list.json')
        if os.path.exists(master_list_path):
            with open(master_list_path, 'r') as f:
                identifiers = json.load(f)
            
            # Check for conflicts
            pubchem_to_internal = {}
            chembl_to_internal = {}
            conflicts = []
            
            for internal_id, data in identifiers.items():
                if 'pubchem_cid' in data and data['pubchem_cid']:
                    pubchem_cid = str(data['pubchem_cid'])
                    if pubchem_cid in pubchem_to_internal and pubchem_to_internal[pubchem_cid] != internal_id:
                        conflicts.append({
                            "type": "pubchem_cid",
                            "value": pubchem_cid,
                            "internal_id_1": pubchem_to_internal[pubchem_cid],
                            "internal_id_2": internal_id
                        })
                    else:
                        pubchem_to_internal[pubchem_cid] = internal_id
                
                if 'chembl_id' in data and data['chembl_id']:
                    chembl_id = data['chembl_id']
                    if chembl_id in chembl_to_internal and chembl_to_internal[chembl_id] != internal_id:
                        conflicts.append({
                            "type": "chembl_id",
                            "value": chembl_id,
                            "internal_id_1": chembl_to_internal[chembl_id],
                            "internal_id_2": internal_id
                        })
                    else:
                        chembl_to_internal[chembl_id] = internal_id
            
            results["cross_references"]["identifier_conflicts"] = {
                "total_molecules": len(identifiers),
                "conflicts_found": len(conflicts),
                "conflicts": conflicts
            }
            
            # If there are conflicts, set assessment to FAILURE
            if conflicts:
                results["assessments"]["cross_references"] = STATUS_FAILURE
                results["overall_assessment"] = STATUS_FAILURE
                logger.warning(f"Found {len(conflicts)} identifier conflicts")
        else:
            logger.warning(f"Master list file not found at {master_list_path}")
            results["cross_references"]["identifier_conflicts"] = {
                "status": "SKIPPED",
                "reason": f"Master list file not found at {master_list_path}"
            }
    except Exception as e:
        logger.error(f"Error checking identifier conflicts: {str(e)}")
        results["cross_references"]["identifier_conflicts"] = {
            "status": "ERROR",
            "reason": str(e)
        }
    
    # Write report to file
    report_file = os.path.join("reports", "cross_reference_validation.json")
    with open(report_file, "w") as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"Verification report saved to {report_file}")
    
    # Print summary to console
    print("\n" + "=" * 60)
    print("Cross-Reference Validation Summary")
    print("=" * 60)
    print(f"Overall Assessment: {results['overall_assessment']}")
    
    if "cross_references" in results and "source_overlap" in results["cross_references"]:
        overlap = results["cross_references"]["source_overlap"]
        print(f"Total Molecules: {overlap.get('total_molecules', 'N/A')}")
        print(f"Molecules with both PubChem and ChEMBL IDs: {overlap.get('molecules_with_both_sources', 'N/A')}")
        print(f"Molecules with PubChem ID only: {overlap.get('pubchem_only', 'N/A')}")
        print(f"Molecules with ChEMBL ID only: {overlap.get('chembl_only', 'N/A')}")
        
        if "both_percent" in overlap:
            print(f"Percentage with both IDs: {overlap['both_percent']:.1f}%")
    
    if "cross_references" in results and "identifier_conflicts" in results["cross_references"]:
        conflicts = results["cross_references"]["identifier_conflicts"]
        if "conflicts_found" in conflicts:
            print(f"\nIdentifier Conflicts: {conflicts['conflicts_found']}")
            
            if conflicts['conflicts_found'] > 0 and "conflicts" in conflicts:
                print("\nConflict Details:")
                for conflict in conflicts["conflicts"]:
                    print(f"  {conflict['type']} {conflict['value']} is assigned to both {conflict['internal_id_1']} and {conflict['internal_id_2']}")
    
    print("\nAssessments:")
    for check, status in results["assessments"].items():
        print(f"  {check}: {status}")
    
    print("=" * 60)
    
    # Return success if overall assessment is not FAILURE
    return results["overall_assessment"] != STATUS_FAILURE

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)