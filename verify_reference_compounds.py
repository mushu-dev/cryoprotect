#!/usr/bin/env python3
"""
Reference Compound Verification Script

This script verifies that all reference cryoprotectant compounds have been correctly
imported with complete data. It adapts the SQL queries from the validation directive
to match the actual database schema structure.
"""

import os
import sys
import json
import logging
from datetime import datetime
from typing import Dict, List, Any, Optional

# Import the Supabase MCP tools
from supabase_mcp_tools import execute_sql_on_supabase
from chembl.reference_compounds import get_reference_compound_ids, get_extended_reference_compound_ids

# Set up logging
os.makedirs('logs', exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/reference_verification.log')
    ]
)

logger = logging.getLogger(__name__)

# Status constants
STATUS_SUCCESS = "SUCCESS"
STATUS_WARNING = "WARNING"
STATUS_FAILURE = "FAILURE"

# Supabase project ID
# This is the example ID from supabase_mcp_tools.py
# Replace with the actual project ID if different
PROJECT_ID = "tsdlmynydfuypiugmkev"

class ReferenceCompoundVerifier:
    """Verification system for reference compounds in the database"""
    
    def __init__(self, project_id=None, extended=False):
        """Initialize verification system with optional project ID"""
        self.project_id = project_id or PROJECT_ID
        self.extended = extended
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "reference_compounds": {},
            "critical_properties": {},
            "overall_assessment": None
        }
        
    def verify_reference_compounds(self):
        """Verify that all reference compounds were imported"""
        try:
            # Get reference compound IDs
            if self.extended:
                reference_ids = get_extended_reference_compound_ids()
                logger.info(f"Verifying {len(reference_ids)} extended reference compounds")
            else:
                reference_ids = get_reference_compound_ids()
                logger.info(f"Verifying {len(reference_ids)} standard reference compounds")
            
            # Check each reference compound
            present_count = 0
            for chembl_id in reference_ids:
                # Check if compound exists
                # Adapted query to match actual schema structure
                query = f"SELECT id, name, molecular_weight, inchikey, smiles FROM molecules WHERE chembl_id = '{chembl_id}'"
                
                try:
                    result = execute_sql_on_supabase(self.project_id, query)
                    
                    is_present = len(result) > 0
                    compound_data = {}
                    
                    if is_present:
                        present_count += 1
                        compound_data = result[0]
                    
                    # Store result
                    self.results["reference_compounds"][chembl_id] = {
                        "present": is_present,
                        "data": compound_data if is_present else None
                    }
                    
                    if is_present:
                        logger.info(f"Found reference compound {chembl_id}: {compound_data.get('name', 'Unknown')}")
                    else:
                        logger.warning(f"Missing reference compound: {chembl_id}")
                        
                except Exception as e:
                    logger.error(f"Error checking reference compound {chembl_id}: {e}")
                    self.results["reference_compounds"][chembl_id] = {
                        "present": False,
                        "error": str(e)
                    }
            
            # Calculate overall success rate
            total = len(reference_ids)
            success_rate = present_count / total if total > 0 else 0
            
            self.results["reference_compounds"]["summary"] = {
                "total": total,
                "present": present_count,
                "missing": total - present_count,
                "success_rate": success_rate
            }
            
            # Add assessment
            if success_rate >= 0.9:  # At least 90% present
                self.results["reference_compounds"]["assessment"] = STATUS_SUCCESS
                logger.info(f"Reference compound verification: {present_count}/{total} present ({success_rate*100:.1f}%)")
            elif success_rate >= 0.7:  # At least 70% present
                self.results["reference_compounds"]["assessment"] = STATUS_WARNING
                logger.warning(f"Reference compound verification: {present_count}/{total} present ({success_rate*100:.1f}%)")
            else:
                self.results["reference_compounds"]["assessment"] = STATUS_FAILURE
                logger.warning(f"Reference compound verification: {present_count}/{total} present ({success_rate*100:.1f}%)")
            
            # Log missing compounds
            if present_count < total:
                missing = [cid for cid in reference_ids
                         if not self.results["reference_compounds"].get(cid, {}).get("present", False)]
                logger.warning(f"Missing reference compounds: {', '.join(missing)}")
            
            return self.results["reference_compounds"]
            
        except Exception as e:
            logger.error(f"Error verifying reference compounds: {e}")
            self.results["reference_compounds"]["assessment"] = "ERROR"
            self.results["reference_compounds"]["error"] = str(e)
            return {}
    
    def check_critical_properties(self):
        """Check for missing or null values in critical properties"""
        try:
            # Adapted query to match actual schema structure
            query = """
            SELECT 
              COUNT(*) as molecules_with_incomplete_data
            FROM molecules
            WHERE 
              molecular_weight IS NULL OR
              inchikey IS NULL OR
              smiles IS NULL;
            """
            
            result = execute_sql_on_supabase(self.project_id, query)
            
            if result and len(result) > 0:
                incomplete_count = result[0].get("molecules_with_incomplete_data", 0)
                
                # Get total count for percentage calculation
                total_query = "SELECT COUNT(*) as total FROM molecules;"
                total_result = execute_sql_on_supabase(self.project_id, total_query)
                total_count = total_result[0].get("total", 0) if total_result and len(total_result) > 0 else 0
                
                if total_count > 0:
                    incomplete_percent = (incomplete_count / total_count) * 100
                else:
                    incomplete_percent = 0
                
                self.results["critical_properties"] = {
                    "total_molecules": total_count,
                    "molecules_with_incomplete_data": incomplete_count,
                    "incomplete_percent": incomplete_percent
                }
                
                # Add assessment
                if incomplete_percent <= 5:  # Less than 5% incomplete
                    self.results["critical_properties"]["assessment"] = STATUS_SUCCESS
                    logger.info(f"Critical properties check: {incomplete_count}/{total_count} incomplete ({incomplete_percent:.1f}%)")
                elif incomplete_percent <= 15:  # Less than 15% incomplete
                    self.results["critical_properties"]["assessment"] = STATUS_WARNING
                    logger.warning(f"Critical properties check: {incomplete_count}/{total_count} incomplete ({incomplete_percent:.1f}%)")
                else:
                    self.results["critical_properties"]["assessment"] = STATUS_FAILURE
                    logger.warning(f"Critical properties check: {incomplete_count}/{total_count} incomplete ({incomplete_percent:.1f}%)")
                
                # Check specific reference compounds for complete data
                if self.results["reference_compounds"]:
                    # Get the key reference compounds (first 5)
                    key_reference_ids = list(get_reference_compound_ids())[:5]
                    
                    # Adapted query to match actual schema structure
                    key_query = f"""
                    SELECT 
                      id, 
                      name, 
                      molecular_weight,
                      inchikey,
                      smiles,
                      logp as logP,
                      h_bond_donor_count as h_donors,
                      h_bond_acceptor_count as h_acceptors
                    FROM molecules
                    WHERE chembl_id IN ('{"', '".join(key_reference_ids)}');
                    """
                    
                    key_result = execute_sql_on_supabase(self.project_id, key_query)
                    
                    if key_result and len(key_result) > 0:
                        self.results["critical_properties"]["key_reference_data"] = key_result
                        
                        # Check completeness of key reference compounds
                        complete_key_refs = 0
                        for compound in key_result:
                            is_complete = all([
                                compound.get("molecular_weight") is not None,
                                compound.get("inchikey") is not None,
                                compound.get("smiles") is not None,
                                compound.get("logP") is not None,
                                compound.get("h_donors") is not None,
                                compound.get("h_acceptors") is not None
                            ])
                            
                            if is_complete:
                                complete_key_refs += 1
                        
                        key_completeness = complete_key_refs / len(key_result) if len(key_result) > 0 else 0
                        self.results["critical_properties"]["key_reference_completeness"] = {
                            "total": len(key_result),
                            "complete": complete_key_refs,
                            "completeness_rate": key_completeness
                        }
                        
                        # Update assessment if key references are incomplete
                        if key_completeness < 0.9 and self.results["critical_properties"]["assessment"] != STATUS_FAILURE:
                            self.results["critical_properties"]["assessment"] = STATUS_WARNING
                            logger.warning(f"Key reference compounds completeness: {complete_key_refs}/{len(key_result)} complete ({key_completeness*100:.1f}%)")
                
                return self.results["critical_properties"]
            
            else:
                logger.error("Failed to execute critical properties check query")
                self.results["critical_properties"]["assessment"] = "ERROR"
                return {}
                
        except Exception as e:
            logger.error(f"Error checking critical properties: {e}")
            self.results["critical_properties"]["assessment"] = "ERROR"
            self.results["critical_properties"]["error"] = str(e)
            return {}
    
    def run_verification(self):
        """Run all verification checks"""
        logger.info("Running reference compound verification")
        
        # Run verification checks
        self.verify_reference_compounds()
        self.check_critical_properties()
        
        # Determine overall assessment
        assessments = [
            self.results["reference_compounds"].get("assessment"),
            self.results["critical_properties"].get("assessment")
        ]
        
        if STATUS_FAILURE in assessments:
            self.results["overall_assessment"] = STATUS_FAILURE
        elif STATUS_WARNING in assessments:
            self.results["overall_assessment"] = STATUS_WARNING
        elif all(status == STATUS_SUCCESS for status in assessments if status):
            self.results["overall_assessment"] = STATUS_SUCCESS
        else:
            self.results["overall_assessment"] = "INCOMPLETE"
        
        return self.results
    
    def generate_report(self, output_dir="reports"):
        """Generate a report of verification results"""
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate report filename
        report_file = os.path.join(output_dir, "reference_compounds_verification.json")
        
        # Write report to file
        with open(report_file, "w") as f:
            json.dump(self.results, f, indent=2)
        
        logger.info(f"Verification report saved to {report_file}")
        
        # Print summary to console
        print("\n" + "=" * 60)
        print("Reference Compound Verification Summary")
        print("=" * 60)
        print(f"Overall Assessment: {self.results['overall_assessment']}")
        
        ref_summary = self.results["reference_compounds"].get("summary", {})
        print(f"Reference Compounds: {ref_summary.get('present', 0)}/{ref_summary.get('total', 0)} present ({ref_summary.get('success_rate', 0)*100:.1f}%)")
        
        crit_props = self.results["critical_properties"]
        if "molecules_with_incomplete_data" in crit_props and "total_molecules" in crit_props:
            incomplete = crit_props["molecules_with_incomplete_data"]
            total = crit_props["total_molecules"]
            print(f"Critical Properties: {total - incomplete}/{total} complete ({(1 - incomplete/total if total > 0 else 0)*100:.1f}%)")
        
        print("\nAssessments:")
        print(f"  Reference Compounds: {self.results['reference_compounds'].get('assessment', 'N/A')}")
        print(f"  Critical Properties: {self.results['critical_properties'].get('assessment', 'N/A')}")
        
        print("=" * 60)
        
        return report_file

def main():
    """Main function"""
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Verify reference compounds in the database")
    parser.add_argument("--project-id", help="Supabase project ID")
    parser.add_argument("--extended", action="store_true", help="Check extended reference compounds")
    parser.add_argument("--output", default="reports", help="Directory for report output")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize verification system
    verifier = ReferenceCompoundVerifier(project_id=args.project_id, extended=args.extended)
    
    # Run verification
    verifier.run_verification()
    
    # Generate report
    report_file = verifier.generate_report(args.output)
    
    # Return success if overall assessment is not FAILURE
    return verifier.results["overall_assessment"] != STATUS_FAILURE

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)