#!/usr/bin/env python3
"""
PubChem Data Verification Script

This script verifies the integrity and completeness of PubChem data
imported into the CryoProtect database.
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime
from typing import Dict, List, Any, Tuple, Optional
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Set up logging
os.makedirs('logs', exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/pubchem_verification.log')
    ]
)

logger = logging.getLogger(__name__)

# Status constants
STATUS_SUCCESS = "SUCCESS"
STATUS_WARNING = "WARNING"
STATUS_FAILURE = "FAILURE"

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Verify PubChem data in CryoProtect database")
    parser.add_argument("--project-id", help="Supabase project ID")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--output", default="reports", help="Directory for report output")
    return parser.parse_args()

class PubChemDataVerification:
    """Verification system for PubChem data in the database"""
    
    def __init__(self, project_id=None):
        """Initialize verification system with optional project ID"""
        self.project_id = project_id
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "data_source": "pubchem",
            "counts": {},
            "null_checks": {},
            "consistency_checks": {},
            "property_completeness": {},
            "reference_compounds": {},
            "assessments": {}
        }
        
        # Print environment variables for debugging
        logger.info(f"SUPABASE_DB_HOST: {os.environ.get('SUPABASE_DB_HOST')}")
        logger.info(f"SUPABASE_DB_PORT: {os.environ.get('SUPABASE_DB_PORT')}")
        logger.info(f"SUPABASE_DB_NAME: {os.environ.get('SUPABASE_DB_NAME')}")
        logger.info(f"SUPABASE_DB_USER: {os.environ.get('SUPABASE_DB_USER')}")
        logger.info(f"SUPABASE_DB_PASSWORD: {'*****' if os.environ.get('SUPABASE_DB_PASSWORD') else 'Not set'}")
    
    def verify_counts(self):
        """Verify counts of PubChem data in the database"""
        try:
            # Import necessary tools
            try:
                from supabase_direct import SupabaseDirectConnection
                db = SupabaseDirectConnection.get_instance()
                execute_query = db.execute_query
            except ImportError:
                # Try using MCP tool if direct connection is not available
                logger.info("Using MCP tool for database queries")
                try:
                    from use_mcp_tool import execute_sql
                    execute_query = lambda query: execute_sql(query, self.project_id)
                except ImportError:
                    logger.error("Failed to import database connection modules")
                    self.results["assessments"]["counts"] = "ERROR"
                    return {}
            
            # Query to count PubChem molecules
            query = """
            SELECT COUNT(*) as pubchem_count
            FROM molecules
            WHERE source = 'pubchem' OR pubchem_cid IS NOT NULL
            """
            
            result = execute_query(query)
            
            if not result:
                logger.error("No result returned from count query")
                self.results["assessments"]["counts"] = "ERROR"
                return {}
            
            pubchem_count = result[0]["pubchem_count"]
            self.results["counts"]["pubchem_molecules"] = pubchem_count
            
            # Query to count PubChem properties
            query = """
            SELECT COUNT(*) as property_count
            FROM molecular_properties
            WHERE molecule_id IN (
                SELECT id FROM molecules
                WHERE source = 'pubchem' OR pubchem_cid IS NOT NULL
            )
            """
            
            result = execute_query(query)
            
            if not result:
                logger.error("No result returned from property count query")
                self.results["assessments"]["counts"] = "ERROR"
                return {}
            
            property_count = result[0]["property_count"]
            self.results["counts"]["pubchem_properties"] = property_count
            
            # Determine assessment based on counts
            if pubchem_count > 0:
                if property_count > 0:
                    self.results["assessments"]["counts"] = STATUS_SUCCESS
                    logger.info(f"Found {pubchem_count} PubChem molecules with {property_count} properties")
                else:
                    self.results["assessments"]["counts"] = STATUS_WARNING
                    logger.warning(f"Found {pubchem_count} PubChem molecules but no properties")
            else:
                self.results["assessments"]["counts"] = STATUS_FAILURE
                logger.warning("No PubChem molecules found in the database")
            
            return self.results["counts"]
            
        except Exception as e:
            logger.error(f"Error verifying counts: {e}")
            self.results["assessments"]["counts"] = "ERROR"
            return {}
    
    def verify_null_values(self):
        """Verify null values in PubChem data"""
        try:
            # Import necessary tools
            try:
                from supabase_direct import SupabaseDirectConnection
                db = SupabaseDirectConnection.get_instance()
                execute_query = db.execute_query
            except ImportError:
                # Try using MCP tool if direct connection is not available
                logger.info("Using MCP tool for database queries")
                try:
                    from use_mcp_tool import execute_sql
                    execute_query = lambda query: execute_sql(query, self.project_id)
                except ImportError:
                    logger.error("Failed to import database connection modules")
                    self.results["assessments"]["null_checks"] = "ERROR"
                    return {}
            
            # Query to check for null values in PubChem molecules
            query = """
            SELECT
                COUNT(*) as total,
                SUM(CASE WHEN name IS NULL THEN 1 ELSE 0 END) as null_name,
                SUM(CASE WHEN smiles IS NULL THEN 1 ELSE 0 END) as null_smiles,
                SUM(CASE WHEN pubchem_cid IS NULL THEN 1 ELSE 0 END) as null_pubchem_cid
            FROM molecules
            WHERE source = 'pubchem' OR pubchem_cid IS NOT NULL
            """
            
            result = execute_query(query)
            
            if not result:
                logger.error("No result returned from null check query")
                self.results["assessments"]["null_checks"] = "ERROR"
                return {}
            
            total = result[0]["total"]
            null_name = result[0]["null_name"]
            null_smiles = result[0]["null_smiles"]
            null_pubchem_cid = result[0]["null_pubchem_cid"]
            
            self.results["null_checks"]["molecules"] = {
                "total": total,
                "null_name": null_name,
                "null_smiles": null_smiles,
                "null_pubchem_cid": null_pubchem_cid,
                "null_name_percent": (null_name / total * 100) if total > 0 else 0,
                "null_smiles_percent": (null_smiles / total * 100) if total > 0 else 0,
                "null_pubchem_cid_percent": (null_pubchem_cid / total * 100) if total > 0 else 0
            }
            
            # Query to check for null values in PubChem properties
            query = """
            SELECT
                COUNT(*) as total,
                SUM(CASE WHEN value IS NULL THEN 1 ELSE 0 END) as null_value,
                SUM(CASE WHEN property_type IS NULL THEN 1 ELSE 0 END) as null_property_type
            FROM molecular_properties
            WHERE molecule_id IN (
                SELECT id FROM molecules
                WHERE source = 'pubchem' OR pubchem_cid IS NOT NULL
            )
            """
            
            result = execute_query(query)
            
            if not result:
                logger.error("No result returned from property null check query")
                self.results["assessments"]["null_checks"] = "ERROR"
                return {}
            
            total = result[0]["total"]
            null_value = result[0]["null_value"]
            null_property_type = result[0]["null_property_type"]
            
            self.results["null_checks"]["properties"] = {
                "total": total,
                "null_value": null_value,
                "null_property_type": null_property_type,
                "null_value_percent": (null_value / total * 100) if total > 0 else 0,
                "null_property_type_percent": (null_property_type / total * 100) if total > 0 else 0
            }
            
            # Determine assessment based on null percentages
            molecules_null_percent = self.results["null_checks"]["molecules"]["null_name_percent"]
            properties_null_percent = self.results["null_checks"]["properties"]["null_value_percent"]
            
            if molecules_null_percent > 10 or properties_null_percent > 10:
                self.results["assessments"]["null_checks"] = STATUS_FAILURE
                logger.warning(f"High percentage of null values: {molecules_null_percent:.1f}% null names, {properties_null_percent:.1f}% null property values")
            elif molecules_null_percent > 5 or properties_null_percent > 5:
                self.results["assessments"]["null_checks"] = STATUS_WARNING
                logger.warning(f"Moderate percentage of null values: {molecules_null_percent:.1f}% null names, {properties_null_percent:.1f}% null property values")
            else:
                self.results["assessments"]["null_checks"] = STATUS_SUCCESS
                logger.info(f"Low percentage of null values: {molecules_null_percent:.1f}% null names, {properties_null_percent:.1f}% null property values")
            
            return self.results["null_checks"]
            
        except Exception as e:
            logger.error(f"Error verifying null values: {e}")
            self.results["assessments"]["null_checks"] = "ERROR"
            return {}
    
    def verify_property_completeness(self):
        """Verify completeness of PubChem properties"""
        try:
            # Import necessary tools
            try:
                from supabase_direct import SupabaseDirectConnection
                db = SupabaseDirectConnection.get_instance()
                execute_query = db.execute_query
            except ImportError:
                # Try using MCP tool if direct connection is not available
                logger.info("Using MCP tool for database queries")
                try:
                    from use_mcp_tool import execute_sql
                    execute_query = lambda query: execute_sql(query, self.project_id)
                except ImportError:
                    logger.error("Failed to import database connection modules")
                    self.results["assessments"]["property_completeness"] = "ERROR"
                    return {}
            
            # Query to check for PubChem property completeness
            query = """
            SELECT
                COUNT(*) as total_pubchem_compounds,
                SUM(CASE WHEN EXISTS (
                    SELECT 1 FROM molecular_properties mp
                    WHERE mp.molecule_id = m.id AND mp.property_type = 'logP'
                ) THEN 1 ELSE 0 END) as logp_count,
                SUM(CASE WHEN EXISTS (
                    SELECT 1 FROM molecular_properties mp
                    WHERE mp.molecule_id = m.id AND mp.property_type = 'h_bond_donor_count'
                ) THEN 1 ELSE 0 END) as hbond_donor_count,
                SUM(CASE WHEN EXISTS (
                    SELECT 1 FROM molecular_properties mp
                    WHERE mp.molecule_id = m.id AND mp.property_type = 'h_bond_acceptor_count'
                ) THEN 1 ELSE 0 END) as hbond_acceptor_count
            FROM molecules m
            WHERE m.source = 'pubchem' OR m.pubchem_cid IS NOT NULL
            """
            
            result = execute_query(query)
            
            if not result:
                logger.error("No result returned from property completeness query")
                self.results["assessments"]["property_completeness"] = "ERROR"
                return {}
            
            total_pubchem_compounds = result[0]["total_pubchem_compounds"]
            logp_count = result[0]["logp_count"]
            hbond_donor_count = result[0]["hbond_donor_count"]
            hbond_acceptor_count = result[0]["hbond_acceptor_count"]
            
            property_coverage = {}
            
            if total_pubchem_compounds > 0:
                property_coverage["logP"] = {
                    "count": logp_count,
                    "coverage_percent": (logp_count / total_pubchem_compounds * 100)
                }
                
                property_coverage["h_bond_donor_count"] = {
                    "count": hbond_donor_count,
                    "coverage_percent": (hbond_donor_count / total_pubchem_compounds * 100)
                }
                
                property_coverage["h_bond_acceptor_count"] = {
                    "count": hbond_acceptor_count,
                    "coverage_percent": (hbond_acceptor_count / total_pubchem_compounds * 100)
                }
            
            # Store results
            self.results["property_completeness"] = property_coverage
            
            # Calculate average coverage
            if property_coverage:
                coverage_values = [prop["coverage_percent"] for prop in property_coverage.values()]
                avg_coverage = sum(coverage_values) / len(coverage_values)
                self.results["property_completeness"]["average_coverage"] = avg_coverage
                
                # Determine assessment
                if avg_coverage >= 90:
                    self.results["assessments"]["property_completeness"] = STATUS_SUCCESS
                    logger.info(f"Excellent property coverage: {avg_coverage:.1f}%")
                elif avg_coverage >= 70:
                    self.results["assessments"]["property_completeness"] = STATUS_WARNING
                    logger.warning(f"Moderate property coverage: {avg_coverage:.1f}%")
                else:
                    self.results["assessments"]["property_completeness"] = STATUS_FAILURE
                    logger.warning(f"Poor property coverage: {avg_coverage:.1f}%")
            else:
                self.results["assessments"]["property_completeness"] = STATUS_FAILURE
                logger.warning("No property coverage data available")
            
            return self.results["property_completeness"]
            
        except Exception as e:
            logger.error(f"Error verifying property completeness: {e}")
            self.results["assessments"]["property_completeness"] = "ERROR"
            return {}
    
    def run_all_verifications(self):
        """Run all verification checks"""
        logger.info(f"Running all verification checks for PubChem data")
        
        # Run verification checks
        self.verify_counts()
        self.verify_null_values()
        self.verify_property_completeness()
        
        # Determine overall assessment
        assessments = self.results["assessments"]
        
        if STATUS_FAILURE in assessments.values():
            self.results["overall_assessment"] = STATUS_FAILURE
        elif STATUS_WARNING in assessments.values():
            self.results["overall_assessment"] = STATUS_WARNING
        elif all(status == STATUS_SUCCESS for status in assessments.values()):
            self.results["overall_assessment"] = STATUS_SUCCESS
        else:
            self.results["overall_assessment"] = "INCOMPLETE"
        
        return self.results
    
    def generate_report(self, output_dir="reports"):
        """Generate a report of verification results"""
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate report filename
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = os.path.join(output_dir, f"pubchem_validation_{timestamp}.json")
        
        # Write report to file
        with open(report_file, "w") as f:
            json.dump(self.results, f, indent=2)
        
        logger.info(f"Verification report saved to {report_file}")
        
        # Print summary to console
        print("\n" + "=" * 60)
        print(f"PubChem Data Verification Summary")
        print("=" * 60)
        print(f"Overall Assessment: {self.results['overall_assessment']}")
        print(f"PubChem Molecules: {self.results['counts'].get('pubchem_molecules', 'N/A')}")
        print(f"PubChem Properties: {self.results['counts'].get('pubchem_properties', 'N/A')}")
        
        print("\nAssessments:")
        for check, status in self.results["assessments"].items():
            print(f"  {check}: {status}")
        
        print("=" * 60)
        
        return report_file

def main():
    """Main function"""
    # Parse command line arguments
    args = parse_arguments()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize verification system
    verification = PubChemDataVerification(project_id=args.project_id)
    
    # Run all verifications
    verification.run_all_verifications()
    
    # Generate report
    report_file = verification.generate_report(args.output)
    
    # Return success if overall assessment is not FAILURE
    return verification.results["overall_assessment"] != STATUS_FAILURE

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)