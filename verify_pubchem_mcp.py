#!/usr/bin/env python3
"""
PubChem Data Verification Script using Supabase MCP

This script verifies the integrity and completeness of PubChem data
imported into the CryoProtect database using the Supabase MCP tools.
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime
from typing import Dict, List, Any, Tuple, Optional

# Set up logging
os.makedirs('logs', exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/pubchem_verification_mcp.log')
    ]
)

logger = logging.getLogger(__name__)

# Status constants
STATUS_SUCCESS = "SUCCESS"
STATUS_WARNING = "WARNING"
STATUS_FAILURE = "FAILURE"

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Verify PubChem data in CryoProtect database using MCP")
    parser.add_argument("--project-id", help="Supabase project ID")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--output", default="reports", help="Directory for report output")
    return parser.parse_args()

class PubChemDataVerificationMCP:
    """Verification system for PubChem data in the database using MCP"""
    
    def __init__(self, project_id=None):
        """Initialize verification system with optional project ID"""
        self.project_id = project_id or os.environ.get('SUPABASE_PROJECT_ID')
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "data_source": "pubchem",
            "counts": {},
            "property_completeness": {},
            "assessments": {}
        }
        
        logger.info(f"Using Supabase Project ID: {self.project_id}")
    
    def verify_counts(self):
        """Verify counts of PubChem data in the database"""
        try:
            # Import MCP tool
            from use_mcp_tool import execute_sql
            
            # Query to count PubChem molecules
            query = """
            SELECT COUNT(*) as pubchem_count
            FROM molecules
            WHERE source = 'pubchem' OR pubchem_cid IS NOT NULL
            """
            
            result = execute_sql(query, self.project_id)
            
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
            
            result = execute_sql(query, self.project_id)
            
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
    
    def verify_property_completeness(self):
        """Verify completeness of PubChem properties"""
        try:
            # Import MCP tool
            from use_mcp_tool import execute_sql
            
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
            
            result = execute_sql(query, self.project_id)
            
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
        
        # Generate markdown report
        md_report_file = os.path.join(output_dir, "pubchem_data_quality.md")
        self.generate_markdown_report(md_report_file)
        
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
        
        return [report_file, md_report_file]
    
    def generate_markdown_report(self, output_file):
        """Generate a markdown report of verification results"""
        try:
            # Create report content
            report = []
            
            # Title and timestamp
            report.append("# PubChem Data Quality Report")
            report.append(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            report.append("")
            
            # Overall assessment
            overall_assessment = self.results.get("overall_assessment", "UNKNOWN")
            report.append("## Overall Assessment")
            
            if overall_assessment == STATUS_SUCCESS:
                report.append("**Status: SUCCESS** ✅")
                report.append("The PubChem data in the database is of high quality and completeness.")
            elif overall_assessment == STATUS_WARNING:
                report.append("**Status: WARNING** ⚠️")
                report.append("The PubChem data in the database has some quality or completeness issues that should be addressed.")
            elif overall_assessment == STATUS_FAILURE:
                report.append("**Status: FAILURE** ❌")
                report.append("The PubChem data in the database has significant quality or completeness issues that must be addressed.")
            else:
                report.append(f"**Status: {overall_assessment}**")
                report.append("The assessment of PubChem data quality could not be completed.")
            
            report.append("")
            
            # Summary statistics
            report.append("## Summary Statistics")
            report.append("")
            
            pubchem_molecules = self.results.get("counts", {}).get("pubchem_molecules", "N/A")
            pubchem_properties = self.results.get("counts", {}).get("pubchem_properties", "N/A")
            
            report.append("| Metric | Count |")
            report.append("|--------|-------|")
            report.append(f"| PubChem Molecules | {pubchem_molecules} |")
            report.append(f"| PubChem Properties | {pubchem_properties} |")
            
            report.append("")
            
            # Property completeness
            report.append("## Property Completeness")
            report.append("")
            
            property_completeness = self.results.get("property_completeness", {})
            if property_completeness:
                report.append("| Property | Count | Coverage (%) |")
                report.append("|----------|-------|-------------|")
                
                for prop_name, prop_data in property_completeness.items():
                    if prop_name != "average_coverage":
                        count = prop_data.get("count", "N/A")
                        coverage = prop_data.get("coverage_percent", "N/A")
                        if isinstance(coverage, (int, float)):
                            coverage = f"{coverage:.1f}%"
                        report.append(f"| {prop_name} | {count} | {coverage} |")
                
                avg_coverage = property_completeness.get("average_coverage", "N/A")
                if isinstance(avg_coverage, (int, float)):
                    avg_coverage = f"{avg_coverage:.1f}%"
                report.append(f"| **Average Coverage** | - | **{avg_coverage}** |")
            else:
                report.append("No property completeness data available.")
            
            report.append("")
            
            # Recommendations
            report.append("## Recommendations")
            report.append("")
            
            if overall_assessment == STATUS_SUCCESS:
                report.append("The PubChem data quality is good. Recommended next steps:")
                report.append("1. Continue with regular data quality monitoring")
                report.append("2. Consider enhancing the database with additional properties")
                report.append("3. Implement automated data quality checks in the data pipeline")
            elif overall_assessment == STATUS_WARNING:
                report.append("The PubChem data quality has some issues that should be addressed. Recommended next steps:")
                report.append("1. Address the property completeness issues identified above")
                report.append("2. Re-run the verification after fixes are implemented")
                report.append("3. Implement data quality gates in the import process")
            else:  # FAILURE or other
                report.append("The PubChem data quality has critical issues that must be addressed. Recommended next steps:")
                report.append("1. Immediately address the critical issues identified above")
                report.append("2. Investigate the root causes of data quality issues")
                report.append("3. Re-run the verification after each fix to track progress")
                report.append("4. Review and potentially revise the PubChem data import process")
            
            # Write report to file
            with open(output_file, "w") as f:
                f.write("\n".join(report))
            
            logger.info(f"Markdown report saved to {output_file}")
            
            return output_file
            
        except Exception as e:
            logger.error(f"Error generating markdown report: {e}")
            with open(output_file, "w") as f:
                f.write(f"# Error Generating Report\n\nAn error occurred while generating the report: {e}")
            return output_file

def main():
    """Main function"""
    # Parse command line arguments
    args = parse_arguments()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize verification system
    verification = PubChemDataVerificationMCP(project_id=args.project_id)
    
    # Run all verifications
    verification.run_all_verifications()
    
    # Generate report
    report_files = verification.generate_report(args.output)
    
    # Return success if overall assessment is not FAILURE
    return verification.results["overall_assessment"] != STATUS_FAILURE

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)