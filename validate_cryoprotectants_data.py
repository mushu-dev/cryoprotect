#!/usr/bin/env python3
"""
CryoProtect Analyzer - Data Validation Script

This script validates the integrity of cryoprotectant data imported into the Supabase database.
It performs various checks on the molecule and molecular_property tables to ensure data quality.

Validation checks include:
- Required fields presence and non-null values
- Valid property type references
- Numeric values within reasonable ranges
- No duplicate molecules
- Consistency of relationships between molecules and properties

Prerequisites:
- Python 3.6+ installed
- Supabase project with the CryoProtect schema applied
- supabase-py package installed (pip install supabase)
- python-dotenv package installed (pip install python-dotenv)
"""

import os
import json
import logging
import argparse
from datetime import datetime
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
LOG_FILE = "data_validation.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables from .env file
load_dotenv()

# Supabase connection
supabase_url = os.getenv("SUPABASE_URL")
supabase_key = os.getenv("SUPABASE_KEY")

if not supabase_url or not supabase_key:
    raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")

supabase: Client = create_client(supabase_url, supabase_key)

# Validation parameters
VALIDATION_PARAMS = {
    "molecular_weight_range": (0, 2000),  # Reasonable range for molecular weights
    "logp_range": (-10, 10),              # Reasonable range for LogP values
    "tpsa_range": (0, 500),               # Reasonable range for TPSA values
    "hbond_donors_range": (0, 20),        # Reasonable range for H-Bond donors
    "hbond_acceptors_range": (0, 20),     # Reasonable range for H-Bond acceptors
    "score_range": (0, 200)               # Range for total score (based on WEIGHTS sum)
}

class ValidationReport:
    """Class to track validation results and generate a report."""
    
    def __init__(self):
        self.validation_time = datetime.now()
        self.total_molecules = 0
        self.total_properties = 0
        self.issues = []
        self.summary = {
            "molecules_checked": 0,
            "properties_checked": 0,
            "missing_required_fields": 0,
            "invalid_numeric_values": 0,
            "duplicate_molecules": 0,
            "orphaned_properties": 0,
            "invalid_property_references": 0,
            "total_issues": 0
        }
    
    def add_issue(self, issue_type, description, entity_id=None, severity="ERROR"):
        """Add an issue to the report."""
        self.issues.append({
            "type": issue_type,
            "description": description,
            "entity_id": entity_id,
            "severity": severity
        })
        
        # Update summary counts
        self.summary["total_issues"] += 1
        if issue_type in self.summary:
            self.summary[issue_type] += 1
    
    def generate_report(self):
        """Generate a formatted validation report."""
        report = {
            "validation_time": self.validation_time.isoformat(),
            "summary": self.summary,
            "issues": self.issues
        }
        return report
    
    def save_report(self, filename="validation_report.json"):
        """Save the validation report to a file."""
        report = self.generate_report()
        with open(filename, "w") as f:
            json.dump(report, f, indent=2)
        logger.info(f"Validation report saved to {filename}")
    
    def print_summary(self):
        """Print a summary of the validation results."""
        print("\n" + "="*50)
        print("VALIDATION SUMMARY")
        print("="*50)
        print(f"Validation completed at: {self.validation_time}")
        print(f"Molecules checked: {self.summary['molecules_checked']}")
        print(f"Properties checked: {self.summary['properties_checked']}")
        print(f"Total issues found: {self.summary['total_issues']}")
        print("\nIssue breakdown:")
        for issue_type, count in self.summary.items():
            if issue_type not in ["molecules_checked", "properties_checked", "total_issues"]:
                print(f"  - {issue_type}: {count}")
        print("="*50)


def authenticate_supabase():
    """Authenticate with Supabase if credentials are provided."""
    supabase_user = os.getenv("SUPABASE_USER")
    supabase_password = os.getenv("SUPABASE_PASSWORD")
    
    if supabase_user and supabase_password:
        try:
            response = supabase.auth.sign_in_with_password({
                "email": supabase_user,
                "password": supabase_password
            })
            if hasattr(response, 'error') and response.error:
                logger.warning(f"Authentication error: {response.error}")
                logger.warning("Continuing without authentication. Some database operations may fail.")
            else:
                logger.info(f"Authenticated as {supabase_user}")
                return True
        except Exception as e:
            logger.warning(f"Authentication error: {str(e)}")
            logger.warning("Continuing without authentication. Some database operations may fail.")
    else:
        logger.warning("No authentication credentials provided. Continuing without authentication.")
        logger.warning("Some database operations may fail due to Row Level Security (RLS) policies.")
    
    return False


def validate_molecules(report, project_id=None):
    """Validate molecules in the database."""
    logger.info("Starting molecule validation...")
    
    # Query to get all molecules, optionally filtered by project_id
    query = supabase.table("molecule").select("*")
    if project_id:
        query = query.eq("project_id", project_id)
    
    response = query.execute()
    
    if hasattr(response, 'error') and response.error:
        logger.error(f"Error fetching molecules: {response.error}")
        return
    
    molecules = response.data
    report.summary["molecules_checked"] = len(molecules)
    report.total_molecules = len(molecules)
    
    logger.info(f"Found {len(molecules)} molecules to validate")
    
    # Check for duplicate InChIKeys
    inchikeys = {}
    for molecule in molecules:
        inchikey = molecule.get("inchikey")
        if inchikey in inchikeys:
            report.add_issue(
                "duplicate_molecules",
                f"Duplicate molecule with InChIKey {inchikey}. IDs: {inchikeys[inchikey]}, {molecule['id']}",
                molecule['id']
            )
        else:
            inchikeys[inchikey] = molecule['id']
    
    # Validate each molecule
    for molecule in molecules:
        # Check required fields
        required_fields = ["name", "smiles", "inchi", "inchikey", "formula", "molecular_weight"]
        for field in required_fields:
            if field not in molecule or molecule[field] is None:
                report.add_issue(
                    "missing_required_fields",
                    f"Molecule {molecule['id']} is missing required field: {field}",
                    molecule['id']
                )
        
        # Validate molecular weight range
        if "molecular_weight" in molecule and molecule["molecular_weight"] is not None:
            try:
                mw = float(molecule["molecular_weight"])
                min_mw, max_mw = VALIDATION_PARAMS["molecular_weight_range"]
                if not (min_mw <= mw <= max_mw):
                    report.add_issue(
                        "invalid_numeric_values",
                        f"Molecule {molecule['id']} has molecular_weight {mw} outside valid range ({min_mw}-{max_mw})",
                        molecule['id'],
                        severity="WARNING"
                    )
            except (ValueError, TypeError):
                report.add_issue(
                    "invalid_numeric_values",
                    f"Molecule {molecule['id']} has invalid molecular_weight value: {molecule['molecular_weight']}",
                    molecule['id']
                )
    
    logger.info("Molecule validation completed")


def validate_molecular_properties(report, project_id=None):
    """Validate molecular properties in the database."""
    logger.info("Starting molecular property validation...")
    
    # First, get all molecules to check for orphaned properties
    molecule_query = supabase.table("molecule").select("id")
    if project_id:
        molecule_query = molecule_query.eq("project_id", project_id)
    
    molecule_response = molecule_query.execute()
    
    if hasattr(molecule_response, 'error') and molecule_response.error:
        logger.error(f"Error fetching molecules: {molecule_response.error}")
        return
    
    molecule_ids = [m["id"] for m in molecule_response.data]
    molecule_id_set = set(molecule_ids)
    
    # Query to get all molecular properties
    property_query = supabase.table("molecular_property").select("*")
    if molecule_ids:
        # Filter properties to only those related to our molecules
        property_query = property_query.in_("molecule_id", molecule_ids)
    
    property_response = property_query.execute()
    
    if hasattr(property_response, 'error') and property_response.error:
        logger.error(f"Error fetching molecular properties: {property_response.error}")
        return
    
    properties = property_response.data
    report.summary["properties_checked"] = len(properties)
    report.total_properties = len(properties)
    
    logger.info(f"Found {len(properties)} molecular properties to validate")
    
    # Validate each property
    for prop in properties:
        # Check if property references a valid molecule
        if "molecule_id" not in prop or prop["molecule_id"] is None:
            report.add_issue(
                "missing_required_fields",
                f"Property {prop['id']} is missing required molecule_id",
                prop['id']
            )
        elif prop["molecule_id"] not in molecule_id_set:
            report.add_issue(
                "orphaned_properties",
                f"Property {prop['id']} references non-existent molecule {prop['molecule_id']}",
                prop['id']
            )
        
        # Check required fields
        if "property_type" not in prop or prop["property_type"] is None:
            report.add_issue(
                "missing_required_fields",
                f"Property {prop['id']} is missing required property_type",
                prop['id']
            )
        
        if "value" not in prop or prop["value"] is None:
            report.add_issue(
                "missing_required_fields",
                f"Property {prop['id']} is missing required value",
                prop['id']
            )
        
        # Validate numeric values based on property type
        if "property_type" in prop and "value" in prop and prop["value"] is not None:
            try:
                value = float(prop["value"])
                
                # Validate based on property type
                if prop["property_type"] == "LogP":
                    min_val, max_val = VALIDATION_PARAMS["logp_range"]
                    if not (min_val <= value <= max_val):
                        report.add_issue(
                            "invalid_numeric_values",
                            f"Property {prop['id']} (LogP) has value {value} outside valid range ({min_val}-{max_val})",
                            prop['id'],
                            severity="WARNING"
                        )
                
                elif prop["property_type"] == "TPSA":
                    min_val, max_val = VALIDATION_PARAMS["tpsa_range"]
                    if not (min_val <= value <= max_val):
                        report.add_issue(
                            "invalid_numeric_values",
                            f"Property {prop['id']} (TPSA) has value {value} outside valid range ({min_val}-{max_val})",
                            prop['id'],
                            severity="WARNING"
                        )
                
                elif prop["property_type"] == "H-Bond Donors":
                    min_val, max_val = VALIDATION_PARAMS["hbond_donors_range"]
                    if not (min_val <= value <= max_val):
                        report.add_issue(
                            "invalid_numeric_values",
                            f"Property {prop['id']} (H-Bond Donors) has value {value} outside valid range ({min_val}-{max_val})",
                            prop['id'],
                            severity="WARNING"
                        )
                
                elif prop["property_type"] == "H-Bond Acceptors":
                    min_val, max_val = VALIDATION_PARAMS["hbond_acceptors_range"]
                    if not (min_val <= value <= max_val):
                        report.add_issue(
                            "invalid_numeric_values",
                            f"Property {prop['id']} (H-Bond Acceptors) has value {value} outside valid range ({min_val}-{max_val})",
                            prop['id'],
                            severity="WARNING"
                        )
                
                elif prop["property_type"] == "Total Score":
                    min_val, max_val = VALIDATION_PARAMS["score_range"]
                    if not (min_val <= value <= max_val):
                        report.add_issue(
                            "invalid_numeric_values",
                            f"Property {prop['id']} (Total Score) has value {value} outside valid range ({min_val}-{max_val})",
                            prop['id'],
                            severity="WARNING"
                        )
                
            except (ValueError, TypeError):
                report.add_issue(
                    "invalid_numeric_values",
                    f"Property {prop['id']} has invalid numeric value: {prop['value']}",
                    prop['id']
                )
    
    logger.info("Molecular property validation completed")


def validate_property_consistency(report, project_id=None):
    """Validate consistency between molecules and their properties."""
    logger.info("Starting property consistency validation...")
    
    # Get all molecules with their properties
    molecule_query = supabase.table("molecule").select("id, name, inchikey")
    if project_id:
        molecule_query = molecule_query.eq("project_id", project_id)
    
    molecule_response = molecule_query.execute()
    
    if hasattr(molecule_response, 'error') and molecule_response.error:
        logger.error(f"Error fetching molecules: {molecule_response.error}")
        return
    
    molecules = molecule_response.data
    
    for molecule in molecules:
        # Get properties for this molecule
        property_response = supabase.table("molecular_property").select("*").eq("molecule_id", molecule["id"]).execute()
        
        if hasattr(property_response, 'error') and property_response.error:
            logger.error(f"Error fetching properties for molecule {molecule['id']}: {property_response.error}")
            continue
        
        properties = property_response.data
        
        # Check if molecule has expected properties
        expected_properties = ["LogP", "TPSA", "H-Bond Donors", "H-Bond Acceptors", "Total Score", "PubChem CID"]
        found_properties = [p["property_type"] for p in properties]
        
        for expected in expected_properties:
            if expected not in found_properties:
                report.add_issue(
                    "missing_required_fields",
                    f"Molecule {molecule['id']} ({molecule['name']}) is missing expected property: {expected}",
                    molecule['id'],
                    severity="WARNING"
                )
    
    logger.info("Property consistency validation completed")


def main():
    """Main function to run the validation script."""
    parser = argparse.ArgumentParser(description="CryoProtect Data Validation Script")
    parser.add_argument("--project-id", type=str, help="Validate data for a specific project ID")
    parser.add_argument("--report", type=str, default="validation_report.json", help="Path to save the validation report")
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect data validation...")
    
    # Authenticate with Supabase
    authenticated = authenticate_supabase()
    if not authenticated:
        logger.warning("Proceeding without authentication. Some validations may fail due to RLS policies.")
    
    # Create validation report
    report = ValidationReport()
    
    # Run validation checks
    validate_molecules(report, args.project_id)
    validate_molecular_properties(report, args.project_id)
    validate_property_consistency(report, args.project_id)
    
    # Generate and save report
    report.save_report(args.report)
    report.print_summary()
    
    logger.info(f"Validation completed. Found {report.summary['total_issues']} issues.")
    
    return report.summary["total_issues"] == 0


if __name__ == "__main__":
    try:
        success = main()
        exit(0 if success else 1)
    except Exception as e:
        logger.error(f"ERROR: Fatal error during validation: {str(e)}")
        exit(1)