#!/usr/bin/env python3
"""
Test script for validating ChEMBL data imports.

This script performs a comprehensive validation of ChEMBL data that has been
imported into the CryoProtect database. It checks for data quality, completeness,
and consistency, generating a detailed validation report.
"""

import os
import sys
import json
import argparse
import logging
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("chembl_validation.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

try:
    # Try to import database modules
    from dotenv import load_dotenv
    import psycopg2
    from psycopg2.extras import RealDictCursor
    
    # Optional imports for enhanced functionality
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski
        HAS_RDKIT = True
    except ImportError:
        logger.warning("RDKit not available. Some validation checks will be skipped.")
        HAS_RDKIT = False
        
    OFFLINE_MODE = False
except ImportError:
    logger.warning("Database modules not available. Running in offline simulation mode.")
    OFFLINE_MODE = True


def get_db_connection():
    """Get a connection to the database."""
    if OFFLINE_MODE:
        return None
        
    try:
        load_dotenv()
        
        # Get connection parameters from environment variables
        db_host = os.getenv("SUPABASE_DB_HOST")
        db_port = os.getenv("SUPABASE_DB_PORT", "5432")
        db_name = os.getenv("SUPABASE_DB_NAME")
        db_user = os.getenv("SUPABASE_DB_USER")
        db_password = os.getenv("SUPABASE_DB_PASSWORD")
        
        # Validate connection parameters
        for param_name, param_value in [
            ("SUPABASE_DB_HOST", db_host),
            ("SUPABASE_DB_NAME", db_name),
            ("SUPABASE_DB_USER", db_user),
            ("SUPABASE_DB_PASSWORD", db_password)
        ]:
            if not param_value:
                raise ValueError(f"Environment variable {param_name} is not set")
        
        # Connect to the database
        conn = psycopg2.connect(
            host=db_host,
            port=db_port,
            dbname=db_name,
            user=db_user,
            password=db_password
        )
        
        return conn
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        raise


def get_imported_chembl_compounds(conn) -> List[Dict[str, Any]]:
    """
    Get all compounds imported from ChEMBL.
    
    Args:
        conn: Database connection
        
    Returns:
        List of compounds with their properties
    """
    if OFFLINE_MODE:
        # Return simulated data for testing
        return [
            {"id": "1", "name": "Glycerol", "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N", "source": "ChEMBL ID: CHEMBL135"},
            {"id": "2", "name": "DMSO", "inchikey": "IAZDPXIOMUYVGZ-UHFFFAOYSA-N", "source": "ChEMBL ID: CHEMBL679"},
            {"id": "3", "name": "Trehalose", "inchikey": "OHCIJQIZZHCUEZ-SGMDXRHPSA-N", "source": "ChEMBL ID: CHEMBL449"}
        ]
    
    try:
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        cursor.execute("""
            SELECT m.id, m.name, m.smiles, m.inchi, m.inchikey, m.formula, 
                   m.molecular_weight, m.pubchem_cid, m.data_source
            FROM molecules m
            WHERE m.data_source LIKE 'ChEMBL ID:%'
            ORDER BY m.id;
        """)
        compounds = cursor.fetchall()
        
        # Get properties for each compound
        for compound in compounds:
            cursor.execute("""
                SELECT mp.id, mp.property_type_id, pt.name as property_name, 
                       mp.numeric_value, mp.text_value, mp.boolean_value,
                       mp.unit, mp.data_source
                FROM molecular_properties mp
                JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE mp.molecule_id = %s;
            """, (compound['id'],))
            compound['properties'] = cursor.fetchall()
            
        cursor.close()
        return compounds
    except Exception as e:
        logger.error(f"Error getting ChEMBL compounds: {str(e)}")
        if not OFFLINE_MODE:
            conn.rollback()
        raise


def validate_molecule_structure(compound: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate the molecular structure of a compound.
    
    Args:
        compound: Compound data
        
    Returns:
        Validation results
    """
    results = {
        "molecule_id": compound.get("id"),
        "name": compound.get("name"),
        "structure_valid": False,
        "has_smiles": bool(compound.get("smiles")),
        "has_inchi": bool(compound.get("inchi")),
        "has_inchikey": bool(compound.get("inchikey")),
        "has_formula": bool(compound.get("formula")),
        "has_molecular_weight": bool(compound.get("molecular_weight")),
        "issues": []
    }
    
    # Basic validation
    if not results["has_smiles"] and not results["has_inchi"]:
        results["issues"].append("Missing both SMILES and InChI")
    
    if not results["has_inchikey"]:
        results["issues"].append("Missing InChIKey")
    
    # RDKit validation
    if HAS_RDKIT and compound.get("smiles"):
        try:
            mol = Chem.MolFromSmiles(compound["smiles"])
            if mol:
                results["structure_valid"] = True
                
                # Check if formula matches
                if results["has_formula"]:
                    rdkit_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                    if rdkit_formula != compound["formula"]:
                        results["issues"].append(f"Formula mismatch: {compound['formula']} vs {rdkit_formula}")
                
                # Check if molecular weight matches
                if results["has_molecular_weight"]:
                    rdkit_mw = round(Descriptors.MolWt(mol), 2)
                    db_mw = round(float(compound["molecular_weight"]), 2)
                    if abs(rdkit_mw - db_mw) > 0.1:  # Allow small floating point differences
                        results["issues"].append(f"Molecular weight mismatch: {db_mw} vs {rdkit_mw}")
            else:
                results["issues"].append("Invalid SMILES structure")
        except Exception as e:
            results["issues"].append(f"Error validating structure: {str(e)}")
    
    return results


def validate_properties(compound: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate the properties of a compound.
    
    Args:
        compound: Compound data
        
    Returns:
        Validation results
    """
    results = {
        "molecule_id": compound.get("id"),
        "name": compound.get("name"),
        "property_count": len(compound.get("properties", [])),
        "has_required_properties": False,
        "property_issues": []
    }
    
    # Define required properties
    required_properties = [
        "molecular weight",
        "logp",
        "hydrogen bond acceptor count",
        "hydrogen bond donor count"
    ]
    
    # Check if required properties are present
    property_names = [p.get("property_name", "").lower() for p in compound.get("properties", [])]
    missing_properties = [p for p in required_properties if p not in property_names]
    
    if missing_properties:
        results["property_issues"].append(f"Missing required properties: {', '.join(missing_properties)}")
    else:
        results["has_required_properties"] = True
    
    # Check property values
    for prop in compound.get("properties", []):
        prop_name = prop.get("property_name", "").lower()
        
        # Check for NULL values
        if prop.get("numeric_value") is None and prop.get("text_value") is None and prop.get("boolean_value") is None:
            results["property_issues"].append(f"Property '{prop_name}' has NULL value")
            continue
        
        # Check numeric properties
        if prop_name == "molecular weight" and prop.get("numeric_value") is not None:
            if float(prop["numeric_value"]) <= 0:
                results["property_issues"].append(f"Invalid molecular weight: {prop['numeric_value']}")
        
        if prop_name == "logp" and prop.get("numeric_value") is not None:
            if float(prop["numeric_value"]) < -10 or float(prop["numeric_value"]) > 10:
                results["property_issues"].append(f"Suspicious LogP value: {prop['numeric_value']}")
    
    return results


def validate_cross_references(compound: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate cross-references to other databases.
    
    Args:
        compound: Compound data
        
    Returns:
        Validation results
    """
    results = {
        "molecule_id": compound.get("id"),
        "name": compound.get("name"),
        "has_pubchem_cid": bool(compound.get("pubchem_cid")),
        "has_chembl_id": False,
        "cross_ref_issues": []
    }
    
    # Extract ChEMBL ID from data_source
    data_source = compound.get("data_source", "")
    if "ChEMBL ID:" in data_source:
        chembl_id = data_source.split("ChEMBL ID:")[1].strip()
        results["has_chembl_id"] = bool(chembl_id)
        results["chembl_id"] = chembl_id
    else:
        results["cross_ref_issues"].append("ChEMBL ID not found in data_source")
    
    # Check if PubChem CID is present
    if not results["has_pubchem_cid"]:
        results["cross_ref_issues"].append("Missing PubChem CID")
    
    return results


def generate_validation_report(validation_results: List[Dict[str, Any]], report_file: str = None) -> Dict[str, Any]:
    """
    Generate a validation report from validation results.
    
    Args:
        validation_results: List of validation results
        report_file: Path to save the report to
        
    Returns:
        Summary statistics
    """
    # Compile statistics
    stats = {
        "total_compounds": len(validation_results),
        "valid_structures": 0,
        "missing_structures": 0,
        "missing_required_properties": 0,
        "missing_pubchem_cid": 0,
        "missing_chembl_id": 0,
        "total_structure_issues": 0,
        "total_property_issues": 0,
        "total_cross_ref_issues": 0,
        "compounds_with_issues": 0,
        "compounds_without_issues": 0
    }
    
    for result in validation_results:
        # Count valid structures
        if result.get("structure", {}).get("structure_valid", False):
            stats["valid_structures"] += 1
        
        # Count issues
        structure_issues = len(result.get("structure", {}).get("issues", []))
        property_issues = len(result.get("properties", {}).get("property_issues", []))
        cross_ref_issues = len(result.get("cross_references", {}).get("cross_ref_issues", []))
        
        if structure_issues + property_issues + cross_ref_issues > 0:
            stats["compounds_with_issues"] += 1
        else:
            stats["compounds_without_issues"] += 1
        
        stats["total_structure_issues"] += structure_issues
        stats["total_property_issues"] += property_issues
        stats["total_cross_ref_issues"] += cross_ref_issues
        
        # Count missing data
        if structure_issues > 0 and not result.get("structure", {}).get("has_smiles", False):
            stats["missing_structures"] += 1
        
        if not result.get("properties", {}).get("has_required_properties", False):
            stats["missing_required_properties"] += 1
        
        if not result.get("cross_references", {}).get("has_pubchem_cid", False):
            stats["missing_pubchem_cid"] += 1
        
        if not result.get("cross_references", {}).get("has_chembl_id", False):
            stats["missing_chembl_id"] += 1
    
    # Calculate percentages
    if stats["total_compounds"] > 0:
        stats["valid_structures_pct"] = round(stats["valid_structures"] / stats["total_compounds"] * 100, 2)
        stats["compounds_with_issues_pct"] = round(stats["compounds_with_issues"] / stats["total_compounds"] * 100, 2)
        stats["missing_properties_pct"] = round(stats["missing_required_properties"] / stats["total_compounds"] * 100, 2)
        stats["missing_pubchem_cid_pct"] = round(stats["missing_pubchem_cid"] / stats["total_compounds"] * 100, 2)
    
    # Prepare full report
    report = {
        "timestamp": datetime.now().isoformat(),
        "statistics": stats,
        "validation_results": validation_results
    }
    
    # Save report to file if specified
    if report_file:
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        # Also save a markdown summary
        md_file = report_file.replace('.json', '.md')
        with open(md_file, 'w') as f:
            f.write("# ChEMBL Import Validation Report\n\n")
            f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Summary Statistics\n\n")
            f.write(f"- **Total compounds:** {stats['total_compounds']}\n")
            f.write(f"- **Valid structures:** {stats['valid_structures']} ({stats.get('valid_structures_pct', 0)}%)\n")
            f.write(f"- **Compounds with issues:** {stats['compounds_with_issues']} ({stats.get('compounds_with_issues_pct', 0)}%)\n")
            f.write(f"- **Missing required properties:** {stats['missing_required_properties']} ({stats.get('missing_properties_pct', 0)}%)\n")
            f.write(f"- **Missing PubChem CID:** {stats['missing_pubchem_cid']} ({stats.get('missing_pubchem_cid_pct', 0)}%)\n\n")
            
            f.write("## Issue Breakdown\n\n")
            f.write(f"- **Structure issues:** {stats['total_structure_issues']}\n")
            f.write(f"- **Property issues:** {stats['total_property_issues']}\n")
            f.write(f"- **Cross-reference issues:** {stats['total_cross_ref_issues']}\n\n")
            
            f.write("## Recommendations\n\n")
            if stats['missing_pubchem_cid'] > 0:
                f.write("- Run the PubChem cross-reference resolution script to fill in missing CIDs\n")
            if stats['missing_required_properties'] > 0:
                f.write("- Use RDKit to calculate missing required properties\n")
            if stats['missing_structures'] > 0:
                f.write("- Review compounds with missing structural information\n")
            
            f.write("\n*See the accompanying JSON file for detailed validation results.*\n")
    
    return stats


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Validate ChEMBL imported data")
    parser.add_argument("--output", type=str, default=f"chembl_validation_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                        help="Output file for validation report")
    parser.add_argument("--offline", action="store_true", help="Run in offline simulation mode")
    args = parser.parse_args()
    
    global OFFLINE_MODE
    if args.offline:
        OFFLINE_MODE = True
        logger.info("Running in offline simulation mode")
    
    try:
        # Connect to database
        conn = get_db_connection() if not OFFLINE_MODE else None
        
        # Get imported compounds
        logger.info("Fetching ChEMBL compounds from database...")
        compounds = get_imported_chembl_compounds(conn)
        logger.info(f"Found {len(compounds)} compounds imported from ChEMBL")
        
        # Perform validation
        logger.info("Validating compounds...")
        validation_results = []
        
        for i, compound in enumerate(compounds):
            if i % 10 == 0:
                logger.info(f"Validating compound {i+1}/{len(compounds)}")
            
            result = {
                "molecule_id": compound.get("id"),
                "name": compound.get("name"),
                "inchikey": compound.get("inchikey")
            }
            
            # Validate structure
            result["structure"] = validate_molecule_structure(compound)
            
            # Validate properties
            result["properties"] = validate_properties(compound)
            
            # Validate cross-references
            result["cross_references"] = validate_cross_references(compound)
            
            validation_results.append(result)
        
        # Generate report
        logger.info("Generating validation report...")
        stats = generate_validation_report(validation_results, args.output)
        
        # Print summary
        logger.info(f"Validation complete. Results saved to {args.output}")
        logger.info(f"Total compounds: {stats['total_compounds']}")
        logger.info(f"Valid structures: {stats['valid_structures']} ({stats.get('valid_structures_pct', 0)}%)")
        logger.info(f"Compounds with issues: {stats['compounds_with_issues']} ({stats.get('compounds_with_issues_pct', 0)}%)")
        
        # Close database connection
        if conn and not OFFLINE_MODE:
            conn.close()
        
    except Exception as e:
        logger.error(f"Error in validation: {str(e)}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())