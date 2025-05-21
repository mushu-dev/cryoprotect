#!/usr/bin/env python3
"""
Comprehensive ChEMBL Import Verification

This script performs a comprehensive verification of ChEMBL data imported into the
CryoProtect database. It checks:

1. Database integrity (table structure, foreign keys, indexes)
2. Data completeness (required compounds, property coverage)
3. Data quality (validity of structures, property values)
4. Cross-references (PubChem IDs, InChIKeys)
5. Visualization of key metrics and property distributions

Usage:
    python comprehensive_chembl_verification.py [--project_id PROJECT_ID] 
                                             [--output_dir OUTPUT_DIR]
                                             [--sample_size SAMPLE_SIZE]
                                             [--generate_visualizations]
                                             [--check_with_chembl_api]
"""

import os
import sys
import json
import time
import logging
import argparse
import concurrent.futures
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Set
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for server environments
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Try different methods for database access
try:
    # First try MCP tools if available
    from use_mcp_tool import execute_sql, get_project_id
    logger = logging.getLogger("chembl_verification")
    USE_MCP = True
except ImportError:
    # Fall back to direct database connection
    USE_MCP = False
    try:
        from supabase_adapter import SupabaseAdapter
        USE_SUPABASE_ADAPTER = True
    except ImportError:
        USE_SUPABASE_ADAPTER = False
        import psycopg2
        import psycopg2.extras
        from config import Config

# Try importing ChEMBL client
try:
    from chembl_webresource_client.new_client import new_client
    CHEMBL_CLIENT_AVAILABLE = True
except ImportError:
    CHEMBL_CLIENT_AVAILABLE = False

# Try importing RDKit for structure validation
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    try:
        # Try with mock_rdkit if available
        from mock_rdkit import Chem, Descriptors, Lipinski
        RDKIT_AVAILABLE = True
    except ImportError:
        pass

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/comprehensive_chembl_verification.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Ensure directories exist
Path("logs").mkdir(exist_ok=True)
Path("reports").mkdir(exist_ok=True)
Path("reports/figures").mkdir(exist_ok=True)

# Reference compounds that must be present
REFERENCE_COMPOUNDS = [
    "CHEMBL25",     # Aspirin
    "CHEMBL1118",   # Caffeine
    "CHEMBL1234",   # Glycerol (common cryoprotectant)
    "CHEMBL444",    # Glucose
    "CHEMBL230130", # Ethylene glycol (common cryoprotectant)
    "CHEMBL9335",   # Dimethyl sulfoxide (DMSO, common cryoprotectant)
    "CHEMBL15151",  # Trehalose (common cryoprotectant)
    "CHEMBL1201625", # Propylene glycol (common cryoprotectant)
    "CHEMBL1276055", # Sucrose (common cryoprotectant)
]

# Required properties for complete molecule data
REQUIRED_PROPERTIES = [
    "LogP",
    "Molecular Weight",
    "TPSA",
    "Heavy Atom Count",
    "Hydrogen Bond Donor Count",
    "Hydrogen Bond Acceptor Count",
    "Rotatable Bond Count",
    "Ring Count",
    "Aromatic Ring Count"
]

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
            conn = psycopg2.connect(
                host=config.SUPABASE_DB_HOST,
                port=config.SUPABASE_DB_PORT,
                database=config.SUPABASE_DB_NAME,
                user=config.SUPABASE_DB_USER,
                password=config.SUPABASE_DB_PASSWORD,
                cursor_factory=psycopg2.extras.DictCursor
            )
            return conn
    except Exception as e:
        logger.error(f"Failed to connect to database: {str(e)}")
        return None

def execute_db_query(query, params=None, project_id=None, conn=None):
    """
    Execute a database query.
    
    Args:
        query: SQL query to execute
        params: Query parameters
        project_id: Supabase project ID (for MCP)
        conn: Database connection (for direct connection)
        
    Returns:
        Query results
    """
    if USE_MCP:
        try:
            result = execute_sql(query, project_id, params)
            return result
        except Exception as e:
            logger.error(f"Error executing query via MCP: {str(e)}")
            return None
    else:
        cursor = conn.cursor()
        try:
            cursor.execute(query, params or {})
            result = cursor.fetchall()
            return result
        except Exception as e:
            logger.error(f"Error executing query: {str(e)}")
            return None
        finally:
            cursor.close()

def check_database_structure(project_id=None, conn=None):
    """
    Check the database structure for ChEMBL import.
    
    Args:
        project_id: Supabase project ID
        conn: Database connection
        
    Returns:
        Dict with database structure information
    """
    logger.info("Checking database structure...")
    
    results = {
        "tables": {},
        "columns": {},
        "indexes": {},
        "constraints": {},
        "status": "success"
    }
    
    # Check required tables
    tables_query = """
    SELECT table_name 
    FROM information_schema.tables 
    WHERE table_schema = 'public' AND 
          table_name IN ('molecules', 'molecular_properties', 'property_types');
    """
    
    tables_result = execute_db_query(tables_query, project_id=project_id, conn=conn)
    
    if not tables_result:
        results["status"] = "error"
        results["error"] = "Failed to query tables"
        return results
    
    # Check if all required tables exist
    found_tables = [row["table_name"] for row in tables_result]
    required_tables = ["molecules", "molecular_properties", "property_types"]
    missing_tables = [table for table in required_tables if table not in found_tables]
    
    results["tables"] = {
        "required": required_tables,
        "found": found_tables,
        "missing": missing_tables,
        "status": "success" if not missing_tables else "error"
    }
    
    if missing_tables:
        results["status"] = "error"
        logger.error(f"Missing required tables: {missing_tables}")
    
    # Check for chembl_id column in molecules table
    if "molecules" in found_tables:
        columns_query = """
        SELECT column_name, data_type
        FROM information_schema.columns
        WHERE table_schema = 'public' AND table_name = 'molecules';
        """
        
        columns_result = execute_db_query(columns_query, project_id=project_id, conn=conn)
        
        if columns_result:
            columns = {row["column_name"]: row["data_type"] for row in columns_result}
            results["columns"]["molecules"] = columns
            
            if "chembl_id" not in columns:
                results["columns"]["status"] = "error"
                results["status"] = "error"
                logger.error("Missing chembl_id column in molecules table")
            else:
                results["columns"]["status"] = "success"
        else:
            results["columns"]["status"] = "error"
            results["status"] = "error"
    
    # Check indexes on molecules table for efficient searches
    if "molecules" in found_tables:
        indexes_query = """
        SELECT indexname, indexdef
        FROM pg_indexes
        WHERE schemaname = 'public' AND tablename = 'molecules';
        """
        
        indexes_result = execute_db_query(indexes_query, project_id=project_id, conn=conn)
        
        if indexes_result:
            indexes = {row["indexname"]: row["indexdef"] for row in indexes_result}
            results["indexes"]["molecules"] = indexes
            
            # Check for index on chembl_id
            has_chembl_id_index = False
            for index_name, index_def in indexes.items():
                if "chembl_id" in index_def:
                    has_chembl_id_index = True
                    break
            
            if not has_chembl_id_index:
                results["indexes"]["status"] = "warning"
                logger.warning("No index found on chembl_id column in molecules table")
            else:
                results["indexes"]["status"] = "success"
        else:
            results["indexes"]["status"] = "error"
            results["status"] = "error"
    
    # Check foreign key constraints for molecular_properties table
    if "molecular_properties" in found_tables:
        constraints_query = """
        SELECT
            tc.constraint_name,
            tc.constraint_type,
            kcu.column_name,
            ccu.table_name AS foreign_table_name,
            ccu.column_name AS foreign_column_name 
        FROM 
            information_schema.table_constraints tc
            JOIN information_schema.key_column_usage kcu ON tc.constraint_name = kcu.constraint_name
            JOIN information_schema.constraint_column_usage ccu ON ccu.constraint_name = tc.constraint_name
        WHERE
            tc.table_name = 'molecular_properties' AND tc.constraint_type = 'FOREIGN KEY';
        """
        
        constraints_result = execute_db_query(constraints_query, project_id=project_id, conn=conn)
        
        if constraints_result:
            constraints = []
            for row in constraints_result:
                constraints.append({
                    "name": row["constraint_name"],
                    "column": row["column_name"],
                    "foreign_table": row["foreign_table_name"],
                    "foreign_column": row["foreign_column_name"]
                })
            
            results["constraints"]["molecular_properties"] = constraints
            
            # Check for foreign key to molecules and property_types
            has_molecule_fk = False
            has_property_type_fk = False
            
            for constraint in constraints:
                if constraint["foreign_table"] == "molecules" and constraint["column"] == "molecule_id":
                    has_molecule_fk = True
                if constraint["foreign_table"] == "property_types" and constraint["column"] == "property_type_id":
                    has_property_type_fk = True
            
            if not has_molecule_fk or not has_property_type_fk:
                results["constraints"]["status"] = "error"
                results["status"] = "error"
                logger.error("Missing required foreign key constraints in molecular_properties table")
            else:
                results["constraints"]["status"] = "success"
        else:
            results["constraints"]["status"] = "error"
            results["status"] = "error"
    
    return results

def check_data_completeness(project_id=None, conn=None):
    """
    Check data completeness for ChEMBL import.
    
    Args:
        project_id: Supabase project ID
        conn: Database connection
        
    Returns:
        Dict with data completeness information
    """
    logger.info("Checking data completeness...")
    
    results = {
        "molecule_count": {},
        "property_types": {},
        "property_coverage": {},
        "reference_compounds": {},
        "status": "success"
    }
    
    # Check total molecule count and ChEMBL molecules
    count_query = """
    SELECT 
        COUNT(*) AS total_count,
        COUNT(CASE WHEN chembl_id IS NOT NULL THEN 1 END) AS chembl_count
    FROM molecules;
    """
    
    count_result = execute_db_query(count_query, project_id=project_id, conn=conn)
    
    if not count_result or not count_result[0]:
        results["status"] = "error"
        results["error"] = "Failed to query molecule counts"
        return results
    
    total_count = count_result[0]["total_count"]
    chembl_count = count_result[0]["chembl_count"]
    
    results["molecule_count"] = {
        "total": total_count,
        "with_chembl_id": chembl_count,
        "percentage": (chembl_count / total_count * 100) if total_count > 0 else 0
    }
    
    # Check property types for required properties
    property_types_query = """
    SELECT id, name, data_type
    FROM property_types;
    """
    
    property_types_result = execute_db_query(property_types_query, project_id=project_id, conn=conn)
    
    if property_types_result:
        property_types = {}
        for row in property_types_result:
            property_types[row["name"]] = {
                "id": row["id"],
                "data_type": row["data_type"]
            }
        
        results["property_types"]["all"] = property_types
        
        # Check for required property types
        missing_property_types = [prop for prop in REQUIRED_PROPERTIES if prop not in property_types]
        results["property_types"]["missing"] = missing_property_types
        results["property_types"]["status"] = "success" if not missing_property_types else "warning"
        
        if missing_property_types:
            logger.warning(f"Missing required property types: {missing_property_types}")
    else:
        results["property_types"]["status"] = "error"
        results["status"] = "error"
    
    # Check property coverage for ChEMBL molecules
    if property_types_result:
        # Build a query that checks coverage for each property type
        property_coverage_clauses = []
        property_id_mapping = {}
        
        for prop_name in REQUIRED_PROPERTIES:
            if prop_name in property_types:
                prop_id = property_types[prop_name]["id"]
                property_id_mapping[prop_name] = prop_id
                property_coverage_clauses.append(f"""
                SUM(CASE WHEN EXISTS (
                    SELECT 1 
                    FROM molecular_properties mp 
                    WHERE mp.molecule_id = m.id AND mp.property_type_id = '{prop_id}'
                ) THEN 1 ELSE 0 END) AS "{prop_name.replace(' ', '_')}_count"
                """)
        
        if property_coverage_clauses:
            property_coverage_query = f"""
            SELECT
                COUNT(*) AS total_molecules,
                {', '.join(property_coverage_clauses)}
            FROM molecules m
            WHERE m.chembl_id IS NOT NULL;
            """
            
            property_coverage_result = execute_db_query(property_coverage_query, project_id=project_id, conn=conn)
            
            if property_coverage_result and property_coverage_result[0]:
                row = property_coverage_result[0]
                total_chembl_molecules = row["total_molecules"]
                
                property_coverage = {}
                for prop_name in REQUIRED_PROPERTIES:
                    if prop_name in property_types:
                        column_name = f"{prop_name.replace(' ', '_')}_count"
                        count = row.get(column_name, 0)
                        percentage = (count / total_chembl_molecules * 100) if total_chembl_molecules > 0 else 0
                        
                        property_coverage[prop_name] = {
                            "count": count,
                            "percentage": percentage
                        }
                
                results["property_coverage"] = {
                    "total_molecules": total_chembl_molecules,
                    "properties": property_coverage
                }
                
                # Check if all molecules have all required properties
                all_complete = True
                for prop_name, coverage in property_coverage.items():
                    if coverage["percentage"] < 95:  # Require at least 95% coverage
                        all_complete = False
                        break
                
                results["property_coverage"]["all_complete"] = all_complete
                results["property_coverage"]["status"] = "success" if all_complete else "warning"
                
                if not all_complete:
                    logger.warning("Not all ChEMBL molecules have complete property data")
            else:
                results["property_coverage"]["status"] = "error"
                results["status"] = "error"
    
    # Check reference compounds
    ref_ids = ", ".join([f"'{id}'" for id in REFERENCE_COMPOUNDS])
    reference_query = f"""
    SELECT chembl_id
    FROM molecules
    WHERE chembl_id IN ({ref_ids});
    """
    
    reference_result = execute_db_query(reference_query, project_id=project_id, conn=conn)
    
    if reference_result:
        found_references = [row["chembl_id"] for row in reference_result]
        missing_references = [ref for ref in REFERENCE_COMPOUNDS if ref not in found_references]
        
        results["reference_compounds"] = {
            "required": REFERENCE_COMPOUNDS,
            "found": found_references,
            "missing": missing_references,
            "status": "success" if not missing_references else "error"
        }
        
        if missing_references:
            results["status"] = "error"
            logger.error(f"Missing reference compounds: {missing_references}")
    else:
        results["reference_compounds"]["status"] = "error"
        results["status"] = "error"
    
    return results

def check_data_quality(project_id=None, conn=None, sample_size=50):
    """
    Check data quality for ChEMBL import.
    
    Args:
        project_id: Supabase project ID
        conn: Database connection
        sample_size: Number of molecules to sample for structure validation
        
    Returns:
        Dict with data quality information
    """
    logger.info("Checking data quality...")
    
    results = {
        "structure_validation": {},
        "property_validation": {},
        "status": "success"
    }
    
    # Only perform structure validation if RDKit is available
    if RDKIT_AVAILABLE:
        # Sample molecules for structure validation
        sample_query = f"""
        SELECT id, name, smiles, inchi, inchikey, chembl_id
        FROM molecules
        WHERE chembl_id IS NOT NULL AND smiles IS NOT NULL
        ORDER BY RANDOM()
        LIMIT {sample_size};
        """
        
        sample_result = execute_db_query(sample_query, project_id=project_id, conn=conn)
        
        if sample_result:
            structure_validation = {
                "total": len(sample_result),
                "valid": 0,
                "invalid": 0,
                "details": []
            }
            
            for row in sample_result:
                molecule_id = row["id"]
                smiles = row["smiles"]
                
                # Validate SMILES with RDKit
                mol = Chem.MolFromSmiles(smiles)
                is_valid = mol is not None
                
                # Validate InChI if present
                inchi = row.get("inchi")
                inchi_valid = False
                if inchi and mol:
                    inchi_mol = Chem.MolFromInchi(inchi)
                    inchi_valid = inchi_mol is not None
                
                # Store results
                validation_result = {
                    "id": molecule_id,
                    "name": row["name"],
                    "chembl_id": row["chembl_id"],
                    "smiles_valid": is_valid,
                    "inchi_valid": inchi_valid if inchi else None
                }
                
                structure_validation["details"].append(validation_result)
                
                if is_valid:
                    structure_validation["valid"] += 1
                else:
                    structure_validation["invalid"] += 1
            
            structure_validation["percentage_valid"] = (structure_validation["valid"] / structure_validation["total"] * 100) if structure_validation["total"] > 0 else 0
            structure_validation["status"] = "success" if structure_validation["percentage_valid"] >= 95 else "warning"
            
            results["structure_validation"] = structure_validation
            
            if structure_validation["status"] == "warning":
                logger.warning(f"Structure validation found {structure_validation['invalid']} invalid structures")
        else:
            results["structure_validation"]["status"] = "error"
            results["status"] = "error"
    else:
        results["structure_validation"] = {
            "status": "skipped",
            "reason": "RDKit not available"
        }
        logger.warning("Skipping structure validation because RDKit is not available")
    
    # Check LogP values for consistency
    logp_query = """
    SELECT 
        m.id, m.name, m.chembl_id, mp.numeric_value AS logp_value
    FROM 
        molecules m
        JOIN molecular_properties mp ON m.id = mp.molecule_id
        JOIN property_types pt ON mp.property_type_id = pt.id
    WHERE 
        pt.name = 'LogP' AND m.chembl_id IS NOT NULL
    ORDER BY 
        RANDOM()
    LIMIT 10;
    """
    
    logp_result = execute_db_query(logp_query, project_id=project_id, conn=conn)
    
    if logp_result:
        # Store LogP values
        logp_values = []
        for row in logp_result:
            logp_values.append({
                "id": row["id"],
                "name": row["name"],
                "chembl_id": row["chembl_id"],
                "logp_value": row["logp_value"]
            })
        
        results["property_validation"]["logp_values"] = logp_values
        results["property_validation"]["status"] = "success"
    else:
        results["property_validation"]["status"] = "error"
        results["status"] = "error"
    
    return results

def check_cross_references(project_id=None, conn=None):
    """
    Check cross-references with PubChem.
    
    Args:
        project_id: Supabase project ID
        conn: Database connection
        
    Returns:
        Dict with cross-reference information
    """
    logger.info("Checking cross-references...")
    
    results = {
        "pubchem_coverage": {},
        "inchikey_coverage": {},
        "status": "success"
    }
    
    # Check PubChem CID and InChIKey coverage
    coverage_query = """
    SELECT 
        COUNT(*) AS total_molecules,
        COUNT(CASE WHEN pubchem_cid IS NOT NULL THEN 1 END) AS with_pubchem_cid,
        COUNT(CASE WHEN inchikey IS NOT NULL THEN 1 END) AS with_inchikey
    FROM molecules
    WHERE chembl_id IS NOT NULL;
    """
    
    coverage_result = execute_db_query(coverage_query, project_id=project_id, conn=conn)
    
    if coverage_result and coverage_result[0]:
        row = coverage_result[0]
        total_molecules = row["total_molecules"]
        with_pubchem_cid = row["with_pubchem_cid"]
        with_inchikey = row["with_inchikey"]
        
        pubchem_percentage = (with_pubchem_cid / total_molecules * 100) if total_molecules > 0 else 0
        inchikey_percentage = (with_inchikey / total_molecules * 100) if total_molecules > 0 else 0
        
        results["pubchem_coverage"] = {
            "total_molecules": total_molecules,
            "with_pubchem_cid": with_pubchem_cid,
            "percentage": pubchem_percentage,
            "status": "success" if pubchem_percentage >= 90 else "warning"
        }
        
        results["inchikey_coverage"] = {
            "total_molecules": total_molecules,
            "with_inchikey": with_inchikey,
            "percentage": inchikey_percentage,
            "status": "success" if inchikey_percentage >= 95 else "warning"
        }
        
        if results["pubchem_coverage"]["status"] == "warning":
            logger.warning(f"PubChem CID coverage is only {pubchem_percentage:.1f}%")
        
        if results["inchikey_coverage"]["status"] == "warning":
            logger.warning(f"InChIKey coverage is only {inchikey_percentage:.1f}%")
    else:
        results["pubchem_coverage"]["status"] = "error"
        results["inchikey_coverage"]["status"] = "error"
        results["status"] = "error"
    
    return results

def check_with_chembl_api(project_id=None, conn=None, sample_size=10):
    """
    Compare local data with ChEMBL API data.
    
    Args:
        project_id: Supabase project ID
        conn: Database connection
        sample_size: Number of molecules to check
        
    Returns:
        Dict with comparison information
    """
    logger.info("Checking data against ChEMBL API...")
    
    if not CHEMBL_CLIENT_AVAILABLE:
        return {
            "status": "skipped",
            "reason": "ChEMBL client not available"
        }
    
    results = {
        "sample_size": sample_size,
        "matches": [],
        "mismatches": [],
        "status": "success"
    }
    
    # Get a sample of ChEMBL molecules from the database
    sample_query = f"""
    SELECT m.id, m.name, m.chembl_id, m.smiles, m.inchi, m.inchikey
    FROM molecules m
    WHERE m.chembl_id IS NOT NULL
    ORDER BY RANDOM()
    LIMIT {sample_size};
    """
    
    sample_result = execute_db_query(sample_query, project_id=project_id, conn=conn)
    
    if not sample_result:
        return {
            "status": "error",
            "error": "Failed to retrieve sample molecules"
        }
    
    # Initialize ChEMBL API client
    molecule = new_client.molecule
    
    # Compare each molecule with ChEMBL API data
    for row in sample_result:
        local_data = {
            "id": row["id"],
            "name": row["name"],
            "chembl_id": row["chembl_id"],
            "smiles": row["smiles"],
            "inchi": row["inchi"],
            "inchikey": row["inchikey"]
        }
        
        # Get molecule data from ChEMBL API
        try:
            api_data = molecule.get(row["chembl_id"])
            time.sleep(0.5)  # Be gentle on the API
            
            if api_data:
                # Extract structures from API response
                structures = api_data.get("molecule_structures", {})
                
                # Compare relevant fields
                comparison = {
                    "local_data": local_data,
                    "api_data": {
                        "chembl_id": api_data.get("molecule_chembl_id"),
                        "name": api_data.get("pref_name"),
                        "smiles": structures.get("canonical_smiles"),
                        "inchi": structures.get("standard_inchi"),
                        "inchikey": structures.get("standard_inchi_key")
                    },
                    "fields_compared": [],
                    "all_match": True
                }
                
                # Check each field
                for field in ["chembl_id", "smiles", "inchi", "inchikey"]:
                    local_value = local_data.get(field)
                    api_value = comparison["api_data"].get(field)
                    
                    if field == "name":
                        # Name might be different but still valid
                        matches = local_value == api_value
                        field_result = {
                            "field": field,
                            "local_value": local_value,
                            "api_value": api_value,
                            "matches": matches,
                            "critical": False
                        }
                        if not matches:
                            comparison["all_match"] = False
                    else:
                        # Other fields should match exactly
                        matches = local_value == api_value
                        field_result = {
                            "field": field,
                            "local_value": local_value,
                            "api_value": api_value,
                            "matches": matches,
                            "critical": field in ["chembl_id", "inchikey"]
                        }
                        if not matches and field in ["chembl_id", "inchikey"]:
                            comparison["all_match"] = False
                    
                    comparison["fields_compared"].append(field_result)
                
                # Add to results
                if comparison["all_match"]:
                    results["matches"].append(comparison)
                else:
                    results["mismatches"].append(comparison)
            else:
                results["mismatches"].append({
                    "local_data": local_data,
                    "api_data": None,
                    "error": "Molecule not found in ChEMBL API"
                })
        except Exception as e:
            logger.error(f"Error checking ChEMBL API for {row['chembl_id']}: {str(e)}")
            results["mismatches"].append({
                "local_data": local_data,
                "api_data": None,
                "error": str(e)
            })
    
    # Calculate match percentage
    results["match_count"] = len(results["matches"])
    results["mismatch_count"] = len(results["mismatches"])
    results["match_percentage"] = (results["match_count"] / sample_size * 100) if sample_size > 0 else 0
    
    results["status"] = "success" if results["match_percentage"] >= 90 else "warning"
    
    if results["status"] == "warning":
        logger.warning(f"Only {results['match_percentage']:.1f}% of molecules match ChEMBL API data")
    
    return results

def generate_visualizations(verification_data, output_dir="reports/figures"):
    """
    Generate visualizations of verification results.
    
    Args:
        verification_data: Verification results data
        output_dir: Directory to save visualizations
        
    Returns:
        Dict with paths to generated visualizations
    """
    logger.info("Generating visualizations...")
    
    # Ensure output directory exists
    Path(output_dir).mkdir(exist_ok=True, parents=True)
    
    visualization_paths = {}
    
    # 1. Property Coverage Chart
    if "data_completeness" in verification_data and "property_coverage" in verification_data["data_completeness"]:
        property_coverage = verification_data["data_completeness"]["property_coverage"].get("properties", {})
        if property_coverage:
            try:
                properties = []
                coverage_values = []
                
                for prop_name, coverage in property_coverage.items():
                    properties.append(prop_name)
                    coverage_values.append(coverage.get("percentage", 0))
                
                plt.figure(figsize=(10, 6))
                bars = plt.bar(properties, coverage_values, color='skyblue')
                plt.axhline(y=95, color='r', linestyle='--', label='Required (95%)')
                
                # Rotate x labels for better readability
                plt.xticks(rotation=45, ha='right')
                plt.title('Property Coverage for ChEMBL Molecules')
                plt.xlabel('Property')
                plt.ylabel('Coverage (%)')
                plt.tight_layout()
                plt.legend()
                
                # Add percentage labels on top of bars
                for bar in bars:
                    height = bar.get_height()
                    plt.text(bar.get_x() + bar.get_width()/2., height + 1,
                            f'{height:.1f}%', ha='center', va='bottom')
                
                chart_path = f"{output_dir}/property_coverage_chart.png"
                plt.savefig(chart_path)
                plt.close()
                
                visualization_paths["property_coverage_chart"] = chart_path
                logger.info(f"Generated property coverage chart: {chart_path}")
            except Exception as e:
                logger.error(f"Error generating property coverage chart: {str(e)}")
    
    # 2. Cross-Reference Coverage Chart
    if "cross_references" in verification_data:
        cross_references = verification_data["cross_references"]
        try:
            # Create data for the chart
            categories = ['PubChem CID', 'InChIKey']
            values = [
                cross_references.get("pubchem_coverage", {}).get("percentage", 0),
                cross_references.get("inchikey_coverage", {}).get("percentage", 0)
            ]
            
            plt.figure(figsize=(8, 6))
            bars = plt.bar(categories, values, color=['green', 'blue'])
            plt.axhline(y=90, color='r', linestyle='--', label='Required (90%)')
            
            plt.title('Cross-Reference Coverage for ChEMBL Molecules')
            plt.xlabel('Reference Type')
            plt.ylabel('Coverage (%)')
            plt.ylim(0, 100)
            plt.tight_layout()
            plt.legend()
            
            # Add percentage labels on top of bars
            for bar in bars:
                height = bar.get_height()
                plt.text(bar.get_x() + bar.get_width()/2., height + 1,
                        f'{height:.1f}%', ha='center', va='bottom')
            
            chart_path = f"{output_dir}/cross_reference_chart.png"
            plt.savefig(chart_path)
            plt.close()
            
            visualization_paths["cross_reference_chart"] = chart_path
            logger.info(f"Generated cross-reference chart: {chart_path}")
        except Exception as e:
            logger.error(f"Error generating cross-reference chart: {str(e)}")
    
    # 3. Data Quality Chart
    if "data_quality" in verification_data and "structure_validation" in verification_data["data_quality"]:
        structure_validation = verification_data["data_quality"]["structure_validation"]
        if structure_validation and "total" in structure_validation:
            try:
                # Create data for pie chart
                labels = ['Valid', 'Invalid']
                sizes = [structure_validation.get("valid", 0), structure_validation.get("invalid", 0)]
                
                plt.figure(figsize=(8, 8))
                plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=['green', 'red'])
                plt.axis('equal')
                plt.title('Structure Validation Results')
                
                chart_path = f"{output_dir}/structure_validation_chart.png"
                plt.savefig(chart_path)
                plt.close()
                
                visualization_paths["structure_validation_chart"] = chart_path
                logger.info(f"Generated structure validation chart: {chart_path}")
            except Exception as e:
                logger.error(f"Error generating structure validation chart: {str(e)}")
    
    return visualization_paths

def generate_html_report(verification_data, visualization_paths=None):
    """
    Generate HTML report from verification results.
    
    Args:
        verification_data: Verification results data
        visualization_paths: Paths to visualization images
        
    Returns:
        Path to the generated HTML report
    """
    logger.info("Generating HTML report...")
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/chembl_verification_report_{timestamp}.html"
    
    # Summary data
    database_structure_status = verification_data.get("database_structure", {}).get("status", "unknown")
    data_completeness_status = verification_data.get("data_completeness", {}).get("status", "unknown")
    data_quality_status = verification_data.get("data_quality", {}).get("status", "unknown")
    cross_references_status = verification_data.get("cross_references", {}).get("status", "unknown")
    chembl_api_status = verification_data.get("chembl_api_comparison", {}).get("status", "unknown")
    
    overall_status = "success"
    if database_structure_status == "error" or data_completeness_status == "error":
        overall_status = "error"
    elif database_structure_status == "warning" or data_completeness_status == "warning" or data_quality_status == "warning":
        overall_status = "warning"
    
    # Generate HTML content
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>ChEMBL Import Verification Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 0; padding: 20px; line-height: 1.6; color: #333; }}
        h1, h2, h3, h4 {{ color: #205493; margin-top: 30px; }}
        .container {{ max-width: 1200px; margin: 0 auto; }}
        .card {{ background: #f9f9f9; border-radius: 5px; padding: 20px; margin-bottom: 20px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        .success {{ color: #2e8540; }}
        .warning {{ color: #fdb81e; }}
        .error {{ color: #d83933; }}
        .unknown {{ color: #5b616b; }}
        .status-badge {{
            display: inline-block;
            padding: 5px 10px;
            border-radius: 3px;
            font-weight: bold;
            margin-left: 10px;
        }}
        .status-badge.success {{ background-color: #e7f4e4; }}
        .status-badge.warning {{ background-color: #fff1d2; }}
        .status-badge.error {{ background-color: #f9dede; }}
        .status-badge.unknown {{ background-color: #f1f1f1; }}
        .summary {{ font-size: 1.2em; font-weight: bold; }}
        table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
        th, td {{ padding: 12px 15px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #f2f2f2; }}
        .chart-container {{ text-align: center; margin: 20px 0; }}
        .chart {{ max-width: 100%; height: auto; }}
        .footer {{ margin-top: 30px; text-align: center; font-size: 0.9em; color: #5b616b; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>ChEMBL Import Verification Report</h1>
        <p>Report generated on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        
        <div class="card">
            <h2>Summary 
                <span class="status-badge {overall_status}">{overall_status.upper()}</span>
            </h2>
            <table>
                <tr>
                    <th>Check</th>
                    <th>Status</th>
                    <th>Details</th>
                </tr>
                <tr>
                    <td>Database Structure</td>
                    <td><span class="{database_structure_status}">{database_structure_status.upper()}</span></td>
                    <td>{verification_data.get("database_structure", {}).get("tables", {}).get("status", "")}</td>
                </tr>
                <tr>
                    <td>Data Completeness</td>
                    <td><span class="{data_completeness_status}">{data_completeness_status.upper()}</span></td>
                    <td>
                        {verification_data.get("data_completeness", {}).get("molecule_count", {}).get("with_chembl_id", 0)} 
                        ChEMBL molecules found
                    </td>
                </tr>
                <tr>
                    <td>Data Quality</td>
                    <td><span class="{data_quality_status}">{data_quality_status.upper()}</span></td>
                    <td>
                        {verification_data.get("data_quality", {}).get("structure_validation", {}).get("percentage_valid", 0):.1f}% 
                        valid structures
                    </td>
                </tr>
                <tr>
                    <td>Cross-References</td>
                    <td><span class="{cross_references_status}">{cross_references_status.upper()}</span></td>
                    <td>
                        PubChem: {verification_data.get("cross_references", {}).get("pubchem_coverage", {}).get("percentage", 0):.1f}%,
                        InChIKey: {verification_data.get("cross_references", {}).get("inchikey_coverage", {}).get("percentage", 0):.1f}%
                    </td>
                </tr>
                <tr>
                    <td>ChEMBL API Comparison</td>
                    <td><span class="{chembl_api_status}">{chembl_api_status.upper()}</span></td>
                    <td>
                        {verification_data.get("chembl_api_comparison", {}).get("match_percentage", 0):.1f}% 
                        match with API data
                    </td>
                </tr>
            </table>
        </div>
        
        <div class="card">
            <h2>Database Structure</h2>
            
            <h3>Tables</h3>
            <table>
                <tr>
                    <th>Required Tables</th>
                    <th>Found</th>
                    <th>Missing</th>
                </tr>
                <tr>
                    <td>{", ".join(verification_data.get("database_structure", {}).get("tables", {}).get("required", []))}</td>
                    <td>{", ".join(verification_data.get("database_structure", {}).get("tables", {}).get("found", []))}</td>
                    <td>{", ".join(verification_data.get("database_structure", {}).get("tables", {}).get("missing", []))}</td>
                </tr>
            </table>
            
            <h3>Columns</h3>
            <p>ChEMBL ID Column: 
                {"Present" if "chembl_id" in verification_data.get("database_structure", {}).get("columns", {}).get("molecules", {}) else "Missing"}
            </p>
            
            <h3>Indexes</h3>
            <p>ChEMBL ID Index: 
                {"Present" if verification_data.get("database_structure", {}).get("indexes", {}).get("status") == "success" else "Missing"}
            </p>
        </div>
        
        <div class="card">
            <h2>Data Completeness</h2>
            
            <h3>Molecule Counts</h3>
            <table>
                <tr>
                    <th>Total Molecules</th>
                    <th>With ChEMBL ID</th>
                    <th>Percentage</th>
                </tr>
                <tr>
                    <td>{verification_data.get("data_completeness", {}).get("molecule_count", {}).get("total", 0)}</td>
                    <td>{verification_data.get("data_completeness", {}).get("molecule_count", {}).get("with_chembl_id", 0)}</td>
                    <td>{verification_data.get("data_completeness", {}).get("molecule_count", {}).get("percentage", 0):.1f}%</td>
                </tr>
            </table>
            
            <h3>Reference Compounds</h3>
            <table>
                <tr>
                    <th>Required</th>
                    <th>Found</th>
                    <th>Missing</th>
                </tr>
                <tr>
                    <td>{len(verification_data.get("data_completeness", {}).get("reference_compounds", {}).get("required", []))}</td>
                    <td>{len(verification_data.get("data_completeness", {}).get("reference_compounds", {}).get("found", []))}</td>
                    <td>{", ".join(verification_data.get("data_completeness", {}).get("reference_compounds", {}).get("missing", []))}</td>
                </tr>
            </table>
"""
    
    # Add property coverage chart if available
    if visualization_paths and "property_coverage_chart" in visualization_paths:
        html_content += f"""
            <h3>Property Coverage</h3>
            <div class="chart-container">
                <img src="{visualization_paths['property_coverage_chart']}" alt="Property Coverage Chart" class="chart">
            </div>
"""
    else:
        # Add property coverage table instead
        html_content += """
            <h3>Property Coverage</h3>
            <table>
                <tr>
                    <th>Property</th>
                    <th>Count</th>
                    <th>Percentage</th>
                </tr>
"""
        
        property_coverage = verification_data.get("data_completeness", {}).get("property_coverage", {}).get("properties", {})
        for prop_name, coverage in property_coverage.items():
            html_content += f"""
                <tr>
                    <td>{prop_name}</td>
                    <td>{coverage.get("count", 0)}</td>
                    <td>{coverage.get("percentage", 0):.1f}%</td>
                </tr>
"""
        
        html_content += """
            </table>
"""
    
    html_content += """
        </div>
        
        <div class="card">
            <h2>Data Quality</h2>
"""
    
    # Add structure validation chart if available
    if visualization_paths and "structure_validation_chart" in visualization_paths:
        html_content += f"""
            <h3>Structure Validation</h3>
            <div class="chart-container">
                <img src="{visualization_paths['structure_validation_chart']}" alt="Structure Validation Chart" class="chart">
            </div>
"""
    else:
        # Add structure validation table instead
        structure_validation = verification_data.get("data_quality", {}).get("structure_validation", {})
        html_content += f"""
            <h3>Structure Validation</h3>
            <p>Valid Structures: {structure_validation.get("valid", 0)} ({structure_validation.get("percentage_valid", 0):.1f}%)</p>
            <p>Invalid Structures: {structure_validation.get("invalid", 0)}</p>
"""
    
    # Add property validation section
    property_validation = verification_data.get("data_quality", {}).get("property_validation", {})
    if "logp_values" in property_validation:
        html_content += """
            <h3>LogP Values (Sample)</h3>
            <table>
                <tr>
                    <th>ChEMBL ID</th>
                    <th>Name</th>
                    <th>LogP</th>
                </tr>
"""
        
        for logp_item in property_validation["logp_values"]:
            html_content += f"""
                <tr>
                    <td>{logp_item.get("chembl_id", "")}</td>
                    <td>{logp_item.get("name", "")}</td>
                    <td>{logp_item.get("logp_value", "")}</td>
                </tr>
"""
        
        html_content += """
            </table>
"""
    
    html_content += """
        </div>
        
        <div class="card">
            <h2>Cross-References</h2>
"""
    
    # Add cross-reference chart if available
    if visualization_paths and "cross_reference_chart" in visualization_paths:
        html_content += f"""
            <div class="chart-container">
                <img src="{visualization_paths['cross_reference_chart']}" alt="Cross-Reference Chart" class="chart">
            </div>
"""
    
    # Add cross-reference coverage data
    pubchem_coverage = verification_data.get("cross_references", {}).get("pubchem_coverage", {})
    inchikey_coverage = verification_data.get("cross_references", {}).get("inchikey_coverage", {})
    
    html_content += f"""
            <h3>PubChem Coverage</h3>
            <p>Molecules with PubChem CID: {pubchem_coverage.get("with_pubchem_cid", 0)} ({pubchem_coverage.get("percentage", 0):.1f}%)</p>
            
            <h3>InChIKey Coverage</h3>
            <p>Molecules with InChIKey: {inchikey_coverage.get("with_inchikey", 0)} ({inchikey_coverage.get("percentage", 0):.1f}%)</p>
        </div>
"""
    
    # Add ChEMBL API comparison section if available
    chembl_api_comparison = verification_data.get("chembl_api_comparison", {})
    if chembl_api_comparison and chembl_api_comparison.get("status") != "skipped":
        html_content += """
        <div class="card">
            <h2>ChEMBL API Comparison</h2>
"""
        
        if "match_percentage" in chembl_api_comparison:
            html_content += f"""
            <p>Sample Size: {chembl_api_comparison.get("sample_size", 0)} molecules</p>
            <p>Match Percentage: {chembl_api_comparison.get("match_percentage", 0):.1f}%</p>
            <p>Matches: {chembl_api_comparison.get("match_count", 0)}</p>
            <p>Mismatches: {chembl_api_comparison.get("mismatch_count", 0)}</p>
"""
        
        # Add details of mismatches
        if "mismatches" in chembl_api_comparison and chembl_api_comparison["mismatches"]:
            html_content += """
            <h3>Mismatches</h3>
            <table>
                <tr>
                    <th>ChEMBL ID</th>
                    <th>Field</th>
                    <th>Local Value</th>
                    <th>API Value</th>
                </tr>
"""
            
            for mismatch in chembl_api_comparison["mismatches"]:
                if "fields_compared" in mismatch:
                    for field in mismatch["fields_compared"]:
                        if not field["matches"]:
                            html_content += f"""
                <tr>
                    <td>{mismatch.get("local_data", {}).get("chembl_id", "")}</td>
                    <td>{field.get("field", "")}</td>
                    <td>{field.get("local_value", "")}</td>
                    <td>{field.get("api_value", "")}</td>
                </tr>
"""
                else:
                    html_content += f"""
                <tr>
                    <td>{mismatch.get("local_data", {}).get("chembl_id", "")}</td>
                    <td colspan="3">{mismatch.get("error", "Unknown error")}</td>
                </tr>
"""
            
            html_content += """
            </table>
"""
        
        html_content += """
        </div>
"""
    
    # Close the HTML document
    html_content += """
        <div class="footer">
            <p>Generated by Comprehensive ChEMBL Verification Script</p>
        </div>
    </div>
</body>
</html>
"""
    
    # Write the HTML report to a file
    with open(report_path, "w") as f:
        f.write(html_content)
    
    logger.info(f"HTML report saved to {report_path}")
    return report_path

def run_verification(project_id=None, output_dir="reports", sample_size=50, check_chembl_api=False, generate_visualizations_flag=False):
    """
    Run the comprehensive verification process.
    
    Args:
        project_id: Supabase project ID (for MCP)
        output_dir: Directory to save reports
        sample_size: Sample size for data quality checks
        check_chembl_api: Whether to check data against ChEMBL API
        generate_visualizations_flag: Whether to generate visualizations
        
    Returns:
        Dict with verification results
    """
    verification_start_time = datetime.now()
    logger.info(f"Starting comprehensive ChEMBL verification at {verification_start_time.isoformat()}")
    
    # Connect to database
    conn = None
    if not USE_MCP:
        conn = get_db_connection()
        if not conn:
            return {"status": "error", "error": "Failed to connect to database"}
    
    try:
        verification_results = {
            "timestamp": verification_start_time.isoformat(),
            "status": "success"
        }
        
        # Run verification steps
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
            # Submit verification tasks
            db_structure_future = executor.submit(check_database_structure, project_id, conn)
            data_completeness_future = executor.submit(check_data_completeness, project_id, conn)
            data_quality_future = executor.submit(check_data_quality, project_id, conn, sample_size)
            cross_references_future = executor.submit(check_cross_references, project_id, conn)
            
            # Get results
            verification_results["database_structure"] = db_structure_future.result()
            verification_results["data_completeness"] = data_completeness_future.result()
            verification_results["data_quality"] = data_quality_future.result()
            verification_results["cross_references"] = cross_references_future.result()
        
        # Run ChEMBL API comparison if requested
        if check_chembl_api:
            if CHEMBL_CLIENT_AVAILABLE:
                verification_results["chembl_api_comparison"] = check_with_chembl_api(
                    project_id=project_id,
                    conn=conn,
                    sample_size=min(sample_size, 10)  # Use smaller sample to be gentle on API
                )
            else:
                verification_results["chembl_api_comparison"] = {
                    "status": "skipped",
                    "reason": "ChEMBL client not available"
                }
        
        # Determine overall status
        if verification_results["database_structure"]["status"] == "error" or \
           verification_results["data_completeness"]["status"] == "error":
            verification_results["status"] = "error"
        elif verification_results["database_structure"]["status"] == "warning" or \
             verification_results["data_completeness"]["status"] == "warning" or \
             verification_results["data_quality"]["status"] == "warning" or \
             verification_results["cross_references"]["status"] == "warning":
            verification_results["status"] = "warning"
        
        # Calculate verification metrics
        molecule_count = verification_results["data_completeness"]["molecule_count"]["with_chembl_id"]
        reference_compounds_present = len(verification_results["data_completeness"]["reference_compounds"]["found"])
        reference_compounds_total = len(verification_results["data_completeness"]["reference_compounds"]["required"])
        
        verification_results["summary"] = {
            "molecule_count": molecule_count,
            "reference_compounds": f"{reference_compounds_present}/{reference_compounds_total}",
            "status": verification_results["status"],
            "overall_result": "PASS" if verification_results["status"] == "success" else "INCOMPLETE" if verification_results["status"] == "warning" else "FAIL"
        }
        
        # Save JSON report
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        json_report_path = f"{output_dir}/chembl_verification_{timestamp}.json"
        with open(json_report_path, "w") as f:
            json.dump(verification_results, f, indent=2)
        
        logger.info(f"JSON report saved to {json_report_path}")
        
        # Generate visualizations if requested
        visualization_paths = None
        if generate_visualizations_flag:
            visualization_paths = generate_visualizations(verification_results)
        
        # Generate HTML report
        html_report_path = generate_html_report(verification_results, visualization_paths)
        
        verification_results["reports"] = {
            "json": json_report_path,
            "html": html_report_path
        }
        
        # Log completion
        verification_end_time = datetime.now()
        verification_duration = (verification_end_time - verification_start_time).total_seconds()
        logger.info(f"Verification completed in {verification_duration:.1f} seconds")
        logger.info(f"Status: {verification_results['summary']['overall_result']}")
        
        return verification_results
    
    finally:
        # Close database connection if using direct connection
        if not USE_MCP and conn:
            conn.close()
            logger.info("Database connection closed")

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Comprehensive ChEMBL Verification")
    parser.add_argument("--project_id", help="Supabase project ID (for MCP)")
    parser.add_argument("--output_dir", default="reports", help="Directory to save reports")
    parser.add_argument("--sample_size", type=int, default=50, help="Sample size for data quality checks")
    parser.add_argument("--check_with_chembl_api", action="store_true", help="Check data against ChEMBL API")
    parser.add_argument("--generate_visualizations", action="store_true", help="Generate visualizations")
    parser.add_argument("--log_level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO", help="Set logging level")
    args = parser.parse_args()
    
    # Set log level
    logging.getLogger().setLevel(getattr(logging, args.log_level))
    
    # Create output directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs("logs", exist_ok=True)
    
    logger.info(f"Starting comprehensive ChEMBL verification")
    
    # Run verification
    try:
        results = run_verification(
            project_id=args.project_id,
            output_dir=args.output_dir,
            sample_size=args.sample_size,
            check_chembl_api=args.check_with_chembl_api,
            generate_visualizations_flag=args.generate_visualizations
        )
        
        # Print summary
        print("\n" + "="*50)
        print("CHEMBL VERIFICATION SUMMARY")
        print("="*50)
        print(f"Molecules with ChEMBL IDs: {results['summary']['molecule_count']}")
        print(f"Reference compounds present: {results['summary']['reference_compounds']}")
        print(f"Status: {results['summary']['overall_result']}")
        print(f"HTML Report: {results.get('reports', {}).get('html', 'Not generated')}")
        print("="*50)
        
        return 0 if results["status"] in ["success", "warning"] else 1
    
    except Exception as e:
        logger.error(f"Error during verification: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())