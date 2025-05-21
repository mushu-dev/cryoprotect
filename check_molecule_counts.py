#!/usr/bin/env python3
"""
Check molecule counts and statistics in the Supabase database.

This script provides a basic overview of the molecules in the database,
including counts by source and property coverage.
"""

import os
import sys
import json
import logging
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def check_molecule_counts(project_id="tsdlmynydfuypiugmkev"):
    """Check molecule counts in the database."""
    try:
        from use_mcp_tool import execute_sql
        
        # Get total molecule count
        logger.info("Getting total molecule count...")
        total_count_query = "SELECT COUNT(*) AS count FROM molecules;"
        total_result = execute_sql(project_id, total_count_query)
        total_count = total_result[0]['count'] if total_result else 0
        
        # Get count by data source
        logger.info("Getting molecule counts by data source...")
        source_query = """
        SELECT data_source, COUNT(*) AS count 
        FROM molecules 
        GROUP BY data_source 
        ORDER BY count DESC;
        """
        source_result = execute_sql(project_id, source_query)
        
        # Get molecule counts with ChEMBL and PubChem IDs
        logger.info("Getting molecule counts with ChEMBL and PubChem IDs...")
        chembl_query = "SELECT COUNT(*) AS count FROM molecules WHERE chembl_id IS NOT NULL;"
        pubchem_query = "SELECT COUNT(*) AS count FROM molecules WHERE pubchem_cid IS NOT NULL;"
        both_query = "SELECT COUNT(*) AS count FROM molecules WHERE chembl_id IS NOT NULL AND pubchem_cid IS NOT NULL;"
        
        chembl_result = execute_sql(project_id, chembl_query)
        pubchem_result = execute_sql(project_id, pubchem_query)
        both_result = execute_sql(project_id, both_query)
        
        chembl_count = chembl_result[0]['count'] if chembl_result else 0
        pubchem_count = pubchem_result[0]['count'] if pubchem_result else 0
        both_count = both_result[0]['count'] if both_result else 0
        
        # Get molecular property counts
        logger.info("Getting molecular property counts...")
        property_count_query = """
        SELECT COUNT(DISTINCT molecule_id) AS molecule_count, 
               COUNT(*) AS property_count
        FROM molecular_properties;
        """
        property_result = execute_sql(project_id, property_count_query)
        
        molecules_with_properties = property_result[0]['molecule_count'] if property_result else 0
        total_properties = property_result[0]['property_count'] if property_result else 0
        
        # Get property type counts
        property_types_query = """
        SELECT property_type, COUNT(*) AS count
        FROM molecular_properties
        GROUP BY property_type
        ORDER BY count DESC;
        """
        property_types_result = execute_sql(project_id, property_types_query)
        
        # Get JSONB property counts
        jsonb_count_query = """
        SELECT COUNT(*) AS count
        FROM molecules
        WHERE properties IS NOT NULL AND properties::text != '{}';
        """
        jsonb_result = execute_sql(project_id, jsonb_count_query)
        jsonb_count = jsonb_result[0]['count'] if jsonb_result else 0
        
        # Compile results
        results = {
            "timestamp": datetime.now().isoformat(),
            "total_molecules": total_count,
            "by_source": source_result if source_result else [],
            "with_chembl_id": chembl_count,
            "with_pubchem_cid": pubchem_count,
            "with_both_ids": both_count,
            "with_molecular_properties": molecules_with_properties,
            "total_properties": total_properties,
            "property_types": property_types_result if property_types_result else [],
            "with_jsonb_properties": jsonb_count
        }
        
        # Print summary
        print("\n===== Molecule Database Summary =====")
        print(f"Total molecules: {total_count}")
        print("\nBy Source:")
        for source in (source_result or []):
            print(f"  {source.get('data_source', 'Unknown')}: {source.get('count', 0)}")
        
        print("\nIdentifiers:")
        print(f"  With ChEMBL ID: {chembl_count} ({chembl_count/total_count*100:.1f}% if total > 0 else 0)%)")
        print(f"  With PubChem CID: {pubchem_count} ({pubchem_count/total_count*100:.1f}% if total > 0 else 0)%)")
        print(f"  With both IDs: {both_count} ({both_count/total_count*100:.1f}% if total > 0 else 0)%)")
        
        print("\nProperties:")
        print(f"  Molecules with properties: {molecules_with_properties} ({molecules_with_properties/total_count*100:.1f}% if total > 0 else 0)%)")
        print(f"  Total properties: {total_properties}")
        print(f"  Molecules with JSONB properties: {jsonb_count} ({jsonb_count/total_count*100:.1f}% if total > 0 else 0)%)")
        
        print("\nProperty Types:")
        for prop in (property_types_result or []):
            print(f"  {prop.get('property_type', 'Unknown')}: {prop.get('count', 0)}")
            
        print("=====================================\n")
        
        # Save results to file
        os.makedirs("reports", exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = f"reports/molecule_counts_{timestamp}.json"
        
        with open(report_file, "w") as f:
            json.dump(results, f, indent=2)
            
        logger.info(f"Results saved to {report_file}")
        
        return results
        
    except ImportError:
        logger.error("Failed to import MCP tools. Make sure use_mcp_tool.py exists in the current directory.")
        return None
    except Exception as e:
        logger.error(f"Error checking molecule counts: {str(e)}")
        return None

if __name__ == "__main__":
    check_molecule_counts()