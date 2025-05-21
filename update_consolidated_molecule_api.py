#!/usr/bin/env python3
"""
CryoProtect Analyzer API - Update Consolidated Molecule API

This script updates the existing consolidated molecule API with the improved implementation
that supports consolidated_molecules table with the enhanced data model.

The script performs the following updates:
1. Updates consolidated_utils.py with the improved implementation from consolidated_utils_updated.py
2. Updates consolidated_molecule_resource.py with the enhanced implementation from consolidated_molecule_resource_updated.py
3. Creates a new consolidated_api_updated.py file to register the enhanced resources
4. Updates the API initialization in api/__init__.py to use the new resources
"""

import os
import sys
import shutil
import logging
import datetime
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(f'update_consolidated_api_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger(__name__)

# Setup project paths
API_DIR = Path(os.path.dirname(os.path.abspath(__file__))) / 'api'
BACKUP_DIR = Path(os.path.dirname(os.path.abspath(__file__))) / 'backups' / f'api_backup_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}'

def create_backup():
    """Create backups of files that will be modified."""
    logger.info("Creating backup of API files...")
    
    # Create backup directory if it doesn't exist
    if not BACKUP_DIR.exists():
        BACKUP_DIR.mkdir(parents=True)
    
    # Files to backup
    files_to_backup = [
        'consolidated_utils.py',
        'consolidated_molecule_resource.py',
        'consolidated_api.py',
        '__init__.py'
    ]
    
    for file in files_to_backup:
        source_path = API_DIR / file
        if source_path.exists():
            target_path = BACKUP_DIR / file
            shutil.copy2(source_path, target_path)
            logger.info(f"Backed up {file} to {target_path}")
        else:
            logger.warning(f"File {file} not found, skipping backup")
    
    logger.info(f"Backup complete. Files saved to {BACKUP_DIR}")
    return True

def update_consolidated_utils():
    """
    Update the consolidated_utils.py file with the new version.
    """
    logger.info("Updating consolidated_utils.py...")
    
    source_path = API_DIR / 'consolidated_utils_updated.py'
    target_path = API_DIR / 'consolidated_utils.py'
    
    if not source_path.exists():
        logger.error(f"Source file {source_path} not found")
        return False
    
    try:
        # Copy the updated file over the original
        shutil.copy2(source_path, target_path)
        logger.info(f"Successfully updated {target_path}")
        return True
    except Exception as e:
        logger.error(f"Error updating consolidated_utils.py: {str(e)}")
        return False

def update_consolidated_molecule_resource():
    """
    Update the consolidated_molecule_resource.py file with the new version.
    """
    logger.info("Updating consolidated_molecule_resource.py...")
    
    source_path = API_DIR / 'consolidated_molecule_resource_updated.py'
    target_path = API_DIR / 'consolidated_molecule_resource.py'
    
    if not source_path.exists():
        logger.error(f"Source file {source_path} not found")
        return False
    
    try:
        # Copy the updated file over the original
        shutil.copy2(source_path, target_path)
        logger.info(f"Successfully updated {target_path}")
        return True
    except Exception as e:
        logger.error(f"Error updating consolidated_molecule_resource.py: {str(e)}")
        return False

def create_consolidated_api_updated():
    """
    Create a new consolidated_api_updated.py file to register the enhanced resources.
    """
    logger.info("Creating consolidated_api_updated.py...")
    
    api_updated_content = """
\"\"\"
CryoProtect Analyzer API - Enhanced Consolidated API Resources Registration

This module registers the enhanced API resources for handling consolidated molecules.
It provides functions to register all consolidated molecule resources with
the Flask-RESTful API and the API documentation.
\"\"\"

from flask_restful import Api
from api.api_docs import register_resource

from api.consolidated_molecule_resource import (
    ConsolidatedMoleculeResource,
    ConsolidatedMoleculeBatchResource,
    PrimaryMoleculeResource,
    ConsolidatedMoleculesListResource,
    ConsolidatedMoleculeAuditResource,
    ConsolidatedMoleculePropertyMigrationResource,
    ConsolidatedMoleculeSearchResource
)

from api.differentiation_resources import (
    DifferentiationGroupListResource,
    DifferentiationGroupResource,
    MoleculeDifferentiationResource
)

def register_enhanced_consolidated_resources(api: Api):
    \"\"\"
    Register enhanced consolidated molecule API resources.
    
    Args:
        api: Flask-RESTful API instance
    \"\"\"
    # Register consolidated molecule resources
    api.add_resource(
        ConsolidatedMoleculeResource,
        '/consolidated/molecules/<string:molecule_id>',
        endpoint='consolidated_molecule_enhanced'
    )
    
    api.add_resource(
        ConsolidatedMoleculeBatchResource,
        '/consolidated/batch',
        endpoint='consolidated_batch_enhanced'
    )
    
    api.add_resource(
        PrimaryMoleculeResource,
        '/molecules/<string:molecule_id>/primary',
        endpoint='primary_molecule_enhanced'
    )
    
    api.add_resource(
        ConsolidatedMoleculesListResource,
        '/consolidated',
        endpoint='consolidated_list_enhanced'
    )
    
    # Register new enhanced resources
    api.add_resource(
        ConsolidatedMoleculeAuditResource,
        '/consolidated/molecules/<string:molecule_id>/audit',
        endpoint='consolidated_molecule_audit'
    )
    
    api.add_resource(
        ConsolidatedMoleculePropertyMigrationResource,
        '/consolidated/molecules/<string:molecule_id>/migrate-properties',
        endpoint='consolidated_molecule_property_migration'
    )
    
    api.add_resource(
        ConsolidatedMoleculeSearchResource,
        '/consolidated/search',
        endpoint='consolidated_molecule_search'
    )
    
    # Register differentiation resources
    api.add_resource(
        DifferentiationGroupListResource,
        '/differentiation/groups',
        endpoint='differentiation_groups_enhanced'
    )
    
    api.add_resource(
        DifferentiationGroupResource,
        '/differentiation/groups/<string:group_id>',
        endpoint='differentiation_group_enhanced'
    )
    
    api.add_resource(
        MoleculeDifferentiationResource,
        '/molecules/<string:molecule_id>/differentiation',
        endpoint='molecule_differentiation_enhanced'
    )

def register_enhanced_consolidated_docs(docs):
    \"\"\"
    Register enhanced consolidated molecule API resources for documentation.
    
    Args:
        docs: FlaskApiSpec instance
    \"\"\"
    # Register consolidated molecule resources
    register_resource(docs, ConsolidatedMoleculeResource, 'consolidated_molecule_enhanced')
    register_resource(docs, ConsolidatedMoleculeBatchResource, 'consolidated_batch_enhanced')
    register_resource(docs, PrimaryMoleculeResource, 'primary_molecule_enhanced')
    register_resource(docs, ConsolidatedMoleculesListResource, 'consolidated_list_enhanced')
    
    # Register new enhanced resources
    register_resource(docs, ConsolidatedMoleculeAuditResource, 'consolidated_molecule_audit')
    register_resource(docs, ConsolidatedMoleculePropertyMigrationResource, 'consolidated_molecule_property_migration')
    register_resource(docs, ConsolidatedMoleculeSearchResource, 'consolidated_molecule_search')
    
    # Register differentiation resources
    register_resource(docs, DifferentiationGroupListResource, 'differentiation_groups_enhanced')
    register_resource(docs, DifferentiationGroupResource, 'differentiation_group_enhanced')
    register_resource(docs, MoleculeDifferentiationResource, 'molecule_differentiation_enhanced')
"""
    
    target_path = API_DIR / 'consolidated_api_updated.py'
    
    try:
        with open(target_path, 'w') as f:
            f.write(api_updated_content.strip())
        
        logger.info(f"Successfully created {target_path}")
        return True
    except Exception as e:
        logger.error(f"Error creating consolidated_api_updated.py: {str(e)}")
        return False

def update_api_init():
    """
    Update the API initialization in api/__init__.py to use the enhanced resources.
    """
    logger.info("Updating api/__init__.py...")
    
    init_path = API_DIR / '__init__.py'
    
    if not init_path.exists():
        logger.error(f"Init file {init_path} not found")
        return False
    
    try:
        with open(init_path, 'r') as f:
            init_content = f.read()
        
        # Import the enhanced consolidated API
        import_statement = "# Import and register enhanced consolidated molecule resources\nfrom api.consolidated_api_updated import register_enhanced_consolidated_resources, register_enhanced_consolidated_docs"
        
        # Find where to insert the import
        import_insert_pos = init_content.find("# Import RBAC routes")
        if import_insert_pos == -1:
            logger.error("Could not find insertion point for import statement")
            return False
        
        # Insert the import statement
        init_content = init_content[:import_insert_pos] + import_statement + "\n\n" + init_content[import_insert_pos:]
        
        # Register the enhanced resources
        register_statement = "# Register enhanced consolidated molecule resources\nregister_enhanced_consolidated_resources(api)"
        
        # Find where to insert the registration
        register_insert_pos = init_content.find("# Register Lab Verification resources")
        if register_insert_pos == -1:
            logger.error("Could not find insertion point for registration statement")
            return False
        
        # Insert the registration statement
        init_content = init_content[:register_insert_pos] + register_statement + "\n\n" + init_content[register_insert_pos:]
        
        # Register the documentation
        docs_statement = "# Register enhanced consolidated molecule resources documentation\nregister_enhanced_consolidated_docs(docs)"
        
        # Find where to insert the documentation registration
        docs_insert_pos = init_content.find("# Register property explorer resources documentation")
        if docs_insert_pos == -1:
            logger.error("Could not find insertion point for documentation statement")
            return False
        
        # Insert the documentation statement
        init_content = init_content[:docs_insert_pos] + docs_statement + "\n\n" + init_content[docs_insert_pos:]
        
        # Write the updated content
        with open(init_path, 'w') as f:
            f.write(init_content)
        
        logger.info(f"Successfully updated {init_path}")
        return True
    except Exception as e:
        logger.error(f"Error updating api/__init__.py: {str(e)}")
        return False

def create_rdkit_wrapper_consolidated():
    """
    Create the rdkit_wrapper_consolidated.py file to handle consolidated molecules.
    """
    logger.info("Creating rdkit_wrapper_consolidated.py...")
    
    rdkit_wrapper_content = """
#!/usr/bin/env python3
\"\"\"
CryoProtect Analyzer - RDKit Wrapper for Consolidated Molecules

This module extends the RDKit wrapper to support consolidated molecules.
It ensures that property calculations are only performed on primary molecules,
and provides functions for migrating properties between consolidated molecules.
\"\"\"

import os
import sys
import logging
from typing import Dict, List, Optional, Union, Any, Tuple
import json

# Add the parent directory to the path to import the rdkit_wrapper
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from rdkit_wrapper import (
    calculate_properties,
    generate_visualization,
    search_substructure,
    search_similar_molecules,
    calculate_molecular_descriptors,
    check_rdkit_available
)

# Set up logging
logger = logging.getLogger(__name__)

def get_db_connection():
    \"\"\"
    Get a database connection.
    
    Returns:
        Database connection from the current context
    \"\"\"
    # Import here to avoid circular imports
    from database.adapter import get_connection
    return get_connection()

def get_primary_molecule_id(molecule_id: str, conn=None) -> str:
    \"\"\"
    Get the primary molecule ID for a given molecule.
    
    If the molecule is a secondary (consolidated) molecule, return the ID
    of its primary molecule. Otherwise, return the original ID.
    
    Args:
        molecule_id: The molecule ID to check
        conn: Optional database connection
        
    Returns:
        Primary molecule ID if consolidated, otherwise the original ID
    \"\"\"
    close_conn = False
    try:
        if conn is None:
            conn = get_db_connection()
            close_conn = True
        
        cursor = conn.cursor()
        
        # Query to get primary molecule
        query = \"\"\"
        SELECT 
            CASE 
                WHEN molecule_status = 'duplicate' AND primary_molecule_id IS NOT NULL 
                THEN primary_molecule_id
                ELSE id
            END as primary_id
        FROM consolidated_molecules
        WHERE id = %s
        \"\"\"
        
        cursor.execute(query, (molecule_id,))
        result = cursor.fetchone()
        
        return result[0] if result is not None else molecule_id
        
    except Exception as e:
        logger.error(f"Error getting primary molecule: {str(e)}")
        return molecule_id
    finally:
        if close_conn and conn:
            conn.close()

def calculate_properties_for_consolidated(molecule_id: str, properties: List[str] = None) -> Dict[str, Any]:
    \"\"\"
    Calculate properties for a molecule, handling consolidated molecules.
    
    If the molecule is a duplicate (consolidated) molecule, get the primary molecule
    and calculate properties for it instead.
    
    Args:
        molecule_id: The ID of the molecule
        properties: Optional list of properties to calculate
        
    Returns:
        Dictionary of calculated properties
    \"\"\"
    try:
        # Get the primary molecule ID
        primary_id = get_primary_molecule_id(molecule_id)
        
        # If the molecule is a duplicate, use the primary
        if primary_id != molecule_id:
            logger.info(f"Molecule {molecule_id} is a duplicate, using primary molecule {primary_id}")
            
            # Calculate properties for the primary molecule
            return calculate_properties(primary_id, properties)
        
        # Otherwise, calculate properties for the original molecule
        return calculate_properties(molecule_id, properties)
        
    except Exception as e:
        logger.error(f"Error calculating properties for consolidated molecule: {str(e)}")
        return {}

def migrate_properties(source_id: str, target_id: str, property_types: List[str] = None) -> Dict[str, Any]:
    \"\"\"
    Migrate properties from one molecule to another.
    
    This function is used when consolidating molecules to ensure that properties
    are preserved and migrated to the primary molecule.
    
    Args:
        source_id: The source molecule ID
        target_id: The target molecule ID
        property_types: Optional list of property types to migrate
        
    Returns:
        Dictionary with migration results
    \"\"\"
    conn = None
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # If no property types specified, get all property types
        if not property_types:
            query = \"\"\"
            SELECT DISTINCT pt.name
            FROM molecular_properties mp
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE mp.molecule_id = %s
            \"\"\"
            
            cursor.execute(query, (source_id,))
            results = cursor.fetchall()
            property_types = [row[0] for row in results]
        
        migrated_properties = []
        skipped_properties = []
        
        # For each property type, check if it exists in the target molecule
        for prop_type in property_types:
            # Get property type ID
            cursor.execute("SELECT id FROM property_types WHERE name = %s", (prop_type,))
            type_result = cursor.fetchone()
            if not type_result:
                skipped_properties.append({
                    'property_type': prop_type,
                    'reason': 'Property type not found'
                })
                continue
            
            property_type_id = type_result[0]
            
            # Check if property exists in source molecule
            cursor.execute(
                "SELECT id, property_value, calculation_method FROM molecular_properties WHERE molecule_id = %s AND property_type_id = %s",
                (source_id, property_type_id)
            )
            source_result = cursor.fetchone()
            if not source_result:
                skipped_properties.append({
                    'property_type': prop_type,
                    'reason': 'Property not found in source molecule'
                })
                continue
            
            source_property_id, property_value, calculation_method = source_result
            
            # Check if property exists in target molecule
            cursor.execute(
                "SELECT id FROM molecular_properties WHERE molecule_id = %s AND property_type_id = %s",
                (target_id, property_type_id)
            )
            target_result = cursor.fetchone()
            
            if target_result:
                # Property exists in target, update it
                target_property_id = target_result[0]
                cursor.execute(
                    \"\"\"
                    UPDATE molecular_properties 
                    SET property_value = %s, calculation_method = %s, updated_at = NOW()
                    WHERE id = %s
                    \"\"\"
                    ,
                    (property_value, calculation_method, target_property_id)
                )
            else:
                # Property doesn't exist in target, insert it
                cursor.execute(
                    \"\"\"
                    INSERT INTO molecular_properties (molecule_id, property_type_id, property_value, calculation_method, created_at, updated_at)
                    VALUES (%s, %s, %s, %s, NOW(), NOW())
                    \"\"\"
                    ,
                    (target_id, property_type_id, property_value, calculation_method)
                )
            
            migrated_properties.append({
                'property_type': prop_type,
                'value': property_value
            })
        
        # Commit the transaction
        conn.commit()
        
        return {
            'source_id': source_id,
            'target_id': target_id,
            'migrated_properties': migrated_properties,
            'skipped_properties': skipped_properties,
            'property_count': len(migrated_properties)
        }
        
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error migrating properties: {str(e)}")
        return {
            'source_id': source_id,
            'target_id': target_id,
            'error': str(e),
            'migrated_properties': [],
            'skipped_properties': [],
            'property_count': 0
        }
    finally:
        if conn:
            conn.close()

def generate_visualization_for_consolidated(molecule_id: str, options: Dict[str, Any] = None) -> Dict[str, Any]:
    \"\"\"
    Generate visualization for a molecule, handling consolidated molecules.
    
    If the molecule is a duplicate (consolidated) molecule, get the primary molecule
    and generate visualization for it instead.
    
    Args:
        molecule_id: The ID of the molecule
        options: Optional visualization options
        
    Returns:
        Dictionary with visualization data
    \"\"\"
    try:
        # Get the primary molecule ID
        primary_id = get_primary_molecule_id(molecule_id)
        
        # Generate visualization
        result = generate_visualization(primary_id, options)
        
        # If the molecule is a duplicate, add information about the consolidated status
        if primary_id != molecule_id:
            result['is_consolidated'] = True
            result['primary_molecule_id'] = primary_id
            result['consolidated_note'] = f"Visualization generated for the primary molecule {primary_id}"
        
        return result
        
    except Exception as e:
        logger.error(f"Error generating visualization for consolidated molecule: {str(e)}")
        return {
            'molecule_id': molecule_id,
            'error': str(e)
        }

# Other functions that need to be consolidated-molecule aware
def search_substructure_with_consolidated(query: str, options: Dict[str, Any] = None) -> Dict[str, Any]:
    \"\"\"
    Search for substructures, handling consolidated molecules.
    
    When returning search results, indicate whether molecules are consolidated
    and provide their primary molecule ID if applicable.
    
    Args:
        query: The substructure query
        options: Optional search options
        
    Returns:
        Dictionary with search results
    \"\"\"
    try:
        # Perform the substructure search
        results = search_substructure(query, options)
        
        if 'matches' in results:
            # Process the matches to add consolidated molecule information
            conn = get_db_connection()
            cursor = conn.cursor()
            
            for match in results['matches']:
                molecule_id = match.get('molecule_id')
                if molecule_id:
                    # Check if the molecule is consolidated
                    query = \"\"\"
                    SELECT 
                        molecule_status,
                        CASE 
                            WHEN molecule_status = 'duplicate' AND primary_molecule_id IS NOT NULL
                            THEN primary_molecule_id
                            ELSE NULL
                        END as primary_id
                    FROM consolidated_molecules
                    WHERE id = %s
                    \"\"\"
                    
                    cursor.execute(query, (molecule_id,))
                    result = cursor.fetchone()
                    
                    if result:
                        status, primary_id = result
                        match['molecule_status'] = status
                        
                        if primary_id:
                            match['primary_molecule_id'] = primary_id
                            match['is_consolidated'] = True
                        else:
                            match['is_consolidated'] = False
            
            conn.close()
        
        return results
        
    except Exception as e:
        logger.error(f"Error in substructure search with consolidated molecules: {str(e)}")
        return {
            'query': query,
            'error': str(e),
            'matches': []
        }

def search_similar_molecules_with_consolidated(query: str, options: Dict[str, Any] = None) -> Dict[str, Any]:
    \"\"\"
    Search for similar molecules, handling consolidated molecules.
    
    When returning search results, indicate whether molecules are consolidated
    and provide their primary molecule ID if applicable.
    
    Args:
        query: The similarity query
        options: Optional search options
        
    Returns:
        Dictionary with search results
    \"\"\"
    try:
        # Perform the similarity search
        results = search_similar_molecules(query, options)
        
        if 'matches' in results:
            # Process the matches to add consolidated molecule information
            conn = get_db_connection()
            cursor = conn.cursor()
            
            for match in results['matches']:
                molecule_id = match.get('molecule_id')
                if molecule_id:
                    # Check if the molecule is consolidated
                    query = \"\"\"
                    SELECT 
                        molecule_status,
                        CASE 
                            WHEN molecule_status = 'duplicate' AND primary_molecule_id IS NOT NULL
                            THEN primary_molecule_id
                            ELSE NULL
                        END as primary_id
                    FROM consolidated_molecules
                    WHERE id = %s
                    \"\"\"
                    
                    cursor.execute(query, (molecule_id,))
                    result = cursor.fetchone()
                    
                    if result:
                        status, primary_id = result
                        match['molecule_status'] = status
                        
                        if primary_id:
                            match['primary_molecule_id'] = primary_id
                            match['is_consolidated'] = True
                        else:
                            match['is_consolidated'] = False
            
            conn.close()
        
        return results
        
    except Exception as e:
        logger.error(f"Error in similarity search with consolidated molecules: {str(e)}")
        return {
            'query': query,
            'error': str(e),
            'matches': []
        }
"""
    
    target_path = Path(os.path.dirname(os.path.abspath(__file__))) / 'rdkit_wrapper_consolidated.py'
    
    try:
        with open(target_path, 'w') as f:
            f.write(rdkit_wrapper_content.strip())
        
        logger.info(f"Successfully created {target_path}")
        return True
    except Exception as e:
        logger.error(f"Error creating rdkit_wrapper_consolidated.py: {str(e)}")
        return False

def update_rdkit_utils():
    """
    Update rdkit_utils.py to import and use the consolidated RDKit wrapper.
    """
    logger.info("Updating rdkit_utils.py...")
    
    rdkit_utils_path = API_DIR / 'rdkit_utils.py'
    
    if not rdkit_utils_path.exists():
        logger.error(f"RDKit utils file {rdkit_utils_path} not found")
        return False
    
    try:
        with open(rdkit_utils_path, 'r') as f:
            rdkit_utils_content = f.read()
        
        # Add the import for rdkit_wrapper_consolidated
        import_statement = "# Import consolidated molecule RDKit wrapper\ntry:\n    from rdkit_wrapper_consolidated import (\n        get_primary_molecule_id,\n        calculate_properties_for_consolidated,\n        migrate_properties,\n        generate_visualization_for_consolidated,\n        search_substructure_with_consolidated,\n        search_similar_molecules_with_consolidated\n    )\n    CONSOLIDATED_WRAPPER_AVAILABLE = True\nexcept ImportError:\n    logger.warning(\"Consolidated molecule RDKit wrapper not available\")\n    CONSOLIDATED_WRAPPER_AVAILABLE = False"
        
        # Find where to insert the import
        import_insert_pos = rdkit_utils_content.find("# Set up logging")
        if import_insert_pos == -1:
            logger.error("Could not find insertion point for import statement")
            return False
        
        # Find where to actually insert (after the logging setup)
        line_end = rdkit_utils_content.find("\n", import_insert_pos)
        next_line_end = rdkit_utils_content.find("\n", line_end + 1)
        
        # Insert the import statement
        rdkit_utils_content = rdkit_utils_content[:next_line_end + 1] + "\n" + import_statement + "\n" + rdkit_utils_content[next_line_end + 1:]
        
        # Write the updated content
        with open(rdkit_utils_path, 'w') as f:
            f.write(rdkit_utils_content)
        
        logger.info(f"Successfully updated {rdkit_utils_path}")
        return True
    except Exception as e:
        logger.error(f"Error updating rdkit_utils.py: {str(e)}")
        return False

def main():
    """
    Main function to update the API.
    """
    logger.info("Starting update of Consolidated Molecule API...")
    
    # Create backup
    if not create_backup():
        logger.error("Failed to create backup, aborting")
        return
    
    # Update consolidated_utils.py
    if not update_consolidated_utils():
        logger.error("Failed to update consolidated_utils.py, aborting")
        return
    
    # Update consolidated_molecule_resource.py
    if not update_consolidated_molecule_resource():
        logger.error("Failed to update consolidated_molecule_resource.py, aborting")
        return
    
    # Create consolidated_api_updated.py
    if not create_consolidated_api_updated():
        logger.error("Failed to create consolidated_api_updated.py, aborting")
        return
    
    # Update api/__init__.py
    if not update_api_init():
        logger.error("Failed to update api/__init__.py, aborting")
        return
    
    # Create rdkit_wrapper_consolidated.py
    if not create_rdkit_wrapper_consolidated():
        logger.error("Failed to create rdkit_wrapper_consolidated.py, aborting")
        return
    
    # Update rdkit_utils.py
    if not update_rdkit_utils():
        logger.error("Failed to update rdkit_utils.py, aborting")
        return
    
    logger.info("Successfully updated Consolidated Molecule API!")
    logger.info("""
==========================================================================
Update Complete! The following changes have been made:

1. Updated consolidated_utils.py with the enhanced implementation
2. Updated consolidated_molecule_resource.py with the enhanced resources
3. Created consolidated_api_updated.py to register the enhanced resources
4. Updated api/__init__.py to use the enhanced resources
5. Created rdkit_wrapper_consolidated.py to handle consolidated molecules
6. Updated rdkit_utils.py to import and use the consolidated RDKit wrapper

New API endpoints are now available:
- /api/v1/consolidated/molecules/{molecule_id}
- /api/v1/consolidated/batch
- /api/v1/molecules/{molecule_id}/primary
- /api/v1/consolidated
- /api/v1/consolidated/molecules/{molecule_id}/audit
- /api/v1/consolidated/molecules/{molecule_id}/migrate-properties
- /api/v1/consolidated/search

The RDKit wrapper now handles consolidated molecules properly, ensuring
that property calculations are performed on primary molecules.

A backup of the original files has been created in the backups directory.
==========================================================================
""")

if __name__ == "__main__":
    main()