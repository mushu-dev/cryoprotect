#!/usr/bin/env python3
"""
Enhanced molecular property calculator for CryoProtect.

This module provides a comprehensive property calculation system for molecules,
ensuring all required properties are calculated and properly stored in both
the molecular_properties table and the JSONB properties field.
"""

import os
import sys
import uuid
import json
import logging
import time
from typing import Dict, List, Any, Optional, Tuple, Set
import psycopg2
from psycopg2.extras import RealDictCursor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Try importing RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, MolSurf, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    logger.warning("RDKit not available. Using fallback property calculations.")
    RDKIT_AVAILABLE = False

# Database connection parameters (from environment or .env file)
DB_PARAMS = {
    'host': os.getenv('SUPABASE_DB_HOST', 'aws-0-us-east-1.pooler.supabase.com'),
    'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
    'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
    'user': os.getenv('SUPABASE_DB_USER', 'postgres.tsdlmynydfuypiugmkev'),
    'password': os.getenv('SUPABASE_DB_PASSWORD', 'LDHt$rkaM&Gmf3X@LQ37'),
    'sslmode': 'require'
}

# Property definitions
PROPERTY_DEFINITIONS = {
    'LogP': {
        'description': 'Octanol-water partition coefficient',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Descriptors.MolLogP
    },
    'TPSA': {
        'description': 'Topological polar surface area',
        'data_type': 'numeric',
        'unit': 'Å²',
        'rdkit_func': MolSurf.TPSA
    },
    'Molecular Weight': {
        'description': 'The molecular weight of the compound',
        'data_type': 'numeric',
        'unit': 'g/mol',
        'rdkit_func': Descriptors.MolWt
    },
    'Heavy Atom Count': {
        'description': 'Number of non-hydrogen atoms',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Lipinski.HeavyAtomCount
    },
    'Hydrogen Bond Donor Count': {
        'description': 'Number of hydrogen bond donors',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Lipinski.NumHDonors
    },
    'Hydrogen Bond Acceptor Count': {
        'description': 'Number of hydrogen bond acceptors',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Lipinski.NumHAcceptors
    },
    'Rotatable Bond Count': {
        'description': 'Number of rotatable bonds',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Descriptors.NumRotatableBonds
    },
    'Ring Count': {
        'description': 'Number of rings',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Descriptors.RingCount
    },
    'Aromatic Ring Count': {
        'description': 'Number of aromatic rings',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': lambda mol: rdMolDescriptors.CalcNumAromaticRings(mol)
    }
}

class PropertyCalculator:
    """
    Enhanced property calculator for CryoProtect.
    
    This class provides methods to calculate and update molecular properties
    for molecules in the database, ensuring all required properties are
    present and correctly stored.
    """
    
    def __init__(self, db_params: Dict[str, Any]):
        """
        Initialize the property calculator.
        
        Args:
            db_params: Database connection parameters
        """
        self.db_params = db_params
        self.conn = None
        self.property_types = {}
        
        # Connect to database
        self._connect()
        
        # Load property types
        self._load_property_types()
    
    def _connect(self):
        """Connect to the database."""
        try:
            self.conn = psycopg2.connect(**self.db_params)
            logger.info("Connected to database")
        except Exception as e:
            logger.error(f"Database connection error: {str(e)}")
            sys.exit(1)
    
    def _load_property_types(self):
        """Load property types from the database."""
        try:
            query = "SELECT id, name, data_type FROM property_types;"
            with self.conn.cursor(cursor_factory=RealDictCursor) as cursor:
                cursor.execute(query)
                results = cursor.fetchall()
                
                for row in results:
                    self.property_types[row['name']] = {
                        'id': row['id'],
                        'data_type': row['data_type']
                    }
                
                logger.info(f"Loaded {len(self.property_types)} property types from database")
                
                # Check if all properties exist, create missing ones
                self._ensure_property_types_exist()
        except Exception as e:
            logger.error(f"Error loading property types: {str(e)}")
            raise
    
    def _ensure_property_types_exist(self):
        """Ensure all required property types exist in the database."""
        for prop_name, definition in PROPERTY_DEFINITIONS.items():
            if prop_name not in self.property_types:
                logger.info(f"Creating missing property type: {prop_name}")
                self._create_property_type(
                    prop_name, 
                    definition['description'],
                    definition['data_type'],
                    definition['unit']
                )
    
    def _create_property_type(self, name: str, description: str, data_type: str, unit: str):
        """
        Create a new property type in the database.
        
        Args:
            name: Property type name
            description: Property type description
            data_type: Data type (numeric, text, boolean)
            unit: Unit of measurement
        """
        try:
            query = """
            INSERT INTO property_types 
            (id, name, description, data_type, units, created_at, updated_at)
            VALUES (%s, %s, %s, %s, %s, NOW(), NOW())
            RETURNING id, name, data_type;
            """
            
            with self.conn.cursor(cursor_factory=RealDictCursor) as cursor:
                prop_id = str(uuid.uuid4())
                cursor.execute(query, (prop_id, name, description, data_type, unit))
                result = cursor.fetchone()
                self.conn.commit()
                
                # Add to property types dict
                self.property_types[name] = {
                    'id': result['id'],
                    'data_type': result['data_type']
                }
                
                logger.info(f"Created property type: {name} (ID: {prop_id})")
        except Exception as e:
            self.conn.rollback()
            logger.error(f"Error creating property type {name}: {str(e)}")
            raise
    
    def calculate_rdkit_properties(self, smiles: str, inchi: str = None) -> Optional[Dict[str, Any]]:
        """
        Calculate molecular properties using RDKit.
        
        Args:
            smiles: SMILES string
            inchi: Optional InChI string
            
        Returns:
            Dictionary of calculated properties or None if calculation failed
        """
        if not RDKIT_AVAILABLE:
            logger.warning("RDKit not available. Cannot calculate properties.")
            return None
        
        try:
            # Try to get a valid RDKit molecule
            mol = None
            
            # First try from SMILES
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
            
            # If SMILES failed, try from InChI
            if mol is None and inchi:
                mol = Chem.MolFromInchi(inchi)
            
            # If both failed, return None
            if mol is None:
                logger.warning(f"Could not convert to RDKit molecule: SMILES='{smiles}', InChI='{inchi}'")
                return None
            
            # Calculate properties
            properties = {}
            
            # Add standard molecular properties
            for prop_name, definition in PROPERTY_DEFINITIONS.items():
                try:
                    value = definition['rdkit_func'](mol)
                    properties[prop_name] = value
                except Exception as e:
                    logger.warning(f"Error calculating {prop_name}: {str(e)}")
            
            # Add calculated molecular formula
            try:
                properties['Molecular Formula'] = AllChem.CalcMolFormula(mol)
            except Exception as e:
                logger.warning(f"Error calculating molecular formula: {str(e)}")
            
            return properties
        except Exception as e:
            logger.error(f"Error calculating properties: {str(e)}")
            return None
    
    def get_molecules_with_missing_properties(self, batch_size: int = 100) -> List[Dict[str, Any]]:
        """
        Get molecules with missing properties.
        
        Args:
            batch_size: Number of molecules to fetch at once
            
        Returns:
            List of molecules with missing properties
        """
        # Get molecules with missing properties
        query = """
        WITH expected_properties AS (
            SELECT unnest(ARRAY[
                'LogP', 'TPSA', 'Molecular Weight', 'Heavy Atom Count',
                'Hydrogen Bond Donor Count', 'Hydrogen Bond Acceptor Count',
                'Rotatable Bond Count', 'Ring Count', 'Aromatic Ring Count'
            ]) AS property_type
        ),
        molecule_properties AS (
            SELECT m.id, m.name, mp.property_type
            FROM molecules m
            LEFT JOIN molecular_properties mp ON m.id = mp.molecule_id
            WHERE m.smiles IS NOT NULL OR m.inchi IS NOT NULL
        ),
        missing_properties AS (
            SELECT m.id, m.name,
                COUNT(DISTINCT mp.property_type) AS property_count,
                COUNT(DISTINCT ep.property_type) AS expected_count
            FROM molecules m
            CROSS JOIN expected_properties ep
            LEFT JOIN molecular_properties mp ON m.id = mp.molecule_id AND mp.property_type = ep.property_type
            WHERE m.smiles IS NOT NULL OR m.inchi IS NOT NULL
            GROUP BY m.id, m.name
            HAVING COUNT(DISTINCT mp.property_type) < COUNT(DISTINCT ep.property_type)
        )
        SELECT m.id, m.name, m.smiles, m.inchi, 
               m.properties, mp.property_count, mp.expected_count,
               COALESCE(m.properties = '{}' OR m.properties IS NULL, TRUE) AS needs_jsonb_update
        FROM missing_properties mp
        JOIN molecules m ON mp.id = m.id
        ORDER BY m.updated_at DESC
        LIMIT %s;
        """
        
        try:
            with self.conn.cursor(cursor_factory=RealDictCursor) as cursor:
                cursor.execute(query, (batch_size,))
                return cursor.fetchall()
        except Exception as e:
            logger.error(f"Error fetching molecules with missing properties: {str(e)}")
            return []
    
    def get_existing_properties(self, molecule_id: str) -> Dict[str, Any]:
        """
        Get existing properties for a molecule.
        
        Args:
            molecule_id: Molecule ID
            
        Returns:
            Dictionary mapping property types to values
        """
        query = """
        SELECT property_type, numeric_value, text_value, boolean_value
        FROM molecular_properties
        WHERE molecule_id = %s;
        """
        
        try:
            with self.conn.cursor(cursor_factory=RealDictCursor) as cursor:
                cursor.execute(query, (molecule_id,))
                results = cursor.fetchall()
                
                properties = {}
                for row in results:
                    prop_type = row['property_type']
                    
                    # Get the appropriate value based on data type
                    if row['numeric_value'] is not None:
                        properties[prop_type] = row['numeric_value']
                    elif row['text_value'] is not None:
                        properties[prop_type] = row['text_value']
                    elif row['boolean_value'] is not None:
                        properties[prop_type] = row['boolean_value']
                
                return properties
        except Exception as e:
            logger.error(f"Error fetching properties for molecule {molecule_id}: {str(e)}")
            return {}
    
    def update_molecule_properties(self, molecule_id: str, properties: Dict[str, Any], existing_properties: Dict[str, Any] = None) -> bool:
        """
        Update molecular properties for a molecule.
        
        Args:
            molecule_id: Molecule ID
            properties: Dictionary of property name to value
            existing_properties: Optional dictionary of existing properties
            
        Returns:
            True if update was successful, False otherwise
        """
        if not properties:
            logger.warning(f"No properties to update for molecule {molecule_id}")
            return False
        
        # Get existing properties if not provided
        if existing_properties is None:
            existing_properties = self.get_existing_properties(molecule_id)
        
        try:
            # Start transaction
            with self.conn:
                with self.conn.cursor() as cursor:
                    # Update molecular_properties table
                    for prop_name, value in properties.items():
                        # Skip if property already exists
                        if prop_name in existing_properties:
                            continue
                        
                        # Skip if property type not defined
                        if prop_name not in self.property_types:
                            logger.warning(f"Property type {prop_name} not found in database")
                            continue
                        
                        # Get property type details
                        prop_type_id = self.property_types[prop_name]['id']
                        data_type = self.property_types[prop_name]['data_type']
                        
                        # Determine which value field to use
                        numeric_value = None
                        text_value = None
                        boolean_value = None
                        
                        if data_type == 'numeric':
                            numeric_value = value
                        elif data_type == 'text':
                            text_value = str(value)
                        elif data_type == 'boolean':
                            boolean_value = bool(value)
                        
                        # Insert property
                        query = """
                        INSERT INTO molecular_properties 
                        (id, molecule_id, property_type_id, property_type, property_name,
                         numeric_value, text_value, boolean_value, 
                         source, created_at, updated_at)
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW());
                        """
                        
                        cursor.execute(
                            query, 
                            (
                                str(uuid.uuid4()),  # id
                                molecule_id,        # molecule_id
                                prop_type_id,       # property_type_id
                                prop_name,          # property_type
                                prop_name,          # property_name
                                numeric_value,      # numeric_value
                                text_value,         # text_value
                                boolean_value,      # boolean_value
                                'Calculated'        # source
                            )
                        )
                    
                    # Update the JSONB properties field
                    jsonb_properties = {}
                    
                    # Combine existing and new properties
                    for prop_name, value in {**existing_properties, **properties}.items():
                        jsonb_properties[prop_name] = value
                    
                    # Convert to JSONB
                    jsonb_str = json.dumps(jsonb_properties)
                    
                    # Update molecule
                    query = """
                    UPDATE molecules
                    SET properties = %s, updated_at = NOW()
                    WHERE id = %s;
                    """
                    
                    cursor.execute(query, (jsonb_str, molecule_id))
                    
                    # Commit transaction
                    self.conn.commit()
                    
                    logger.info(f"Updated properties for molecule {molecule_id}")
                    return True
        except Exception as e:
            self.conn.rollback()
            logger.error(f"Error updating properties for molecule {molecule_id}: {str(e)}")
            return False
    
    def update_all_missing_properties(self, batch_size: int = 100, limit: int = None) -> Dict[str, int]:
        """
        Update all molecules with missing properties.
        
        Args:
            batch_size: Number of molecules to process at once
            limit: Maximum number of molecules to process (None for all)
            
        Returns:
            Statistics about the update process
        """
        stats = {
            'total_molecules': 0,
            'successful_updates': 0,
            'failed_updates': 0,
            'skipped_no_structure': 0,
            'skipped_error': 0
        }
        
        processed_count = 0
        
        while True:
            # Check if we've reached the limit
            if limit is not None and processed_count >= limit:
                logger.info(f"Reached limit of {limit} molecules")
                break
            
            # Get batch of molecules with missing properties
            molecules = self.get_molecules_with_missing_properties(batch_size)
            
            if not molecules:
                logger.info("No more molecules with missing properties")
                break
            
            logger.info(f"Processing batch of {len(molecules)} molecules")
            stats['total_molecules'] += len(molecules)
            
            # Process each molecule
            for molecule in molecules:
                molecule_id = molecule['id']
                name = molecule['name']
                smiles = molecule['smiles']
                inchi = molecule['inchi']
                
                logger.info(f"Processing molecule: {name} (ID: {molecule_id})")
                
                # Skip if no structure
                if not smiles and not inchi:
                    logger.warning(f"Skipping molecule {name} (ID: {molecule_id}): No structure")
                    stats['skipped_no_structure'] += 1
                    continue
                
                try:
                    # Calculate properties
                    calculated_properties = self.calculate_rdkit_properties(smiles, inchi)
                    
                    if calculated_properties:
                        # Get existing properties
                        existing_properties = self.get_existing_properties(molecule_id)
                        
                        # Update properties
                        success = self.update_molecule_properties(
                            molecule_id, 
                            calculated_properties,
                            existing_properties
                        )
                        
                        if success:
                            stats['successful_updates'] += 1
                        else:
                            stats['failed_updates'] += 1
                    else:
                        logger.warning(f"Could not calculate properties for molecule {name} (ID: {molecule_id})")
                        stats['skipped_error'] += 1
                except Exception as e:
                    logger.error(f"Error processing molecule {name} (ID: {molecule_id}): {str(e)}")
                    stats['skipped_error'] += 1
                
                processed_count += 1
                
                # Check if we've reached the limit
                if limit is not None and processed_count >= limit:
                    break
                
                # Small delay to avoid overloading the database
                time.sleep(0.01)
        
        return stats

def main():
    """Main function."""
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Update molecular properties in the database')
    parser.add_argument('--batch-size', type=int, default=100, help='Batch size')
    parser.add_argument('--limit', type=int, default=None, help='Maximum number of molecules to process')
    parser.add_argument('--molecule-id', type=str, default=None, help='Process a specific molecule ID')
    args = parser.parse_args()
    
    # Initialize property calculator
    calculator = PropertyCalculator(DB_PARAMS)
    
    # Process specific molecule if requested
    if args.molecule_id:
        logger.info(f"Processing molecule ID: {args.molecule_id}")
        
        # Get molecule info
        query = "SELECT name, smiles, inchi FROM molecules WHERE id = %s;"
        with calculator.conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute(query, (args.molecule_id,))
            molecule = cursor.fetchone()
            
            if not molecule:
                logger.error(f"Molecule ID {args.molecule_id} not found")
                return
            
            name = molecule['name']
            smiles = molecule['smiles']
            inchi = molecule['inchi']
            
            logger.info(f"Molecule: {name}")
            
            # Calculate properties
            calculated_properties = calculator.calculate_rdkit_properties(smiles, inchi)
            
            if calculated_properties:
                # Get existing properties
                existing_properties = calculator.get_existing_properties(args.molecule_id)
                
                # Update properties
                success = calculator.update_molecule_properties(
                    args.molecule_id, 
                    calculated_properties,
                    existing_properties
                )
                
                if success:
                    logger.info(f"Successfully updated properties for molecule {name}")
                else:
                    logger.error(f"Failed to update properties for molecule {name}")
            else:
                logger.error(f"Could not calculate properties for molecule {name}")
    else:
        # Process all molecules with missing properties
        logger.info("Processing all molecules with missing properties")
        
        stats = calculator.update_all_missing_properties(
            batch_size=args.batch_size,
            limit=args.limit
        )
        
        # Print statistics
        logger.info("=== Property Update Statistics ===")
        logger.info(f"Total molecules processed: {stats['total_molecules']}")
        logger.info(f"Successful updates: {stats['successful_updates']}")
        logger.info(f"Failed updates: {stats['failed_updates']}")
        logger.info(f"Skipped (no structure): {stats['skipped_no_structure']}")
        logger.info(f"Skipped (error): {stats['skipped_error']}")
        logger.info("=================================")

if __name__ == "__main__":
    main()