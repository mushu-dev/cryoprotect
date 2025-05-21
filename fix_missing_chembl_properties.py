#!/usr/bin/env python3
"""
Fix Missing Properties in ChEMBL Molecules

This script identifies and fixes ChEMBL molecules with missing properties:
1. Calculates missing molecular properties using RDKit
2. Updates both the dedicated molecular_properties table and JSONB properties field
3. Generates a detailed report of fixed properties

Usage:
    python fix_missing_chembl_properties.py [--batch-size BATCH_SIZE] [--limit LIMIT]
                                           [--dry-run] [--verbose]
"""

import os
import sys
import json
import time
import uuid
import logging
import argparse
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Set
import psycopg2
from psycopg2.extras import RealDictCursor, Json

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/fix_chembl_properties.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Try importing RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, MolSurf
    RDKIT_AVAILABLE = True
except ImportError:
    logger.warning("RDKit not available. Trying to import mock_rdkit...")
    try:
        from mock_rdkit import Chem, Descriptors, Lipinski, MolSurf
        RDKIT_AVAILABLE = True
    except ImportError:
        logger.error("Neither RDKit nor mock_rdkit available. This script requires RDKit.")
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

# Required property names and their calculation functions
PROPERTY_CALCULATORS = {
    "LogP": lambda mol: Descriptors.MolLogP(mol),
    "TPSA": lambda mol: Descriptors.TPSA(mol),
    "Molecular Weight": lambda mol: Descriptors.MolWt(mol),
    "Heavy Atom Count": lambda mol: Descriptors.HeavyAtomCount(mol),
    "Hydrogen Bond Donor Count": lambda mol: Descriptors.NumHDonors(mol),
    "Hydrogen Bond Acceptor Count": lambda mol: Descriptors.NumHAcceptors(mol),
    "Rotatable Bond Count": lambda mol: Descriptors.NumRotatableBonds(mol),
    "Ring Count": lambda mol: Descriptors.RingCount(mol),
    "Aromatic Ring Count": lambda mol: Descriptors.NumAromaticRings(mol),
    "Fraction Sp3": lambda mol: Descriptors.FractionCSP3(mol),
    "Num Stereo Centers": lambda mol: Chem.CalcNumAtomStereoCenters(mol)
}

# Property units for database storage
PROPERTY_UNITS = {
    "LogP": "",
    "TPSA": "Å²",
    "Molecular Weight": "g/mol",
    "Heavy Atom Count": "",
    "Hydrogen Bond Donor Count": "",
    "Hydrogen Bond Acceptor Count": "",
    "Rotatable Bond Count": "",
    "Ring Count": "",
    "Aromatic Ring Count": "",
    "Fraction Sp3": "",
    "Num Stereo Centers": ""
}

class MolecularPropertyFixer:
    """Fix missing properties in ChEMBL molecules."""
    
    def __init__(self, db_params: Dict[str, Any], dry_run: bool = False):
        """
        Initialize the property fixer.
        
        Args:
            db_params: Database connection parameters
            dry_run: Whether to run in dry-run mode (no database changes)
        """
        self.db_params = db_params
        self.dry_run = dry_run
        self.conn = None
        self.property_types = {}
        
        # Connect to database
        self._connect()
        
        # Get property types
        self._get_property_types()
    
    def _connect(self):
        """Connect to the database."""
        try:
            self.conn = psycopg2.connect(**self.db_params)
            logger.info("Connected to database")
        except Exception as e:
            logger.error(f"Database connection error: {str(e)}")
            sys.exit(1)
    
    def _get_property_types(self):
        """Get property types from database."""
        try:
            with self.conn.cursor(cursor_factory=RealDictCursor) as cursor:
                cursor.execute("SELECT id, name FROM property_types")
                
                for row in cursor.fetchall():
                    self.property_types[row['name']] = row['id']
                
                logger.info(f"Loaded {len(self.property_types)} property types")
                
                # Check if all required properties are defined
                missing_property_types = []
                for prop_name in PROPERTY_CALCULATORS.keys():
                    if prop_name not in self.property_types:
                        missing_property_types.append(prop_name)
                
                if missing_property_types:
                    logger.warning(f"Missing property types: {', '.join(missing_property_types)}")
                    
                    # Create missing property types
                    if not self.dry_run:
                        for prop_name in missing_property_types:
                            self._create_property_type(prop_name)
        except Exception as e:
            logger.error(f"Error getting property types: {str(e)}")
            sys.exit(1)
    
    def _create_property_type(self, property_name: str):
        """
        Create a new property type.
        
        Args:
            property_name: Name of the property
        """
        try:
            with self.conn.cursor() as cursor:
                property_id = str(uuid.uuid4())
                unit = PROPERTY_UNITS.get(property_name, "")
                
                cursor.execute(
                    """
                    INSERT INTO property_types (id, name, data_type, unit)
                    VALUES (%s, %s, %s, %s)
                    RETURNING id
                    """,
                    (property_id, property_name, "numeric", unit)
                )
                
                self.conn.commit()
                self.property_types[property_name] = property_id
                
                logger.info(f"Created property type: {property_name} (ID: {property_id})")
        except Exception as e:
            self.conn.rollback()
            logger.error(f"Error creating property type {property_name}: {str(e)}")
    
    def get_molecules_with_missing_properties(self, batch_size: int = 50) -> List[Dict[str, Any]]:
        """
        Get ChEMBL molecules with missing properties.
        
        Args:
            batch_size: Number of molecules to fetch at once
            
        Returns:
            List of molecules with missing properties
        """
        # Step 1: Find molecules with missing jsonb properties
        query_jsonb = """
        SELECT m.id, m.name, m.chembl_id, m.smiles, m.properties
        FROM molecules m
        WHERE m.chembl_id IS NOT NULL
        AND (
            m.properties IS NULL
            OR m.properties = '{}'::jsonb
            OR NOT (m.properties ? 'LogP')
            OR NOT (m.properties ? 'Molecular Weight')
            OR NOT (m.properties ? 'TPSA')
        )
        LIMIT %s;
        """
        
        # Step 2: Find molecules with missing properties in the molecular_properties table
        query_properties = """
        SELECT m.id, m.name, m.chembl_id, m.smiles, m.properties
        FROM molecules m
        WHERE m.chembl_id IS NOT NULL
        AND NOT EXISTS (
            SELECT 1 FROM molecular_properties mp
            WHERE mp.molecule_id = m.id
            AND mp.property_type_id = %s
        )
        LIMIT %s;
        """
        
        try:
            result = []
            molecule_ids = set()
            
            # Check for missing JSONB properties
            with self.conn.cursor(cursor_factory=RealDictCursor) as cursor:
                cursor.execute(query_jsonb, (batch_size,))
                for row in cursor.fetchall():
                    if row['id'] not in molecule_ids:
                        molecule_ids.add(row['id'])
                        result.append(row)
            
            # For each important property, check for missing values in molecular_properties
            for property_name in ["LogP", "Molecular Weight", "TPSA"]:
                if property_name in self.property_types:
                    with self.conn.cursor(cursor_factory=RealDictCursor) as cursor:
                        remaining = batch_size - len(result)
                        if remaining <= 0:
                            break
                        
                        cursor.execute(query_properties, (self.property_types[property_name], remaining))
                        for row in cursor.fetchall():
                            if row['id'] not in molecule_ids:
                                molecule_ids.add(row['id'])
                                result.append(row)
            
            logger.info(f"Found {len(result)} molecules with missing properties")
            return result
        except Exception as e:
            logger.error(f"Error fetching molecules with missing properties: {str(e)}")
            return []
    
    def calculate_properties(self, molecule: Dict[str, Any]) -> Dict[str, float]:
        """
        Calculate properties for a molecule using RDKit.
        
        Args:
            molecule: Molecule data with SMILES
            
        Returns:
            Dictionary of calculated properties
        """
        smiles = molecule.get('smiles')
        if not smiles:
            logger.warning(f"Molecule {molecule.get('name', 'Unknown')} (ID: {molecule.get('id')}) has no SMILES, skipping")
            return {}
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                logger.warning(f"Could not parse SMILES for molecule {molecule.get('name', 'Unknown')} (ID: {molecule.get('id')})")
                return {}
            
            properties = {}
            for prop_name, calculator in PROPERTY_CALCULATORS.items():
                try:
                    value = calculator(mol)
                    properties[prop_name] = value
                except Exception as e:
                    logger.warning(f"Error calculating {prop_name} for molecule {molecule.get('name', 'Unknown')}: {str(e)}")
            
            return properties
        except Exception as e:
            logger.error(f"Error calculating properties for molecule {molecule.get('name', 'Unknown')}: {str(e)}")
            return {}
    
    def update_jsonb_properties(self, molecule_id: str, properties: Dict[str, float]) -> bool:
        """
        Update JSONB properties for a molecule.
        
        Args:
            molecule_id: Molecule ID
            properties: Dictionary of properties to update
            
        Returns:
            True if update was successful, False otherwise
        """
        if self.dry_run:
            logger.info(f"[DRY RUN] Would update JSONB properties for molecule {molecule_id}")
            return True
        
        try:
            with self.conn.cursor() as cursor:
                query = """
                UPDATE molecules
                SET properties = COALESCE(properties, '{}'::jsonb) || %s::jsonb,
                    updated_at = NOW()
                WHERE id = %s;
                """
                
                cursor.execute(query, (json.dumps(properties), molecule_id))
                self.conn.commit()
                
                return True
        except Exception as e:
            self.conn.rollback()
            logger.error(f"Error updating JSONB properties for molecule {molecule_id}: {str(e)}")
            return False
    
    def insert_molecular_properties(self, molecule_id: str, properties: Dict[str, float]) -> bool:
        """
        Insert properties into molecular_properties table.
        
        Args:
            molecule_id: Molecule ID
            properties: Dictionary of properties to insert
            
        Returns:
            True if insert was successful, False otherwise
        """
        if self.dry_run:
            logger.info(f"[DRY RUN] Would insert {len(properties)} properties for molecule {molecule_id}")
            return True
        
        try:
            with self.conn.cursor() as cursor:
                # Get existing properties to avoid duplicates
                cursor.execute(
                    """
                    SELECT property_type_id
                    FROM molecular_properties
                    WHERE molecule_id = %s
                    """,
                    (molecule_id,)
                )
                
                existing_property_types = set(row[0] for row in cursor.fetchall())
                
                # Insert each property if it doesn't already exist
                for prop_name, value in properties.items():
                    if prop_name not in self.property_types:
                        logger.warning(f"Property type '{prop_name}' not found in database, skipping")
                        continue
                    
                    property_type_id = self.property_types[prop_name]
                    
                    if property_type_id in existing_property_types:
                        # Update existing property
                        cursor.execute(
                            """
                            UPDATE molecular_properties
                            SET numeric_value = %s,
                                updated_at = NOW()
                            WHERE molecule_id = %s AND property_type_id = %s
                            """,
                            (value, molecule_id, property_type_id)
                        )
                    else:
                        # Insert new property
                        property_id = str(uuid.uuid4())
                        cursor.execute(
                            """
                            INSERT INTO molecular_properties (id, molecule_id, property_type_id, numeric_value)
                            VALUES (%s, %s, %s, %s)
                            """,
                            (property_id, molecule_id, property_type_id, value)
                        )
                
                self.conn.commit()
                return True
        except Exception as e:
            self.conn.rollback()
            logger.error(f"Error inserting properties for molecule {molecule_id}: {str(e)}")
            return False
    
    def fix_molecule_properties(self, molecule: Dict[str, Any]) -> Dict[str, Any]:
        """
        Fix properties for a molecule.
        
        Args:
            molecule: Molecule data with SMILES
            
        Returns:
            Dictionary with results of the fix operation
        """
        molecule_id = molecule['id']
        name = molecule.get('name', 'Unknown')
        chembl_id = molecule.get('chembl_id', 'Unknown')
        
        logger.info(f"Fixing properties for molecule: {name} (ChEMBL ID: {chembl_id})")
        
        result = {
            'molecule_id': molecule_id,
            'name': name,
            'chembl_id': chembl_id,
            'properties_calculated': 0,
            'jsonb_updated': False,
            'properties_inserted': 0,
            'success': False,
            'calculated_properties': {}
        }
        
        try:
            # Calculate properties
            calculated_properties = self.calculate_properties(molecule)
            
            if not calculated_properties:
                logger.warning(f"No properties calculated for molecule {name}")
                return result
            
            result['properties_calculated'] = len(calculated_properties)
            result['calculated_properties'] = calculated_properties
            
            # Update JSONB properties
            jsonb_success = self.update_jsonb_properties(molecule_id, calculated_properties)
            result['jsonb_updated'] = jsonb_success
            
            # Insert into molecular_properties table
            properties_success = self.insert_molecular_properties(molecule_id, calculated_properties)
            result['properties_inserted'] = len(calculated_properties) if properties_success else 0
            
            result['success'] = jsonb_success and properties_success
            
            if result['success']:
                logger.info(f"Successfully fixed properties for molecule {name}")
            else:
                logger.error(f"Failed to fix properties for molecule {name}")
            
            return result
        except Exception as e:
            logger.error(f"Error fixing properties for molecule {name}: {str(e)}")
            return result
    
    def fix_all_missing_properties(self, batch_size: int = 50, limit: int = None) -> Dict[str, Any]:
        """
        Fix all molecules with missing properties.
        
        Args:
            batch_size: Number of molecules to process at once
            limit: Maximum number of molecules to process (None for all)
            
        Returns:
            Statistics about the fix operation
        """
        stats = {
            'total_molecules': 0,
            'molecules_with_missing_properties': 0,
            'successfully_fixed': 0,
            'failed_to_fix': 0,
            'properties_fixed': 0,
            'start_time': datetime.now().isoformat(),
            'end_time': None,
            'duration_seconds': 0,
            'detailed_results': []
        }
        
        processed_count = 0
        
        while True:
            # Check if we've reached the limit
            if limit is not None and processed_count >= limit:
                logger.info(f"Reached limit of {limit} molecules")
                break
            
            # Get batch of molecules with missing properties
            molecules = self.get_molecules_with_missing_properties(
                batch_size=min(batch_size, limit - processed_count if limit is not None else batch_size)
            )
            
            if not molecules:
                logger.info("No more molecules with missing properties")
                break
            
            logger.info(f"Processing batch of {len(molecules)} molecules")
            stats['molecules_with_missing_properties'] += len(molecules)
            
            # Process each molecule
            for molecule in molecules:
                result = self.fix_molecule_properties(molecule)
                stats['total_molecules'] += 1
                
                if result['success']:
                    stats['successfully_fixed'] += 1
                    stats['properties_fixed'] += result['properties_calculated']
                else:
                    stats['failed_to_fix'] += 1
                
                stats['detailed_results'].append(result)
                
                processed_count += 1
                
                # Check if we've reached the limit
                if limit is not None and processed_count >= limit:
                    break
        
        # Calculate duration
        end_time = datetime.now()
        stats['end_time'] = end_time.isoformat()
        stats['duration_seconds'] = (end_time - datetime.fromisoformat(stats['start_time'])).total_seconds()
        
        return stats

def main():
    """Main function."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Fix missing properties in ChEMBL molecules')
    parser.add_argument('--batch-size', type=int, default=50, help='Batch size')
    parser.add_argument('--limit', type=int, default=None, help='Maximum number of molecules to process')
    parser.add_argument('--dry-run', action='store_true', help='Run in dry-run mode (no database changes)')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--output', type=str, default='reports/property_fix_report.json', help='Output report file')
    args = parser.parse_args()
    
    # Set log level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    if not RDKIT_AVAILABLE:
        logger.error("RDKit is required for this script")
        return 1
    
    # Create output directory
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Initialize property fixer
    fixer = MolecularPropertyFixer(DB_PARAMS, dry_run=args.dry_run)
    
    # Fix all missing properties
    logger.info(f"Starting property fix operation (batch size: {args.batch_size}, limit: {args.limit or 'none'}, dry run: {args.dry_run})")
    
    start_time = time.time()
    stats = fixer.fix_all_missing_properties(batch_size=args.batch_size, limit=args.limit)
    elapsed_time = time.time() - start_time
    
    # Generate timestamp for report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = args.output.replace('.json', f'_{timestamp}.json')
    
    # Save report
    with open(report_path, 'w') as f:
        json.dump(stats, f, indent=2)
    
    # Print summary
    logger.info("=== Property Fix Operation Complete ===")
    logger.info(f"Total molecules checked: {stats['total_molecules']}")
    logger.info(f"Molecules with missing properties: {stats['molecules_with_missing_properties']}")
    logger.info(f"Successfully fixed: {stats['successfully_fixed']}")
    logger.info(f"Failed to fix: {stats['failed_to_fix']}")
    logger.info(f"Total properties fixed: {stats['properties_fixed']}")
    logger.info(f"Duration: {elapsed_time:.2f} seconds")
    logger.info(f"Report saved to: {report_path}")
    logger.info("========================================")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())