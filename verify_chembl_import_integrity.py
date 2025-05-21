#!/usr/bin/env python3
"""
ChEMBL Import Integrity Verification Tool

This script performs a comprehensive verification of ChEMBL data imports including:
1. Data integrity checks (missing fields, required properties)
2. Cross-reference validation with PubChem
3. Structure validation with RDKit (when available)
4. Statistical analysis of import completeness

Usage:
    python verify_chembl_import_integrity.py [--config CONFIG] [--output OUTPUT]

Options:
    --config CONFIG   Path to configuration file (default: config.py)
    --output OUTPUT   Path to output report directory (default: ./reports)
"""

import argparse
import json
import os
import sys
import time
from datetime import datetime
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional, Union

try:
    import psycopg2
    from psycopg2.extras import DictCursor
except ImportError:
    print("psycopg2 is required. Install it with 'pip install psycopg2-binary'")
    sys.exit(1)

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    print("Warning: pandas not available. Statistical analysis will be limited.")

# Try to import RDKit for structure validation
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Structure validation will be skipped.")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('chembl_verification')

# Constants
REQUIRED_PROPERTIES = [
    'molecular_weight', 
    'xlogp', 
    'h_bond_donor_count', 
    'h_bond_acceptor_count',
    'rotatable_bond_count',
    'exact_mass',
    'topological_polar_surface_area'
]

ESSENTIAL_FIELDS = [
    'chembl_id',
    'name',
    'smiles',
    'inchi',
    'inchi_key'
]

class ChEMBLVerifier:
    """Verifies ChEMBL import data integrity and cross-references."""
    
    def __init__(self, config_path: str = 'config.py', output_dir: str = './reports'):
        """Initialize the verifier with configuration."""
        self.config = self._load_config(config_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Connection setup
        self.conn = None
        self.conn_params = self._get_connection_params()
        
        # Results containers
        self.results = {
            'general': {
                'timestamp': datetime.now().isoformat(),
                'total_compounds': 0,
                'verified_compounds': 0,
                'missing_properties': 0,
                'missing_pubchem_refs': 0,
                'invalid_structures': 0
            },
            'property_stats': {},
            'pubchem_stats': {},
            'structure_stats': {},
            'specific_issues': []
        }
    
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration from a Python file."""
        config = {}
        try:
            sys.path.insert(0, os.path.dirname(os.path.abspath(config_path)))
            config_module = __import__(os.path.basename(config_path).replace('.py', ''))
            
            # Extract relevant database connection attributes
            for attr in ['SUPABASE_URL', 'SUPABASE_KEY', 'SUPABASE_DB_HOST', 
                         'SUPABASE_DB_PORT', 'SUPABASE_DB_NAME', 
                         'SUPABASE_DB_USER', 'SUPABASE_DB_PASSWORD',
                         'DB_HOST', 'DB_PORT', 'DB_NAME', 'DB_USER', 'DB_PASSWORD']:
                if hasattr(config_module, attr):
                    config[attr] = getattr(config_module, attr)
            
            return config
        except ImportError as e:
            logger.error(f"Failed to import config: {e}")
            sys.exit(1)
    
    def _get_connection_params(self) -> Dict:
        """Determine the correct connection parameters from config."""
        params = {}
        
        # First try Supabase configuration
        if all(k in self.config for k in ['SUPABASE_DB_HOST', 'SUPABASE_DB_PORT', 
                                         'SUPABASE_DB_NAME', 'SUPABASE_DB_USER', 
                                         'SUPABASE_DB_PASSWORD']):
            params = {
                'host': self.config['SUPABASE_DB_HOST'],
                'port': self.config['SUPABASE_DB_PORT'],
                'dbname': self.config['SUPABASE_DB_NAME'],
                'user': self.config['SUPABASE_DB_USER'],
                'password': self.config['SUPABASE_DB_PASSWORD'],
                'sslmode': 'require'
            }
        # Next try direct DB configuration
        elif all(k in self.config for k in ['DB_HOST', 'DB_PORT', 'DB_NAME', 'DB_USER', 'DB_PASSWORD']):
            params = {
                'host': self.config['DB_HOST'],
                'port': self.config['DB_PORT'],
                'dbname': self.config['DB_NAME'],
                'user': self.config['DB_USER'],
                'password': self.config['DB_PASSWORD']
            }
        else:
            logger.error("Database connection parameters not found in config")
            sys.exit(1)
            
        return params
    
    def connect(self) -> None:
        """Establish connection to the database."""
        try:
            logger.info(f"Connecting to database at {self.conn_params['host']}:{self.conn_params['port']}")
            self.conn = psycopg2.connect(**self.conn_params)
            logger.info("Connection established successfully")
        except Exception as e:
            logger.error(f"Failed to connect to database: {str(e)}")
            sys.exit(1)
    
    def verify_data_integrity(self) -> None:
        """Check for data integrity issues in ChEMBL imports."""
        logger.info("Verifying ChEMBL data integrity...")
        
        if not self.conn:
            self.connect()
        
        try:
            with self.conn.cursor(cursor_factory=DictCursor) as cursor:
                # Count total compounds
                cursor.execute("""
                    SELECT COUNT(*) 
                    FROM molecules 
                    WHERE source = 'ChEMBL'
                """)
                total_count = cursor.fetchone()[0]
                self.results['general']['total_compounds'] = total_count
                logger.info(f"Found {total_count} ChEMBL compounds in the database")
                
                # Check for essential fields
                for field in ESSENTIAL_FIELDS:
                    cursor.execute(f"""
                        SELECT COUNT(*) 
                        FROM molecules 
                        WHERE source = 'ChEMBL' AND ({field} IS NULL OR {field} = '')
                    """)
                    missing_count = cursor.fetchone()[0]
                    self.results['general'][f'missing_{field}'] = missing_count
                    logger.info(f"Found {missing_count} compounds with missing {field}")
                    
                    if missing_count > 0:
                        # Get sample of compounds with missing field
                        cursor.execute(f"""
                            SELECT id, name, chembl_id, pubchem_cid
                            FROM molecules
                            WHERE source = 'ChEMBL' AND ({field} IS NULL OR {field} = '')
                            LIMIT 10
                        """)
                        samples = cursor.fetchall()
                        self.results['specific_issues'].append({
                            'issue_type': f'missing_{field}',
                            'count': missing_count,
                            'percentage': round((missing_count / total_count) * 100, 2),
                            'samples': [dict(s) for s in samples]
                        })
                
                # Check for required properties
                cursor.execute("""
                    SELECT m.id, m.name, m.chembl_id, m.pubchem_cid, 
                           COUNT(DISTINCT mp.property_type_id) as property_count
                    FROM molecules m
                    LEFT JOIN molecular_properties mp ON m.id = mp.molecule_id
                    LEFT JOIN molecular_property_types mpt ON mp.property_type_id = mpt.id
                    WHERE m.source = 'ChEMBL' AND mpt.name IN (
                        %s, %s, %s, %s, %s, %s, %s
                    )
                    GROUP BY m.id, m.name, m.chembl_id, m.pubchem_cid
                    HAVING COUNT(DISTINCT mp.property_type_id) < %s
                """, (*REQUIRED_PROPERTIES, len(REQUIRED_PROPERTIES)))
                
                compounds_with_missing_props = cursor.fetchall()
                missing_props_count = len(compounds_with_missing_props)
                self.results['general']['missing_properties'] = missing_props_count
                logger.info(f"Found {missing_props_count} compounds with missing required properties")
                
                if missing_props_count > 0:
                    # Add samples to specific issues
                    self.results['specific_issues'].append({
                        'issue_type': 'incomplete_properties',
                        'count': missing_props_count,
                        'percentage': round((missing_props_count / total_count) * 100, 2),
                        'samples': [dict(s) for s in compounds_with_missing_props[:10]]
                    })
                
                # Property distribution statistics
                for prop in REQUIRED_PROPERTIES:
                    cursor.execute("""
                        SELECT COUNT(mp.value) as count, 
                               AVG(mp.value::float) as avg, 
                               MIN(mp.value::float) as min, 
                               MAX(mp.value::float) as max
                        FROM molecules m
                        JOIN molecular_properties mp ON m.id = mp.molecule_id
                        JOIN molecular_property_types mpt ON mp.property_type_id = mpt.id
                        WHERE m.source = 'ChEMBL' AND mpt.name = %s AND mp.value ~ '^[0-9]+\.?[0-9]*$'
                    """, (prop,))
                    stats = cursor.fetchone()
                    if stats:
                        self.results['property_stats'][prop] = {
                            'count': stats['count'],
                            'coverage': round((stats['count'] / total_count) * 100, 2),
                            'avg': stats['avg'],
                            'min': stats['min'],
                            'max': stats['max']
                        }
                        logger.info(f"Property {prop}: {stats['count']} values (coverage: {self.results['property_stats'][prop]['coverage']}%)")
        
        except Exception as e:
            logger.error(f"Error during data integrity verification: {str(e)}")
            if self.conn:
                self.conn.rollback()
    
    def verify_pubchem_cross_references(self) -> None:
        """Verify PubChem cross-references for ChEMBL compounds."""
        logger.info("Verifying PubChem cross-references...")
        
        if not self.conn:
            self.connect()
        
        try:
            with self.conn.cursor(cursor_factory=DictCursor) as cursor:
                # Count compounds with missing PubChem CIDs
                cursor.execute("""
                    SELECT COUNT(*) 
                    FROM molecules 
                    WHERE source = 'ChEMBL' AND (pubchem_cid IS NULL OR pubchem_cid = '')
                """)
                missing_cid_count = cursor.fetchone()[0]
                self.results['general']['missing_pubchem_refs'] = missing_cid_count
                logger.info(f"Found {missing_cid_count} ChEMBL compounds without PubChem CIDs")
                
                if missing_cid_count > 0:
                    # Get sample of compounds with missing CIDs
                    cursor.execute("""
                        SELECT id, name, chembl_id, inchi_key
                        FROM molecules
                        WHERE source = 'ChEMBL' AND (pubchem_cid IS NULL OR pubchem_cid = '')
                        LIMIT 10
                    """)
                    samples = cursor.fetchall()
                    self.results['specific_issues'].append({
                        'issue_type': 'missing_pubchem_cid',
                        'count': missing_cid_count,
                        'percentage': round((missing_cid_count / self.results['general']['total_compounds']) * 100, 2),
                        'samples': [dict(s) for s in samples]
                    })
                
                # Check if there are any duplicates (same ChEMBL ID but different PubChem CIDs)
                cursor.execute("""
                    SELECT chembl_id, COUNT(DISTINCT pubchem_cid) as cid_count
                    FROM molecules
                    WHERE source = 'ChEMBL' AND chembl_id IS NOT NULL AND pubchem_cid IS NOT NULL
                    GROUP BY chembl_id
                    HAVING COUNT(DISTINCT pubchem_cid) > 1
                """)
                duplicate_mappings = cursor.fetchall()
                duplicate_count = len(duplicate_mappings)
                self.results['pubchem_stats']['duplicate_mappings'] = duplicate_count
                logger.info(f"Found {duplicate_count} ChEMBL IDs with multiple PubChem CIDs")
                
                if duplicate_count > 0:
                    # Get details for duplicates
                    duplicate_examples = []
                    for dup in duplicate_mappings[:5]:  # Limit to 5 examples
                        cursor.execute("""
                            SELECT id, name, chembl_id, pubchem_cid, inchi_key
                            FROM molecules
                            WHERE chembl_id = %s AND pubchem_cid IS NOT NULL
                        """, (dup['chembl_id'],))
                        compounds = cursor.fetchall()
                        duplicate_examples.append({
                            'chembl_id': dup['chembl_id'],
                            'cid_count': dup['cid_count'],
                            'compounds': [dict(c) for c in compounds]
                        })
                    
                    self.results['specific_issues'].append({
                        'issue_type': 'duplicate_pubchem_mappings',
                        'count': duplicate_count,
                        'examples': duplicate_examples
                    })
                
                # Check for consistency between InChI Keys and PubChem CIDs
                cursor.execute("""
                    SELECT m1.id, m1.chembl_id, m1.inchi_key, m1.pubchem_cid, m2.id as conflict_id, m2.pubchem_cid as conflict_cid
                    FROM molecules m1
                    JOIN molecules m2 ON m1.inchi_key = m2.inchi_key AND m1.pubchem_cid != m2.pubchem_cid
                    WHERE m1.source = 'ChEMBL' AND m2.source = 'ChEMBL'
                    LIMIT 100
                """)
                inchikey_conflicts = cursor.fetchall()
                conflict_count = len(inchikey_conflicts)
                self.results['pubchem_stats']['inchikey_cid_conflicts'] = conflict_count
                logger.info(f"Found {conflict_count} InChI Key conflicts with different PubChem CIDs")
                
                if conflict_count > 0:
                    self.results['specific_issues'].append({
                        'issue_type': 'inchikey_cid_conflicts',
                        'count': conflict_count,
                        'samples': [dict(c) for c in inchikey_conflicts[:10]]  # Limit to 10 examples
                    })
        
        except Exception as e:
            logger.error(f"Error during PubChem cross-reference verification: {str(e)}")
            if self.conn:
                self.conn.rollback()
    
    def verify_structures(self) -> None:
        """Verify molecular structures using RDKit when available."""
        if not RDKIT_AVAILABLE:
            logger.info("Skipping structure verification as RDKit is not available")
            self.results['structure_stats']['skipped'] = True
            return
        
        logger.info("Verifying molecular structures with RDKit...")
        
        if not self.conn:
            self.connect()
        
        try:
            with self.conn.cursor(cursor_factory=DictCursor) as cursor:
                # Get a sample of compounds to validate
                cursor.execute("""
                    SELECT id, name, chembl_id, smiles, inchi
                    FROM molecules
                    WHERE source = 'ChEMBL' AND smiles IS NOT NULL
                    ORDER BY RANDOM()
                    LIMIT 1000  -- Sample size for structure validation
                """)
                compounds = cursor.fetchall()
                
                valid_count = 0
                invalid_count = 0
                invalid_samples = []
                
                for compound in compounds:
                    # Validate SMILES
                    mol = None
                    try:
                        mol = Chem.MolFromSmiles(compound['smiles'])
                    except:
                        pass
                    
                    if mol is None:
                        invalid_count += 1
                        invalid_samples.append({
                            'id': compound['id'],
                            'name': compound['name'],
                            'chembl_id': compound['chembl_id'],
                            'smiles': compound['smiles'],
                            'issue': 'Invalid SMILES'
                        })
                    else:
                        valid_count += 1
                
                self.results['structure_stats']['validated_count'] = len(compounds)
                self.results['structure_stats']['valid_count'] = valid_count
                self.results['structure_stats']['invalid_count'] = invalid_count
                self.results['structure_stats']['valid_percentage'] = round((valid_count / len(compounds)) * 100, 2)
                
                logger.info(f"Validated {len(compounds)} structures: {valid_count} valid, {invalid_count} invalid")
                
                if invalid_count > 0:
                    self.results['specific_issues'].append({
                        'issue_type': 'invalid_structures',
                        'count': invalid_count,
                        'percentage': round((invalid_count / len(compounds)) * 100, 2),
                        'samples': invalid_samples[:10]  # Limit to 10 examples
                    })
                
                # Update the general results
                self.results['general']['invalid_structures'] = invalid_count
        
        except Exception as e:
            logger.error(f"Error during structure verification: {str(e)}")
            if self.conn:
                self.conn.rollback()
    
    def generate_recommendations(self) -> List[Dict]:
        """Generate recommendations based on verification results."""
        recommendations = []
        
        total = self.results['general']['total_compounds']
        
        # Missing properties recommendations
        missing_props = self.results['general'].get('missing_properties', 0)
        if missing_props > 0:
            missing_percentage = round((missing_props / total) * 100, 2)
            if missing_percentage > 20:
                recommendations.append({
                    'priority': 'high',
                    'issue': f'High number of missing properties ({missing_percentage}%)',
                    'action': 'Run property calculation scripts to fill in missing properties'
                })
            elif missing_percentage > 5:
                recommendations.append({
                    'priority': 'medium',
                    'issue': f'Moderate number of missing properties ({missing_percentage}%)',
                    'action': 'Consider running property calculation for specific property types'
                })
        
        # Missing PubChem CIDs
        missing_cids = self.results['general'].get('missing_pubchem_refs', 0)
        if missing_cids > 0:
            missing_cid_percentage = round((missing_cids / total) * 100, 2)
            if missing_cid_percentage > 30:
                recommendations.append({
                    'priority': 'high',
                    'issue': f'High number of missing PubChem CIDs ({missing_cid_percentage}%)',
                    'action': 'Run PubChem resolver to establish cross-references'
                })
            elif missing_cid_percentage > 10:
                recommendations.append({
                    'priority': 'medium',
                    'issue': f'Moderate number of missing PubChem CIDs ({missing_cid_percentage}%)',
                    'action': 'Run targeted PubChem resolution for high-priority compounds'
                })
        
        # Structure validation issues
        if RDKIT_AVAILABLE:
            invalid_structures = self.results['structure_stats'].get('invalid_count', 0)
            validated_count = self.results['structure_stats'].get('validated_count', 0)
            if validated_count > 0 and invalid_structures > 0:
                invalid_percentage = round((invalid_structures / validated_count) * 100, 2)
                if invalid_percentage > 10:
                    recommendations.append({
                        'priority': 'high',
                        'issue': f'High number of invalid structures ({invalid_percentage}%)',
                        'action': 'Review and fix SMILES/InChI representations in the database'
                    })
                elif invalid_percentage > 2:
                    recommendations.append({
                        'priority': 'medium',
                        'issue': f'Some invalid structures detected ({invalid_percentage}%)',
                        'action': 'Check and fix specific invalid structures'
                    })
        
        # Add recommendations to results
        self.results['recommendations'] = recommendations
        return recommendations
    
    def run_verification(self) -> Dict:
        """Run all verification steps and return results."""
        start_time = time.time()
        logger.info("Starting ChEMBL import verification")
        
        try:
            # Connect to database
            self.connect()
            
            # Run verification steps
            self.verify_data_integrity()
            self.verify_pubchem_cross_references()
            self.verify_structures()
            
            # Generate recommendations
            self.generate_recommendations()
            
            # Add summary stats
            self.results['general']['verified_compounds'] = self.results['general']['total_compounds']
            self.results['general']['duration_seconds'] = round(time.time() - start_time, 2)
            
            # Save results to file
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_path = self.output_dir / f'chembl_verification_report_{timestamp}.json'
            
            with open(output_path, 'w') as f:
                json.dump(self.results, f, indent=2)
            
            logger.info(f"Verification complete. Results saved to {output_path}")
            
            # Generate human-readable markdown report
            self.generate_markdown_report(output_path.with_suffix('.md'))
            
            return self.results
            
        except Exception as e:
            logger.error(f"Error during verification: {str(e)}")
            raise
        finally:
            if self.conn:
                self.conn.close()
    
    def generate_markdown_report(self, output_path: Path) -> None:
        """Generate a human-readable Markdown report."""
        with open(output_path, 'w') as f:
            f.write("# ChEMBL Import Verification Report\n\n")
            f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # General statistics
            f.write("## General Statistics\n\n")
            f.write(f"- Total ChEMBL compounds: {self.results['general']['total_compounds']}\n")
            f.write(f"- Compounds with missing properties: {self.results['general']['missing_properties']} ")
            f.write(f"({round((self.results['general']['missing_properties'] / self.results['general']['total_compounds']) * 100, 2)}%)\n")
            f.write(f"- Compounds with missing PubChem CIDs: {self.results['general']['missing_pubchem_refs']} ")
            f.write(f"({round((self.results['general']['missing_pubchem_refs'] / self.results['general']['total_compounds']) * 100, 2)}%)\n")
            
            if RDKIT_AVAILABLE:
                f.write(f"- Invalid structures (from sample): {self.results['structure_stats']['invalid_count']} ")
                f.write(f"({self.results['structure_stats']['valid_percentage']}% valid)\n")
            
            f.write(f"- Verification duration: {self.results['general']['duration_seconds']} seconds\n\n")
            
            # Property statistics
            f.write("## Property Coverage\n\n")
            f.write("| Property | Count | Coverage (%) | Average | Min | Max |\n")
            f.write("|----------|-------|-------------|---------|-----|-----|\n")
            
            for prop, stats in self.results['property_stats'].items():
                f.write(f"| {prop} | {stats['count']} | {stats['coverage']}% | {stats['avg']:.2f} | {stats['min']:.2f} | {stats['max']:.2f} |\n")
            
            f.write("\n")
            
            # PubChem cross-reference statistics
            f.write("## PubChem Cross-References\n\n")
            f.write(f"- ChEMBL compounds with PubChem CIDs: {self.results['general']['total_compounds'] - self.results['general']['missing_pubchem_refs']} ")
            f.write(f"({round(((self.results['general']['total_compounds'] - self.results['general']['missing_pubchem_refs']) / self.results['general']['total_compounds']) * 100, 2)}%)\n")
            
            if 'duplicate_mappings' in self.results['pubchem_stats']:
                f.write(f"- ChEMBL IDs with multiple PubChem CIDs: {self.results['pubchem_stats']['duplicate_mappings']}\n")
            
            if 'inchikey_cid_conflicts' in self.results['pubchem_stats']:
                f.write(f"- InChI Key conflicts with different PubChem CIDs: {self.results['pubchem_stats']['inchikey_cid_conflicts']}\n")
            
            f.write("\n")
            
            # Recommendations
            f.write("## Recommendations\n\n")
            
            if not self.results.get('recommendations'):
                f.write("No specific recommendations at this time.\n\n")
            else:
                for rec in self.results['recommendations']:
                    f.write(f"**{rec['priority'].upper()}**: {rec['issue']}\n")
                    f.write(f"- Action: {rec['action']}\n\n")
            
            # Specific issues
            f.write("## Specific Issues\n\n")
            
            if not self.results.get('specific_issues'):
                f.write("No specific issues identified.\n")
            else:
                for issue in self.results['specific_issues']:
                    f.write(f"### {issue['issue_type'].replace('_', ' ').title()}\n\n")
                    f.write(f"- Count: {issue['count']}\n")
                    if 'percentage' in issue:
                        f.write(f"- Percentage: {issue['percentage']}%\n")
                    
                    if 'samples' in issue and issue['samples']:
                        f.write("\nSample issues:\n\n")
                        f.write("```\n")
                        for sample in issue['samples'][:5]:  # Limit to 5 samples
                            f.write(f"{sample}\n")
                        f.write("```\n\n")
            
            logger.info(f"Markdown report generated at {output_path}")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="ChEMBL Import Verification Tool")
    parser.add_argument('--config', default='config.py', help='Path to configuration file')
    parser.add_argument('--output', default='./reports', help='Path to output report directory')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    
    try:
        verifier = ChEMBLVerifier(config_path=args.config, output_dir=args.output)
        results = verifier.run_verification()
        
        # Print summary to console
        print("\nVerification Summary:")
        print(f"Total ChEMBL compounds: {results['general']['total_compounds']}")
        print(f"Missing properties: {results['general']['missing_properties']} ({round((results['general']['missing_properties'] / results['general']['total_compounds']) * 100, 2)}%)")
        print(f"Missing PubChem CIDs: {results['general']['missing_pubchem_refs']} ({round((results['general']['missing_pubchem_refs'] / results['general']['total_compounds']) * 100, 2)}%)")
        
        if results.get('recommendations'):
            print("\nKey Recommendations:")
            for rec in results['recommendations']:
                if rec['priority'] == 'high':
                    print(f"- {rec['action']}")
        
        sys.exit(0)
    except Exception as e:
        logger.error(f"Verification failed: {str(e)}")
        sys.exit(1)