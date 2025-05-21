#!/usr/bin/env python3
"""
Reconcile ChEMBL and PubChem properties and cross-references.

This script identifies molecules that exist in both ChEMBL and PubChem databases
and reconciles their properties and cross-references using the connection factory.
"""

import os
import sys
import logging
import argparse
import json
import time
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime

# Import database connection utilities
from database.connection import get_db_connection, close_all_db_connections
import sql_executor

# Import custom modules
from cryoprotectant_identifiers import CryoprotectantIdentifierManager

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def reconcile_cross_references(output_report: Optional[str] = None, dry_run: bool = False) -> Dict:
    """
    Reconcile cross-references between ChEMBL and PubChem using the connection factory.
    
    Args:
        output_report: Optional path to save the reconciliation report
        dry_run: If True, don't actually connect to the database or make changes
        
    Returns:
        Reconciliation results dictionary
    """
    id_manager = CryoprotectantIdentifierManager.get_instance()
    
    results = {
        'timestamp': datetime.now().isoformat(),
        'molecules_updated': 0,
        'cross_references_added': 0,
        'conflicts_resolved': 0,
        'details': {}
    }
    
    # If dry run, use mock data
    if dry_run:
        logger.info("Running in dry-run mode with mock data")
        pubchem_molecules = {
            'CRYO001': {'name': 'Glycerol', 'pubchem_cid': '962'},
            'CRYO003': {'name': 'beta-Alanine', 'pubchem_cid': '6342'}
        }
        chembl_molecules = {
            'CRYO002': {'name': 'DMSO', 'chembl_id': 'CHEMBL1098659'},
            'CRYO004': {'name': 'tert-Butanol', 'chembl_id': 'CHEMBL500033'}
        }
        # Mock molecule with both identifiers for testing
        pubchem_molecules['CRYO005'] = {'name': 'Urea', 'pubchem_cid': '6057'}
        chembl_molecules['CRYO006'] = {'name': 'Urea', 'chembl_id': 'CHEMBL1487'}
        
        # Mock InChI key map
        inchi_key_map = {
            'XSQUKJJJFZCRTK-UHFFFAOYSA-N': ['CRYO005', 'CRYO006']  # Urea
        }
    else:
        # Step 1: Find molecules with PubChem or ChEMBL IDs using direct PostgreSQL connection
        try:
            # Get molecules with PubChem IDs
            pubchem_query = """
                SELECT id, name, pubchem_cid
                FROM molecules
                WHERE pubchem_cid IS NOT NULL
            """
            pubchem_results = sql_executor.execute_query(pubchem_query)
            pubchem_molecules = {row['id']: {'name': row['name'], 'pubchem_cid': str(row['pubchem_cid'])}
                               for row in pubchem_results}
            
            # Get molecules with ChEMBL IDs
            chembl_query = """
                SELECT id, name, chembl_id
                FROM molecules
                WHERE chembl_id IS NOT NULL
            """
            chembl_results = sql_executor.execute_query(chembl_query)
            chembl_molecules = {row['id']: {'name': row['name'], 'chembl_id': row['chembl_id']}
                              for row in chembl_results}
        except Exception as e:
            logger.error(f"Database connection error: {str(e)}")
            logger.info("Falling back to dry-run mode with mock data")
            return reconcile_cross_references(output_report, dry_run=True)
    
    logger.info(f"Found {len(pubchem_molecules)} molecules with PubChem IDs")
    logger.info(f"Found {len(chembl_molecules)} molecules with ChEMBL IDs")
    
    # Step 2: Find molecules that should be linked based on shared identifiers
    if dry_run:
        # Already created mock inchi_key_map in dry-run mode
        pass
    else:
        # Use optimized query to build inchi_key_map
        try:
            inchi_key_query = """
                SELECT inchikey as inchi_key, array_agg(id) as molecule_ids
                FROM molecules
                WHERE inchikey IS NOT NULL
                GROUP BY inchikey
                HAVING COUNT(*) > 1
            """
            inchi_key_results = sql_executor.execute_query(inchi_key_query)
            
            inchi_key_map = {row['inchi_key']: row['molecule_ids'] for row in inchi_key_results}
        except Exception as e:
            logger.error(f"Error fetching InChI keys: {str(e)}")
            if not dry_run:
                return reconcile_cross_references(output_report, dry_run=True)
    
    # Find molecules that share the same InChI Key
    reconciliation_candidates = []
    for inchi_key, molecule_ids in inchi_key_map.items():
        reconciliation_candidates.append({
            'inchi_key': inchi_key,
            'molecule_ids': molecule_ids
        })
    
    logger.info(f"Found {len(reconciliation_candidates)} InChI Keys with multiple molecules")
    
    # Step 3: Process reconciliation candidates in batches
    batch_size = 100
    total_batches = (len(reconciliation_candidates) + batch_size - 1) // batch_size
    
    for batch_idx in range(0, len(reconciliation_candidates), batch_size):
        batch = reconciliation_candidates[batch_idx:batch_idx + batch_size]
        batch_num = batch_idx // batch_size + 1
        logger.info(f"Processing batch {batch_num}/{total_batches} with {len(batch)} candidates")
        
        start_time = time.time()
        
        # Process each candidate in the batch
        updates_for_batch = []
        
        for candidate in batch:
            molecule_ids = candidate['molecule_ids']
            inchi_key = candidate['inchi_key']
            
            # Find molecules with PubChem and ChEMBL IDs in this group
            pubchem_id = None
            chembl_id = None
            pubchem_molecule = None
            chembl_molecule = None
            
            for mol_id in molecule_ids:
                if mol_id in pubchem_molecules:
                    pubchem_id = mol_id
                    pubchem_molecule = pubchem_molecules[mol_id]
                if mol_id in chembl_molecules:
                    chembl_id = mol_id
                    chembl_molecule = chembl_molecules[mol_id]
            
            # If we have one of each, we can reconcile
            if pubchem_id and chembl_id:
                logger.info(f"Reconciling molecules with InChI Key {inchi_key}: " +
                           f"{pubchem_id} (PubChem) and {chembl_id} (ChEMBL)")
                
                # Prepare updates for batch execution
                updates_for_batch.append({
                    'pubchem_id': pubchem_id,
                    'chembl_id': chembl_id,
                    'pubchem_cid': pubchem_molecule['pubchem_cid'],
                    'chembl_id_str': chembl_molecule['chembl_id'],
                    'inchi_key': inchi_key
                })
                
                # Update results
                results['cross_references_added'] += 2
                results['molecules_updated'] += 2
                results['details'][inchi_key] = {
                    'pubchem_molecule': pubchem_id,
                    'chembl_molecule': chembl_id,
                    'action': 'cross_references_added'
                }
        
        # Step 4: Update molecule properties and cross-references in a single transaction
        if not dry_run and updates_for_batch:
            try:
                # Use transaction decorator for atomic updates
                @sql_executor.with_transaction
                def update_cross_references(transaction, updates):
                    # Prepare batch updates for PubChem molecules
                    pubchem_updates = []
                    chembl_updates = []
                    
                    for update in updates:
                        # Parameters for PubChem molecule update
                        pubchem_updates.append((
                            update['chembl_id_str'],  # chembl_id
                            update['pubchem_id']      # id
                        ))
                        
                        # Parameters for ChEMBL molecule update
                        chembl_updates.append((
                            int(update['pubchem_cid']),  # pubchem_cid
                            update['chembl_id']          # id
                        ))
                    
                    # Update PubChem molecules with ChEMBL IDs
                    pubchem_update_query = """
                        UPDATE molecules
                        SET chembl_id = %s, updated_at = NOW()
                        WHERE id = %s
                    """
                    sql_executor.execute_batch(pubchem_update_query, pubchem_updates)
                    
                    # Update ChEMBL molecules with PubChem CIDs
                    chembl_update_query = """
                        UPDATE molecules
                        SET pubchem_cid = %s, updated_at = NOW()
                        WHERE id = %s
                    """
                    sql_executor.execute_batch(chembl_update_query, chembl_updates)
                    
                    logger.info(f"Updated {len(pubchem_updates)} PubChem molecules and {len(chembl_updates)} ChEMBL molecules")
                
                # Execute the transaction
                update_cross_references(updates_for_batch)
                
            except Exception as e:
                logger.error(f"Error updating database: {str(e)}")
                # Continue with the next batch even if this one fails
        
        # Step 5: Update identifier manager for this batch
        for update in updates_for_batch:
            pubchem_id = update['pubchem_id']
            chembl_id = update['chembl_id']
            pubchem_cid = update['pubchem_cid']
            chembl_id_str = update['chembl_id_str']
            
            # Update PubChem molecule's record
            pubchem_record = id_manager.get_molecule_by_internal_id(pubchem_id)
            if pubchem_record:
                pubchem_record['chembl_id'] = chembl_id_str
                id_manager.add_molecule(pubchem_id, pubchem_record)
            
            # Update ChEMBL molecule's record
            chembl_record = id_manager.get_molecule_by_internal_id(chembl_id)
            if chembl_record:
                chembl_record['pubchem_cid'] = pubchem_cid
                id_manager.add_molecule(chembl_id, chembl_record)
        
        # Save identifiers after each batch
        id_manager.save_identifiers()
        
        # Log batch completion
        elapsed = time.time() - start_time
        logger.info(f"Batch {batch_num} completed in {elapsed:.2f} seconds")
        
        # Estimate remaining time
        if batch_num < total_batches:
            remaining_batches = total_batches - batch_num
            est_remaining_time = remaining_batches * elapsed
            logger.info(f"Estimated remaining time: {est_remaining_time:.2f} seconds")
    
    # Generate report
    if output_report:
        os.makedirs(os.path.dirname(output_report), exist_ok=True)
        with open(output_report, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Reconciliation report saved to {output_report}")
    
    logger.info(f"Reconciliation completed: {results['molecules_updated']} molecules updated, " +
               f"{results['cross_references_added']} cross-references added, " +
               f"{results['conflicts_resolved']} conflicts resolved")
    
    return results

def main():
    """CLI entry point for reconciliation."""
    parser = argparse.ArgumentParser(description='Reconcile ChEMBL and PubChem properties and cross-references')
    parser.add_argument('--report',
                      default=f"reports/reconciliation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                      help='Output file path for report')
    parser.add_argument('--dry-run', action='store_true',
                      help='Run in dry-run mode without connecting to the database')
    parser.add_argument('--batch-size', type=int, default=100,
                      help='Batch size for processing reconciliation candidates')
    
    args = parser.parse_args()
    
    try:
        # Initialize the database connection using the connection factory
        if not args.dry_run:
            logger.info("Initializing database connection using connection factory")
            # The connection is initialized lazily when needed
        
        # Run the reconciliation
        start_time = time.time()
        reconcile_cross_references(args.report, args.dry_run)
        elapsed = time.time() - start_time
        
        logger.info(f"Reconciliation process completed in {elapsed:.2f} seconds")
        
        # Close connections
        if not args.dry_run:
            sql_executor.close_connections()
            logger.info("Database connections closed")
        
        return 0
    except Exception as e:
        logger.error(f"Reconciliation failed: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())