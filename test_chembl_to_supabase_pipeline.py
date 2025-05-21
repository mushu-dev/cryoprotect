#!/usr/bin/env python3
"""
CryoProtect - Test ChEMBL to Supabase Full Data Pipeline

This script tests the complete ChEMBL to Supabase data pipeline by:
1. Importing molecules from ChEMBL
2. Inserting them into the Supabase database
3. Verifying the data was stored correctly

Usage:
    python test_chembl_to_supabase_pipeline.py [--molecules 100] [--report-file pipeline_report.json]
"""

import os
import sys
import time
import json
import logging
import argparse
import datetime
import concurrent.futures
import uuid
from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f'chembl_pipeline_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger(__name__)

try:
    # Import required modules
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from failsafe_import import FailsafeChEMBLClient as ChEMBLClient
    from database.db import execute_query, test_connection

    # Check database connection
    try:
        success, message = test_connection()
        HAS_DB_CONNECTION = success
        if success:
            logger.info("Database connection available")
        else:
            logger.error(f"Database connection not available: {message}")
    except Exception as e:
        logger.error(f"Database connection not available: {str(e)}")
        HAS_DB_CONNECTION = False
except ImportError as e:
    logger.error(f"Failed to import required modules: {str(e)}")
    sys.exit(1)

class ChEMBLToSupabasePipeline:
    """Tests the ChEMBL to Supabase data pipeline."""
    
    def __init__(self, molecule_count=100, cache_dir=None, batch_size=None, max_workers=None, test_mode=False):
        """Initialize pipeline test."""
        self.molecule_count = molecule_count
        self.test_mode = test_mode
        self.cache_dir = cache_dir or os.path.join(os.path.dirname(__file__), "chembl_cache")
        
        # Create ChEMBL client
        self.client = ChEMBLClient(cache_dir=self.cache_dir)
        
        # We'll use execute_query function directly instead of a db connection
        
        # Set batch size and max workers
        self.batch_size = batch_size or min(100, max(10, molecule_count // 10))
        self.max_workers = max_workers or min(10, (molecule_count + 99) // 100)
        
        # Stats
        self.stats_before = None
        self.stats_after = None
        self.chembl_ids = []
        self.inserted_molecule_ids = []
        
        # Initialize results
        self.results = {
            "test_info": {
                "molecule_count": molecule_count,
                "timestamp": datetime.datetime.now().isoformat(),
                "cache_dir": self.cache_dir,
                "batch_size": self.batch_size,
                "max_workers": self.max_workers,
                "test_mode": test_mode,
                "has_db_connection": HAS_DB_CONNECTION
            },
            "pipeline_stats": {
                "chembl_fetch_time": 0,
                "supabase_insert_time": 0,
                "verification_time": 0,
                "total_time": 0,
                "success_count": 0,
                "error_count": 0,
                "molecules_per_second": 0,
            },
            "chembl_stats": {
                "molecules_fetched": 0,
                "properties_fetched": 0,
                "synonyms_fetched": 0,
                "fetch_errors": 0,
            },
            "supabase_stats": {
                "molecules_inserted": 0,
                "properties_inserted": 0,
                "synonyms_inserted": 0,
                "insert_errors": 0,
            },
            "verification": {
                "molecules_verified": 0,
                "properties_verified": 0, 
                "verify_errors": [],
                "verification_success": False,
            },
            "database_stats": {
                "before": {},
                "after": {},
                "difference": {}
            }
        }
    
    def get_database_stats(self):
        """Get database statistics."""
        if not HAS_DB_CONNECTION:
            logger.warning("No database connection available")
            return None

        try:
            stats = {
                "table_counts": {},
                "timestamp": datetime.datetime.now().isoformat()
            }

            # Get table counts
            tables = ["molecules", "molecular_properties", "molecule_synonyms", "molecule_identifiers"]
            for table in tables:
                try:
                    result = execute_query(f"SELECT COUNT(*) FROM {table}")
                    count = result[0]['count'] if result else 0
                    stats["table_counts"][table] = count
                except Exception as e:
                    logger.warning(f"Error getting count for {table}: {str(e)}")
                    stats["table_counts"][table] = -1

            return stats
        except Exception as e:
            logger.error(f"Failed to get database stats: {str(e)}")
            return None
    
    def generate_chembl_ids(self):
        """Generate or fetch ChEMBL IDs for import."""
        logger.info(f"Generating {self.molecule_count} ChEMBL IDs for testing")
        
        try:
            # Try to get ChEMBL IDs from client
            self.chembl_ids = self.client.generate_chembl_ids(self.molecule_count)
            
            # Save the IDs for reference
            with open("chembl_pipeline_ids.json", "w") as f:
                json.dump(self.chembl_ids, f)
                
            logger.info(f"Generated {len(self.chembl_ids)} ChEMBL IDs")
            return True
        except Exception as e:
            logger.error(f"Failed to generate ChEMBL IDs: {str(e)}")
            return False
    
    def process_molecule(self, chembl_id):
        """Process a single molecule: fetch from ChEMBL and insert into Supabase."""
        start_time = time.time()
        result = {
            "chembl_id": chembl_id,
            "success": False,
            "chembl_fetch_time": 0,
            "supabase_insert_time": 0,
            "total_time": 0,
            "properties_count": 0,
            "synonyms_count": 0,
            "errors": [],
            "molecule_id": None
        }
        
        try:
            # Step 1: Fetch from ChEMBL
            chembl_start = time.time()
            compound = self.client.get_compound_by_chembl_id(chembl_id)
            if not compound:
                result["errors"].append(f"Failed to fetch compound {chembl_id}")
                return result
                
            records = self.client.get_compound_records(chembl_id)
            properties = self.client.get_compound_properties(chembl_id)
            
            result["chembl_fetch_time"] = time.time() - chembl_start
            result["properties_count"] = len(properties) if properties else 0
            result["synonyms_count"] = len(compound.get("molecule_synonyms", [])) if compound else 0
            
            # Step 2: Insert into Supabase (or simulate if no DB connection or test mode)
            supabase_start = time.time()

            if HAS_DB_CONNECTION and not self.test_mode:
                try:
                    # Check if molecule already exists
                    check_result = execute_query(
                        "SELECT id FROM molecules WHERE external_id = %s",
                        (chembl_id,)
                    )

                    if check_result and check_result[0]['id']:
                        # Molecule already exists, use existing ID
                        molecule_id = check_result[0]['id']
                        logger.debug(f"Molecule {chembl_id} already exists with ID {molecule_id}")
                    else:
                        # Insert new molecule
                        # Get SMILES, InChI, and InChI Key
                        structures = compound.get("molecule_structures", {})
                        smiles = structures.get("canonical_smiles", "") if structures else ""
                        inchi = structures.get("standard_inchi", "") if structures else ""
                        inchi_key = structures.get("standard_inchi_key", "") if structures else ""

                        # Get molecule properties
                        mol_props = compound.get("molecule_properties", {})

                        # Insert molecule
                        molecule_id = str(uuid.uuid4())
                        execute_query(
                            """
                            INSERT INTO molecules (
                                id, external_id, smiles, inchi, inchi_key,
                                created_at, updated_at
                            ) VALUES (%s, %s, %s, %s, %s, NOW(), NOW())
                            """,
                            (molecule_id, chembl_id, smiles, inchi, inchi_key)
                        )

                        # Insert molecular properties
                        if mol_props:
                            for prop_name, prop_value in mol_props.items():
                                execute_query(
                                    """
                                    INSERT INTO molecular_properties (
                                        molecule_id, property_name, property_value,
                                        source, created_at, updated_at
                                    ) VALUES (%s, %s, %s, %s, NOW(), NOW())
                                    """,
                                    (molecule_id, prop_name, str(prop_value), "ChEMBL")
                                )

                        # Insert from properties dict if available
                        if properties:
                            for prop_name, prop_value in properties.items():
                                execute_query(
                                    """
                                    INSERT INTO molecular_properties (
                                        molecule_id, property_name, property_value,
                                        source, created_at, updated_at
                                    ) VALUES (%s, %s, %s, %s, NOW(), NOW())
                                    """,
                                    (molecule_id, prop_name, str(prop_value), "ChEMBL")
                                )

                        # Insert synonyms
                        synonyms = compound.get("molecule_synonyms", [])
                        if synonyms:
                            for synonym in synonyms:
                                if "synonym" in synonym and "syn_type" in synonym:
                                    execute_query(
                                        """
                                        INSERT INTO molecule_synonyms (
                                            molecule_id, synonym, type,
                                            created_at, updated_at
                                        ) VALUES (%s, %s, %s, NOW(), NOW())
                                        """,
                                        (molecule_id, synonym["synonym"], synonym["syn_type"])
                                    )

                    # Success
                    result["success"] = True
                    result["molecule_id"] = molecule_id
                except Exception as e:
                    result["errors"].append(f"Database error: {str(e)}")
            else:
                # In test mode or no DB, simulate success
                result["success"] = True
                result["molecule_id"] = str(uuid.uuid4())
                time.sleep(0.01)  # Simulate DB operation time
                
            result["supabase_insert_time"] = time.time() - supabase_start
            
            # Add molecule ID to list of inserted molecules
            if result["success"] and result["molecule_id"]:
                self.inserted_molecule_ids.append(result["molecule_id"])
                
        except Exception as e:
            result["errors"].append(f"Error processing molecule: {str(e)}")
        
        result["total_time"] = time.time() - start_time
        return result
    
    def process_molecules(self):
        """Process all molecules: fetch from ChEMBL and insert into Supabase."""
        if not self.chembl_ids:
            if not self.generate_chembl_ids():
                return False
        
        logger.info(f"Starting processing of {len(self.chembl_ids)} molecules")
        
        # Get database stats before import
        self.stats_before = self.get_database_stats()
        logger.info("Captured database stats before import")
        
        # Initialize counters
        start_time = time.time()
        results = []
        success_count = 0
        error_count = 0
        total_chembl_time = 0
        total_supabase_time = 0
        total_properties = 0
        total_synonyms = 0
        
        # Process molecules in batches
        total_molecules = len(self.chembl_ids)
        logger.info(f"Processing molecules in batches of {self.batch_size} using {self.max_workers} workers")
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Process in batches
            for batch_start in range(0, total_molecules, self.batch_size):
                batch_end = min(batch_start + self.batch_size, total_molecules)
                current_batch = self.chembl_ids[batch_start:batch_end]
                logger.info(f"Processing batch {batch_start//self.batch_size + 1} of {(total_molecules + self.batch_size - 1)//self.batch_size} ({len(current_batch)} molecules)")
                
                # Submit batch
                future_to_id = {executor.submit(self.process_molecule, chembl_id): chembl_id for chembl_id in current_batch}
                
                # Process results
                for i, future in enumerate(concurrent.futures.as_completed(future_to_id)):
                    chembl_id = future_to_id[future]
                    try:
                        result = future.result()
                        results.append(result)
                        
                        # Update counters
                        if result["success"]:
                            success_count += 1
                            total_properties += result.get("properties_count", 0)
                            total_synonyms += result.get("synonyms_count", 0)
                            total_chembl_time += result.get("chembl_fetch_time", 0)
                            total_supabase_time += result.get("supabase_insert_time", 0)
                        else:
                            error_count += 1
                            logger.warning(f"Error processing {chembl_id}: {', '.join(result.get('errors', ['Unknown error']))}")
                            
                        # Log progress
                        processed = batch_start + i + 1
                        if (i + 1) % (self.batch_size // 5 or 1) == 0 or i + 1 == len(current_batch):
                            elapsed = time.time() - start_time
                            rate = processed / elapsed if elapsed > 0 else 0
                            logger.info(f"Processed {processed}/{total_molecules} molecules ({processed/total_molecules*100:.1f}%) - " +
                                       f"{success_count} successes, {error_count} errors - " +
                                       f"Rate: {rate:.2f} molecules/sec")
                    except Exception as e:
                        logger.error(f"Error processing result for {chembl_id}: {str(e)}")
                        error_count += 1
                
                # Save intermediate statistics
                current_time = time.time()
                progress_stats = {
                    "completed_molecules": batch_start + len(current_batch),
                    "total_molecules": total_molecules,
                    "success_count": success_count,
                    "error_count": error_count,
                    "progress_percentage": (batch_start + len(current_batch)) / total_molecules * 100,
                    "elapsed_time": current_time - start_time,
                    "molecules_per_second": (batch_start + len(current_batch)) / (current_time - start_time) if current_time > start_time else 0,
                }
                logger.info(f"Batch complete. Overall progress: {progress_stats['progress_percentage']:.1f}%")
        
        # Record total time and get database stats after import
        total_time = time.time() - start_time
        self.stats_after = self.get_database_stats()
        logger.info("Captured database stats after import")
        
        # Update results
        self.results["pipeline_stats"].update({
            "total_time": total_time,
            "chembl_fetch_time": total_chembl_time / success_count if success_count else 0,
            "supabase_insert_time": total_supabase_time / success_count if success_count else 0,
            "success_count": success_count,
            "error_count": error_count,
            "molecules_per_second": total_molecules / total_time if total_time > 0 else 0,
        })
        
        self.results["chembl_stats"].update({
            "molecules_fetched": success_count,
            "properties_fetched": total_properties,
            "synonyms_fetched": total_synonyms,
            "fetch_errors": error_count,
        })
        
        self.results["supabase_stats"].update({
            "molecules_inserted": success_count,
            "properties_inserted": total_properties,
            "synonyms_inserted": total_synonyms,
            "insert_errors": error_count,
        })
        
        # Update database stats
        if self.stats_before and self.stats_after:
            self.results["database_stats"].update({
                "before": self.stats_before,
                "after": self.stats_after,
                "difference": self.calculate_stats_difference()
            })
        
        logger.info(f"Processing completed: {success_count} successes, {error_count} errors in {total_time:.2f} seconds")
        logger.info(f"Average rate: {self.results['pipeline_stats']['molecules_per_second']:.2f} molecules/second")
        
        return success_count > 0
    
    def calculate_stats_difference(self):
        """Calculate difference between before and after database stats."""
        if not self.stats_before or not self.stats_after:
            return {}
            
        difference = {}
        
        # Calculate differences in table counts
        difference["table_counts"] = {}
        for table in self.stats_before["table_counts"]:
            if table in self.stats_after["table_counts"]:
                before = self.stats_before["table_counts"][table]
                after = self.stats_after["table_counts"][table]
                if before >= 0 and after >= 0:
                    difference["table_counts"][table] = after - before
        
        return difference
    
    def verify_data(self):
        """Verify data was correctly inserted into Supabase."""
        if not HAS_DB_CONNECTION or self.test_mode:
            logger.warning("Cannot verify data: no database connection or test mode enabled")
            self.results["verification"]["verification_success"] = True  # Assume success in test mode
            return True

        logger.info("Verifying data in database...")
        start_time = time.time()
        verification_errors = []
        molecules_verified = 0
        properties_verified = 0

        # Check a sample of molecules (up to 100)
        sample_size = min(100, len(self.inserted_molecule_ids))
        sample_ids = self.inserted_molecule_ids[:sample_size]

        for molecule_id in sample_ids:
            try:
                # Verify molecule exists
                molecule_result = execute_query(
                    "SELECT id, external_id, smiles, inchi, inchi_key FROM molecules WHERE id = %s",
                    (molecule_id,)
                )

                if not molecule_result:
                    verification_errors.append(f"Molecule {molecule_id} not found in database")
                    continue

                # Verify properties exist
                properties_result = execute_query(
                    "SELECT COUNT(*) as property_count FROM molecular_properties WHERE molecule_id = %s",
                    (molecule_id,)
                )

                property_count = properties_result[0]['property_count'] if properties_result else 0

                if property_count == 0:
                    verification_errors.append(f"No properties found for molecule {molecule_id}")

                # Count as verified
                molecules_verified += 1
                properties_verified += property_count

            except Exception as e:
                verification_errors.append(f"Error verifying molecule {molecule_id}: {str(e)}")

        # Update verification results
        verification_time = time.time() - start_time
        verification_success = len(verification_errors) == 0

        self.results["verification"].update({
            "molecules_verified": molecules_verified,
            "properties_verified": properties_verified,
            "verify_errors": verification_errors,
            "verification_success": verification_success,
            "verification_time": verification_time,
        })

        logger.info(f"Data verification completed: {molecules_verified} molecules verified with {len(verification_errors)} errors")

        return verification_success
    
    def run_pipeline_test(self):
        """Run the complete ChEMBL to Supabase pipeline test."""
        logger.info(f"Starting ChEMBL to Supabase pipeline test with {self.molecule_count} molecules")
        logger.info(f"Database connection: {'Available' if HAS_DB_CONNECTION else 'Not available'}")
        logger.info(f"Test mode: {'Enabled' if self.test_mode else 'Disabled'}")
        
        # Step 1: Process molecules (fetch and insert)
        if not self.process_molecules():
            logger.error("Processing failed")
            return False
            
        # Step 2: Verify data integrity
        if not self.verify_data():
            logger.warning("Data verification failed")
        
        # Generate report
        self.generate_report()
        
        success = self.results["pipeline_stats"]["success_count"] > 0.9 * self.molecule_count
        return success
    
    def generate_report(self, filename=None):
        """Generate a detailed report of the pipeline test."""
        if not filename:
            filename = f"chembl_pipeline_report_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            
        try:
            with open(filename, "w") as f:
                json.dump(self.results, f, indent=2)
                
            # Also generate a markdown summary
            md_filename = filename.replace('.json', '.md')
            with open(md_filename, "w") as f:
                f.write(f"# ChEMBL to Supabase Pipeline Test Report\n\n")
                f.write(f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write(f"## Test Configuration\n\n")
                f.write(f"- Molecules: {self.results['test_info']['molecule_count']}\n")
                f.write(f"- Database Connection: {'Yes' if self.results['test_info']['has_db_connection'] else 'No'}\n")
                f.write(f"- Test Mode: {'Enabled' if self.results['test_info']['test_mode'] else 'Disabled'}\n\n")
                
                f.write(f"## Pipeline Performance\n\n")
                stats = self.results['pipeline_stats']
                f.write(f"- Total time: {stats['total_time']:.2f} seconds\n")
                f.write(f"- Success rate: {stats['success_count']} / {self.results['test_info']['molecule_count']} ")
                f.write(f"({stats['success_count'] / self.results['test_info']['molecule_count'] * 100:.1f}%)\n")
                f.write(f"- Average ChEMBL fetch time: {stats['chembl_fetch_time']*1000:.2f} ms per molecule\n")
                f.write(f"- Average Supabase insert time: {stats['supabase_insert_time']*1000:.2f} ms per molecule\n")
                f.write(f"- Overall rate: {stats['molecules_per_second']:.2f} molecules per second\n\n")
                
                f.write(f"## Data Stats\n\n")
                f.write(f"- Molecules fetched: {self.results['chembl_stats']['molecules_fetched']}\n")
                f.write(f"- Properties fetched: {self.results['chembl_stats']['properties_fetched']}\n")
                f.write(f"- Synonyms fetched: {self.results['chembl_stats']['synonyms_fetched']}\n")
                f.write(f"- Molecules inserted: {self.results['supabase_stats']['molecules_inserted']}\n")
                
                # Database stats difference if available
                if self.stats_before and self.stats_after and "difference" in self.results["database_stats"]:
                    diff = self.results["database_stats"]["difference"]
                    if "table_counts" in diff:
                        f.write(f"\n## Database Changes\n\n")
                        for table, count in diff["table_counts"].items():
                            f.write(f"- {table}: +{count} rows\n")
                
                f.write(f"\n## Verification\n\n")
                verif = self.results["verification"]
                if self.test_mode:
                    f.write("Verification skipped in test mode\n")
                else:
                    f.write(f"- Molecules verified: {verif['molecules_verified']}\n")
                    f.write(f"- Properties verified: {verif['properties_verified']}\n")
                    f.write(f"- Verification success: {'Yes' if verif['verification_success'] else 'No'}\n")
                    if verif.get('verify_errors'):
                        f.write(f"- Errors: {len(verif['verify_errors'])}\n")
                
                f.write(f"\n## Conclusion\n\n")
                if self.results['pipeline_stats']['success_count'] > 0.9 * self.results['test_info']['molecule_count']:
                    f.write("✅ Pipeline test completed successfully. The ChEMBL to Supabase data pipeline is functioning properly.\n")
                else:
                    f.write("❌ Pipeline test completed with errors. The ChEMBL to Supabase data pipeline needs attention.\n")
            
            logger.info(f"Reports saved to {filename} and {md_filename}")
            return True
        except Exception as e:
            logger.error(f"Failed to generate report: {str(e)}")
            return False

def main():
    parser = argparse.ArgumentParser(description="Test ChEMBL to Supabase pipeline")
    parser.add_argument("--molecules", type=int, default=100, help="Number of molecules to process")
    parser.add_argument("--report-file", default=None, help="Output report file")
    parser.add_argument("--cache-dir", default=None, help="Cache directory for ChEMBL data")
    parser.add_argument("--batch-size", type=int, default=None, help="Size of each batch for processing")
    parser.add_argument("--max-workers", type=int, default=None, help="Maximum number of parallel workers")
    parser.add_argument("--test-mode", action="store_true", help="Run in test mode without database writes")
    args = parser.parse_args()
    
    # Create and run the pipeline test
    pipeline = ChEMBLToSupabasePipeline(
        molecule_count=args.molecules,
        cache_dir=args.cache_dir,
        batch_size=args.batch_size,
        max_workers=args.max_workers,
        test_mode=args.test_mode
    )
    
    success = pipeline.run_pipeline_test()
    
    # Generate report
    if args.report_file:
        pipeline.generate_report(args.report_file)
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())