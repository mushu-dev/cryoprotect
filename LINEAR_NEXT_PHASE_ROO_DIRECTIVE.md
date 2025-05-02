# ROO DIRECTIVE: DATABASE POPULATION & VERIFICATION

## PHASE OVERVIEW

With the ChEMBL integration framework now complete, we're entering the **Database Population & Verification Phase**. This phase will follow a linear execution model, where each task builds upon the previous one and reports back to the main orchestrator.

## TECHNICAL CONTEXT

The CryoProtect v2 application now has:
- ✅ Complete ChEMBL integration framework with error handling
- ✅ Efficient worker pool implementation for parallel processing
- ✅ Checkpoint system for resumable operations
- ✅ Reference compound identification system

Our primary focus is to populate the database with real scientific data and verify its quality.

## LINEAR TASK SEQUENCE

### STAGE 1: FRAMEWORK TESTING

**TASK 1: Basic Framework Unit Tests**
```
TASK: TEST-FRAMEWORK-1: Error Handler Unit Tests

FILE: tests/test_chembl_error_handler.py
REFERENCE: tests/test_pubchem/test_rate_limiter.py:10-50

IMPLEMENTATION:
1. Create unit tests for ErrorCategory enum
2. Create unit tests for classify_error function
3. Create unit tests for get_recovery_strategy function
4. Test all error categories and recovery strategies

INTERFACE:
```python
import unittest
from chembl.error_handler import ErrorCategory, classify_error, get_recovery_strategy

class TestErrorHandler(unittest.TestCase):
    def test_error_categories(self):
        """Test that all expected error categories exist"""
        self.assertTrue(hasattr(ErrorCategory, 'API_RATE_LIMIT'))
        self.assertTrue(hasattr(ErrorCategory, 'CONNECTION_ERROR'))
        # Additional categories...
        
    def test_classify_error(self):
        """Test error classification logic"""
        # Test HTTP errors
        http_error = requests.exceptions.HTTPError()
        http_error.response = type('obj', (object,), {'status_code': 429})
        category, _ = classify_error(http_error)
        self.assertEqual(category, ErrorCategory.API_RATE_LIMIT)
        
        # Additional test cases...
        
    def test_recovery_strategy(self):
        """Test recovery strategy determination"""
        strategy = get_recovery_strategy(ErrorCategory.API_RATE_LIMIT)
        self.assertEqual(strategy, "RETRY")
        
        # Additional test cases...
```

VERIFICATION:
- All tests pass successfully
- Coverage is >90% for error_handler.py
- All error categories are tested
- All recovery strategies are tested
```

**TASK 2: Checkpoint System Unit Tests**
```
TASK: TEST-FRAMEWORK-2: Checkpoint System Unit Tests

FILE: tests/test_chembl_checkpoint.py
REFERENCE: tests/test_pubchem/test_cache.py:5-45

IMPLEMENTATION:
1. Create unit tests for CheckpointManager initialization
2. Create unit tests for save_checkpoint method
3. Create unit tests for load_checkpoint method
4. Create unit tests for progress tracking methods

INTERFACE:
```python
import unittest
import os
import shutil
import tempfile
from chembl.checkpoint import CheckpointManager

class TestCheckpointManager(unittest.TestCase):
    def setUp(self):
        """Set up a temporary directory for checkpoint files"""
        self.temp_dir = tempfile.mkdtemp()
        self.checkpoint_manager = CheckpointManager(self.temp_dir)
        
    def tearDown(self):
        """Remove temporary directory after tests"""
        shutil.rmtree(self.temp_dir)
        
    def test_initialization(self):
        """Test checkpoint manager initialization"""
        self.assertEqual(self.checkpoint_manager.checkpoint_dir, self.temp_dir)
        self.assertIn("processed_compounds", self.checkpoint_manager.state)
        
    def test_save_load_checkpoint(self):
        """Test saving and loading checkpoints"""
        # Update state
        self.checkpoint_manager.state["total_processed"] = 10
        
        # Save checkpoint
        checkpoint_path = self.checkpoint_manager.save_checkpoint()
        self.assertTrue(os.path.exists(checkpoint_path))
        
        # Create new manager and load checkpoint
        new_manager = CheckpointManager(self.temp_dir)
        success = new_manager.load_checkpoint()
        
        # Verify state was loaded
        self.assertTrue(success)
        self.assertEqual(new_manager.state["total_processed"], 10)
        
    # Additional test methods...
```

VERIFICATION:
- All tests pass successfully
- Coverage is >90% for checkpoint.py
- Checkpoint files are correctly saved and loaded
- State is properly maintained between saves and loads
```

**TASK 3: Worker Implementation Tests**
```
TASK: TEST-FRAMEWORK-3: Worker Implementation Tests

FILE: tests/test_chembl_worker.py
REFERENCE: tests/test_pubchem/test_client.py:50-90

IMPLEMENTATION:
1. Create unit tests for ChEMBLWorker initialization
2. Create unit tests for task processing
3. Create unit tests for error handling
4. Create unit tests for result reporting

INTERFACE:
```python
import unittest
from queue import Queue
from chembl.worker import ChEMBLWorker

class TestChEMBLWorker(unittest.TestCase):
    def setUp(self):
        """Set up queues and worker for testing"""
        self.task_queue = Queue()
        self.result_queue = Queue()
        self.worker = ChEMBLWorker(1, self.task_queue, self.result_queue)
        
    def test_initialization(self):
        """Test worker initialization"""
        self.assertEqual(self.worker.worker_id, 1)
        self.assertEqual(self.worker.task_queue, self.task_queue)
        self.assertEqual(self.worker.result_queue, self.result_queue)
        self.assertFalse(self.worker.running)
        
    def test_process_task(self):
        """Test processing a single task"""
        # Mock the worker's methods to avoid actual API calls
        self.worker.fetch_compound_data = lambda cid: {"molecule_chembl_id": cid}
        self.worker.transform_compound_data = lambda data: {"molecule": {"name": "Test"}, "properties": []}
        self.worker.store_compound_data = lambda data: "test_id"
        
        # Process a task
        result = self.worker.process_task({"compound_id": "CHEMBL25"})
        
        # Verify result
        self.assertEqual(result["status"], "success")
        self.assertEqual(result["compound_id"], "CHEMBL25")
        
    # Additional test methods...
```

VERIFICATION:
- All tests pass successfully
- Coverage is >90% for worker.py
- All worker methods are tested
- Error handling is verified
```

**TASK 4: Integration Tests**
```
TASK: TEST-FRAMEWORK-4: Integration Tests

FILE: tests/test_chembl_integration.py
REFERENCE: tests/test_api_integration.py:10-60

IMPLEMENTATION:
1. Create integration test with mock database
2. Test end-to-end import process with small dataset
3. Test checkpoint creation and resumption
4. Test error handling and recovery

INTERFACE:
```python
import unittest
import os
import tempfile
import shutil
from queue import Queue
from chembl.worker import ChEMBLWorker
from chembl.checkpoint import CheckpointManager

class TestChEMBLIntegration(unittest.TestCase):
    def setUp(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp()
        self.checkpoint_manager = CheckpointManager(self.temp_dir)
        
        # Mock database connection
        self.db_mock = MockDatabase()
        
        # Test compounds
        self.test_compounds = ["CHEMBL25", "CHEMBL1118"]
        
    def tearDown(self):
        """Clean up test environment"""
        shutil.rmtree(self.temp_dir)
        
    def test_end_to_end_import(self):
        """Test end-to-end import process"""
        # Set up worker
        task_queue = Queue()
        result_queue = Queue()
        
        # Create worker with mocked database
        worker = ChEMBLWorker(1, task_queue, result_queue)
        worker.store_compound_data = self.db_mock.store_compound
        
        # Start worker
        worker.start()
        
        # Add tasks
        for compound_id in self.test_compounds:
            task_queue.put({"compound_id": compound_id})
        
        # Process results
        results = []
        for _ in range(len(self.test_compounds)):
            results.append(result_queue.get(timeout=10))
            
        # Stop worker
        worker.stop()
        
        # Verify results
        self.assertEqual(len(results), len(self.test_compounds))
        success_count = sum(1 for r in results if r["status"] == "success")
        self.assertEqual(success_count, len(self.test_compounds))
        
        # Verify database state
        self.assertEqual(len(self.db_mock.compounds), len(self.test_compounds))
        
    # Additional test methods...
    
class MockDatabase:
    """Mock database for testing"""
    def __init__(self):
        self.compounds = {}
        self.properties = {}
        
    def store_compound(self, data):
        """Store compound data and return ID"""
        molecule = data.get("molecule", {})
        molecule_id = molecule.get("chembl_id", f"test_{len(self.compounds)}")
        self.compounds[molecule_id] = molecule
        
        # Store properties
        properties = data.get("properties", [])
        if molecule_id not in self.properties:
            self.properties[molecule_id] = []
        self.properties[molecule_id].extend(properties)
        
        return molecule_id
```

VERIFICATION:
- All integration tests pass successfully
- Test covers complete import workflow
- Mock database correctly stores compounds and properties
- Worker correctly processes all test compounds
```

### STAGE 2: DATA POPULATION

**TASK 5: Reference Compounds Import**
```
TASK: DATA-POP-1: Reference Compounds Import

FILE: import_reference_compounds.py
REFERENCE: run_chembl_import.py:1-100

IMPLEMENTATION:
1. Create script for importing reference compounds
2. Add detailed logging of import process
3. Add verification of imported compounds
4. Generate summary report of import results

INTERFACE:
```python
#!/usr/bin/env python3
"""
Reference Compounds Import Script

This script imports all reference compounds defined in chembl/reference_compounds.py
to establish a baseline of well-known compounds in the database.
"""

import os
import sys
import logging
import argparse
from datetime import datetime
from queue import Queue

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/reference_import.log')
    ]
)

logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Import reference compounds into CryoProtect database")
    parser.add_argument("--workers", type=int, default=2,
                        help="Number of worker threads")
    parser.add_argument("--checkpoint-dir", default="./checkpoints",
                        help="Directory for checkpoint files")
    parser.add_argument("--mode", choices=["full", "resume", "dry-run"], default="resume",
                        help="Execution mode")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose logging")
    return parser.parse_args()

def main():
    """Main execution function"""
    args = parse_arguments()
    
    # Set log level based on verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logger.info(f"Starting reference compounds import in {args.mode} mode")
    
    # Import necessary modules
    try:
        from chembl.worker import ChEMBLWorker
        from chembl.checkpoint import CheckpointManager
        from chembl.reference_compounds import get_reference_compound_ids
    except ImportError as e:
        logger.error(f"Failed to import required modules: {e}")
        return 1
    
    # Get reference compounds
    reference_compounds = get_reference_compound_ids()
    logger.info(f"Found {len(reference_compounds)} reference compounds to import")
    
    # Create checkpoint manager
    checkpoint_manager = CheckpointManager(args.checkpoint_dir, prefix="reference_import")
    
    # Load checkpoint if in resume mode
    if args.mode == "resume":
        checkpoint_loaded = checkpoint_manager.load_checkpoint()
        if checkpoint_loaded:
            logger.info("Resuming from checkpoint")
            # Skip already processed compounds
            processed = checkpoint_manager.state.get("processed_compounds", [])
            reference_compounds = [cid for cid in reference_compounds if cid not in processed]
            logger.info(f"Skipping {len(processed)} already processed compounds")
        else:
            logger.info("No checkpoint found, starting fresh")
    
    # Continued in the next tasks...
```

VERIFICATION:
- Script correctly imports reference compounds
- Checkpointing works properly
- Logging provides clear progress information
- All reference compounds are correctly imported
```

**TASK 6: Reference Compounds Import Execution**
```
TASK: DATA-POP-2: Reference Compounds Import Execution

FILE: import_reference_compounds.py:100-200
REFERENCE: run_chembl_import.py:101-200

IMPLEMENTATION:
1. Set up worker pool for processing reference compounds
2. Implement results collection and processing
3. Add summary report generation
4. Add verification of imported compounds

INTERFACE:
```python
    # ...continuing from previous task
    
    # Set up worker pool
    task_queue = Queue()
    result_queue = Queue()
    
    # Create workers
    workers = []
    for i in range(args.workers):
        worker = ChEMBLWorker(i, task_queue, result_queue)
        workers.append(worker)
    
    # Start workers
    for worker in workers:
        worker.start()
        
    logger.info(f"Started {len(workers)} workers")
    
    # Add dry run flag for dry-run mode
    dry_run = (args.mode == "dry-run")
    
    # Add tasks to queue
    for compound_id in reference_compounds:
        task_queue.put({
            "compound_id": compound_id,
            "dry_run": dry_run,
            "reference": True
        })
    
    logger.info(f"Added {len(reference_compounds)} compounds to task queue")
    
    # Statistics
    stats = {
        "total": len(reference_compounds),
        "processed": 0,
        "success": 0,
        "error": 0,
        "start_time": datetime.now(),
        "errors": {}
    }
    
    # Process results
    try:
        for _ in range(len(reference_compounds)):
            try:
                result = result_queue.get(timeout=300)  # 5 minute timeout per compound
                
                # Process the result
                stats["processed"] += 1
                
                if result.get("status") == "success":
                    stats["success"] += 1
                    
                    # Update checkpoint if available
                    checkpoint_manager.update_progress(
                        result["compound_id"], 
                        success=True
                    )
                else:
                    stats["error"] += 1
                    error_category = result.get("error_category", "UNKNOWN")
                    
                    # Track error by category
                    if error_category not in stats["errors"]:
                        stats["errors"][error_category] = 0
                    stats["errors"][error_category] += 1
                    
                    # Update checkpoint
                    checkpoint_manager.update_progress(
                        result["compound_id"],
                        success=False,
                        error=result.get("error")
                    )
                
                # Display progress
                display_progress(stats)
                
                # Save checkpoint periodically
                if stats["processed"] % 5 == 0:
                    checkpoint_manager.save_checkpoint()
                    
                # Mark task as done
                result_queue.task_done()
                
            except Queue.Empty:
                logger.warning("Timeout waiting for result")
                break
                
    except KeyboardInterrupt:
        logger.info("Import interrupted by user")
    finally:
        # Stop workers
        for worker in workers:
            worker.stop()
        
        # Final progress display
        display_progress(stats, final=True)
        
        # Save final checkpoint
        checkpoint_manager.save_checkpoint()
        
        # Generate summary report
        generate_summary_report(stats, args.mode)
    
    return 0

def display_progress(stats, final=False):
    """Display progress information"""
    elapsed = (datetime.now() - stats["start_time"]).total_seconds()
    processed = stats["processed"]
    total = stats["total"]
    
    if processed == 0:
        return
    
    percent = (processed / total) * 100 if total > 0 else 0
    rate = processed / elapsed if elapsed > 0 else 0
    
    # Calculate ETA
    if rate > 0 and not final:
        eta = (total - processed) / rate
        eta_str = f"ETA: {eta:.1f}s"
    else:
        eta_str = "Complete"
    
    # Create progress bar
    bar_length = 30
    filled_length = int(bar_length * percent / 100)
    bar = '█' * filled_length + '░' * (bar_length - filled_length)
    
    # Print progress
    sys.stdout.write(f"\r[{bar}] {percent:.1f}% | {processed}/{total} | " +
                    f"Success: {stats['success']} | Errors: {stats['error']} | " +
                    f"Rate: {rate:.1f}/s | {eta_str}")
    sys.stdout.flush()
    
    if final:
        sys.stdout.write("\n")
        
        # Print error summary if there are errors
        if stats["error"] > 0:
            sys.stdout.write("\nError summary:\n")
            for category, count in stats["errors"].items():
                sys.stdout.write(f"  {category}: {count}\n")

def generate_summary_report(stats, mode):
    """Generate a summary report of the import process"""
    # Create report directory if it doesn't exist
    os.makedirs("reports", exist_ok=True)
    
    # Generate report filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = f"reports/reference_import_{timestamp}.json"
    
    # Add additional report data
    report_data = {
        "timestamp": datetime.now().isoformat(),
        "mode": mode,
        "stats": stats,
        "elapsed_seconds": (datetime.now() - stats["start_time"]).total_seconds(),
        "success_rate": stats["success"] / stats["total"] if stats["total"] > 0 else 0
    }
    
    # Save report
    with open(report_file, 'w') as f:
        json.dump(report_data, f, indent=2, default=str)
    
    logger.info(f"Summary report saved to {report_file}")
    
    return report_file

if __name__ == "__main__":
    sys.exit(main())
```

VERIFICATION:
- Script correctly processes all reference compounds
- Progress is displayed in real-time
- Summary report contains accurate statistics
- All reference compounds are correctly imported
```

**TASK 7: Data Verification Implementation**
```
TASK: DATA-VER-1: Data Verification Implementation

FILE: verify_chembl_data.py:1-100
REFERENCE: verify_database_integrity.py:50-150

IMPLEMENTATION:
1. Create script for verifying ChEMBL data in database
2. Add checks for reference compounds
3. Add verification of property completeness
4. Generate verification report

INTERFACE:
```python
#!/usr/bin/env python3
"""
ChEMBL Data Verification Script

This script verifies the integrity and completeness of ChEMBL data
imported into the CryoProtect database.
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/chembl_verification.log')
    ]
)

logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Verify ChEMBL data in CryoProtect database")
    parser.add_argument("--project-id", help="Supabase project ID")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--output", default="reports", help="Directory for report output")
    return parser.parse_args()

class ChEMBLVerification:
    """Verification system for ChEMBL data in the database"""
    
    def __init__(self, project_id=None):
        """Initialize verification system with optional project ID"""
        self.project_id = project_id
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "counts": {},
            "reference_compounds": {},
            "property_completeness": {},
            "assessments": {}
        }
        
    def verify_counts(self):
        """Verify total counts of imported data"""
        try:
            # Import necessary tools
            try:
                from use_mcp_tool import execute_sql
            except ImportError:
                from supabase_direct import execute_query as execute_sql
            
            # Query the database for counts
            sql = "SELECT COUNT(*) FROM molecules WHERE data_source = 'ChEMBL'"
            result = execute_sql(sql, self.project_id)
            
            self.results["counts"]["molecules"] = result[0]["count"]
            
            sql = "SELECT COUNT(*) FROM molecular_properties WHERE source = 'ChEMBL'"
            result = execute_sql(sql, self.project_id)
            
            self.results["counts"]["properties"] = result[0]["count"]
            
            # Add assessment
            min_expected = 10  # At least reference compounds
            if self.results["counts"]["molecules"] >= min_expected:
                self.results["assessments"]["counts"] = "SUCCESS"
            else:
                self.results["assessments"]["counts"] = "FAILURE"
                
            logger.info(f"Found {self.results['counts']['molecules']} ChEMBL molecules")
            logger.info(f"Found {self.results['counts']['properties']} ChEMBL properties")
            
            return self.results["counts"]
            
        except Exception as e:
            logger.error(f"Error verifying counts: {e}")
            self.results["assessments"]["counts"] = "ERROR"
            return {}
```

VERIFICATION:
- Script correctly verifies ChEMBL data in database
- Successfully identifies missing reference compounds
- Produces accurate verification report
- Reports property completeness percentages
```

**TASK 8: Reference Compound Verification**
```
TASK: DATA-VER-2: Reference Compound Verification

FILE: verify_chembl_data.py:101-200
REFERENCE: verify_database_integrity.py:151-250

IMPLEMENTATION:
1. Add verification of reference compounds
2. Check that all reference compounds are present
3. Verify property completeness for reference compounds
4. Generate reference compound verification report

INTERFACE:
```python
    def verify_reference_compounds(self):
        """Verify that all reference compounds were imported"""
        try:
            # Import necessary tools
            try:
                from use_mcp_tool import execute_sql
            except ImportError:
                from supabase_direct import execute_query as execute_sql
            
            # Import reference compound list
            from chembl.reference_compounds import get_reference_compounds
            reference_compounds = get_reference_compounds()
            
            logger.info(f"Verifying {len(reference_compounds)} reference compounds")
            
            # Check each reference compound
            present_count = 0
            for compound in reference_compounds:
                chembl_id = compound["chembl_id"]
                
                # Check if compound exists
                sql = f"SELECT id, name FROM molecules WHERE chembl_id = '{chembl_id}'"
                result = execute_sql(sql, self.project_id)
                
                is_present = len(result) > 0
                name_matches = False
                db_name = None
                
                if is_present:
                    present_count += 1
                    db_name = result[0]["name"]
                    # Check if name matches (case insensitive)
                    expected_name = compound["name"]
                    name_matches = expected_name.lower() in db_name.lower()
                
                # Store result
                self.results["reference_compounds"][chembl_id] = {
                    "present": is_present,
                    "expected_name": compound["name"],
                    "actual_name": db_name,
                    "name_matches": name_matches,
                    "role": compound.get("role", "unknown")
                }
                
                # If present, check properties
                if is_present:
                    molecule_id = result[0]["id"]
                    sql = f"SELECT COUNT(*) FROM molecular_properties WHERE molecule_id = '{molecule_id}'"
                    prop_result = execute_sql(sql, self.project_id)
                    property_count = prop_result[0]["count"]
                    
                    self.results["reference_compounds"][chembl_id]["property_count"] = property_count
            
            # Calculate overall success rate
            total = len(reference_compounds)
            success_rate = present_count / total if total > 0 else 0
            
            self.results["reference_compounds"]["summary"] = {
                "total": total,
                "present": present_count,
                "missing": total - present_count,
                "success_rate": success_rate
            }
            
            # Add assessment
            if success_rate >= 0.9:  # At least 90% present
                self.results["assessments"]["reference_compounds"] = "SUCCESS"
            elif success_rate >= 0.7:  # At least 70% present
                self.results["assessments"]["reference_compounds"] = "WARNING"
            else:
                self.results["assessments"]["reference_compounds"] = "FAILURE"
                
            logger.info(f"Reference compound verification: {present_count}/{total} present ({success_rate*100:.1f}%)")
            
            # Log missing compounds
            if present_count < total:
                missing = [c["chembl_id"] for c in reference_compounds 
                         if not self.results["reference_compounds"].get(c["chembl_id"], {}).get("present", False)]
                logger.warning(f"Missing reference compounds: {', '.join(missing)}")
            
            return self.results["reference_compounds"]
            
        except Exception as e:
            logger.error(f"Error verifying reference compounds: {e}")
            self.results["assessments"]["reference_compounds"] = "ERROR"
            return {}
```

VERIFICATION:
- Successfully verifies all reference compounds
- Correctly identifies missing compounds
- Checks property counts for each compound
- Produces accurate verification report
```

**TASK 9: Full ChEMBL Import**
```
TASK: DATA-POP-3: Full ChEMBL Import

FILE: run_full_chembl_import.py
REFERENCE: run_chembl_import.py:1-200

IMPLEMENTATION:
1. Create script for full ChEMBL import
2. Enhance with search query for finding cryoprotectants
3. Add optimization for large-scale import
4. Implement enhanced reporting

INTERFACE:
```python
#!/usr/bin/env python3
"""
Full ChEMBL Import Script

This script imports cryoprotectant-related compounds from ChEMBL
into the CryoProtect database using the worker pool architecture.
"""

import os
import sys
import json
import logging
import argparse
import time
from datetime import datetime
from queue import Queue

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/chembl_import.log')
    ]
)

logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Import ChEMBL data into CryoProtect database")
    parser.add_argument("--limit", type=int, default=1000, 
                        help="Maximum number of compounds to import")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Number of compounds per batch")
    parser.add_argument("--workers", type=int, default=4,
                        help="Number of worker threads")
    parser.add_argument("--checkpoint-dir", default="./checkpoints",
                        help="Directory for checkpoint files")
    parser.add_argument("--mode", choices=["full", "resume", "dry-run"], default="resume",
                        help="Execution mode")
    parser.add_argument("--search-terms", nargs="+", 
                        default=["cryoprotectant", "antifreeze", "cryopreservation", "cell preservation"],
                        help="Search terms for finding compounds")
    parser.add_argument("--include-refs", action="store_true",
                        help="Include reference compounds")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose logging")
    return parser.parse_args()

def main():
    """Main execution function"""
    args = parse_arguments()
    
    # Set log level based on verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logger.info(f"Starting full ChEMBL import in {args.mode} mode")
    logger.info(f"Settings: limit={args.limit}, batch_size={args.batch_size}, workers={args.workers}")
    logger.info(f"Search terms: {args.search_terms}")
    
    # Import necessary modules
    try:
        from chembl.worker import ChEMBLWorker
        from chembl.checkpoint import CheckpointManager
        from chembl.client import ChEMBLClient
        from chembl.reference_compounds import get_reference_compound_ids
    except ImportError as e:
        logger.error(f"Failed to import required modules: {e}")
        return 1
    
    # Create checkpoint manager
    checkpoint_manager = CheckpointManager(args.checkpoint_dir)
    
    # Load checkpoint if in resume mode
    if args.mode == "resume":
        checkpoint_loaded = checkpoint_manager.load_checkpoint()
        if checkpoint_loaded:
            logger.info("Resuming from checkpoint")
        else:
            logger.info("No checkpoint found, starting fresh")
    
    # Initialize client for fetching compound IDs
    client = ChEMBLClient()
    
    # Get compound IDs to process
    compound_ids = []
    
    # Add reference compounds if requested
    if args.include_refs:
        reference_compounds = get_reference_compound_ids()
        logger.info(f"Including {len(reference_compounds)} reference compounds")
        compound_ids.extend(reference_compounds)
    
    # Search for compounds by terms
    for term in args.search_terms:
        logger.info(f"Searching for compounds with term: {term}")
        results = client.search_compounds(term)
        logger.info(f"Found {len(results)} compounds for term '{term}'")
        compound_ids.extend(results)
    
    # Deduplicate
    compound_ids = list(set(compound_ids))
    logger.info(f"Total unique compounds after deduplication: {len(compound_ids)}")
    
    # Apply limit
    if args.limit and len(compound_ids) > args.limit:
        logger.info(f"Limiting to {args.limit} compounds")
        compound_ids = compound_ids[:args.limit]
    
    # Skip already processed compounds if resuming
    if args.mode == "resume" and checkpoint_manager.state.get("processed_compounds"):
        processed = checkpoint_manager.state.get("processed_compounds", [])
        skipped = 0
        new_compound_ids = []
        for cid in compound_ids:
            if cid in processed:
                skipped += 1
            else:
                new_compound_ids.append(cid)
        compound_ids = new_compound_ids
        logger.info(f"Skipping {skipped} already processed compounds")
    
    # Proceed with import...
    # (The rest of the implementation will follow in the next task)
```

VERIFICATION:
- Script correctly searches for cryoprotectant-related compounds
- Properly handles reference compounds inclusion
- Implements checkpointing for resumable operation
- Correctly deduplicates compound IDs
```

**TASK 10: Full ChEMBL Import Execution**
```
TASK: DATA-POP-4: Full ChEMBL Import Execution

FILE: run_full_chembl_import.py:150-300
REFERENCE: run_chembl_import.py:201-400

IMPLEMENTATION:
1. Complete the import script with worker pool setup
2. Add result processing and statistics
3. Implement enhanced progress reporting
4. Add summary report generation

INTERFACE:
```python
    # ...continuing from previous task
    
    # Set up worker pool
    task_queue = Queue()
    result_queue = Queue()
    
    # Create workers
    workers = []
    for i in range(args.workers):
        worker = ChEMBLWorker(i, task_queue, result_queue)
        workers.append(worker)
    
    # Start workers
    for worker in workers:
        worker.start()
        
    logger.info(f"Started {len(workers)} workers")
    
    # Add dry run flag for dry-run mode
    dry_run = (args.mode == "dry-run")
    
    # Add tasks to queue
    for compound_id in compound_ids:
        task_queue.put({
            "compound_id": compound_id,
            "dry_run": dry_run
        })
    
    logger.info(f"Added {len(compound_ids)} compounds to task queue")
    
    # Statistics
    stats = {
        "total": len(compound_ids),
        "processed": 0,
        "success": 0,
        "error": 0,
        "start_time": datetime.now(),
        "errors": {},
        "batch_times": [],
        "search_terms": args.search_terms
    }
    
    # Process results in batches
    batch_start_time = time.time()
    batch_count = 0
    batch_size = args.batch_size
    
    try:
        for i in range(len(compound_ids)):
            try:
                result = result_queue.get(timeout=300)  # 5 minute timeout per compound
                
                # Process the result
                stats["processed"] += 1
                
                if result.get("status") == "success":
                    stats["success"] += 1
                    
                    # Update checkpoint if available
                    checkpoint_manager.update_progress(
                        result["compound_id"], 
                        success=True
                    )
                else:
                    stats["error"] += 1
                    error_category = result.get("error_category", "UNKNOWN")
                    
                    # Track error by category
                    if error_category not in stats["errors"]:
                        stats["errors"][error_category] = 0
                    stats["errors"][error_category] += 1
                    
                    # Update checkpoint
                    checkpoint_manager.update_progress(
                        result["compound_id"],
                        success=False,
                        error=result.get("error")
                    )
                
                # Display progress
                display_progress(stats)
                
                # Batch processing
                if (i + 1) % batch_size == 0 or i == len(compound_ids) - 1:
                    batch_end_time = time.time()
                    batch_duration = batch_end_time - batch_start_time
                    stats["batch_times"].append(batch_duration)
                    
                    # Calculate batch statistics
                    compounds_per_second = batch_size / batch_duration if batch_duration > 0 else 0
                    logger.info(f"Batch {batch_count} completed: {batch_size} compounds in {batch_duration:.2f}s " +
                              f"({compounds_per_second:.2f} compounds/s)")
                    
                    # Save checkpoint after each batch
                    checkpoint_manager.save_checkpoint()
                    
                    # Reset for next batch
                    batch_start_time = time.time()
                    batch_count += 1
                    
                # Mark task as done
                result_queue.task_done()
                
            except Queue.Empty:
                logger.warning("Timeout waiting for result")
                break
                
    except KeyboardInterrupt:
        logger.info("Import interrupted by user")
    finally:
        # Stop workers
        for worker in workers:
            worker.stop()
        
        # Final progress display
        display_progress(stats, final=True)
        
        # Save final checkpoint
        checkpoint_manager.save_checkpoint()
        
        # Generate summary report
        report_file = generate_summary_report(stats, args.mode)
        
        # Run verification if not in dry-run mode
        if args.mode != "dry-run" and stats["success"] > 0:
            try:
                from verify_chembl_data import verify_chembl_import
                logger.info("Running verification...")
                verify_chembl_import()
            except ImportError:
                logger.warning("Verification module not available")
    
    return 0

def display_progress(stats, final=False):
    """Display progress information"""
    # Implementation same as previous task

def generate_summary_report(stats, mode):
    """Generate a summary report of the import process"""
    # Implementation same as previous task with additional metrics
    
    # Create report directory if it doesn't exist
    os.makedirs("reports", exist_ok=True)
    
    # Generate report filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = f"reports/chembl_import_{timestamp}.json"
    
    # Calculate additional metrics
    elapsed_seconds = (datetime.now() - stats["start_time"]).total_seconds()
    compounds_per_second = stats["processed"] / elapsed_seconds if elapsed_seconds > 0 else 0
    average_batch_time = sum(stats["batch_times"]) / len(stats["batch_times"]) if stats["batch_times"] else 0
    
    # Add additional report data
    report_data = {
        "timestamp": datetime.now().isoformat(),
        "mode": mode,
        "stats": stats,
        "elapsed_seconds": elapsed_seconds,
        "compounds_per_second": compounds_per_second,
        "average_batch_time": average_batch_time,
        "success_rate": stats["success"] / stats["total"] if stats["total"] > 0 else 0,
        "search_terms": stats["search_terms"]
    }
    
    # Save report
    with open(report_file, 'w') as f:
        json.dump(report_data, f, indent=2, default=str)
    
    logger.info(f"Summary report saved to {report_file}")
    
    return report_file

if __name__ == "__main__":
    sys.exit(main())
```

VERIFICATION:
- Script successfully processes all compounds
- Checkpointing correctly allows resuming interrupted operations
- Progress reporting provides real-time status
- Summary report contains all relevant statistics
```

### STAGE 3: UI OPTIMIZATION

**TASK 11: Implement Pagination**
```
TASK: UI-OPT-1: API Pagination Implementation

FILE: api/resources.py:300-400
REFERENCE: api/resources.py:100-200

IMPLEMENTATION:
1. Add pagination support to molecule listing endpoint
2. Implement limit and offset parameters
3. Add total count to response
4. Add next/prev page links

INTERFACE:
```python
@app.route('/api/molecules', methods=['GET'])
@jwt_required
def get_molecules():
    """
    Get paginated list of molecules.
    
    Query parameters:
    - limit: Maximum number of molecules to return (default: 20, max: 100)
    - offset: Number of molecules to skip (default: 0)
    - search: Optional search term for filtering
    - sort: Field to sort by (default: 'name')
    - order: Sort order ('asc' or 'desc', default: 'asc')
    """
    try:
        # Parse pagination parameters
        limit = min(int(request.args.get('limit', 20)), 100)
        offset = int(request.args.get('offset', 0))
        search = request.args.get('search', '')
        sort = request.args.get('sort', 'name')
        order = request.args.get('order', 'asc')
        
        # Validate sort field
        valid_sort_fields = ['name', 'molecular_weight', 'formula', 'created_at']
        if sort not in valid_sort_fields:
            sort = 'name'
            
        # Validate order direction
        if order not in ['asc', 'desc']:
            order = 'asc'
        
        # Build SQL query with pagination
        base_query = "FROM molecules"
        where_clause = ""
        
        # Add search filter if provided
        if search:
            where_clause = f" WHERE name ILIKE '%{search}%' OR formula ILIKE '%{search}%'"
        
        # Count total matches
        count_sql = f"SELECT COUNT(*) {base_query}{where_clause}"
        count_result = execute_sql(count_sql, get_project_id())
        total_count = count_result[0]['count']
        
        # Get paginated results
        order_direction = "ASC" if order == 'asc' else "DESC"
        result_sql = f"""
            SELECT id, name, formula, molecular_weight, smiles, inchi_key, data_source, created_at
            {base_query}{where_clause}
            ORDER BY {sort} {order_direction}
            LIMIT {limit} OFFSET {offset}
        """
        results = execute_sql(result_sql, get_project_id())
        
        # Calculate pagination metadata
        page = offset // limit + 1
        total_pages = (total_count + limit - 1) // limit
        
        # Add next/prev page links
        base_url = request.base_url
        links = {}
        
        if offset + limit < total_count:
            next_offset = offset + limit
            links['next'] = f"{base_url}?limit={limit}&offset={next_offset}&search={search}&sort={sort}&order={order}"
            
        if offset > 0:
            prev_offset = max(0, offset - limit)
            links['prev'] = f"{base_url}?limit={limit}&offset={prev_offset}&search={search}&sort={sort}&order={order}"
        
        # Return paginated response
        return jsonify({
            'data': results,
            'pagination': {
                'count': total_count,
                'page': page,
                'total_pages': total_pages,
                'limit': limit,
                'offset': offset
            },
            'links': links
        })
        
    except Exception as e:
        logger.error(f"Error in get_molecules: {str(e)}")
        return jsonify({'error': str(e)}), 500
```

VERIFICATION:
- Endpoint correctly implements pagination
- Search filtering works as expected
- Sorting works for all valid fields
- Response includes proper pagination metadata
```

**TASK 12: UI Pagination Integration**
```
TASK: UI-OPT-2: Frontend Pagination Integration

FILE: static/js/pagination.js
REFERENCE: static/js/app.js:50-150

IMPLEMENTATION:
1. Create client-side pagination component
2. Add pagination controls to molecule listing
3. Implement page navigation functionality
4. Add loading indicators for pagination

INTERFACE:
```javascript
/**
 * Pagination component for molecule listings
 * 
 * Provides pagination controls and manages API requests for paginated data.
 */

class MoleculePagination {
  /**
   * Initialize pagination component
   * @param {string} containerId - ID of container element
   * @param {string} apiEndpoint - API endpoint for fetching data
   * @param {Function} renderCallback - Function to render items
   * @param {Object} options - Additional options
   */
  constructor(containerId, apiEndpoint, renderCallback, options = {}) {
    this.container = document.getElementById(containerId);
    if (!this.container) {
      console.error(`Container element #${containerId} not found`);
      return;
    }
    
    this.apiEndpoint = apiEndpoint;
    this.renderCallback = renderCallback;
    
    // Default options
    this.options = {
      limit: 20,
      initialPage: 1,
      searchEnabled: true,
      sortEnabled: true,
      defaultSort: 'name',
      defaultOrder: 'asc',
      ...options
    };
    
    // State
    this.currentPage = this.options.initialPage;
    this.totalPages = 1;
    this.isLoading = false;
    this.searchTerm = '';
    this.sortField = this.options.defaultSort;
    this.sortOrder = this.options.defaultOrder;
    
    // Initialize component
    this.initialize();
  }
  
  /**
   * Initialize pagination component
   */
  initialize() {
    // Create component elements
    this.createElements();
    
    // Add event listeners
    this.setupEventListeners();
    
    // Load initial data
    this.loadPage(this.currentPage);
  }
  
  /**
   * Create pagination UI elements
   */
  createElements() {
    // Create container structure
    this.container.innerHTML = `
      <div class="pagination-controls">
        <div class="search-container ${!this.options.searchEnabled ? 'hidden' : ''}">
          <input type="text" id="search-input" placeholder="Search..." class="search-input">
          <button id="search-button" class="search-button">Search</button>
        </div>
        
        <div class="sort-container ${!this.options.sortEnabled ? 'hidden' : ''}">
          <select id="sort-select" class="sort-select">
            <option value="name">Name</option>
            <option value="molecular_weight">Molecular Weight</option>
            <option value="formula">Formula</option>
            <option value="created_at">Date Added</option>
          </select>
          <button id="sort-order" class="sort-order" data-order="${this.sortOrder}">
            ${this.sortOrder === 'asc' ? '↑' : '↓'}
          </button>
        </div>
      </div>
      
      <div id="content-container" class="content-container">
        <div id="loading-indicator" class="loading-indicator hidden">
          <div class="spinner"></div>
          <p>Loading...</p>
        </div>
        <div id="content-list" class="content-list"></div>
      </div>
      
      <div class="pagination-nav">
        <button id="prev-page" class="pagination-button" disabled>&lt; Previous</button>
        <span id="page-indicator" class="page-indicator">Page 1 of 1</span>
        <button id="next-page" class="pagination-button" disabled>Next &gt;</button>
      </div>
    `;
    
    // Get references to elements
    this.searchInput = document.getElementById('search-input');
    this.searchButton = document.getElementById('search-button');
    this.sortSelect = document.getElementById('sort-select');
    this.sortOrderButton = document.getElementById('sort-order');
    this.contentList = document.getElementById('content-list');
    this.loadingIndicator = document.getElementById('loading-indicator');
    this.prevButton = document.getElementById('prev-page');
    this.nextButton = document.getElementById('next-page');
    this.pageIndicator = document.getElementById('page-indicator');
    
    // Set initial sort field
    if (this.sortSelect) {
      this.sortSelect.value = this.sortField;
    }
  }
  
  /**
   * Set up event listeners for pagination controls
   */
  setupEventListeners() {
    // Search functionality
    if (this.searchButton) {
      this.searchButton.addEventListener('click', () => {
        this.searchTerm = this.searchInput.value;
        this.currentPage = 1;
        this.loadPage(this.currentPage);
      });
    }
    
    if (this.searchInput) {
      this.searchInput.addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
          this.searchTerm = this.searchInput.value;
          this.currentPage = 1;
          this.loadPage(this.currentPage);
        }
      });
    }
    
    // Sort functionality
    if (this.sortSelect) {
      this.sortSelect.addEventListener('change', () => {
        this.sortField = this.sortSelect.value;
        this.loadPage(this.currentPage);
      });
    }
    
    if (this.sortOrderButton) {
      this.sortOrderButton.addEventListener('click', () => {
        this.sortOrder = this.sortOrder === 'asc' ? 'desc' : 'asc';
        this.sortOrderButton.setAttribute('data-order', this.sortOrder);
        this.sortOrderButton.textContent = this.sortOrder === 'asc' ? '↑' : '↓';
        this.loadPage(this.currentPage);
      });
    }
    
    // Pagination navigation
    if (this.prevButton) {
      this.prevButton.addEventListener('click', () => {
        if (this.currentPage > 1) {
          this.loadPage(this.currentPage - 1);
        }
      });
    }
    
    if (this.nextButton) {
      this.nextButton.addEventListener('click', () => {
        if (this.currentPage < this.totalPages) {
          this.loadPage(this.currentPage + 1);
        }
      });
    }
  }
  
  /**
   * Load data for specified page
   * @param {number} page - Page number to load
   */
  loadPage(page) {
    // Prevent loading while already in progress
    if (this.isLoading) {
      return;
    }
    
    this.isLoading = true;
    this.showLoading(true);
    
    // Calculate offset
    const offset = (page - 1) * this.options.limit;
    
    // Build API URL with parameters
    let url = `${this.apiEndpoint}?limit=${this.options.limit}&offset=${offset}`;
    
    // Add search parameter if provided
    if (this.searchTerm) {
      url += `&search=${encodeURIComponent(this.searchTerm)}`;
    }
    
    // Add sort parameters
    url += `&sort=${this.sortField}&order=${this.sortOrder}`;
    
    // Fetch data from API
    fetch(url, {
      method: 'GET',
      headers: {
        'Authorization': `Bearer ${getToken()}`,
        'Content-Type': 'application/json'
      }
    })
    .then(response => {
      if (!response.ok) {
        throw new Error('Network response was not ok');
      }
      return response.json();
    })
    .then(data => {
      // Update state
      this.currentPage = page;
      this.totalPages = data.pagination.total_pages;
      
      // Update UI
      this.updatePageIndicator();
      this.updateNavigationButtons();
      
      // Render data
      if (typeof this.renderCallback === 'function') {
        this.renderCallback(data.data, this.contentList);
      } else {
        console.error('No render callback provided');
      }
    })
    .catch(error => {
      console.error('Error fetching data:', error);
      this.contentList.innerHTML = `
        <div class="error-message">
          <p>Failed to load data. Please try again.</p>
        </div>
      `;
    })
    .finally(() => {
      this.isLoading = false;
      this.showLoading(false);
    });
  }
  
  /**
   * Show or hide loading indicator
   * @param {boolean} show - Whether to show loading indicator
   */
  showLoading(show) {
    if (this.loadingIndicator) {
      if (show) {
        this.loadingIndicator.classList.remove('hidden');
      } else {
        this.loadingIndicator.classList.add('hidden');
      }
    }
  }
  
  /**
   * Update page indicator text
   */
  updatePageIndicator() {
    if (this.pageIndicator) {
      this.pageIndicator.textContent = `Page ${this.currentPage} of ${this.totalPages}`;
    }
  }
  
  /**
   * Update navigation button states
   */
  updateNavigationButtons() {
    if (this.prevButton) {
      this.prevButton.disabled = this.currentPage <= 1;
    }
    
    if (this.nextButton) {
      this.nextButton.disabled = this.currentPage >= this.totalPages;
    }
  }
}

// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
  module.exports = { MoleculePagination };
}
```

VERIFICATION:
- Component correctly implements pagination UI
- Search and sorting functionality works as expected
- Loading indicators show during API requests
- Navigation buttons update based on current page
```

### STAGE 4: VERIFICATION & IMPROVEMENT

**TASK 13: Data Quality Report Generator**
```
TASK: DATA-VER-3: Data Quality Report Generator

FILE: generate_data_quality_report.py
REFERENCE: verify_chembl_data.py:1-100

IMPLEMENTATION:
1. Create script for generating comprehensive data quality report
2. Add property distribution analysis
3. Add property completeness checks
4. Generate visualizations for property distributions

INTERFACE:
```python
#!/usr/bin/env python3
"""
Data Quality Report Generator

This script generates a comprehensive report on the quality and completeness
of data in the CryoProtect database, with a focus on imported ChEMBL data.
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/data_quality_report.log')
    ]
)

logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Generate data quality report")
    parser.add_argument("--project-id", help="Supabase project ID")
    parser.add_argument("--output", default="reports", help="Output directory for report")
    parser.add_argument("--format", choices=["json", "html", "both"], default="both",
                       help="Output format for report")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

class DataQualityReport:
    """Data quality report generator"""
    
    def __init__(self, project_id=None, output_dir="reports"):
        """Initialize report generator"""
        self.project_id = project_id
        self.output_dir = output_dir
        self.report_data = {
            "timestamp": datetime.now().isoformat(),
            "database_stats": {},
            "property_distributions": {},
            "property_completeness": {},
            "assessments": {}
        }
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Create plots directory if it doesn't exist
        self.plots_dir = os.path.join(output_dir, "plots")
        os.makedirs(self.plots_dir, exist_ok=True)
        
    def collect_database_stats(self):
        """Collect basic database statistics"""
        try:
            # Import database tools
            try:
                from use_mcp_tool import execute_sql
            except ImportError:
                from supabase_direct import execute_query as execute_sql
            
            # Get molecule counts by source
            sql = """
                SELECT data_source, COUNT(*) as count 
                FROM molecules 
                GROUP BY data_source
            """
            result = execute_sql(sql, self.project_id)
            
            # Process results
            source_counts = {}
            total_molecules = 0
            
            for row in result:
                source = row["data_source"] or "Unknown"
                count = row["count"]
                source_counts[source] = count
                total_molecules += count
            
            # Get property counts by source
            sql = """
                SELECT source, COUNT(*) as count
                FROM molecular_properties
                GROUP BY source
            """
            result = execute_sql(sql, self.project_id)
            
            # Process results
            property_counts = {}
            total_properties = 0
            
            for row in result:
                source = row["source"] or "Unknown"
                count = row["count"]
                property_counts[source] = count
                total_properties += count
            
            # Get average properties per molecule
            if total_molecules > 0:
                avg_properties = total_properties / total_molecules
            else:
                avg_properties = 0
                
            # Store in report
            self.report_data["database_stats"] = {
                "total_molecules": total_molecules,
                "total_properties": total_properties,
                "average_properties_per_molecule": avg_properties,
                "molecule_counts_by_source": source_counts,
                "property_counts_by_source": property_counts
            }
            
            # Add assessment
            min_expected = 100  # At least this many molecules expected
            if total_molecules >= min_expected:
                self.report_data["assessments"]["database_size"] = "SUCCESS"
            elif total_molecules >= 10:
                self.report_data["assessments"]["database_size"] = "WARNING"
            else:
                self.report_data["assessments"]["database_size"] = "FAILURE"
                
            logger.info(f"Database contains {total_molecules} molecules and {total_properties} properties")
            
            return self.report_data["database_stats"]
            
        except Exception as e:
            logger.error(f"Error collecting database stats: {e}")
            self.report_data["assessments"]["database_size"] = "ERROR"
            return {}
```

VERIFICATION:
- Script successfully collects database statistics
- Generates quality report with property distributions
- Creates visualizations for data quality assessment
- Provides accurate data completeness metrics
```

**TASK 14: MoleculeJS Visualization Enhancement**
```
TASK: UI-VIS-1: Molecule Viewer Enhancement

FILE: static/js/molecular-viewer.js
REFERENCE: static/js/rdkit.js:50-150

IMPLEMENTATION:
1. Enhance molecule visualization component
2. Add rotation and zooming controls
3. Implement measurement tool for bond distances
4. Add property visualization overlay

INTERFACE:
```javascript
/**
 * Enhanced Molecular Viewer Component
 * 
 * Provides interactive 3D visualization of molecules with rotation, zooming,
 * measurement tools, and property visualization.
 */

class EnhancedMoleculeViewer {
  /**
   * Initialize the molecule viewer
   * @param {string} containerId - ID of container element
   * @param {Object} options - Configuration options
   */
  constructor(containerId, options = {}) {
    this.container = document.getElementById(containerId);
    if (!this.container) {
      console.error(`Container element #${containerId} not found`);
      return;
    }
    
    // Default options
    this.options = {
      width: 400,
      height: 300,
      backgroundColor: 0xf0f0f0,
      atomLabelsEnabled: false,
      measurementEnabled: false,
      defaultRepresentation: 'balls-and-sticks',
      ...options
    };
    
    // State
    this.molecule = null;
    this.measurePoints = [];
    this.renderer = null;
    this.scene = null;
    this.camera = null;
    this.controls = null;
    this.raycaster = null;
    this.hoveredAtom = null;
    
    // Initialize viewer
    this.initialize();
  }
  
  /**
   * Initialize the 3D viewer
   */
  initialize() {
    // Create viewer elements
    this.createElements();
    
    // Set up Three.js components
    this.setupThreeJS();
    
    // Add event listeners
    this.setupEventListeners();
    
    // Start animation loop
    this.animate();
  }
  
  /**
   * Create UI elements for the viewer
   */
  createElements() {
    // Create container structure
    this.container.innerHTML = `
      <div class="molecule-viewer-container" style="position: relative; width: ${this.options.width}px; height: ${this.options.height}px;">
        <div class="canvas-container" style="width: 100%; height: 100%;"></div>
        
        <div class="viewer-controls" style="position: absolute; bottom: 10px; left: 10px;">
          <button id="rotate-btn" class="control-btn" title="Toggle Rotation">🔄</button>
          <button id="labels-btn" class="control-btn" title="Toggle Atom Labels">Aa</button>
          <button id="measure-btn" class="control-btn" title="Toggle Measurement Tool">📏</button>
          <button id="reset-btn" class="control-btn" title="Reset View">⟲</button>
        </div>
        
        <div class="viewer-info" style="position: absolute; top: 10px; right: 10px; background: rgba(255,255,255,0.7); padding: 5px; border-radius: 3px;">
          <div id="atom-info" class="atom-info"></div>
          <div id="measure-info" class="measure-info"></div>
        </div>
        
        <div class="representation-controls" style="position: absolute; top: 10px; left: 10px;">
          <select id="representation-select" class="representation-select">
            <option value="balls-and-sticks">Balls and Sticks</option>
            <option value="sticks">Sticks</option>
            <option value="spacefill">Spacefill</option>
            <option value="wireframe">Wireframe</option>
          </select>
        </div>
      </div>
    `;
    
    // Get references to elements
    this.canvasContainer = this.container.querySelector('.canvas-container');
    this.rotateBtnElement = document.getElementById('rotate-btn');
    this.labelsBtnElement = document.getElementById('labels-btn');
    this.measureBtnElement = document.getElementById('measure-btn');
    this.resetBtnElement = document.getElementById('reset-btn');
    this.atomInfoElement = document.getElementById('atom-info');
    this.measureInfoElement = document.getElementById('measure-info');
    this.representationSelect = document.getElementById('representation-select');
    
    // Set initial representation
    if (this.representationSelect) {
      this.representationSelect.value = this.options.defaultRepresentation;
    }
  }
  
  /**
   * Set up Three.js components for 3D rendering
   */
  setupThreeJS() {
    // Create scene
    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(this.options.backgroundColor);
    
    // Create camera
    this.camera = new THREE.PerspectiveCamera(
      75, this.options.width / this.options.height, 0.1, 1000
    );
    this.camera.position.z = 5;
    
    // Create renderer
    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setSize(this.options.width, this.options.height);
    this.canvasContainer.appendChild(this.renderer.domElement);
    
    // Create controls
    this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.25;
    
    // Create raycaster for interaction
    this.raycaster = new THREE.Raycaster();
    this.mouse = new THREE.Vector2();
    
    // Add ambient light
    const ambientLight = new THREE.AmbientLight(0x404040, 1.5);
    this.scene.add(ambientLight);
    
    // Add directional light
    const directionalLight = new THREE.DirectionalLight(0xffffff, 1);
    directionalLight.position.set(1, 1, 1).normalize();
    this.scene.add(directionalLight);
  }
  
  /**
   * Set up event listeners for controls
   */
  setupEventListeners() {
    // Window resize
    window.addEventListener('resize', () => {
      this.updateSize();
    });
    
    // Mouse events for atom picking and measurement
    this.renderer.domElement.addEventListener('mousemove', (event) => {
      this.handleMouseMove(event);
    });
    
    this.renderer.domElement.addEventListener('click', (event) => {
      this.handleMouseClick(event);
    });
    
    // Control buttons
    if (this.rotateBtnElement) {
      this.rotateBtnElement.addEventListener('click', () => {
        this.toggleAutoRotation();
      });
    }
    
    if (this.labelsBtnElement) {
      this.labelsBtnElement.addEventListener('click', () => {
        this.toggleAtomLabels();
      });
    }
    
    if (this.measureBtnElement) {
      this.measureBtnElement.addEventListener('click', () => {
        this.toggleMeasurementMode();
      });
    }
    
    if (this.resetBtnElement) {
      this.resetBtnElement.addEventListener('click', () => {
        this.resetView();
      });
    }
    
    // Representation select
    if (this.representationSelect) {
      this.representationSelect.addEventListener('change', () => {
        this.changeRepresentation(this.representationSelect.value);
      });
    }
  }
  
  /**
   * Load molecule from SMILES string
   * @param {string} smiles - SMILES representation of molecule
   */
  loadMoleculeFromSmiles(smiles) {
    try {
      // Clear any existing molecule
      this.clearMolecule();
      
      // Parse SMILES using RDKit.js
      const mol = RDKit.Molecule.fromSmiles(smiles);
      if (!mol) {
        throw new Error("Failed to parse SMILES string");
      }
      
      this.molecule = mol;
      
      // Add 3D coordinates
      mol.addHs();
      mol.embed3D();
      
      // Create 3D representation
      this.createMoleculeRepresentation(this.options.defaultRepresentation);
      
      // Center and zoom to fit
      this.centerAndZoomToFit();
      
      return true;
    } catch (error) {
      console.error("Error loading molecule:", error);
      return false;
    }
  }
  
  /**
   * Create 3D representation of molecule
   * @param {string} style - Representation style
   */
  createMoleculeRepresentation(style) {
    // Implementation of different molecule representation styles
    // This would use Three.js to create the actual 3D objects
    
    // Clear any existing representation
    this.clearMoleculeRepresentation();
    
    // Create group for molecule
    this.moleculeGroup = new THREE.Group();
    this.scene.add(this.moleculeGroup);
    
    // Create atom and bond representations based on style
    switch (style) {
      case 'balls-and-sticks':
        this.createBallsAndSticksRepresentation();
        break;
      case 'sticks':
        this.createSticksRepresentation();
        break;
      case 'spacefill':
        this.createSpacefillRepresentation();
        break;
      case 'wireframe':
        this.createWireframeRepresentation();
        break;
      default:
        this.createBallsAndSticksRepresentation();
    }
    
    // Add atom labels if enabled
    if (this.options.atomLabelsEnabled) {
      this.addAtomLabels();
    }
  }
  
  /**
   * Animation loop
   */
  animate() {
    requestAnimationFrame(() => this.animate());
    
    // Update controls
    if (this.controls) {
      this.controls.update();
    }
    
    // Render scene
    if (this.renderer && this.scene && this.camera) {
      this.renderer.render(this.scene, this.camera);
    }
  }
  
  // Additional methods for molecule manipulation and interaction
  // (Implementation details omitted for brevity)
}

// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
  module.exports = { EnhancedMoleculeViewer };
}
```

VERIFICATION:
- Component successfully visualizes molecular structures
- Interactive controls work for rotation and zooming
- Measurement tool accurately displays bond distances
- Different representation styles render correctly
```

## LINEAR EXECUTION PATH

Execute the tasks in this exact order:

1. **TEST-FRAMEWORK-1**: Error Handler Unit Tests
2. **TEST-FRAMEWORK-2**: Checkpoint System Unit Tests
3. **TEST-FRAMEWORK-3**: Worker Implementation Tests
4. **TEST-FRAMEWORK-4**: Integration Tests
5. **DATA-POP-1**: Reference Compounds Import
6. **DATA-POP-2**: Reference Compounds Import Execution
7. **DATA-VER-1**: Data Verification Implementation
8. **DATA-VER-2**: Reference Compound Verification
9. **DATA-POP-3**: Full ChEMBL Import
10. **DATA-POP-4**: Full ChEMBL Import Execution
11. **UI-OPT-1**: API Pagination Implementation
12. **UI-OPT-2**: Frontend Pagination Integration
13. **DATA-VER-3**: Data Quality Report Generator
14. **UI-VIS-1**: Molecule Viewer Enhancement

## SUCCESS CRITERIA

The linear execution will be considered successful when:

1. All tests pass for the ChEMBL integration framework
2. Reference compounds are successfully imported and verified
3. At least 1,000 ChEMBL compounds are imported into the database
4. UI enhancements improve user experience with pagination and visualization
5. Data quality reports confirm the integrity and completeness of the data

Each task builds upon the previous ones, creating a cohesive implementation path that maximizes efficiency while maintaining clear dependencies and boundaries.