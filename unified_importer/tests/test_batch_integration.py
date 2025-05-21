"""
Tests for the batch integration module.

This module tests the integration of enhanced batch processing
with the unified molecular importer, including database operations
and transformation pipelines.
"""

import unittest
import asyncio
import logging
import time
from unittest.mock import MagicMock, patch
from typing import Dict, List, Any

# Import modules to test
from unified_importer.core.batch_integration import (
    BatchProcessingManager,
    batch_molecule_transform,
    batch_property_transform,
    batch_db_insert
)

# Create a mock transformer
class MockTransformer:
    """Mock transformer for testing."""
    
    def __init__(self):
        """Initialize the mock transformer."""
        self.logger = logging.getLogger(__name__)
        self.transform_count = 0
        self.error_rate = 0.05  # 5% error rate
    
    def transform(self, item: Dict[str, Any]) -> Dict[str, Any]:
        """Transform a single item."""
        self.transform_count += 1
        
        # Simulate occasional failures
        if self.transform_count % 20 == 0:
            raise ValueError("Simulated transform error")
        
        # Transform the item
        result = item.copy()
        result["transformed"] = True
        result["transform_id"] = self.transform_count
        
        return result


# Create a mock database
class MockDatabase:
    """Mock database for testing."""
    
    def __init__(self):
        """Initialize the mock database."""
        self.logger = logging.getLogger(__name__)
        self.transactions = {}
        self.data = {}
        self.insert_count = 0
        self.error_rate = 0.05  # 5% error rate
    
    @contextlib.asynccontextmanager
    async def transaction_async(self):
        """Mock async transaction context manager."""
        tx_id = f"tx-{time.time()}"
        self.transactions[tx_id] = {
            "started_at": time.time(),
            "operations": 0
        }
        try:
            yield tx_id
            # Commit on exit
            self.transactions[tx_id]["committed"] = True
        except Exception as e:
            # Rollback on error
            self.transactions[tx_id]["rolled_back"] = True
            self.logger.error(f"Transaction {tx_id} rolled back: {str(e)}")
            raise
        finally:
            self.transactions[tx_id]["ended_at"] = time.time()
    
    async def batch_insert(
        self,
        table: str,
        records: List[Dict[str, Any]],
        transaction_id: str = None
    ) -> List[str]:
        """Mock batch insert operation."""
        # Track operation count
        if transaction_id in self.transactions:
            self.transactions[transaction_id]["operations"] += 1
        
        # Initialize table if needed
        if table not in self.data:
            self.data[table] = {}
        
        # Insert records
        inserted_ids = []
        for record in records:
            self.insert_count += 1
            
            # Simulate occasional failures
            if self.insert_count % 20 == 0:
                continue
            
            # Generate ID if not present
            record_id = record.get("id", f"id-{self.insert_count}")
            
            # Store record
            self.data[table][record_id] = record
            inserted_ids.append(record_id)
        
        return inserted_ids


# Import contextlib for async context manager mock
import contextlib


class TestBatchProcessingManager(unittest.TestCase):
    """Tests for the BatchProcessingManager class."""
    
    def setUp(self):
        """Set up the test case."""
        # Configure logging
        logging.basicConfig(level=logging.DEBUG)
        
        # Create mock database
        self.db = MockDatabase()
        
        # Create batch processing manager
        self.manager = BatchProcessingManager(
            db=self.db,
            max_workers=2,
            adaptive_batch_sizing=True,
            parallel_processing=True,
            logger=logging.getLogger(__name__)
        )
    
    def test_get_processor(self):
        """Test getting a batch processor."""
        # Create a function
        def process_func(items):
            return items
        
        # Get a processor
        processor = self.manager.get_processor(
            operation="test",
            data_source="mock",
            process_func=process_func,
            initial_batch_size=50
        )
        
        # Should create and cache the processor
        self.assertIsNotNone(processor)
        self.assertEqual(self.manager.processors["test_mock"], processor)
        
        # Getting again should return the same processor
        processor2 = self.manager.get_processor(
            operation="test",
            data_source="mock",
            process_func=process_func
        )
        self.assertIs(processor, processor2)
    
    def test_process_items(self):
        """Test processing a list of items."""
        # Create items
        items = [{"id": f"item-{i}", "value": i} for i in range(100)]
        
        # Create a processing function
        def process_func(batch):
            return [{"id": item["id"], "processed": True} for item in batch]
        
        # Process items
        results = self.manager.process_items(
            items=items,
            operation="test",
            data_source="mock",
            process_func=process_func,
            initial_batch_size=25
        )
        
        # Check results
        self.assertEqual(len(results), 100)
        self.assertTrue(all(item["processed"] for item in results))
        
        # Check stats
        stats = self.manager.get_stats()
        self.assertEqual(stats["total_items_processed"], 100)
        self.assertGreater(stats["total_batches"], 1)
        self.assertIn("test", stats["operations"])
        self.assertEqual(stats["operations"]["test"]["items_processed"], 100)
    
    def test_async_processing(self):
        """Test asynchronous processing."""
        
        async def run_async_test():
            # Create items
            items = [{"id": f"item-{i}", "value": i} for i in range(100)]
            
            # Create a processing function
            def process_func(batch):
                return [{"id": item["id"], "processed": True} for item in batch]
            
            # Process items asynchronously
            results = await self.manager.process_items_async(
                items=items,
                operation="async_test",
                data_source="mock",
                process_func=process_func,
                initial_batch_size=25
            )
            
            # Check results
            self.assertEqual(len(results), 100)
            self.assertTrue(all(item["processed"] for item in results))
            
            # Check stats
            stats = self.manager.get_stats()
            self.assertEqual(stats["total_items_processed"], 100)
            self.assertGreater(stats["total_batches"], 1)
            self.assertIn("async_test", stats["operations"])
            self.assertEqual(stats["operations"]["async_test"]["items_processed"], 100)
        
        # Run the async test
        asyncio.run(run_async_test())


class TestBatchTransformFunctions(unittest.TestCase):
    """Tests for the batch transformation functions."""
    
    def setUp(self):
        """Set up the test case."""
        # Configure logging
        logging.basicConfig(level=logging.DEBUG)
        
        # Create mock transformer
        self.transformer = MockTransformer()
    
    def test_batch_molecule_transform(self):
        """Test batch molecule transformation."""
        # Create mock molecules
        molecules = [{"id": f"mol-{i}", "smiles": f"C{i}"} for i in range(100)]
        
        # Transform batch
        results = batch_molecule_transform(self.transformer, molecules)
        
        # Check results (some may be filtered due to errors)
        self.assertGreaterEqual(len(results), 95)
        self.assertTrue(all(mol["transformed"] for mol in results))
    
    def test_batch_property_transform(self):
        """Test batch property transformation."""
        # Create mock properties
        properties = [{"id": f"prop-{i}", "name": f"Property {i}", "value": i} for i in range(100)]
        
        # Transform batch
        results = batch_property_transform(self.transformer, properties)
        
        # Check results (some may be filtered due to errors)
        self.assertGreaterEqual(len(results), 95)
        self.assertTrue(all(prop["transformed"] for prop in results))


class TestBatchDatabaseFunctions(unittest.TestCase):
    """Tests for the batch database functions."""
    
    def setUp(self):
        """Set up the test case."""
        # Configure logging
        logging.basicConfig(level=logging.DEBUG)
        
        # Create mock database
        self.db = MockDatabase()
    
    def test_batch_db_insert(self):
        """Test batch database insertion."""
        
        async def run_async_test():
            # Create mock records
            records = [{"id": f"rec-{i}", "value": i} for i in range(100)]
            
            # Insert batch
            inserted_ids = await batch_db_insert(self.db, "test_table", records)
            
            # Check results (some may be filtered due to errors)
            self.assertGreaterEqual(len(inserted_ids), 95)
            self.assertEqual(len(self.db.data["test_table"]), len(inserted_ids))
        
        # Run the async test
        asyncio.run(run_async_test())


class TestEndToEndBatchProcessing(unittest.TestCase):
    """End-to-end tests for batch processing with actual modules."""
    
    @patch('unified_importer.core.enhanced_database.EnhancedDatabaseOperations')
    def test_integration_with_database(self, mock_db_class):
        """Test integration with a mock database."""
        # Configure mock database
        mock_db = MagicMock()
        mock_db.batch_insert.return_value = ["id-1", "id-2", "id-3"]
        mock_db_class.return_value = mock_db
        
        # Create a batch processing manager
        manager = BatchProcessingManager(
            db=mock_db,
            max_workers=2,
            adaptive_batch_sizing=True,
            parallel_processing=True
        )
        
        # Create test items
        items = [{"id": f"item-{i}", "value": i} for i in range(10)]
        
        # Define a processing function that uses the database
        def process_with_db(batch):
            # In a real scenario, this might do more processing
            return batch
        
        # Process items
        results = manager.process_items(
            items=items,
            operation="db_test",
            data_source="mock",
            process_func=process_with_db,
            initial_batch_size=5
        )
        
        # Verify results
        self.assertEqual(len(results), 10)
        
        # Verify manager stats
        stats = manager.get_stats()
        self.assertEqual(stats["total_items_processed"], 10)
        self.assertGreaterEqual(stats["total_batches"], 2)  # At least 2 batches
        
        # Clean up
        manager.close()


if __name__ == "__main__":
    unittest.main()