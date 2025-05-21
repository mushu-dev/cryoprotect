"""
Tests for the enhanced batch processing module.

This module contains tests for the adaptive batch processing
features including batch sizing and parallel processing.
"""

import unittest
import time
import random
import asyncio
import logging
from typing import List, Dict, Any

# Import the module to test
from unified_importer.core.enhanced_batch_processing import (
    BatchStats,
    BatchProcessingStrategy,
    AdaptiveBatchProcessor,
    ChemicalDataBatchStrategy
)


# Mock for CPU-intensive processing
def process_items(items: List[int]) -> List[Dict[str, Any]]:
    """Process a list of integers with simulated delays."""
    results = []
    
    for item in items:
        # Simulate processing
        time.sleep(0.01)  # 10ms per item
        
        # Randomly fail some items
        if random.random() < 0.05:  # 5% failure rate
            continue
        
        # Create result
        result = {
            "id": item,
            "value": item * 2,
            "processed": True
        }
        results.append(result)
    
    return results


# Mock for memory-intensive processing
def process_memory_intensive(items: List[int]) -> List[Dict[str, Any]]:
    """Process a list of integers with memory usage proportional to the input."""
    results = []
    # Hold some data in memory
    memory_sink = []
    
    for item in items:
        # Simulate processing
        time.sleep(0.01)  # 10ms per item
        
        # Simulate memory usage (proportional to the item value)
        # This is just for testing - don't actually waste this much memory
        memory_usage = item % 10
        memory_sink.append([0] * (memory_usage * 100))
        
        # Create result
        result = {
            "id": item,
            "value": item * 2,
            "processed": True
        }
        results.append(result)
    
    return results


class TestBatchStats(unittest.TestCase):
    """Tests for the BatchStats class."""
    
    def test_batch_stats(self):
        """Test batch statistics calculation."""
        # Create batch stats
        stats = BatchStats(batch_size=100)
        
        # Add processing times
        stats.processing_times = [0.1, 0.2, 0.3]
        stats.items_processed = 100
        stats.success_count = 95
        stats.error_count = 5
        
        # Set end time 1 second after start
        stats.start_time = time.time() - 1
        stats.complete()
        
        # Check derived metrics
        self.assertAlmostEqual(stats.duration, 1.0, delta=0.1)
        self.assertIsNotNone(stats.memory_usage_delta)
        self.assertAlmostEqual(stats.items_per_second, 100.0, delta=10.0)
        self.assertAlmostEqual(stats.avg_processing_time, 0.2, delta=0.01)


class TestBatchProcessingStrategy(unittest.TestCase):
    """Tests for the BatchProcessingStrategy class."""
    
    def setUp(self):
        """Set up the test case."""
        # Configure logging
        logging.basicConfig(level=logging.DEBUG)
        
        # Create a strategy
        self.strategy = BatchProcessingStrategy(
            initial_batch_size=100,
            min_batch_size=10,
            max_batch_size=1000,
            target_duration=1.0
        )
    
    def test_initial_state(self):
        """Test the initial state of the strategy."""
        self.assertEqual(self.strategy.current_batch_size, 100)
        self.assertEqual(self.strategy.total_items_processed, 0)
        self.assertEqual(len(self.strategy.batch_stats), 0)
    
    def test_get_next_batch_size_without_stats(self):
        """Test getting the next batch size without stats."""
        # Without stats, should return the current size
        next_size = self.strategy.get_next_batch_size()
        self.assertEqual(next_size, 100)
    
    def test_record_batch_stats(self):
        """Test recording batch statistics."""
        # Create batch stats
        stats = BatchStats(batch_size=100)
        stats.items_processed = 100
        stats.success_count = 95
        stats.error_count = 5
        stats.start_time = time.time() - 1
        stats.complete()
        
        # Record stats
        self.strategy.record_batch_stats(stats)
        
        # Check totals
        self.assertEqual(self.strategy.total_items_processed, 100)
        self.assertEqual(self.strategy.total_success_count, 95)
        self.assertEqual(self.strategy.total_error_count, 5)
        self.assertAlmostEqual(self.strategy.total_duration, 1.0, delta=0.1)
        self.assertEqual(len(self.strategy.batch_stats), 1)
    
    def test_get_overall_stats(self):
        """Test getting overall statistics."""
        # Record some stats
        for i in range(3):
            stats = BatchStats(batch_size=100)
            stats.items_processed = 100
            stats.success_count = 95
            stats.error_count = 5
            stats.start_time = time.time() - 1
            stats.complete()
            self.strategy.record_batch_stats(stats)
        
        # Get overall stats
        overall = self.strategy.get_overall_stats()
        
        # Check stats
        self.assertEqual(overall["total_items_processed"], 300)
        self.assertEqual(overall["total_success_count"], 285)
        self.assertEqual(overall["total_error_count"], 15)
        self.assertAlmostEqual(overall["avg_batch_size"], 100.0)
        self.assertEqual(overall["batch_count"], 3)
    
    def test_adaptive_sizing(self):
        """Test adaptive batch sizing."""
        # Record stats for a batch that's too slow
        for i in range(3):
            stats = BatchStats(batch_size=100)
            stats.items_processed = 100
            stats.success_count = 95
            stats.error_count = 5
            stats.start_time = time.time() - 2  # Twice the target duration
            stats.complete()
            self.strategy.record_batch_stats(stats)
        
        # Get next batch size - should be smaller
        next_size = self.strategy.get_next_batch_size()
        self.assertLess(next_size, 100)
        
        # Record stats for a batch that's too fast
        for i in range(3):
            stats = BatchStats(batch_size=next_size)
            stats.items_processed = next_size
            stats.success_count = next_size - 2
            stats.error_count = 2
            stats.start_time = time.time() - 0.5  # Half the target duration
            stats.complete()
            self.strategy.record_batch_stats(stats)
        
        # Get next batch size - should be larger
        new_size = self.strategy.get_next_batch_size()
        self.assertGreater(new_size, next_size)


class TestAdaptiveBatchProcessor(unittest.TestCase):
    """Tests for the AdaptiveBatchProcessor class."""
    
    def setUp(self):
        """Set up the test case."""
        # Configure logging
        logging.basicConfig(level=logging.DEBUG)
        
        # Create a processor
        self.processor = AdaptiveBatchProcessor(
            process_func=process_items,
            initial_batch_size=50,
            min_batch_size=10,
            max_batch_size=200,
            max_workers=4
        )
    
    def test_process_batch(self):
        """Test processing a single batch."""
        # Create a batch of 50 items
        items = list(range(50))
        
        # Process the batch
        results, stats = self.processor.process_batch(items)
        
        # Check results
        self.assertEqual(len(results), 47)  # Approximately 5% failure rate
        self.assertEqual(stats.batch_size, 50)
        self.assertEqual(stats.items_processed, 50)
        self.assertEqual(stats.success_count, 47)
        self.assertEqual(stats.error_count, 3)
    
    def test_process_items(self):
        """Test processing a list of items."""
        # Create a list of 120 items
        items = list(range(120))
        
        # Process the items
        results = self.processor.process_items(items)
        
        # Check results (should be processed in 3 batches)
        self.assertGreaterEqual(len(results), 110)  # Approximately 5% failure rate
        
        # Check that stats were recorded
        stats = self.processor.get_stats()
        self.assertEqual(stats["total_items_processed"], 120)
        self.assertGreaterEqual(stats["batch_count"], 3)  # 3 batches
    
    def test_process_items_parallel(self):
        """Test processing items in parallel."""
        # Create a list of 200 items
        items = list(range(200))
        
        # Process the items in parallel
        results = self.processor.process_items_parallel(items)
        
        # Check results
        self.assertGreaterEqual(len(results), 180)  # Approximately 5% failure rate
        
        # Check that stats were recorded
        stats = self.processor.get_stats()
        self.assertEqual(stats["total_items_processed"], 200)
        self.assertGreater(stats["batch_count"], 1)  # Multiple batches
    
    def test_batch_size_adaptation(self):
        """Test that batch size adapts over time."""
        # Initial batch size
        initial_size = self.processor.strategy.current_batch_size
        
        # Process several batches to trigger adaptation
        for _ in range(5):
            items = list(range(100))
            _ = self.processor.process_items(items)
        
        # Check that batch size has changed
        final_size = self.processor.strategy.current_batch_size
        self.assertNotEqual(final_size, initial_size)
    
    def test_async_processing(self):
        """Test asynchronous processing."""
        async def run_async_test():
            # Create a list of 150 items
            items = list(range(150))
            
            # Process the items asynchronously
            results = await self.processor.process_items_async(items)
            
            # Check results
            self.assertGreaterEqual(len(results), 140)  # Approximately 5% failure rate
            
            # Process items in parallel asynchronously
            items = list(range(200))
            results = await self.processor.process_items_parallel_async(items)
            
            # Check results
            self.assertGreaterEqual(len(results), 180)  # Approximately 5% failure rate
        
        # Run the async test
        asyncio.run(run_async_test())
    
    def test_memory_intensive_processing(self):
        """Test memory-intensive processing."""
        # Create a processor for memory-intensive operations
        processor = AdaptiveBatchProcessor(
            process_func=process_memory_intensive,
            initial_batch_size=50,
            min_batch_size=10,
            max_batch_size=200
        )
        
        # Process several batches to trigger memory adaptation
        for _ in range(5):
            items = list(range(100))
            _ = processor.process_items(items)
        
        # Check that batch size adjusts based on memory usage
        stats = processor.get_stats()
        self.assertGreater(stats["batch_count"], 1)


class TestChemicalDataBatchStrategy(unittest.TestCase):
    """Tests for the ChemicalDataBatchStrategy class."""
    
    def test_chembl_strategy(self):
        """Test batch strategy for ChEMBL data."""
        # Create a strategy for ChEMBL data
        strategy = ChemicalDataBatchStrategy(
            data_source="chembl",
            initial_batch_size=100,
            min_batch_size=10,
            max_batch_size=500
        )
        
        # Check that initial batch size was adjusted
        self.assertEqual(strategy.initial_batch_size, 50)
        self.assertEqual(strategy.current_batch_size, 50)
        self.assertGreaterEqual(strategy.target_duration, 10.0)
    
    def test_pubchem_strategy(self):
        """Test batch strategy for PubChem data."""
        # Create a strategy for PubChem data
        strategy = ChemicalDataBatchStrategy(
            data_source="pubchem",
            initial_batch_size=100,
            min_batch_size=10,
            max_batch_size=500
        )
        
        # Check that initial values are appropriate
        self.assertEqual(strategy.initial_batch_size, 100)
        self.assertEqual(strategy.current_batch_size, 100)
        self.assertGreaterEqual(strategy.target_duration, 8.0)
    
    def test_chemical_data_adaptation(self):
        """Test adaptation specific to chemical data."""
        # Create a strategy
        strategy = ChemicalDataBatchStrategy(
            data_source="chembl",
            initial_batch_size=50,
            min_batch_size=10,
            max_batch_size=500
        )
        
        # Simulate batches with unstable processing times
        for i in range(5):
            stats = BatchStats(batch_size=50)
            stats.items_processed = 50
            stats.success_count = 48
            stats.error_count = 2
            stats.start_time = time.time() - (1.0 + i)  # Increasingly slower
            stats.processing_times = [0.02 * (1 + i)] * 50
            stats.memory_usage_start = 1000000
            stats.memory_usage_end = 1000000 + 50 * 1000000  # 1MB per item
            stats.complete()
            strategy.record_batch_stats(stats)
        
        # Get next batch size - should be reduced due to instability and memory usage
        next_size = strategy.get_next_batch_size()
        self.assertLess(next_size, 50)


if __name__ == "__main__":
    unittest.main()