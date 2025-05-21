"""
Integration of enhanced batch processing with the unified molecular importer.

This module provides the connection between the enhanced batch processing system
and the unified molecular importer, allowing the importer to take advantage of
adaptive batch sizing, parallel processing, and other optimizations.
"""

import os
import time
import logging
import asyncio
from typing import Dict, List, Any, Optional, Tuple, Union, Set, Callable

# Import batch processing modules
from .enhanced_batch_processing import (
    AdaptiveBatchProcessor,
    ChemicalDataBatchStrategy,
    BatchStats
)

# Import database operations
from .enhanced_database import EnhancedDatabaseOperations


class BatchProcessingManager:
    """
    Manager for batch processing in the unified molecular importer.
    
    This class provides:
    - Integration with data sources (ChEMBL, PubChem)
    - Optimal batch sizing for different operations
    - Management of parallel processing resources
    - Monitoring and reporting of batch processing metrics
    """
    
    def __init__(
        self,
        db: EnhancedDatabaseOperations,
        max_workers: Optional[int] = None,
        adaptive_batch_sizing: bool = True,
        parallel_processing: bool = True,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the batch processing manager.
        
        Args:
            db: Database operations instance
            max_workers: Maximum worker threads/processes
            adaptive_batch_sizing: Whether to use adaptive batch sizing
            parallel_processing: Whether to use parallel processing
            logger: Logger instance
        """
        self.db = db
        self.logger = logger or logging.getLogger(__name__)
        self.adaptive_batch_sizing = adaptive_batch_sizing
        self.parallel_processing = parallel_processing
        
        # Determine optimal number of workers if not specified
        if max_workers is None:
            # Use CPU count, but cap to avoid excessive threads
            max_workers = min(os.cpu_count() or 4, 8)
        
        self.max_workers = max_workers
        
        # Track processors for different operations
        self.processors = {}
        
        # Track processing statistics
        self.stats = {
            "total_items_processed": 0,
            "total_batches": 0,
            "total_duration": 0,
            "operations": {}
        }
        
        self.logger.info(
            f"Initialized BatchProcessingManager with {max_workers} workers, "
            f"adaptive_batch_sizing={adaptive_batch_sizing}, "
            f"parallel_processing={parallel_processing}"
        )
    
    def get_processor(
        self,
        operation: str,
        data_source: str,
        process_func: Callable,
        initial_batch_size: int = 100
    ) -> AdaptiveBatchProcessor:
        """
        Get or create a batch processor for a specific operation.
        
        Args:
            operation: Name of the operation (e.g., 'molecule_transform', 'insert')
            data_source: Source of the data (e.g., 'chembl', 'pubchem')
            process_func: Function to process a batch of items
            initial_batch_size: Initial batch size
            
        Returns:
            AdaptiveBatchProcessor instance
        """
        processor_key = f"{operation}_{data_source}"
        
        # Return existing processor if available
        if processor_key in self.processors:
            return self.processors[processor_key]
        
        # Create a new processor
        if self.adaptive_batch_sizing:
            # Use chemical data specific strategy
            strategy = ChemicalDataBatchStrategy(
                data_source=data_source,
                initial_batch_size=initial_batch_size,
                min_batch_size=max(10, initial_batch_size // 10),
                max_batch_size=min(1000, initial_batch_size * 10),
                logger=self.logger
            )
        else:
            # Use fixed batch size
            from .enhanced_batch_processing import BatchProcessingStrategy
            strategy = BatchProcessingStrategy(
                initial_batch_size=initial_batch_size,
                min_batch_size=initial_batch_size,
                max_batch_size=initial_batch_size,
                logger=self.logger
            )
        
        # Create processor
        processor = AdaptiveBatchProcessor(
            process_func=process_func,
            max_workers=self.max_workers if self.parallel_processing else 1,
            strategy=strategy,
            logger=self.logger
        )
        
        # Store for reuse
        self.processors[processor_key] = processor
        
        # Initialize operation stats
        if operation not in self.stats["operations"]:
            self.stats["operations"][operation] = {
                "items_processed": 0,
                "batches": 0,
                "duration": 0,
                "data_sources": {}
            }
        
        if data_source not in self.stats["operations"][operation]["data_sources"]:
            self.stats["operations"][operation]["data_sources"][data_source] = {
                "items_processed": 0,
                "batches": 0,
                "duration": 0
            }
        
        return processor
    
    def process_items(
        self,
        items: List[Any],
        operation: str,
        data_source: str,
        process_func: Callable,
        initial_batch_size: int = 100
    ) -> List[Any]:
        """
        Process a list of items using batch processing.
        
        Args:
            items: List of items to process
            operation: Name of the operation
            data_source: Source of the data
            process_func: Function to process a batch of items
            initial_batch_size: Initial batch size
            
        Returns:
            List of processed results
        """
        if not items:
            return []
        
        # Get processor
        processor = self.get_processor(
            operation=operation,
            data_source=data_source,
            process_func=process_func,
            initial_batch_size=initial_batch_size
        )
        
        # Process items
        start_time = time.time()
        
        if self.parallel_processing:
            results = processor.process_items_parallel(items)
        else:
            results = processor.process_items(items)
        
        # Update statistics
        duration = time.time() - start_time
        self.stats["total_items_processed"] += len(items)
        self.stats["total_batches"] += 1
        self.stats["total_duration"] += duration
        
        # Update operation stats
        self.stats["operations"][operation]["items_processed"] += len(items)
        self.stats["operations"][operation]["batches"] += 1
        self.stats["operations"][operation]["duration"] += duration
        
        # Update data source stats
        self.stats["operations"][operation]["data_sources"][data_source]["items_processed"] += len(items)
        self.stats["operations"][operation]["data_sources"][data_source]["batches"] += 1
        self.stats["operations"][operation]["data_sources"][data_source]["duration"] += duration
        
        return results
    
    async def process_items_async(
        self,
        items: List[Any],
        operation: str,
        data_source: str,
        process_func: Callable,
        initial_batch_size: int = 100
    ) -> List[Any]:
        """
        Process a list of items asynchronously using batch processing.
        
        Args:
            items: List of items to process
            operation: Name of the operation
            data_source: Source of the data
            process_func: Function to process a batch of items
            initial_batch_size: Initial batch size
            
        Returns:
            List of processed results
        """
        if not items:
            return []
        
        # Get processor
        processor = self.get_processor(
            operation=operation,
            data_source=data_source,
            process_func=process_func,
            initial_batch_size=initial_batch_size
        )
        
        # Process items
        start_time = time.time()
        
        if self.parallel_processing:
            results = await processor.process_items_parallel_async(items)
        else:
            results = await processor.process_items_async(items)
        
        # Update statistics
        duration = time.time() - start_time
        self.stats["total_items_processed"] += len(items)
        self.stats["total_batches"] += 1
        self.stats["total_duration"] += duration
        
        # Update operation stats
        self.stats["operations"][operation]["items_processed"] += len(items)
        self.stats["operations"][operation]["batches"] += 1
        self.stats["operations"][operation]["duration"] += duration
        
        # Update data source stats
        self.stats["operations"][operation]["data_sources"][data_source]["items_processed"] += len(items)
        self.stats["operations"][operation]["data_sources"][data_source]["batches"] += 1
        self.stats["operations"][operation]["data_sources"][data_source]["duration"] += duration
        
        return results
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get detailed batch processing statistics.
        
        Returns:
            Dictionary with processing statistics
        """
        # Clone stats to avoid modification during read
        stats = {
            "total_items_processed": self.stats["total_items_processed"],
            "total_batches": self.stats["total_batches"],
            "total_duration": self.stats["total_duration"],
            "operations": {}
        }
        
        # Add derived metrics
        if stats["total_duration"] > 0:
            stats["items_per_second"] = stats["total_items_processed"] / stats["total_duration"]
        else:
            stats["items_per_second"] = 0
        
        # Add processor stats
        for processor_key, processor in self.processors.items():
            operation, data_source = processor_key.split("_", 1)
            
            if operation not in stats["operations"]:
                stats["operations"][operation] = {
                    "items_processed": self.stats["operations"][operation]["items_processed"],
                    "batches": self.stats["operations"][operation]["batches"],
                    "duration": self.stats["operations"][operation]["duration"],
                    "data_sources": {}
                }
            
            # Add derived metrics for operation
            op_stats = stats["operations"][operation]
            if op_stats["duration"] > 0:
                op_stats["items_per_second"] = op_stats["items_processed"] / op_stats["duration"]
            else:
                op_stats["items_per_second"] = 0
            
            # Add data source stats
            ds_stats = self.stats["operations"][operation]["data_sources"].get(data_source, {})
            stats["operations"][operation]["data_sources"][data_source] = {
                "items_processed": ds_stats.get("items_processed", 0),
                "batches": ds_stats.get("batches", 0),
                "duration": ds_stats.get("duration", 0)
            }
            
            # Add derived metrics for data source
            ds_stats = stats["operations"][operation]["data_sources"][data_source]
            if ds_stats["duration"] > 0:
                ds_stats["items_per_second"] = ds_stats["items_processed"] / ds_stats["duration"]
            else:
                ds_stats["items_per_second"] = 0
            
            # Add processor-specific stats
            processor_stats = processor.get_stats()
            stats["operations"][operation]["data_sources"][data_source]["processor"] = processor_stats
        
        # Add configuration info
        stats["config"] = {
            "max_workers": self.max_workers,
            "adaptive_batch_sizing": self.adaptive_batch_sizing,
            "parallel_processing": self.parallel_processing
        }
        
        return stats
    
    def close(self) -> None:
        """Clean up resources."""
        for processor in self.processors.values():
            processor.close()


# Specialized batch processing for chemical data import

def batch_molecule_transform(transformer, molecules: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Process a batch of molecules through a transformer.
    
    Args:
        transformer: Molecule transformer instance
        molecules: List of molecule dictionaries
        
    Returns:
        List of transformed molecules
    """
    transformed = []
    
    for molecule in molecules:
        try:
            # Apply transformations
            result = transformer.transform(molecule)
            if result:
                transformed.append(result)
        except Exception as e:
            # Log error and continue with the next molecule
            transformer.logger.error(f"Error transforming molecule {molecule.get('id', 'unknown')}: {str(e)}")
    
    return transformed


def batch_property_transform(transformer, properties: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Process a batch of properties through a transformer.
    
    Args:
        transformer: Property transformer instance
        properties: List of property dictionaries
        
    Returns:
        List of transformed properties
    """
    transformed = []
    
    for prop in properties:
        try:
            # Apply transformations
            result = transformer.transform(prop)
            if result:
                transformed.append(result)
        except Exception as e:
            # Log error and continue with the next property
            transformer.logger.error(f"Error transforming property {prop.get('id', 'unknown')}: {str(e)}")
    
    return transformed


async def batch_db_insert(db, table: str, records: List[Dict[str, Any]]) -> List[str]:
    """
    Insert a batch of records into the database.
    
    Args:
        db: Database operations instance
        table: Table name
        records: List of records to insert
        
    Returns:
        List of inserted record IDs
    """
    try:
        # Begin transaction
        async with db.transaction_async() as tx_id:
            # Insert all records in one batch
            inserted_ids = await db.batch_insert(
                table=table,
                records=records,
                transaction_id=tx_id
            )
            return inserted_ids
    except Exception as e:
        db.logger.error(f"Error inserting batch into {table}: {str(e)}")
        return []