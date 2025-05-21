"""
Base class for molecular data sources.

This module defines the interface that all molecular data sources must implement.
"""

import abc
import asyncio
import time
import logging
from typing import Dict, List, Any, Optional, Union, Set, Tuple, AsyncIterator

from ..core.database import DatabaseOperations
from ..core.checkpoint import CheckpointManager
from ..core.progress import ProgressTracker
from ..core.validation import MoleculeValidator
from ..core.batch_integration import BatchProcessingManager


class MolecularDataSource(abc.ABC):
    """
    Abstract base class for molecular data sources.
    
    This class defines the interface and common functionality for all 
    molecular data sources (e.g., PubChem, ChEMBL).
    """
    
    def __init__(
        self,
        db_operations: DatabaseOperations,
        checkpoint_manager: Optional[CheckpointManager] = None,
        progress_tracker: Optional[ProgressTracker] = None,
        validator: Optional[MoleculeValidator] = None,
        batch_manager: Optional[BatchProcessingManager] = None,
        config: Optional[Dict[str, Any]] = None,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the molecular data source.

        Args:
            db_operations: Database operations instance
            checkpoint_manager: Checkpoint manager for resumable imports
            progress_tracker: Progress tracker for monitoring import progress
            validator: Validator for molecular data
            batch_manager: Batch processing manager for optimized batch operations
            config: Configuration dictionary
            logger: Logger instance
        """
        self.db = db_operations
        self.checkpoint_manager = checkpoint_manager
        self.progress_tracker = progress_tracker
        self.validator = validator or MoleculeValidator()
        self.batch_manager = batch_manager
        self.config = config or {}
        self.logger = logger or logging.getLogger(__name__)
        
        # Extract common configuration
        self.batch_size = self.config.get('batch_size', 50)
        self.max_retries = self.config.get('max_retries', 3)
        self.retry_delay = self.config.get('retry_delay', 2.0)
        self.api_delay = self.config.get('api_delay', 0.5)
        self.dry_run = self.config.get('dry_run', False)
    
    @abc.abstractmethod
    async def fetch_compound(self, identifier: str) -> Optional[Dict[str, Any]]:
        """
        Fetch a single compound by identifier.
        
        Args:
            identifier: Compound identifier (e.g., PubChem CID, ChEMBL ID)
            
        Returns:
            Compound data dictionary or None if not found
        """
        pass
    
    @abc.abstractmethod
    async def search_compounds(
        self,
        query: str,
        max_results: Optional[int] = None
    ) -> List[str]:
        """
        Search for compounds matching a query.
        
        Args:
            query: Search query string
            max_results: Maximum number of results to return
            
        Returns:
            List of compound identifiers matching the query
        """
        pass
    
    @abc.abstractmethod
    async def get_compound_batch(
        self,
        identifiers: List[str]
    ) -> List[Dict[str, Any]]:
        """
        Fetch multiple compounds by identifier.
        
        Args:
            identifiers: List of compound identifiers
            
        Returns:
            List of compound data dictionaries
        """
        pass
    
    @abc.abstractmethod
    async def get_property_data(
        self,
        compound_data: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """
        Fetch property data for a compound.
        
        Args:
            compound_data: Compound data dictionary
            
        Returns:
            List of property data dictionaries
        """
        pass
    
    @abc.abstractmethod
    async def get_compound_count(self, query: Optional[str] = None) -> int:
        """
        Get the total count of compounds matching a query.
        
        Args:
            query: Optional search query to filter compounds
            
        Returns:
            Total count of matched compounds
        """
        pass
    
    @abc.abstractmethod
    def get_source_name(self) -> str:
        """
        Get the name of this data source.
        
        Returns:
            Source name (e.g., "PubChem", "ChEMBL")
        """
        pass
    
    @abc.abstractmethod
    async def transform_compound_data(
        self,
        raw_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Transform raw compound data to the unified format.
        
        Args:
            raw_data: Raw compound data dictionary
            
        Returns:
            Transformed compound data in unified format
        """
        pass
    
    @abc.abstractmethod
    async def transform_property_data(
        self,
        raw_properties: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Transform raw property data to the unified format.
        
        Args:
            raw_properties: List of raw property dictionaries
            
        Returns:
            List of transformed property dictionaries in unified format
        """
        pass
    
    async def import_compound(
        self,
        identifier: str
    ) -> Tuple[bool, Optional[str], Optional[Dict[str, Any]]]:
        """
        Import a single compound by identifier.
        
        Args:
            identifier: Compound identifier
            
        Returns:
            Tuple of (success, error_message, imported_data)
        """
        try:
            # Check if already processed
            if self.checkpoint_manager and self.checkpoint_manager.is_processed(identifier):
                self.logger.debug(f"Compound {identifier} already processed, skipping")
                
                if self.progress_tracker:
                    self.progress_tracker.update(skipped=1)
                    
                return True, None, None
            
            # Fetch the compound data
            compound_data = await self.fetch_compound(identifier)
            
            if not compound_data:
                error_message = f"No data found for compound {identifier}"
                self.logger.warning(error_message)
                
                if self.checkpoint_manager:
                    self.checkpoint_manager.mark_skipped(identifier)
                    
                if self.progress_tracker:
                    self.progress_tracker.update(skipped=1)
                    
                return False, error_message, None
            
            # Transform to unified format
            transformed_data = await self.transform_compound_data(compound_data)
            
            # Validate the transformed data
            is_valid, error = self.validator.validate_molecule(transformed_data)
            
            if not is_valid:
                if self.checkpoint_manager:
                    self.checkpoint_manager.mark_failed(identifier)
                    
                if self.progress_tracker:
                    self.progress_tracker.update(failed=1)
                    
                return False, error, None
            
            # Fetch and transform properties
            property_data = await self.get_property_data(compound_data)
            transformed_properties = await self.transform_property_data(property_data)
            
            # Insert into database (or simulate in dry run mode)
            if self.dry_run:
                self.logger.info(
                    f"[DRY RUN] Would insert compound {identifier}",
                    extra={'structured_data': {
                        'molecule': transformed_data,
                        'properties': transformed_properties
                    }}
                )
                molecule_id = "dry-run-id"
            else:
                # Check if molecule already exists in database
                exists, existing_id = await self.db.molecule_exists(identifier)
                
                if exists:
                    self.logger.info(f"Compound {identifier} already exists in database with ID {existing_id}")
                    molecule_id = existing_id
                else:
                    # Insert the molecule
                    try:
                        async with self.db.transaction():
                            molecule_id = await self.db.insert_molecule(transformed_data)
                            
                            # Insert properties
                            if transformed_properties:
                                await self.db.insert_properties(molecule_id, transformed_properties)
                                
                            # Insert synonyms if present
                            if 'synonyms' in transformed_data and transformed_data['synonyms']:
                                await self.db.insert_synonyms(molecule_id, transformed_data['synonyms'])
                    except Exception as e:
                        error_message = f"Database error inserting compound {identifier}: {str(e)}"
                        self.logger.error(error_message)
                        
                        if self.checkpoint_manager:
                            self.checkpoint_manager.mark_failed(identifier)
                            
                        if self.progress_tracker:
                            self.progress_tracker.update(failed=1)
                            
                        return False, error_message, None
            
            # Mark as processed
            if self.checkpoint_manager:
                self.checkpoint_manager.mark_processed(identifier)
                
            if self.progress_tracker:
                self.progress_tracker.update(successful=1)
                
            self.logger.info(f"Successfully imported compound {identifier}")
            
            # Return the result
            return True, None, {
                'molecule_id': molecule_id,
                'identifier': identifier,
                'name': transformed_data.get('name', '')
            }
            
        except Exception as e:
            error_message = f"Error importing compound {identifier}: {str(e)}"
            self.logger.exception(error_message)
            
            if self.checkpoint_manager:
                self.checkpoint_manager.mark_failed(identifier)
                
            if self.progress_tracker:
                self.progress_tracker.update(failed=1)
                
            return False, error_message, None
    
    async def import_compounds(
        self,
        identifiers: List[str]
    ) -> Tuple[int, int, List[Tuple[str, str]]]:
        """
        Import multiple compounds by identifier.
        
        Args:
            identifiers: List of compound identifiers
            
        Returns:
            Tuple of (success_count, failure_count, failures_with_reasons)
        """
        success_count = 0
        failure_count = 0
        failures: List[Tuple[str, str]] = []
        
        # Process compounds in batches
        for i in range(0, len(identifiers), self.batch_size):
            batch = identifiers[i:i + self.batch_size]
            
            if self.checkpoint_manager:
                self.checkpoint_manager.update_batch(batch)
                
            # Log batch progress
            self.logger.info(f"Processing batch {i//self.batch_size + 1}/{(len(identifiers) + self.batch_size - 1)//self.batch_size}")
            
            # Create tasks for batch processing
            tasks = []
            for identifier in batch:
                # Skip already processed compounds
                if self.checkpoint_manager and self.checkpoint_manager.is_processed(identifier):
                    continue
                    
                tasks.append(self.import_compound(identifier))
            
            # Process batch
            batch_results = await asyncio.gather(*tasks, return_exceptions=True)
            
            # Process results
            for j, result in enumerate(batch_results):
                if isinstance(result, Exception):
                    # Handle exceptions
                    error_message = f"Exception during import: {str(result)}"
                    self.logger.error(error_message)
                    
                    if j < len(batch):
                        failures.append((batch[j], error_message))
                        failure_count += 1
                else:
                    # Process normal results
                    success, error, _ = result
                    
                    if success:
                        success_count += 1
                    else:
                        failure_count += 1
                        if j < len(batch) and error:
                            failures.append((batch[j], error))
            
            # Save checkpoint after each batch
            if self.checkpoint_manager:
                self.checkpoint_manager.save()
                
            # Update progress
            if self.progress_tracker:
                self.progress_tracker.log_progress()
                
            # Apply API delay
            if i + self.batch_size < len(identifiers) and self.api_delay > 0:
                await asyncio.sleep(self.api_delay)
        
        return success_count, failure_count, failures
    
    async def search_and_import(
        self,
        query: str,
        max_results: Optional[int] = None
    ) -> Tuple[int, int, List[Tuple[str, str]]]:
        """
        Search for compounds and import the results.
        
        Args:
            query: Search query string
            max_results: Maximum number of results to import
            
        Returns:
            Tuple of (success_count, failure_count, failures_with_reasons)
        """
        # Search for compounds
        self.logger.info(f"Searching for compounds matching query: {query}")
        identifiers = await self.search_compounds(query, max_results)
        
        if not identifiers:
            self.logger.warning(f"No compounds found matching query: {query}")
            return 0, 0, []
        
        self.logger.info(f"Found {len(identifiers)} compounds matching query")
        
        # Set total items in progress tracker
        if self.progress_tracker:
            self.progress_tracker.total_items = len(identifiers)
            
        # Import the found compounds
        return await self.import_compounds(identifiers)
    
    async def import_all(
        self,
        max_results: Optional[int] = None
    ) -> Tuple[int, int, List[Tuple[str, str]]]:
        """
        Import all compounds from this data source.
        
        Warning: This can be a very large operation!
        
        Args:
            max_results: Maximum number of results to import
            
        Returns:
            Tuple of (success_count, failure_count, failures_with_reasons)
        """
        # This is a placeholder implementation
        # Subclasses should override with more efficient implementations
        self.logger.warning(
            "Using default import_all implementation which may be inefficient. "
            "Consider overriding this method in the data source subclass."
        )
        
        return await self.search_and_import("*", max_results)
    
    async def stream_compound_identifiers(
        self,
        query: Optional[str] = None,
        batch_size: int = 1000
    ) -> AsyncIterator[List[str]]:
        """
        Stream compound identifiers in batches.
        
        This is a base implementation that subclasses can override
        with more efficient source-specific implementations.
        
        Args:
            query: Optional query string to filter compounds
            batch_size: Number of identifiers to yield in each batch
            
        Yields:
            Batches of compound identifiers
        """
        # Subclasses should override with more efficient implementations
        self.logger.warning(
            "Using default stream_compound_identifiers implementation which may be inefficient. "
            "Consider overriding this method in the data source subclass."
        )
        
        # Get all identifiers
        all_ids = await self.search_compounds(query or "*")
        
        # Yield in batches
        for i in range(0, len(all_ids), batch_size):
            yield all_ids[i:i + batch_size]