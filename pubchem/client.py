"""
Resilient PubChem API client for CryoProtect v2.

This module provides a robust client for interacting with the PubChem API,
featuring adaptive rate limiting, caching, retry logic, and circuit breaking.
"""

import os
import time
import logging
import json
import requests
from typing import Dict, Any, Optional, List, Union, Callable
from datetime import datetime
from pathlib import Path

from .cache import PubChemCache
from .rate_limiter import AdaptiveRateLimiter
from .scheduler import WeekendJobScheduler
from .utils import (
    CircuitBreaker, 
    retry_with_backoff, 
    create_request_key,
    CircuitBreakerError
)

logger = logging.getLogger(__name__)

class ResilientPubChemClient:
    """
    Resilient client for PubChem API.
    
    Features:
    - Adaptive rate limiting (slower on weekdays, faster on weekends)
    - Multi-level caching (in-memory and disk-based)
    - Exponential backoff retry logic
    - Circuit breaker to prevent repeated failures
    - Fallback to pre-cached data when API is unavailable
    - Weekend job scheduler for bulk operations
    """
    
    def __init__(
        self,
        cache_dir: str = "cache/pubchem",
        weekday_requests_per_second: float = 2.0,
        weekend_requests_per_second: float = 5.0,
        max_retries: int = 5,
        failure_threshold: int = 3,
        recovery_timeout: int = 60,
        cache_ttl: int = 86400 * 30,  # 30 days
        memory_cache_size: int = 1000,
        enable_scheduler: bool = True
    ):
        """
        Initialize the client.
        
        Args:
            cache_dir: Directory to store cache files
            weekday_requests_per_second: Maximum requests per second on weekdays
            weekend_requests_per_second: Maximum requests per second on weekends
            max_retries: Maximum number of retries for failed requests
            failure_threshold: Number of failures before opening the circuit
            recovery_timeout: Time in seconds to wait before trying again
            cache_ttl: Time-to-live for cache entries in seconds
            memory_cache_size: Maximum number of items to keep in memory cache
            enable_scheduler: Whether to enable the weekend job scheduler
        """
        # Create cache directory if it doesn't exist
        os.makedirs(cache_dir, exist_ok=True)
        
        # Initialize components
        self.cache = PubChemCache(
            cache_dir=cache_dir,
            memory_size=memory_cache_size,
            ttl=cache_ttl
        )
        
        self.rate_limiter = AdaptiveRateLimiter(
            weekday_requests_per_second=weekday_requests_per_second,
            weekend_requests_per_second=weekend_requests_per_second
        )
        
        self.scheduler = WeekendJobScheduler(
            storage_dir=os.path.join(cache_dir, "scheduler")
        )
        
        # Circuit breaker for API requests
        self.circuit_breaker = CircuitBreaker(
            failure_threshold=failure_threshold,
            recovery_timeout=recovery_timeout,
            expected_exceptions=(requests.RequestException,),
            name="pubchem_api"
        )
        
        # Configuration
        self.max_retries = max_retries
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        
        # Start scheduler if enabled
        if enable_scheduler:
            self.scheduler.start()
        
        logger.info("ResilientPubChemClient initialized")
    
    def __del__(self):
        """Clean up resources when the client is destroyed."""
        try:
            self.scheduler.stop()
        except:
            pass
    
    @retry_with_backoff(expected_exceptions=(requests.RequestException, CircuitBreakerError))
    def _make_request(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Make a request to the PubChem API with circuit breaker protection.
        
        Args:
            endpoint: API endpoint
            params: Request parameters
            
        Returns:
            JSON response from the API
            
        Raises:
            requests.RequestException: If the request fails
            CircuitBreakerError: If the circuit is open
        """
        url = f"{self.base_url}/{endpoint}"
        
        # Check circuit state before making request
        circuit_stats = self.get_circuit_breaker_stats()
        if circuit_stats["state"] == "open":
            recovery_time = self.recovery_timeout - (time.time() - circuit_stats["last_failure_time"])
            if recovery_time > 0:
                logger.warning(
                    f"Circuit {circuit_stats['name']} is OPEN. Request rejected. "
                    f"Retry after {recovery_time:.1f} seconds."
                )
                raise CircuitBreakerError(
                    f"Circuit {circuit_stats['name']} is OPEN. Request rejected. "
                    f"Retry after {recovery_time:.1f} seconds."
                )
        
        # Wait for rate limiter
        self.rate_limiter.wait()
        
        # Make request with circuit breaker
        @self.circuit_breaker
        def do_request():
            response = requests.get(url, params=params)
            response.raise_for_status()
            return response.json()
        
        return do_request()
    
    def get_molecule_properties(
        self, 
        cid: Union[str, int],
        use_cache: bool = True,
        fallback_to_cache: bool = True,
        progress_tracker: Optional[Any] = None
    ) -> Dict[str, Any]:
        """
        Fetch molecular properties and names from PubChem.
        
        Args:
            cid: PubChem Compound ID
            use_cache: Whether to use cached data if available
            fallback_to_cache: Whether to fallback to cached data if API fails
            progress_tracker: Optional progress tracker object
            
        Returns:
            Dictionary with molecular properties
        """
        # Convert CID to string
        cid_str = str(cid)
        
        # Create cache key
        endpoint = f"compound/cid/{cid_str}/property/MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,IsomericSMILES,InChI,InChIKey,IUPACName,Title/JSON"
        cache_key = create_request_key(endpoint)
        
        # Check cache if enabled
        if use_cache:
            cached_data = self.cache.get(cache_key)
            if cached_data:
                logger.debug(f"Cache hit for CID {cid_str}")
                return cached_data
        
        try:
            # Make request to PubChem API
            response = self._make_request(endpoint)
            
            # Extract properties
            properties = response["PropertyTable"]["Properties"][0]
            
            # Format result
            result = {
                "CID": cid_str,
                "Molecular Formula": properties.get("MolecularFormula"),
                "Molecular Weight": properties.get("MolecularWeight"),
                "LogP": properties.get("XLogP"),
                "TPSA": properties.get("TPSA"),
                "H-Bond Donors": properties.get("HBondDonorCount"),
                "H-Bond Acceptors": properties.get("HBondAcceptorCount"),
                "SMILES": properties.get("IsomericSMILES"),
                "InChI": properties.get("InChI"),
                "InChIKey": properties.get("InChIKey"),
                "IUPACName": properties.get("IUPACName"),
                "Title": properties.get("Title"),
                "PubChem Link": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid_str}"
            }
            
            # Cache result
            self.cache.set(cache_key, result)
            
            return result
            
        except Exception as e:
            error_msg = f"Error fetching properties for CID {cid_str}: {str(e)}"
            logger.warning(error_msg)
            
            if progress_tracker:
                progress_tracker.add_error(error_msg)
            
            # Fallback to cache if enabled
            if fallback_to_cache:
                cached_data = self.cache.get(cache_key)
                if cached_data:
                    logger.info(f"Falling back to cached data for CID {cid_str}")
                    return cached_data
            
            # Return error if no cached data available
            return {
                "CID": cid_str,
                "Error": str(e)
            }
    
    def get_compound_synonyms(
        self, 
        cid: Union[str, int],
        use_cache: bool = True,
        fallback_to_cache: bool = True
    ) -> Dict[str, Any]:
        """
        Fetch compound synonyms from PubChem.
        
        Args:
            cid: PubChem Compound ID
            use_cache: Whether to use cached data if available
            fallback_to_cache: Whether to fallback to cached data if API fails
            
        Returns:
            Dictionary with compound synonyms
        """
        # Convert CID to string
        cid_str = str(cid)
        
        # Create cache key
        endpoint = f"compound/cid/{cid_str}/synonyms/JSON"
        cache_key = create_request_key(endpoint)
        
        # Check cache if enabled
        if use_cache:
            cached_data = self.cache.get(cache_key)
            if cached_data:
                logger.debug(f"Cache hit for synonyms of CID {cid_str}")
                return cached_data
        
        try:
            # Make request to PubChem API
            response = self._make_request(endpoint)
            
            # Extract synonyms
            information = response.get("InformationList", {}).get("Information", [])
            
            if information and "Synonym" in information[0]:
                synonyms = information[0]["Synonym"]
            else:
                synonyms = []
            
            # Format result
            result = {
                "CID": cid_str,
                "Synonyms": synonyms
            }
            
            # Cache result
            self.cache.set(cache_key, result)
            
            return result
            
        except Exception as e:
            error_msg = f"Error fetching synonyms for CID {cid_str}: {str(e)}"
            logger.warning(error_msg)
            
            # Fallback to cache if enabled
            if fallback_to_cache:
                cached_data = self.cache.get(cache_key)
                if cached_data:
                    logger.info(f"Falling back to cached data for synonyms of CID {cid_str}")
                    return cached_data
            
            # Return error if no cached data available
            return {
                "CID": cid_str,
                "Error": str(e)
            }
    
    def search_compounds(
        self, 
        query: str,
        use_cache: bool = True,
        fallback_to_cache: bool = True
    ) -> Dict[str, Any]:
        """
        Search for compounds by name or other identifiers.
        
        Args:
            query: Search query
            use_cache: Whether to use cached data if available
            fallback_to_cache: Whether to fallback to cached data if API fails
            
        Returns:
            Dictionary with search results
        """
        # Create cache key
        endpoint = "compound/name/JSON"
        cache_key = create_request_key(endpoint, {"name": query})
        
        # Check cache if enabled
        if use_cache:
            cached_data = self.cache.get(cache_key)
            if cached_data:
                logger.debug(f"Cache hit for search query '{query}'")
                return cached_data
        
        try:
            # Make request to PubChem API
            response = self._make_request(endpoint, {"name": query})
            
            # Extract CIDs
            if "IdentifierList" in response and "CID" in response["IdentifierList"]:
                cids = response["IdentifierList"]["CID"]
            else:
                cids = []
            
            # Format result
            result = {
                "Query": query,
                "CIDs": cids
            }
            
            # Cache result
            self.cache.set(cache_key, result)
            
            return result
            
        except Exception as e:
            error_msg = f"Error searching for '{query}': {str(e)}"
            logger.warning(error_msg)
            
            # Fallback to cache if enabled
            if fallback_to_cache:
                cached_data = self.cache.get(cache_key)
                if cached_data:
                    logger.info(f"Falling back to cached data for search query '{query}'")
                    return cached_data
            
            # Return error if no cached data available
            return {
                "Query": query,
                "Error": str(e)
            }
    
    def prefetch_molecule_properties(self, cids: List[Union[str, int]]) -> str:
        """
        Schedule a job to prefetch molecule properties for a list of CIDs.
        
        Args:
            cids: List of PubChem Compound IDs
            
        Returns:
            Job ID
        """
        # Convert CIDs to strings
        cid_strs = [str(cid) for cid in cids]
        
        # Schedule job
        job_id = self.scheduler.schedule_job(
            job_type="prefetch",
            params={"cids": cid_strs}
        )
        
        return job_id
    
    def batch_update_cache(self, cids: List[Union[str, int]]) -> str:
        """
        Schedule a job to update cache for a list of CIDs.
        
        Args:
            cids: List of PubChem Compound IDs
            
        Returns:
            Job ID
        """
        # Convert CIDs to strings
        cid_strs = [str(cid) for cid in cids]
        
        # Schedule job
        job_id = self.scheduler.schedule_job(
            job_type="batch_update",
            params={"cids": cid_strs}
        )
        
        return job_id
    
    def get_job_status(self, job_id: str) -> Optional[Dict[str, Any]]:
        """
        Get the status of a scheduled job.
        
        Args:
            job_id: Job ID
            
        Returns:
            Job status or None if job not found
        """
        return self.scheduler.get_job_status(job_id)
    
    def clear_cache(self) -> None:
        """Clear the cache."""
        self.cache.clear()
    
    def get_cache_stats(self) -> Dict[str, int]:
        """
        Get cache statistics.
        
        Returns:
            Dictionary with cache statistics
        """
        return self.cache.get_stats()
    
    def get_rate_limiter_stats(self) -> Dict[str, Any]:
        """
        Get rate limiter statistics.
        
        Returns:
            Dictionary with rate limiter statistics
        """
        return self.rate_limiter.get_stats()
    
    def get_circuit_breaker_stats(self) -> Dict[str, Any]:
        """
        Get circuit breaker statistics.
        
        Returns:
            Dictionary with circuit breaker statistics
        """
        stats = self.circuit_breaker.get_stats()
        
        # Add additional information about circuit state
        if stats["state"] == "open":
            recovery_time = self.recovery_timeout - (time.time() - stats["last_failure_time"])
            stats["recovery_time_remaining"] = max(0, recovery_time)
            stats["recovery_percentage"] = min(100, (recovery_time / self.recovery_timeout) * 100)
        else:
            stats["recovery_time_remaining"] = 0
            stats["recovery_percentage"] = 0
            
        return stats
    
    def reset_circuit_breaker(self) -> None:
        """Reset the circuit breaker to its initial state."""
        self.circuit_breaker.reset()
        logger.info("Circuit breaker has been manually reset")