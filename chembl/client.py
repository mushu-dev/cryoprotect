"""
Resilient ChEMBL API client for CryoProtect v2.

This module provides a robust client for interacting with the ChEMBL API,
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

from .cache import ChEMBLCache
from .rate_limiter import AdaptiveRateLimiter
from .utils import (
    CircuitBreaker, 
    retry_with_backoff, 
    create_request_key,
    normalize_molecule_data,
    CircuitBreakerError
)

logger = logging.getLogger(__name__)

class ResilientChEMBLClient:
    """
    Resilient client for ChEMBL API.
    
    Features:
    - Adaptive rate limiting
    - Multi-level caching (in-memory and disk-based)
    - Exponential backoff retry logic
    - Circuit breaker to prevent repeated failures
    - Fallback to pre-cached data when API is unavailable
    """
    
    def __init__(
        self,
        cache_dir: str = "cache/chembl",
        weekday_requests_per_second: float = 3.0,
        weekend_requests_per_second: float = 4.5,
        monday_requests_per_second: float = 2.0,  # Stricter limit for Mondays
        max_retries: int = 5,
        failure_threshold: int = 3,
        recovery_timeout: int = 60,
        cache_ttl: int = 86400 * 30,  # 30 days
        error_cache_ttl: int = 3600,  # 1 hour for error responses
        memory_cache_size: int = 1000
    ):
        """
        Initialize the client.
        
        Args:
            cache_dir: Directory to store cache files
            weekday_requests_per_second: Maximum requests per second on weekdays (Tue-Fri)
            weekend_requests_per_second: Maximum requests per second on weekends
            monday_requests_per_second: Maximum requests per second on Mondays (stricter)
            max_retries: Maximum number of retries for failed requests
            failure_threshold: Number of failures before opening the circuit
            recovery_timeout: Time in seconds to wait before trying again
            cache_ttl: Time-to-live for cache entries in seconds
            error_cache_ttl: Time-to-live for error responses in seconds (shorter)
            memory_cache_size: Maximum number of items to keep in memory cache
        """
        # Create cache directory if it doesn't exist
        os.makedirs(cache_dir, exist_ok=True)
        
        # Initialize components
        self.cache = ChEMBLCache(
            cache_dir=cache_dir,
            memory_size=memory_cache_size,
            ttl=cache_ttl,
            error_ttl=error_cache_ttl,
            max_retries=3
        )
        
        self.rate_limiter = AdaptiveRateLimiter(
            weekday_requests_per_second=weekday_requests_per_second,
            weekend_requests_per_second=weekend_requests_per_second,
            monday_requests_per_second=monday_requests_per_second
        )
        
        # Circuit breaker for API requests
        self.circuit_breaker = CircuitBreaker(
            failure_threshold=failure_threshold,
            recovery_timeout=recovery_timeout,
            expected_exceptions=(requests.RequestException, json.JSONDecodeError),
            name="chembl_api"
        )
        
        # Track consecutive errors for adaptive backoff
        self.consecutive_errors = 0
        self.last_error_time = 0
        
        # Configuration
        self.max_retries = max_retries
        self.base_url = "https://www.ebi.ac.uk/chembl/api/data"
        
        logger.info("ResilientChEMBLClient initialized")
    
    @retry_with_backoff(
        max_retries=5,
        initial_backoff=1.0,
        backoff_factor=2.0,
        max_backoff=60.0,
        expected_exceptions=(requests.RequestException, CircuitBreakerError, json.JSONDecodeError)
    )
    def _make_request(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Make a request to the ChEMBL API with circuit breaker protection.
        
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
        
        # Wait for rate limiter
        self.rate_limiter.wait()
        
        # Make request with circuit breaker
        @self.circuit_breaker
        def do_request():
            try:
                logger.debug(f"Making request to {url} with params {params}")
                # Set Accept header to request JSON format
                headers = {
                    "Accept": "application/json"
                }
                response = requests.get(url, params=params, headers=headers)
                response.raise_for_status()
                
                # Check if response is valid JSON
                try:
                    result = response.json()
                    return result
                except json.JSONDecodeError as e:
                    logger.error(f"Invalid JSON response: {response.text[:200]}...")
                    logger.error(f"JSON decode error: {str(e)}")
                    raise
            except requests.RequestException as e:
                logger.error(f"Request error: {str(e)}")
                # Track consecutive errors for adaptive rate limiting
                self.consecutive_errors += 1
                self.last_error_time = time.time()
                # Report rate limit error to rate limiter if it's a 429 Too Many Requests
                if hasattr(e, 'response') and e.response is not None and e.response.status_code == 429:
                    self.rate_limiter.report_rate_limit_error()
                    logger.warning("Rate limit exceeded, backing off")
                raise
        
        try:
            result = do_request()
            
            # Report success to rate limiter
            self.rate_limiter.report_success()
            
            # Reset consecutive errors on success
            if self.consecutive_errors > 0:
                self.consecutive_errors = 0
                
            return result
            
        except Exception as e:
            # If this is a circuit breaker error, we should check if we can use cached data
            if isinstance(e, CircuitBreakerError):
                logger.warning(f"Circuit breaker open: {str(e)}")
                # We'll let the calling method handle fallback to cache
            
            # Re-raise the exception
            raise
    
    def get_molecule_by_chembl_id(
        self, 
        chembl_id: str,
        use_cache: bool = True,
        fallback_to_cache: bool = True
    ) -> Dict[str, Any]:
        """
        Fetch molecule data from ChEMBL by ChEMBL ID.
        
        Args:
            chembl_id: ChEMBL ID (e.g., "CHEMBL25")
            use_cache: Whether to use cached data if available
            fallback_to_cache: Whether to fallback to cached data if API fails
            
        Returns:
            Dictionary with molecule data
        """
        # Create cache key
        endpoint = f"molecule/{chembl_id}"
        cache_key = create_request_key(endpoint)
        
        # Check cache if enabled
        if use_cache:
            cached_data = self.cache.get(cache_key)
            if cached_data:
                logger.debug(f"Cache hit for ChEMBL ID {chembl_id}")
                # Check if circuit is open and we're using cached data
                if self.circuit_breaker.state == "open":
                    logger.info(f"Using cached data for {chembl_id} due to open circuit")
                return cached_data
        
        try:
            # Make request to ChEMBL API
            response = self._make_request(endpoint)
            
            # Normalize the data
            result = normalize_molecule_data(response)
            
            # Cache result
            self.cache.set(cache_key, result)
            
            return result
            
        except Exception as e:
            error_msg = f"Error fetching molecule for ChEMBL ID {chembl_id}: {str(e)}"
            logger.warning(error_msg)
            
            # Fallback to cache if enabled
            if fallback_to_cache:
                cached_data = self.cache.get(cache_key)
                if cached_data:
                    logger.info(f"Falling back to cached data for ChEMBL ID {chembl_id}")
                    return cached_data
            
            # Create error response
            error_response = {
                "ChEMBL ID": chembl_id,
                "Error": str(e),
                "ErrorType": type(e).__name__,
                "Timestamp": datetime.now().isoformat()
            }
            
            # Cache the error response with a shorter TTL
            if use_cache:
                self.cache.set(cache_key, error_response)
                
            return error_response
    
    def search_molecules(
        self,
        query: str,
        limit: int = 10,
        use_cache: bool = True,
        fallback_to_cache: bool = True
    ) -> Dict[str, Any]:
        """
        Search for molecules by name, SMILES, InChI, etc.
        
        Args:
            query: Search query
            limit: Maximum number of results to return
            use_cache: Whether to use cached data if available
            fallback_to_cache: Whether to fallback to cached data if API fails
            
        Returns:
            Dictionary with search results
        """
        # Create cache key
        endpoint = "molecule"
        params = {
            "pref_name__icontains": query,
            "limit": limit
        }
        cache_key = create_request_key(endpoint, params)
        
        # Check cache if enabled
        if use_cache:
            cached_data = self.cache.get(cache_key)
            if cached_data:
                logger.debug(f"Cache hit for search query '{query}'")
                return cached_data
        
        try:
            # Make request to ChEMBL API
            response = self._make_request(endpoint, params)
            
            # Extract molecules
            molecules = []
            for mol in response.get("molecules", []):
                molecules.append(normalize_molecule_data(mol))
            
            # Format result
            result = {
                "Query": query,
                "Count": len(molecules),
                "Molecules": molecules
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
            
            # Create error response
            error_response = {
                "Query": query,
                "Error": str(e),
                "ErrorType": type(e).__name__,
                "Timestamp": datetime.now().isoformat(),
                "Count": 0,
                "Molecules": []
            }
            
            # Cache the error response with a shorter TTL
            if use_cache:
                self.cache.set(cache_key, error_response)
                
            return error_response
    
    def get_molecule_properties(
        self,
        chembl_id: str,
        use_cache: bool = True,
        fallback_to_cache: bool = True
    ) -> Dict[str, Any]:
        """
        Fetch molecular properties from ChEMBL.
        
        Args:
            chembl_id: ChEMBL ID
            use_cache: Whether to use cached data if available
            fallback_to_cache: Whether to fallback to cached data if API fails
            
        Returns:
            Dictionary with molecular properties
        """
        # Get the molecule data - all properties are already included in the molecule data
        molecule_data = self.get_molecule_by_chembl_id(
            chembl_id=chembl_id,
            use_cache=use_cache,
            fallback_to_cache=fallback_to_cache
        )
        
        # Check for error
        if "Error" in molecule_data:
            return molecule_data
        
        # Extract additional properties if available
        properties = {}
        
        # Add any additional properties you want to extract here
        # For now, we'll just return the molecule data as is
        
        return molecule_data
    
    def get_molecule_by_inchikey(
        self, 
        inchikey: str,
        use_cache: bool = True,
        fallback_to_cache: bool = True
    ) -> Dict[str, Any]:
        """
        Fetch molecule data from ChEMBL by InChIKey.
        
        Args:
            inchikey: InChIKey
            use_cache: Whether to use cached data if available
            fallback_to_cache: Whether to fallback to cached data if API fails
            
        Returns:
            Dictionary with molecule data
        """
        # Create cache key
        endpoint = "molecule"
        params = {
            "molecule_structures__standard_inchi_key": inchikey,
            "limit": 1
        }
        cache_key = create_request_key(endpoint, params)
        
        # Check cache if enabled
        if use_cache:
            cached_data = self.cache.get(cache_key)
            if cached_data:
                logger.debug(f"Cache hit for InChIKey {inchikey}")
                return cached_data
        
        try:
            # Make request to ChEMBL API
            response = self._make_request(endpoint, params)
            
            # Extract first molecule
            molecules = response.get("molecules", [])
            if not molecules:
                return {
                    "InChIKey": inchikey,
                    "Error": "No molecule found with this InChIKey"
                }
            
            # Normalize the data
            result = normalize_molecule_data(molecules[0])
            
            # Cache result
            self.cache.set(cache_key, result)
            
            return result
            
        except Exception as e:
            error_msg = f"Error fetching molecule for InChIKey {inchikey}: {str(e)}"
            logger.warning(error_msg)
            
            # Fallback to cache if enabled
            if fallback_to_cache:
                cached_data = self.cache.get(cache_key)
                if cached_data:
                    logger.info(f"Falling back to cached data for InChIKey {inchikey}")
                    return cached_data
            
            # Create error response
            error_response = {
                "InChIKey": inchikey,
                "Error": str(e),
                "ErrorType": type(e).__name__,
                "Timestamp": datetime.now().isoformat()
            }
            
            # Cache the error response with a shorter TTL
            if use_cache:
                self.cache.set(cache_key, error_response)
                
            return error_response
    
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
        return self.circuit_breaker.get_stats()
    
    def reset_circuit_breaker(self) -> None:
        """
        Manually reset the circuit breaker to closed state.
        This can be useful for testing or after fixing issues.
        """
        if self.circuit_breaker.state != "closed":
            logger.info(f"Manually resetting circuit breaker from {self.circuit_breaker.state} to closed")
            self.circuit_breaker.state = "closed"
            self.circuit_breaker.failures = 0
            self.consecutive_errors = 0
            
    def get_similar_compounds(
        self,
        chembl_id: str,
        similarity: int = 70,
        limit: int = 20,
        use_cache: bool = True,
        fallback_to_cache: bool = True
    ) -> List[str]:
        """
        Fetch similar compounds to a given ChEMBL ID based on structural similarity.
        
        Args:
            chembl_id: ChEMBL ID of the reference compound
            similarity: Minimum similarity percentage (0-100)
            limit: Maximum number of similar compounds to return
            use_cache: Whether to use cached data if available
            fallback_to_cache: Whether to fallback to cached data if API fails
            
        Returns:
            List of ChEMBL IDs of similar compounds
        """
        # Create cache key
        endpoint = f"similarity/{chembl_id}/{similarity}"
        params = {
            "limit": limit
        }
        cache_key = create_request_key(endpoint, params)
        
        # Check cache if enabled
        if use_cache:
            cached_data = self.cache.get(cache_key)
            if cached_data:
                logger.debug(f"Cache hit for similar compounds to {chembl_id}")
                return cached_data
        
        try:
            # Make request to ChEMBL API
            response = self._make_request(endpoint, params)
            
            # Extract ChEMBL IDs of similar molecules
            similar_ids = []
            for mol in response.get("molecules", []):
                mol_id = mol.get("molecule_chembl_id")
                if mol_id and mol_id != chembl_id:  # Exclude the reference compound itself
                    similar_ids.append(mol_id)
            
            # Cache result
            self.cache.set(cache_key, similar_ids)
            
            return similar_ids
            
        except Exception as e:
            error_msg = f"Error fetching similar compounds for {chembl_id}: {str(e)}"
            logger.warning(error_msg)
            
            # Fallback to cache if enabled
            if fallback_to_cache:
                cached_data = self.cache.get(cache_key)
                if cached_data:
                    logger.info(f"Falling back to cached data for similar compounds to {chembl_id}")
                    return cached_data
            
            # Return empty list on error
            return []

# Simplified ChEMBL client for backward compatibility
class ChEMBLClient:
    """
    Simplified ChEMBL client that uses the ResilientChEMBLClient internally.
    This class provides backward compatibility with existing code.
    """
    
    def __init__(self, base_url: str = "https://www.ebi.ac.uk/chembl/api/data"):
        """
        Initialize the ChEMBL client.
        
        Args:
            base_url: Base URL for the ChEMBL API
        """
        self.resilient_client = ResilientChEMBLClient()
        self.base_url = base_url
        
    def get_compound(self, chembl_id: str) -> Dict[str, Any]:
        """
        Get compound data by ChEMBL ID.
        
        Args:
            chembl_id: ChEMBL ID
            
        Returns:
            Dictionary with compound data in the format expected by property filtering functions:
            {
                'molecule_structures': {
                    'canonical_smiles': '...'
                },
                'molecule_chembl_id': 'CHEMBLXXXX',
                ...
            }
        """
        try:
            # Get molecule data from resilient client
            molecule_data = self.resilient_client.get_molecule_by_chembl_id(chembl_id)
            
            # Check if there was an error
            if "Error" in molecule_data:
                logger.warning(f"Error retrieving compound {chembl_id}: {molecule_data.get('Error')}")
                return {
                    "molecule_chembl_id": chembl_id,
                    "error": molecule_data.get("Error", "Unknown error")
                }
            
            # Ensure the data has the expected structure
            if "molecule_structures" not in molecule_data:
                # Create the expected structure
                molecule_data["molecule_structures"] = {}
            
            # Make sure canonical_smiles is present
            if "canonical_smiles" not in molecule_data.get("molecule_structures", {}):
                # Try to get SMILES from other fields if available
                smiles = None
                
                # Check if we have a SMILES string in another location
                if "SMILES" in molecule_data:
                    smiles = molecule_data["SMILES"]
                elif "smiles" in molecule_data:
                    smiles = molecule_data["smiles"]
                elif "canonical_smiles" in molecule_data:
                    smiles = molecule_data["canonical_smiles"]
                
                # Add it to the expected location
                if smiles:
                    molecule_data["molecule_structures"]["canonical_smiles"] = smiles
            
            # Ensure molecule_chembl_id is present
            if "molecule_chembl_id" not in molecule_data and "ChEMBL ID" in molecule_data:
                molecule_data["molecule_chembl_id"] = molecule_data["ChEMBL ID"]
            
            return molecule_data
            
        except Exception as e:
            logger.error(f"Error in get_compound for {chembl_id}: {str(e)}")
            # Return a minimal structure that won't cause errors in property filtering
            return {
                "molecule_chembl_id": chembl_id,
                "error": str(e),
                "molecule_structures": {}
            }
        
    def search_compounds(self, term: str, limit: int = 100) -> List[str]:
        """
        Search for compounds by term and return their ChEMBL IDs.
        
        Args:
            term: Search term
            limit: Maximum number of results to return
            
        Returns:
            List of ChEMBL IDs
        """
        results = self.resilient_client.search_molecules(term, limit)
        chembl_ids = []
        
        for mol in results.get("Molecules", []):
            if "molecule_chembl_id" in mol:
                chembl_ids.append(mol["molecule_chembl_id"])
                
        return chembl_ids
        
    def get_similar_compounds(self, chembl_id: str, similarity: int = 70, limit: int = 20) -> List[str]:
        """
        Get similar compounds to a given ChEMBL ID.
        
        Args:
            chembl_id: ChEMBL ID
            similarity: Minimum similarity percentage (0-100)
            limit: Maximum number of similar compounds to return
            
        Returns:
            List of ChEMBL IDs of similar compounds
        """
        return self.resilient_client.get_similar_compounds(chembl_id, similarity, limit)