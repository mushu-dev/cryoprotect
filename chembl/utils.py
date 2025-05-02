"""
Utility functions for ChEMBL API client.

This module provides utility functions for the ChEMBL API client,
including retry logic, circuit breaking, and request key generation.
"""

import time
import logging
import functools
import hashlib
import json
from typing import Dict, Any, Optional, List, Union, Callable, Type, Tuple, Set

logger = logging.getLogger(__name__)

class CircuitBreakerError(Exception):
    """Exception raised when the circuit breaker is open."""
    pass

class CircuitBreaker:
    """
    Circuit breaker pattern implementation.
    
    Prevents repeated calls to a failing service by "opening the circuit"
    after a threshold of failures is reached.
    """
    
    def __init__(
        self,
        failure_threshold: int = 3,
        recovery_timeout: int = 60,
        expected_exceptions: Tuple[Type[Exception], ...] = (Exception,),
        name: str = "default"
    ):
        """
        Initialize the circuit breaker.
        
        Args:
            failure_threshold: Number of failures before opening the circuit
            recovery_timeout: Time in seconds to wait before trying again
            expected_exceptions: Exceptions that count as failures
            name: Name of the circuit breaker for logging
        """
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.expected_exceptions = expected_exceptions
        self.name = name
        
        self.failures = 0
        self.state = "closed"  # closed, open, half-open
        self.last_failure_time = 0
        self.total_failures = 0
        self.successful_calls = 0
        
        logger.info(f"CircuitBreaker '{name}' initialized (threshold: {failure_threshold}, timeout: {recovery_timeout}s)")
    
    def __call__(self, func: Callable) -> Callable:
        """
        Decorator to apply circuit breaker to a function.
        
        Args:
            func: Function to wrap
            
        Returns:
            Wrapped function
        """
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if self.state == "open":
                # Check if recovery timeout has elapsed
                if time.time() - self.last_failure_time >= self.recovery_timeout:
                    logger.info(f"CircuitBreaker '{self.name}' transitioning from open to half-open")
                    self.state = "half-open"
                else:
                    # Circuit is open and timeout hasn't elapsed
                    raise CircuitBreakerError(f"Circuit '{self.name}' is open")
            
            try:
                # Call the function
                result = func(*args, **kwargs)
                
                # If successful and in half-open state, close the circuit
                if self.state == "half-open":
                    logger.info(f"CircuitBreaker '{self.name}' transitioning from half-open to closed")
                    self.state = "closed"
                    self.failures = 0
                
                # Reset failures on success if in closed state
                if self.state == "closed" and self.failures > 0:
                    self.failures = 0
                
                self.successful_calls += 1
                return result
                
            except self.expected_exceptions as e:
                # Count as a failure
                self.failures += 1
                self.total_failures += 1
                self.last_failure_time = time.time()
                
                # If failures exceed threshold, open the circuit
                if self.failures >= self.failure_threshold:
                    if self.state != "open":
                        logger.warning(f"CircuitBreaker '{self.name}' transitioning to open after {self.failures} failures")
                        self.state = "open"
                
                # Re-raise the exception
                raise
        
        return wrapper
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get circuit breaker statistics.
        
        Returns:
            Dictionary with circuit breaker statistics
        """
        return {
            "name": self.name,
            "state": self.state,
            "failures": self.failures,
            "total_failures": self.total_failures,
            "successful_calls": self.successful_calls,
            "last_failure_time": self.last_failure_time,
            "seconds_since_last_failure": time.time() - self.last_failure_time if self.last_failure_time > 0 else None
        }

def retry_with_backoff(
    max_retries: int = 5,
    initial_backoff: float = 1.0,
    backoff_factor: float = 2.0,
    max_backoff: float = 60.0,
    expected_exceptions: Tuple[Type[Exception], ...] = (Exception,)
):
    """
    Decorator for retrying a function with exponential backoff.
    
    Args:
        max_retries: Maximum number of retries
        initial_backoff: Initial backoff time in seconds
        backoff_factor: Factor to increase backoff by on each retry
        max_backoff: Maximum backoff time in seconds
        expected_exceptions: Exceptions to retry on
        
    Returns:
        Decorated function
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            retries = 0
            backoff = initial_backoff
            
            while True:
                try:
                    return func(*args, **kwargs)
                except expected_exceptions as e:
                    retries += 1
                    if retries > max_retries:
                        logger.warning(f"Max retries ({max_retries}) exceeded for {func.__name__}")
                        raise
                    
                    # Calculate backoff time
                    backoff = min(backoff * backoff_factor, max_backoff)
                    logger.info(f"Retry {retries}/{max_retries} for {func.__name__} after {backoff:.2f}s: {str(e)}")
                    
                    # Wait before retrying
                    time.sleep(backoff)
        
        return wrapper
    
    return decorator

def create_request_key(endpoint: str, params: Optional[Dict[str, Any]] = None) -> str:
    """
    Create a unique key for a request to use in caching.
    
    Args:
        endpoint: API endpoint
        params: Request parameters
        
    Returns:
        Unique key for the request
    """
    # Create a string representation of the request
    request_str = endpoint
    if params:
        # Sort params to ensure consistent keys
        sorted_params = sorted(params.items())
        request_str += "?" + "&".join(f"{k}={v}" for k, v in sorted_params)
    
    # Create a hash of the request string
    hash_obj = hashlib.md5(request_str.encode())
    return hash_obj.hexdigest()

def normalize_molecule_data(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize molecule data from ChEMBL to a consistent format.
    
    Args:
        data: Raw molecule data from ChEMBL
        
    Returns:
        Normalized molecule data
    """
    # Extract relevant fields and normalize names
    normalized = {}
    
    # Map ChEMBL fields to our standard fields
    field_mapping = {
        "molecule_chembl_id": "ChEMBL ID",
        "pref_name": "Name",
    }
    
    # Extract top-level fields
    for chembl_field, our_field in field_mapping.items():
        normalized[our_field] = data.get(chembl_field)
    
    # Extract nested fields from molecule_structures
    if "molecule_structures" in data and data["molecule_structures"]:
        structures = data["molecule_structures"]
        normalized["SMILES"] = structures.get("canonical_smiles")
        normalized["InChI"] = structures.get("standard_inchi")
        normalized["InChIKey"] = structures.get("standard_inchi_key")
    
    # Extract nested fields from molecule_properties
    if "molecule_properties" in data and data["molecule_properties"]:
        props = data["molecule_properties"]
        normalized["Molecular Formula"] = props.get("full_molformula")
        normalized["Molecular Weight"] = props.get("full_mwt")
        normalized["LogP"] = props.get("alogp")
        normalized["TPSA"] = props.get("psa")
        normalized["H-Bond Acceptors"] = props.get("hba")
        normalized["H-Bond Donors"] = props.get("hbd")
    
    # Add ChEMBL link
    if "ChEMBL ID" in normalized and normalized["ChEMBL ID"]:
        normalized["ChEMBL Link"] = f"https://www.ebi.ac.uk/chembl/compound_report_card/{normalized['ChEMBL ID']}"
    
    return normalized