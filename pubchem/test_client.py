"""
Unit tests for the ResilientPubChemClient.
"""

import os
import json
import shutil
import unittest
from unittest import mock
import requests
from pathlib import Path

from .client import ResilientPubChemClient
from .utils import CircuitBreakerError

class TestResilientPubChemClient(unittest.TestCase):
    """Test cases for ResilientPubChemClient."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_cache_dir = "test_pubchem_client"
        self.client = ResilientPubChemClient(
            cache_dir=self.test_cache_dir,
            weekday_requests_per_second=10.0,  # High RPS to minimize delays in tests
            weekend_requests_per_second=10.0,
            max_retries=2,
            failure_threshold=2,
            recovery_timeout=1,
            enable_scheduler=False  # Disable scheduler for testing
        )
    
    def tearDown(self):
        """Clean up test environment."""
        # Remove test directory
        if os.path.exists(self.test_cache_dir):
            shutil.rmtree(self.test_cache_dir)
    
    def test_get_molecule_properties_success(self):
        """Test successful molecule properties retrieval."""
        # Mock successful API response
        mock_response = mock.Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "PropertyTable": {
                "Properties": [
                    {
                        "CID": 123,
                        "MolecularFormula": "C2H6O",
                        "MolecularWeight": 46.07,
                        "XLogP": -0.3,
                        "TPSA": 20.2,
                        "HBondDonorCount": 1,
                        "HBondAcceptorCount": 1,
                        "IsomericSMILES": "CCO",
                        "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                        "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                        "IUPACName": "ethanol",
                        "Title": "Ethanol"
                    }
                ]
            }
        }
        
        # Mock requests.get to return our mock response
        with mock.patch('requests.get', return_value=mock_response):
            # Mock rate limiter to avoid delays
            with mock.patch.object(self.client.rate_limiter, 'wait'):
                # Get molecule properties
                result = self.client.get_molecule_properties(123)
                
                # Check result
                self.assertEqual(result["CID"], "123")
                self.assertEqual(result["Molecular Formula"], "C2H6O")
                self.assertEqual(result["Molecular Weight"], 46.07)
                self.assertEqual(result["LogP"], -0.3)
                self.assertEqual(result["TPSA"], 20.2)
                self.assertEqual(result["H-Bond Donors"], 1)
                self.assertEqual(result["H-Bond Acceptors"], 1)
                self.assertEqual(result["SMILES"], "CCO")
                self.assertEqual(result["InChI"], "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
                self.assertEqual(result["InChIKey"], "LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
                self.assertEqual(result["IUPACName"], "ethanol")
                self.assertEqual(result["Title"], "Ethanol")
                self.assertEqual(result["PubChem Link"], "https://pubchem.ncbi.nlm.nih.gov/compound/123")
    
    def test_get_molecule_properties_cache(self):
        """Test molecule properties retrieval with caching."""
        # Mock successful API response
        mock_response = mock.Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "PropertyTable": {
                "Properties": [
                    {
                        "CID": 123,
                        "MolecularFormula": "C2H6O",
                        "MolecularWeight": 46.07,
                        "XLogP": -0.3,
                        "TPSA": 20.2,
                        "HBondDonorCount": 1,
                        "HBondAcceptorCount": 1,
                        "IsomericSMILES": "CCO",
                        "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                        "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                        "IUPACName": "ethanol",
                        "Title": "Ethanol"
                    }
                ]
            }
        }
        
        # Mock requests.get to return our mock response
        with mock.patch('requests.get', return_value=mock_response) as mock_get:
            # Mock rate limiter to avoid delays
            with mock.patch.object(self.client.rate_limiter, 'wait'):
                # Get molecule properties (first call)
                result1 = self.client.get_molecule_properties(123)
                
                # Get molecule properties again (should use cache)
                result2 = self.client.get_molecule_properties(123)
                
                # Check that both results are the same
                self.assertEqual(result1, result2)
                
                # Check that requests.get was called only once
                self.assertEqual(mock_get.call_count, 1)
    
    def test_get_molecule_properties_error(self):
        """Test molecule properties retrieval with API error."""
        # Mock failed API response
        with mock.patch('requests.get', side_effect=requests.RequestException("Test error")):
            # Mock rate limiter to avoid delays
            with mock.patch.object(self.client.rate_limiter, 'wait'):
                # Get molecule properties
                result = self.client.get_molecule_properties(123)
                
                # Check result
                self.assertEqual(result["CID"], "123")
                self.assertIn("Error", result)
                self.assertIn("Test error", result["Error"])
    
    def test_get_molecule_properties_fallback(self):
        """Test fallback to cached data when API fails."""
        # Mock successful API response for first call
        mock_success = mock.Mock()
        mock_success.status_code = 200
        mock_success.json.return_value = {
            "PropertyTable": {
                "Properties": [
                    {
                        "CID": 123,
                        "MolecularFormula": "C2H6O",
                        "MolecularWeight": 46.07,
                        "XLogP": -0.3,
                        "TPSA": 20.2,
                        "HBondDonorCount": 1,
                        "HBondAcceptorCount": 1,
                        "IsomericSMILES": "CCO",
                        "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                        "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                        "IUPACName": "ethanol",
                        "Title": "Ethanol"
                    }
                ]
            }
        }
        
        # Set up mock to return success then fail
        side_effects = [mock_success, requests.RequestException("Test error")]
        
        # Mock requests.get to return our mock responses
        with mock.patch('requests.get', side_effect=side_effects) as mock_get:
            # Mock rate limiter to avoid delays
            with mock.patch.object(self.client.rate_limiter, 'wait'):
                # Get molecule properties (first call, should succeed)
                result1 = self.client.get_molecule_properties(123)
                
                # Get molecule properties again (should fail but fallback to cache)
                result2 = self.client.get_molecule_properties(123)
                
                # Check that both results are the same
                self.assertEqual(result1, result2)
                
                # Check that requests.get was called twice
                self.assertEqual(mock_get.call_count, 2)
    
    def test_circuit_breaker(self):
        """Test circuit breaker functionality."""
        # Mock failed API response
        with mock.patch('requests.get', side_effect=requests.RequestException("Test error")):
            # Mock rate limiter to avoid delays
            with mock.patch.object(self.client.rate_limiter, 'wait'):
                # Call API multiple times to trigger circuit breaker
                for _ in range(2):
                    self.client.get_molecule_properties(123)
                
                # Next call should be rejected by circuit breaker
                with mock.patch.object(self.client, '_make_request', side_effect=CircuitBreakerError("Circuit open")):
                    result = self.client.get_molecule_properties(123)
                    
                    # Check result
                    self.assertEqual(result["CID"], "123")
                    self.assertIn("Error", result)
                    self.assertIn("Circuit open", result["Error"])
    
    def test_get_compound_synonyms(self):
        """Test compound synonyms retrieval."""
        # Mock successful API response
        mock_response = mock.Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "InformationList": {
                "Information": [
                    {
                        "CID": 123,
                        "Synonym": ["Ethanol", "Ethyl alcohol", "Alcohol"]
                    }
                ]
            }
        }
        
        # Mock requests.get to return our mock response
        with mock.patch('requests.get', return_value=mock_response):
            # Mock rate limiter to avoid delays
            with mock.patch.object(self.client.rate_limiter, 'wait'):
                # Get compound synonyms
                result = self.client.get_compound_synonyms(123)
                
                # Check result
                self.assertEqual(result["CID"], "123")
                self.assertEqual(result["Synonyms"], ["Ethanol", "Ethyl alcohol", "Alcohol"])
    
    def test_search_compounds(self):
        """Test compound search."""
        # Mock successful API response
        mock_response = mock.Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "IdentifierList": {
                "CID": [123, 456, 789]
            }
        }
        
        # Mock requests.get to return our mock response
        with mock.patch('requests.get', return_value=mock_response):
            # Mock rate limiter to avoid delays
            with mock.patch.object(self.client.rate_limiter, 'wait'):
                # Search compounds
                result = self.client.search_compounds("ethanol")
                
                # Check result
                self.assertEqual(result["Query"], "ethanol")
                self.assertEqual(result["CIDs"], [123, 456, 789])
    
    def test_prefetch_molecule_properties(self):
        """Test prefetching molecule properties."""
        # Mock scheduler
        with mock.patch.object(self.client.scheduler, 'schedule_job', return_value="job123"):
            # Prefetch molecule properties
            job_id = self.client.prefetch_molecule_properties([123, 456, 789])
            
            # Check job ID
            self.assertEqual(job_id, "job123")
    
    def test_batch_update_cache(self):
        """Test batch updating cache."""
        # Mock scheduler
        with mock.patch.object(self.client.scheduler, 'schedule_job', return_value="job123"):
            # Batch update cache
            job_id = self.client.batch_update_cache([123, 456, 789])
            
            # Check job ID
            self.assertEqual(job_id, "job123")
    
    def test_get_job_status(self):
        """Test getting job status."""
        # Mock scheduler
        with mock.patch.object(self.client.scheduler, 'get_job_status', return_value={"status": "pending"}):
            # Get job status
            status = self.client.get_job_status("job123")
            
            # Check status
            self.assertEqual(status["status"], "pending")
    
    def test_clear_cache(self):
        """Test clearing cache."""
        # Mock cache
        with mock.patch.object(self.client.cache, 'clear') as mock_clear:
            # Clear cache
            self.client.clear_cache()
            
            # Check that cache.clear was called
            mock_clear.assert_called_once()
    
    def test_get_stats(self):
        """Test getting statistics."""
        # Mock component stats
        cache_stats = {"memory_hits": 10, "disk_hits": 5, "misses": 2, "writes": 15}
        rate_limiter_stats = {"total_requests": 20, "weekday_requests": 15, "weekend_requests": 5}
        circuit_breaker_stats = {"name": "pubchem_api", "state": "closed", "failure_count": 0}
        
        with mock.patch.object(self.client.cache, 'get_stats', return_value=cache_stats):
            with mock.patch.object(self.client.rate_limiter, 'get_stats', return_value=rate_limiter_stats):
                with mock.patch.object(self.client.circuit_breaker, 'get_stats', return_value=circuit_breaker_stats):
                    # Get stats
                    cache_result = self.client.get_cache_stats()
                    rate_limiter_result = self.client.get_rate_limiter_stats()
                    circuit_breaker_result = self.client.get_circuit_breaker_stats()
                    
                    # Check stats
                    self.assertEqual(cache_result, cache_stats)
                    self.assertEqual(rate_limiter_result, rate_limiter_stats)
                    self.assertEqual(circuit_breaker_result, circuit_breaker_stats)


if __name__ == "__main__":
    unittest.main()