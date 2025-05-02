"""
Toxicity Scoring System for CryoProtect v2.

This module implements a toxicity scoring system based on Tox21 data.
It provides functionality for:
- Calculating toxicity scores for individual assays
- Grouping assays by toxicological endpoints
- Computing endpoint-specific scores
- Calculating weighted aggregate toxicity scores
- Providing confidence and data completeness metrics
- Handling missing data appropriately

The scoring system is designed to be scientifically sound, efficient for large datasets,
and well-documented.
"""

import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any, Union, Set
from datetime import datetime
from supabase import Client, create_client
import os
from pathlib import Path

# Local imports
from config import Config
from chemical_data.toxicity.endpoint_classification import (
    classify_assay,
    get_endpoint_weight,
    get_all_endpoints,
    get_endpoint_description,
    classify_assays_in_database,
    get_assays_by_endpoint
)

# Set up logging
logger = logging.getLogger(__name__)

class ToxicityScorer:
    """
    Class for calculating toxicity scores based on Tox21 data.
    """

    def __init__(self, supabase_client: Optional[Client] = None, config: Optional[Config] = None):
        """
        Initialize the ToxicityScorer.
        
        Args:
            supabase_client: Optional Supabase client
            config: Optional configuration
        """
        self.config = config or Config()
        
        # Initialize Supabase client if not provided
        if supabase_client:
            self.supabase = supabase_client
        else:
            supabase_url = os.environ.get("SUPABASE_URL") or self.config.SUPABASE_URL
            supabase_key = os.environ.get("SUPABASE_KEY") or self.config.SUPABASE_KEY
            self.supabase = create_client(supabase_url, supabase_key)
        
        # Get Tox21 source ID from database
        self.tox21_source_id = self._get_tox21_source_id()
        
        # Get toxicity calculation method ID
        self.calculation_method_id = self._get_toxicity_calculation_method_id()
    
    def _get_tox21_source_id(self) -> str:
        """Get the Tox21 source ID from the database."""
        try:
            response = self.supabase.table("toxicity_data_source").select("id").eq("name", "Tox21").execute()
            if response.data and len(response.data) > 0:
                return response.data[0]["id"]
            else:
                logger.error("Tox21 data source not found in database")
                raise ValueError("Tox21 data source not found in database. Run migrations first.")
        except Exception as e:
            logger.error(f"Error getting Tox21 source ID: {str(e)}")
            raise
    
    def _get_toxicity_calculation_method_id(self) -> str:
        """Get or create the toxicity calculation method ID."""
        method_name = "Tox21 Endpoint-Based Scoring"
        method_description = "Toxicity scoring based on Tox21 assays grouped by toxicological endpoints"
        
        try:
            # Check if method already exists
            response = self.supabase.table("calculation_method").select("id").eq("name", method_name).execute()
            
            if response.data and len(response.data) > 0:
                return response.data[0]["id"]
            
            # Create new method
            response = self.supabase.table("calculation_method").insert({
                "name": method_name,
                "description": method_description,
                "method_type": "toxicity_scoring",
                "version": "1.0",
                "parameters": {
                    "data_source": "Tox21",
                    "endpoint_weights": {endpoint: weight for endpoint, weight in zip(get_all_endpoints(), [get_endpoint_weight(ep) for ep in get_all_endpoints()])}
                }
            }).execute()
            
            if response.data and len(response.data) > 0:
                return response.data[0]["id"]
            else:
                logger.error("Failed to create toxicity calculation method")
                raise ValueError("Failed to create toxicity calculation method")
        except Exception as e:
            logger.error(f"Error getting/creating toxicity calculation method: {str(e)}")
            raise
    
    def calculate_toxicity_scores(self, molecule_ids: Optional[List[str]] = None) -> int:
        """
        Calculate toxicity scores for molecules.
        
        Args:
            molecule_ids: Optional list of molecule IDs to calculate scores for.
                          If None, scores will be calculated for all molecules with Tox21 data.
        
        Returns:
            Number of molecules scored
        """
        try:
            # Ensure assays are classified by endpoint
            classify_assays_in_database(self.supabase, self.tox21_source_id)
            
            # Get molecules with Tox21 data
            if molecule_ids:
                # Get specific molecules
                molecules_query = self.supabase.table("toxicity_data").select("distinct(molecule_id)").eq("source_id", self.tox21_source_id).in_("molecule_id", molecule_ids).execute()
            else:
                # Get all molecules with Tox21 data
                molecules_query = self.supabase.table("toxicity_data").select("distinct(molecule_id)").eq("source_id", self.tox21_source_id).execute()
            
            if not molecules_query.data:
                logger.warning("No molecules found with Tox21 data")
                return 0
            
            # Extract molecule IDs
            molecules_to_score = [item["molecule_id"] for item in molecules_query.data]
            logger.info(f"Calculating toxicity scores for {len(molecules_to_score)} molecules")
            
            # Process molecules in batches
            batch_size = 50
            molecules_scored = 0
            
            for i in range(0, len(molecules_to_score), batch_size):
                batch = molecules_to_score[i:i+batch_size]
                logger.info(f"Processing batch {i//batch_size + 1}/{(len(molecules_to_score) + batch_size - 1)//batch_size}")
                
                for molecule_id in batch:
                    try:
                        # Calculate toxicity score for the molecule
                        score_result = self._calculate_molecule_toxicity_score(molecule_id)
                        
                        if score_result["success"]:
                            molecules_scored += 1
                            
                            if molecules_scored % 10 == 0:
                                logger.info(f"Scored {molecules_scored}/{len(molecules_to_score)} molecules")
                    except Exception as e:
                        logger.error(f"Error calculating toxicity score for molecule {molecule_id}: {str(e)}")
            
            logger.info(f"Successfully calculated toxicity scores for {molecules_scored} molecules")
            return molecules_scored
        
        except Exception as e:
            logger.error(f"Error calculating toxicity scores: {str(e)}")
            raise
    
    def _calculate_molecule_toxicity_score(self, molecule_id: str) -> Dict[str, Any]:
        """
        Calculate toxicity score for a single molecule.
        
        Args:
            molecule_id: Molecule ID
        
        Returns:
            Dictionary with score results
        """
        try:
            # Get all toxicity data for the molecule
            response = self.supabase.table("toxicity_data").select("id, assay_id, value, is_active").eq("molecule_id", molecule_id).eq("source_id", self.tox21_source_id).execute()
            
            if not response.data:
                logger.warning(f"No toxicity data found for molecule {molecule_id}")
                return {"success": False, "reason": "no_data"}
            
            # Get assay information
            assay_ids = [item["assay_id"] for item in response.data]
            assay_response = self.supabase.table("toxicity_assay").select("id, toxicological_endpoint").in_("id", assay_ids).execute()
            
            if not assay_response.data:
                logger.warning(f"No assay information found for molecule {molecule_id}")
                return {"success": False, "reason": "no_assay_info"}
            
            # Create assay endpoint mapping
            assay_endpoints = {assay["id"]: assay.get("toxicological_endpoint", "other") for assay in assay_response.data}
            
            # Group toxicity data by endpoint
            endpoint_data = {}
            for item in response.data:
                endpoint = assay_endpoints.get(item["assay_id"], "other")
                
                if endpoint not in endpoint_data:
                    endpoint_data[endpoint] = []
                
                endpoint_data[endpoint].append(item)
            
            # Calculate endpoint-specific scores
            endpoint_scores = {}
            endpoint_confidences = {}
            endpoint_completeness = {}
            
            for endpoint, data in endpoint_data.items():
                score, confidence, completeness = self._calculate_endpoint_score(data)
                endpoint_scores[endpoint] = score
                endpoint_confidences[endpoint] = confidence
                endpoint_completeness[endpoint] = completeness
            
            # Calculate overall toxicity score
            overall_score, overall_confidence, overall_completeness = self._calculate_overall_toxicity_score(
                endpoint_scores, endpoint_confidences, endpoint_completeness
            )
            
            # Store the results
            score_data = {
                "molecule_id": molecule_id,
                "calculation_method_id": self.calculation_method_id,
                "toxicity_score": overall_score,
                "confidence_score": overall_confidence,
                "data_completeness": overall_completeness,
                "calculation_date": datetime.now().isoformat(),
                "endpoint_scores": {
                    endpoint: {
                        "score": score,
                        "confidence": endpoint_confidences.get(endpoint, 0.0),
                        "completeness": endpoint_completeness.get(endpoint, 0.0)
                    } for endpoint, score in endpoint_scores.items()
                }
            }
            
            # Save to database
            response = self.supabase.table("molecule_toxicity_score").upsert(score_data).execute()
            
            if not response.data:
                logger.warning(f"Failed to save toxicity score for molecule {molecule_id}")
                return {"success": False, "reason": "save_failed"}
            
            return {
                "success": True,
                "molecule_id": molecule_id,
                "overall_score": overall_score,
                "overall_confidence": overall_confidence,
                "overall_completeness": overall_completeness,
                "endpoint_scores": endpoint_scores
            }
        
        except Exception as e:
            logger.error(f"Error calculating toxicity score for molecule {molecule_id}: {str(e)}")
            return {"success": False, "reason": str(e)}
    
    def _calculate_endpoint_score(self, toxicity_data: List[Dict]) -> Tuple[float, float, float]:
        """
        Calculate toxicity score for a specific endpoint.
        
        Args:
            toxicity_data: List of toxicity data points for the endpoint
        
        Returns:
            Tuple of (score, confidence, completeness)
        """
        if not toxicity_data:
            return 0.0, 0.0, 0.0
        
        # Extract active/inactive values
        active_values = []
        for item in toxicity_data:
            if item.get("is_active") is not None:
                active_values.append(1.0 if item["is_active"] else 0.0)
            elif item.get("value") is not None:
                # If binary activity not available, use normalized value
                # Assuming values are normalized between 0 and 1
                active_values.append(float(item["value"]))
        
        if not active_values:
            return 0.0, 0.0, 0.0
        
        # Calculate score as the mean of active values
        score = sum(active_values) / len(active_values)
        
        # Calculate confidence based on number of data points and consistency
        num_points = len(active_values)
        consistency = 1.0 - np.std(active_values) if len(active_values) > 1 else 1.0
        
        # Confidence increases with more data points and higher consistency
        confidence = min(1.0, (0.5 * min(1.0, num_points / 5.0) + 0.5 * consistency))
        
        # Calculate data completeness
        # Assuming a "complete" endpoint has at least 5 assays
        completeness = min(1.0, num_points / 5.0)
        
        return score, confidence, completeness
    
    def _calculate_overall_toxicity_score(
        self,
        endpoint_scores: Dict[str, float],
        endpoint_confidences: Dict[str, float],
        endpoint_completeness: Dict[str, float]
    ) -> Tuple[float, float, float]:
        """
        Calculate overall toxicity score from endpoint-specific scores.
        
        Args:
            endpoint_scores: Dictionary mapping endpoints to scores
            endpoint_confidences: Dictionary mapping endpoints to confidence scores
            endpoint_completeness: Dictionary mapping endpoints to completeness scores
        
        Returns:
            Tuple of (overall_score, overall_confidence, overall_completeness)
        """
        if not endpoint_scores:
            return 0.0, 0.0, 0.0
        
        # Get all possible endpoints
        all_endpoints = get_all_endpoints()
        
        # Calculate weighted sum of endpoint scores
        weighted_scores = []
        weights = []
        
        for endpoint in endpoint_scores:
            if endpoint in all_endpoints:
                weight = get_endpoint_weight(endpoint)
                score = endpoint_scores[endpoint]
                confidence = endpoint_confidences.get(endpoint, 0.0)
                
                # Weight by both endpoint importance and confidence
                adjusted_weight = weight * confidence
                weighted_scores.append(score * adjusted_weight)
                weights.append(adjusted_weight)
        
        if not weighted_scores:
            return 0.0, 0.0, 0.0
        
        # Calculate weighted average
        overall_score = sum(weighted_scores) / sum(weights) if sum(weights) > 0 else 0.0
        
        # Calculate overall confidence
        # Based on average confidence and coverage of endpoints
        avg_confidence = sum(endpoint_confidences.values()) / len(endpoint_confidences) if endpoint_confidences else 0.0
        endpoint_coverage = len(endpoint_scores) / len(all_endpoints)
        overall_confidence = 0.7 * avg_confidence + 0.3 * endpoint_coverage
        
        # Calculate overall completeness
        # Based on average completeness and coverage of endpoints
        avg_completeness = sum(endpoint_completeness.values()) / len(endpoint_completeness) if endpoint_completeness else 0.0
        overall_completeness = 0.7 * avg_completeness + 0.3 * endpoint_coverage
        
        return overall_score, overall_confidence, overall_completeness
    
    def get_molecule_toxicity_score(self, molecule_id: str) -> Dict[str, Any]:
        """
        Get the toxicity score for a molecule.
        
        Args:
            molecule_id: Molecule ID
        
        Returns:
            Dictionary with toxicity score information
        """
        try:
            response = self.supabase.table("molecule_toxicity_score").select("*").eq("molecule_id", molecule_id).eq("calculation_method_id", self.calculation_method_id).execute()
            
            if not response.data:
                logger.warning(f"No toxicity score found for molecule {molecule_id}")
                return None
            
            return response.data[0]
        except Exception as e:
            logger.error(f"Error getting toxicity score for molecule {molecule_id}: {str(e)}")
            return None
    
    def get_toxicity_scores_by_threshold(self, threshold: float, comparison: str = ">=") -> List[Dict[str, Any]]:
        """
        Get molecules with toxicity scores meeting a threshold.
        
        Args:
            threshold: Toxicity score threshold
            comparison: Comparison operator (">", ">=", "<", "<=", "=")
        
        Returns:
            List of dictionaries with molecule IDs and toxicity scores
        """
        try:
            # Construct filter based on comparison
            if comparison == ">":
                response = self.supabase.table("molecule_toxicity_score").select("molecule_id, toxicity_score, confidence_score").eq("calculation_method_id", self.calculation_method_id).gt("toxicity_score", threshold).execute()
            elif comparison == ">=":
                response = self.supabase.table("molecule_toxicity_score").select("molecule_id, toxicity_score, confidence_score").eq("calculation_method_id", self.calculation_method_id).gte("toxicity_score", threshold).execute()
            elif comparison == "<":
                response = self.supabase.table("molecule_toxicity_score").select("molecule_id, toxicity_score, confidence_score").eq("calculation_method_id", self.calculation_method_id).lt("toxicity_score", threshold).execute()
            elif comparison == "<=":
                response = self.supabase.table("molecule_toxicity_score").select("molecule_id, toxicity_score, confidence_score").eq("calculation_method_id", self.calculation_method_id).lte("toxicity_score", threshold).execute()
            elif comparison == "=":
                response = self.supabase.table("molecule_toxicity_score").select("molecule_id, toxicity_score, confidence_score").eq("calculation_method_id", self.calculation_method_id).eq("toxicity_score", threshold).execute()
            else:
                raise ValueError(f"Invalid comparison operator: {comparison}")
            
            return response.data if response.data else []
        except Exception as e:
            logger.error(f"Error getting toxicity scores by threshold: {str(e)}")
            return []
    
    def get_endpoint_specific_scores(self, molecule_id: str) -> Dict[str, Dict[str, float]]:
        """
        Get endpoint-specific toxicity scores for a molecule.
        
        Args:
            molecule_id: Molecule ID
        
        Returns:
            Dictionary mapping endpoints to score information
        """
        try:
            response = self.supabase.table("molecule_toxicity_score").select("endpoint_scores").eq("molecule_id", molecule_id).eq("calculation_method_id", self.calculation_method_id).execute()
            
            if not response.data or not response.data[0].get("endpoint_scores"):
                logger.warning(f"No endpoint scores found for molecule {molecule_id}")
                return {}
            
            return response.data[0]["endpoint_scores"]
        except Exception as e:
            logger.error(f"Error getting endpoint scores for molecule {molecule_id}: {str(e)}")
            return {}