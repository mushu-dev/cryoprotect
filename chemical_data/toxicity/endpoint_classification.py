"""
Endpoint Classification Module for CryoProtect v2.

This module provides functionality for classifying Tox21 assays into toxicological endpoints
and grouping related assays for toxicity scoring. It includes functions for:
- Mapping assays to toxicological endpoints
- Grouping assays by endpoint
- Determining endpoint weights for aggregate scoring
"""

import logging
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Union, Set
from supabase import Client

# Set up logging
logger = logging.getLogger(__name__)

# Define toxicological endpoints and their descriptions
TOXICOLOGICAL_ENDPOINTS = {
    "nuclear_receptor": "Nuclear receptor activity (e.g., androgen, estrogen, glucocorticoid receptors)",
    "stress_response": "Stress response pathways (e.g., oxidative stress, heat shock)",
    "dna_damage": "DNA damage and repair pathways",
    "cell_cycle": "Cell cycle regulation and proliferation",
    "apoptosis": "Programmed cell death pathways",
    "mitochondrial_function": "Mitochondrial function and integrity",
    "enzyme_inhibition": "Enzyme inhibition (e.g., CYP450 enzymes)",
    "cell_signaling": "Cell signaling pathways",
    "immunotoxicity": "Immune system effects",
    "neurotoxicity": "Neurological system effects",
    "developmental_toxicity": "Developmental and reproductive toxicity",
    "other": "Other toxicological endpoints"
}

# Define endpoint weights for aggregate scoring
# Higher weights indicate greater importance in the overall toxicity score
ENDPOINT_WEIGHTS = {
    "nuclear_receptor": 1.0,
    "stress_response": 1.0,
    "dna_damage": 1.2,
    "cell_cycle": 0.9,
    "apoptosis": 1.1,
    "mitochondrial_function": 1.0,
    "enzyme_inhibition": 0.8,
    "cell_signaling": 0.8,
    "immunotoxicity": 1.1,
    "neurotoxicity": 1.2,
    "developmental_toxicity": 1.3,
    "other": 0.7
}

# Tox21 assay to endpoint mapping
# This maps specific Tox21 assay IDs or patterns to toxicological endpoints
TOX21_ASSAY_ENDPOINT_MAPPING = {
    # Nuclear receptor assays
    "tox21-er-bla-agonist": "nuclear_receptor",
    "tox21-er-bla-antagonist": "nuclear_receptor",
    "tox21-ar-bla-agonist": "nuclear_receptor",
    "tox21-ar-bla-antagonist": "nuclear_receptor",
    "tox21-aromatase": "nuclear_receptor",
    "tox21-pparg-bla-agonist": "nuclear_receptor",
    "tox21-ppard-bla-agonist": "nuclear_receptor",
    "tox21-ahr-luc-agonist": "nuclear_receptor",
    
    # Stress response assays
    "tox21-are-bla": "stress_response",
    "tox21-hse-bla": "stress_response",
    "tox21-nrf2-are-luc": "stress_response",
    
    # DNA damage assays
    "tox21-p53-bla": "dna_damage",
    
    # Cell cycle assays
    "tox21-pcna": "cell_cycle",
    
    # Apoptosis assays
    "tox21-casp3": "apoptosis",
    
    # Mitochondrial function assays
    "tox21-mitotox": "mitochondrial_function",
    "tox21-mmp": "mitochondrial_function",
    
    # Default mapping for unclassified assays
    "default": "other"
}

def classify_assay(assay_id: str, assay_name: str, assay_target: str) -> str:
    """
    Classify a Tox21 assay into a toxicological endpoint.
    
    Args:
        assay_id: The assay ID
        assay_name: The assay name
        assay_target: The assay target
        
    Returns:
        The toxicological endpoint for the assay
    """
    # Convert to lowercase for case-insensitive matching
    assay_id_lower = assay_id.lower()
    assay_name_lower = assay_name.lower()
    assay_target_lower = assay_target.lower() if assay_target else ""
    
    # Check for exact matches in the mapping
    if assay_id_lower in TOX21_ASSAY_ENDPOINT_MAPPING:
        return TOX21_ASSAY_ENDPOINT_MAPPING[assay_id_lower]
    
    # Check for partial matches in assay ID
    for pattern, endpoint in TOX21_ASSAY_ENDPOINT_MAPPING.items():
        if pattern != "default" and pattern in assay_id_lower:
            return endpoint
    
    # Check for keywords in assay name and target
    keywords_to_endpoints = {
        # Nuclear receptor keywords
        "estrogen": "nuclear_receptor",
        "androgen": "nuclear_receptor",
        "glucocorticoid": "nuclear_receptor",
        "thyroid": "nuclear_receptor",
        "ppar": "nuclear_receptor",
        "ahr": "nuclear_receptor",
        
        # Stress response keywords
        "oxidative stress": "stress_response",
        "heat shock": "stress_response",
        "nrf2": "stress_response",
        "antioxidant": "stress_response",
        
        # DNA damage keywords
        "dna damage": "dna_damage",
        "genotoxic": "dna_damage",
        "p53": "dna_damage",
        
        # Cell cycle keywords
        "cell cycle": "cell_cycle",
        "proliferation": "cell_cycle",
        
        # Apoptosis keywords
        "apoptosis": "apoptosis",
        "caspase": "apoptosis",
        "cell death": "apoptosis",
        
        # Mitochondrial function keywords
        "mitochondria": "mitochondrial_function",
        "membrane potential": "mitochondrial_function",
        
        # Enzyme inhibition keywords
        "enzyme": "enzyme_inhibition",
        "cyp": "enzyme_inhibition",
        "cytochrome": "enzyme_inhibition",
        
        # Cell signaling keywords
        "signaling": "cell_signaling",
        "kinase": "cell_signaling",
        "phosphatase": "cell_signaling",
        
        # Immunotoxicity keywords
        "immune": "immunotoxicity",
        "cytokine": "immunotoxicity",
        "inflammation": "immunotoxicity",
        
        # Neurotoxicity keywords
        "neuro": "neurotoxicity",
        "brain": "neurotoxicity",
        "neuron": "neurotoxicity",
        
        # Developmental toxicity keywords
        "development": "developmental_toxicity",
        "reproductive": "developmental_toxicity",
        "teratogen": "developmental_toxicity"
    }
    
    for keyword, endpoint in keywords_to_endpoints.items():
        if keyword in assay_name_lower or keyword in assay_target_lower:
            return endpoint
    
    # Default to "other" if no match found
    return TOX21_ASSAY_ENDPOINT_MAPPING["default"]

def get_endpoint_weight(endpoint: str) -> float:
    """
    Get the weight for a toxicological endpoint.
    
    Args:
        endpoint: The toxicological endpoint
        
    Returns:
        The weight for the endpoint
    """
    return ENDPOINT_WEIGHTS.get(endpoint, 1.0)

def get_all_endpoints() -> List[str]:
    """
    Get all defined toxicological endpoints.
    
    Returns:
        List of all endpoint names
    """
    return list(TOXICOLOGICAL_ENDPOINTS.keys())

def get_endpoint_description(endpoint: str) -> str:
    """
    Get the description for a toxicological endpoint.
    
    Args:
        endpoint: The toxicological endpoint
        
    Returns:
        The description for the endpoint
    """
    return TOXICOLOGICAL_ENDPOINTS.get(endpoint, "Unknown endpoint")

def classify_assays_in_database(supabase: Client, source_id: str) -> int:
    """
    Classify all Tox21 assays in the database and update their endpoint classifications.
    
    Args:
        supabase: Supabase client
        source_id: The Tox21 source ID
        
    Returns:
        Number of assays classified
    """
    try:
        # Get all Tox21 assays from the database
        response = supabase.table("toxicity_assay").select("id, assay_id, assay_name, assay_target").eq("source_id", source_id).execute()
        
        if not response.data:
            logger.warning(f"No assays found for source ID {source_id}")
            return 0
        
        # Classify each assay and prepare update data
        assays_classified = 0
        batch_size = 100
        
        for i in range(0, len(response.data), batch_size):
            batch = response.data[i:i+batch_size]
            
            for assay in batch:
                # Classify the assay
                endpoint = classify_assay(
                    assay.get("assay_id", ""),
                    assay.get("assay_name", ""),
                    assay.get("assay_target", "")
                )
                
                # Update the assay with the endpoint classification
                try:
                    update_response = supabase.table("toxicity_assay").update({
                        "toxicological_endpoint": endpoint
                    }).eq("id", assay["id"]).execute()
                    
                    if update_response.data:
                        assays_classified += 1
                except Exception as e:
                    logger.error(f"Error updating assay {assay['id']}: {str(e)}")
        
        logger.info(f"Classified {assays_classified} Tox21 assays")
        return assays_classified
    
    except Exception as e:
        logger.error(f"Error classifying assays: {str(e)}")
        raise

def get_assays_by_endpoint(supabase: Client, source_id: str) -> Dict[str, List[Dict]]:
    """
    Group Tox21 assays by toxicological endpoint.
    
    Args:
        supabase: Supabase client
        source_id: The Tox21 source ID
        
    Returns:
        Dictionary mapping endpoints to lists of assays
    """
    try:
        # Get all Tox21 assays with their endpoint classifications
        response = supabase.table("toxicity_assay").select("id, assay_id, assay_name, toxicological_endpoint").eq("source_id", source_id).execute()
        
        if not response.data:
            logger.warning(f"No assays found for source ID {source_id}")
            return {}
        
        # Group assays by endpoint
        assays_by_endpoint = {}
        
        for assay in response.data:
            endpoint = assay.get("toxicological_endpoint", "other")
            
            if endpoint not in assays_by_endpoint:
                assays_by_endpoint[endpoint] = []
            
            assays_by_endpoint[endpoint].append(assay)
        
        return assays_by_endpoint
    
    except Exception as e:
        logger.error(f"Error getting assays by endpoint: {str(e)}")
        return {}