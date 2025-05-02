"""
Molecule Identifier Mapping Module for CryoProtect v2.

This module provides functionality for mapping external compound identifiers
(from Tox21, etc.) to internal molecule IDs in the CryoProtect database.
It includes functions for:
- Mapping compounds by various identifiers (InChIKey, InChI, SMILES, name, etc.)
- Calculating confidence scores for mappings
- Storing and retrieving mappings
"""

import logging
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Union
from supabase import Client

# Set up logging
logger = logging.getLogger(__name__)

def map_compound_identifiers(
    supabase: Client,
    identifiers: Dict[str, str],
    source_id: str,
    external_id: str,
    external_id_type: str
) -> Dict[str, Any]:
    """
    Map external compound identifiers to internal molecule IDs.
    
    Args:
        supabase: Supabase client
        identifiers: Dictionary of identifiers (e.g., {"inchikey": "...", "smiles": "..."})
        source_id: UUID of the toxicity data source
        external_id: External identifier value
        external_id_type: Type of external identifier (e.g., "DTXSID", "NCGC_ID")
        
    Returns:
        Dict with molecule_id, mapping_id, and created flag
    """
    # Check if mapping already exists
    try:
        response = supabase.table("molecule_identifier_mapping").select("id, molecule_id").eq("source_id", source_id).eq("external_id_type", external_id_type).eq("external_id", external_id).execute()
        
        if response.data and len(response.data) > 0:
            # Mapping exists
            return {
                "molecule_id": response.data[0]["molecule_id"],
                "mapping_id": response.data[0]["id"],
                "created": 0
            }
    except Exception as e:
        logger.error(f"Error checking existing mapping: {str(e)}")
    
    # No mapping exists, try to map using available identifiers
    molecule_id = None
    confidence = 0.0
    mapping_method = None
    
    # Try mapping by InChIKey (most reliable)
    if "inchikey" in identifiers:
        try:
            response = supabase.table("molecule").select("id").eq("inchikey", identifiers["inchikey"]).execute()
            if response.data and len(response.data) > 0:
                molecule_id = response.data[0]["id"]
                confidence = 1.0
                mapping_method = "inchikey_match"
        except Exception as e:
            logger.error(f"Error mapping by InChIKey: {str(e)}")
    
    # Try mapping by InChI if InChIKey failed
    if not molecule_id and "inchi" in identifiers:
        try:
            response = supabase.table("molecule").select("id").eq("inchi", identifiers["inchi"]).execute()
            if response.data and len(response.data) > 0:
                molecule_id = response.data[0]["id"]
                confidence = 0.95
                mapping_method = "inchi_match"
        except Exception as e:
            logger.error(f"Error mapping by InChI: {str(e)}")
    
    # Try mapping by SMILES if InChI failed
    if not molecule_id and "smiles" in identifiers:
        try:
            response = supabase.table("molecule").select("id").eq("smiles", identifiers["smiles"]).execute()
            if response.data and len(response.data) > 0:
                molecule_id = response.data[0]["id"]
                confidence = 0.9
                mapping_method = "smiles_match"
        except Exception as e:
            logger.error(f"Error mapping by SMILES: {str(e)}")
    
    # Try mapping by CAS Registry Number if available
    if not molecule_id and "casrn" in identifiers:
        try:
            # Check if there's a mapping from CAS RN to our molecules
            casrn_mapping_response = supabase.table("molecule_identifier_mapping").select("molecule_id").eq("external_id_type", "CASRN").eq("external_id", identifiers["casrn"]).execute()
            
            if casrn_mapping_response.data and len(casrn_mapping_response.data) > 0:
                molecule_id = casrn_mapping_response.data[0]["molecule_id"]
                confidence = 0.85
                mapping_method = "casrn_match"
        except Exception as e:
            logger.error(f"Error mapping by CAS RN: {str(e)}")
    
    # Try mapping by PubChem CID if available
    if not molecule_id and "pubchem_cid" in identifiers:
        try:
            # Check if there's a mapping from PubChem to our molecules
            pubchem_mapping_response = supabase.table("molecule_identifier_mapping").select("molecule_id").eq("external_id_type", "PUBCHEM_CID").eq("external_id", identifiers["pubchem_cid"]).execute()
            
            if pubchem_mapping_response.data and len(pubchem_mapping_response.data) > 0:
                molecule_id = pubchem_mapping_response.data[0]["molecule_id"]
                confidence = 0.85
                mapping_method = "pubchem_cid_match"
        except Exception as e:
            logger.error(f"Error mapping by PubChem CID: {str(e)}")
    
    # Try mapping by name as a last resort
    if not molecule_id and "name" in identifiers:
        try:
            response = supabase.table("molecule").select("id").ilike("name", f"%{identifiers['name']}%").execute()
            if response.data and len(response.data) > 0:
                molecule_id = response.data[0]["id"]
                confidence = 0.7
                mapping_method = "name_match"
        except Exception as e:
            logger.error(f"Error mapping by name: {str(e)}")
    
    # If mapping successful, create mapping record
    if molecule_id:
        try:
            mapping_data = {
                "molecule_id": molecule_id,
                "source_id": source_id,
                "external_id": external_id,
                "external_id_type": external_id_type,
                "confidence_score": confidence,
                "mapping_method": mapping_method
            }
            
            response = supabase.table("molecule_identifier_mapping").insert(mapping_data).execute()
            
            if response.data and len(response.data) > 0:
                return {
                    "molecule_id": molecule_id,
                    "mapping_id": response.data[0]["id"],
                    "created": 1
                }
            else:
                logger.error(f"Failed to create mapping for {external_id}")
                return {
                    "molecule_id": molecule_id,
                    "mapping_id": None,
                    "created": 0
                }
        except Exception as e:
            logger.error(f"Error creating mapping: {str(e)}")
            return {
                "molecule_id": molecule_id,
                "mapping_id": None,
                "created": 0
            }
    
    # No mapping found
    return {
        "molecule_id": None,
        "mapping_id": None,
        "created": 0
    }

def get_molecule_by_external_id(
    supabase: Client,
    external_id: str,
    external_id_type: str,
    source_id: Optional[str] = None
) -> Optional[str]:
    """
    Get a molecule ID by external identifier.
    
    Args:
        supabase: Supabase client
        external_id: External identifier value
        external_id_type: Type of external identifier (e.g., "DTXSID", "NCGC_ID")
        source_id: Optional UUID of the toxicity data source
        
    Returns:
        Molecule ID if found, None otherwise
    """
    try:
        query = supabase.table("molecule_identifier_mapping").select("molecule_id").eq("external_id", external_id).eq("external_id_type", external_id_type)
        
        if source_id:
            query = query.eq("source_id", source_id)
        
        response = query.execute()
        
        if response.data and len(response.data) > 0:
            return response.data[0]["molecule_id"]
        else:
            return None
    except Exception as e:
        logger.error(f"Error getting molecule by external ID: {str(e)}")
        return None

def get_external_ids_for_molecule(
    supabase: Client,
    molecule_id: str,
    source_id: Optional[str] = None
) -> List[Dict[str, Any]]:
    """
    Get all external identifiers for a molecule.
    
    Args:
        supabase: Supabase client
        molecule_id: Internal molecule ID
        source_id: Optional UUID of the toxicity data source
        
    Returns:
        List of dictionaries with external_id, external_id_type, and source_id
    """
    try:
        query = supabase.table("molecule_identifier_mapping").select("external_id, external_id_type, source_id").eq("molecule_id", molecule_id)
        
        if source_id:
            query = query.eq("source_id", source_id)
        
        response = query.execute()
        
        return response.data if response.data else []
    except Exception as e:
        logger.error(f"Error getting external IDs for molecule: {str(e)}")
        return []

def calculate_mapping_confidence(
    identifiers: Dict[str, str],
    mapping_method: str
) -> float:
    """
    Calculate a confidence score for a mapping based on the method and available identifiers.
    
    Args:
        identifiers: Dictionary of identifiers
        mapping_method: Method used for mapping
        
    Returns:
        Confidence score (0.0-1.0)
    """
    # Base confidence by mapping method
    base_confidence = {
        "inchikey_match": 1.0,
        "inchi_match": 0.95,
        "smiles_match": 0.9,
        "casrn_match": 0.85,
        "pubchem_cid_match": 0.85,
        "name_match": 0.7
    }.get(mapping_method, 0.5)
    
    # Adjust confidence based on available identifiers
    # More identifiers = higher confidence
    identifier_bonus = min(0.1, len(identifiers) * 0.02)
    
    # Cap at 1.0
    return min(1.0, base_confidence + identifier_bonus)