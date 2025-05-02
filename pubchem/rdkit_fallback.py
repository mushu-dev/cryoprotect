"""
RDKit property calculation fallback system for PubChem API.

This module provides a robust fallback system that uses RDKit to calculate
molecular properties when PubChem API fails or returns incomplete data.
"""

import logging
import time
from typing import Dict, Any, Optional, List, Union, Tuple

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. Property calculation fallback will not work.")

# Local imports
from .cache import store_compound, get_compound

logger = logging.getLogger(__name__)

# Property mapping between PubChem and RDKit
PROPERTY_MAPPING = {
    "Molecular Formula": "MolecularFormula",
    "Molecular Weight": "MolecularWeight",
    "LogP": "XLogP",
    "TPSA": "TPSA",
    "H-Bond Donors": "HBondDonorCount",
    "H-Bond Acceptors": "HBondAcceptorCount",
    "SMILES": "IsomericSMILES",
    "InChI": "InChI",
    "InChIKey": "InChIKey",
    "IUPACName": None,  # Not available in RDKit
    "Title": None,      # Not available in RDKit
}

# Additional RDKit properties not in PubChem
ADDITIONAL_PROPERTIES = {
    "Rotatable Bonds": "NumRotatableBonds",
    "Ring Count": "RingCount",
    "Aromatic Ring Count": "NumAromaticRings",
    "Fraction CSP3": "FractionCSP3",
    "Heavy Atom Count": "HeavyAtomCount",
}

def calculate_rdkit_properties(
    molecule_data: str,
    input_format: str = 'smiles',
    include_additional_properties: bool = False
) -> Dict[str, Any]:
    """
    Calculate molecular properties using RDKit.
    
    Args:
        molecule_data: Molecular data (SMILES, InChI, etc.)
        input_format: Format of molecule_data ('smiles', 'inchi', 'inchikey')
        include_additional_properties: Whether to include additional RDKit properties
        
    Returns:
        Dictionary of calculated properties
    """
    if not RDKIT_AVAILABLE:
        return {"Error": "RDKit not available"}
    
    # Initialize result dictionary
    result = {}
    
    try:
        # Create RDKit molecule object based on input format
        mol = None
        if input_format.lower() == 'smiles':
            mol = Chem.MolFromSmiles(molecule_data)
        elif input_format.lower() == 'inchi':
            mol = Chem.MolFromInchi(molecule_data)
        elif input_format.lower() == 'inchikey':
            logger.warning("Cannot create molecule directly from InChIKey. Need SMILES or InChI.")
            return {"Error": "Cannot create molecule directly from InChIKey"}
        else:
            logger.warning(f"Unsupported input format: {input_format}")
            return {"Error": f"Unsupported input format: {input_format}"}
        
        # Check if molecule was successfully created
        if mol is None:
            logger.warning(f"Failed to create molecule from {input_format}: {molecule_data}")
            return {"Error": f"Failed to create molecule from {input_format}"}
        
        # Calculate standard properties that match PubChem
        result["Molecular Formula"] = rdMolDescriptors.CalcMolFormula(mol)
        result["Molecular Weight"] = Descriptors.MolWt(mol)
        result["LogP"] = Descriptors.MolLogP(mol)
        result["TPSA"] = Descriptors.TPSA(mol)
        result["H-Bond Donors"] = Lipinski.NumHDonors(mol)
        result["H-Bond Acceptors"] = Lipinski.NumHAcceptors(mol)
        result["SMILES"] = Chem.MolToSmiles(mol, isomericSmiles=True)
        
        # Try to calculate InChI and InChIKey
        try:
            result["InChI"] = Chem.MolToInchi(mol)
            result["InChIKey"] = Chem.MolToInchiKey(mol)
        except Exception as e:
            logger.warning(f"Failed to calculate InChI/InChIKey: {str(e)}")
            result["InChI"] = None
            result["InChIKey"] = None
        
        # These properties are not available in RDKit
        result["IUPACName"] = None
        result["Title"] = None
        
        # Calculate additional RDKit properties if requested
        if include_additional_properties:
            result["Rotatable Bonds"] = Descriptors.NumRotatableBonds(mol)
            result["Ring Count"] = Descriptors.RingCount(mol)
            result["Aromatic Ring Count"] = Lipinski.NumAromaticRings(mol)
            result["Fraction CSP3"] = Descriptors.FractionCSP3(mol)
            result["Heavy Atom Count"] = Descriptors.HeavyAtomCount(mol)
        
        # Add metadata
        result["Source"] = "rdkit"
        result["Calculation Time"] = time.time()
        
        return result
        
    except Exception as e:
        logger.error(f"Error calculating properties with RDKit: {str(e)}")
        return {
            "Error": f"Error calculating properties with RDKit: {str(e)}",
            "Source": "rdkit",
            "Calculation Time": time.time()
        }

# Default prioritization order for properties
DEFAULT_PRIORITY = {
    "Molecular Formula": ["pubchem", "rdkit"],
    "Molecular Weight": ["pubchem", "rdkit"],
    "LogP": ["pubchem", "rdkit"],
    "TPSA": ["pubchem", "rdkit"],
    "H-Bond Donors": ["pubchem", "rdkit"],
    "H-Bond Acceptors": ["pubchem", "rdkit"],
    "SMILES": ["pubchem", "rdkit"],
    "InChI": ["pubchem", "rdkit"],
    "InChIKey": ["pubchem", "rdkit"],
    "IUPACName": ["pubchem"],
    "Title": ["pubchem"],
    # Additional properties
    "Rotatable Bonds": ["rdkit"],
    "Ring Count": ["rdkit"],
    "Aromatic Ring Count": ["rdkit"],
    "Fraction CSP3": ["rdkit"],
    "Heavy Atom Count": ["rdkit"]
}

def merge_compound_data(
    data_sources: Dict[str, Dict[str, Any]],
    priority_config: Dict[str, List[str]] = None
) -> Dict[str, Any]:
    """
    Merge compound data from multiple sources based on priority rules.
    
    Args:
        data_sources: Dictionary of data sources (key: source name, value: compound data)
        priority_config: Configuration for property source prioritization
        
    Returns:
        Merged compound data
    """
    if not data_sources:
        logger.warning("No data sources provided for merging")
        return {}
    
    # Use default priority config if none provided
    if priority_config is None:
        priority_config = DEFAULT_PRIORITY
    
    # Initialize merged data with metadata
    merged_data = {
        "Source": "merged",
        "Merge Time": time.time(),
        "Merge Sources": list(data_sources.keys()),
        "Merge Details": {}
    }
    
    # Get all available properties from all sources
    all_properties = set()
    for source_data in data_sources.values():
        all_properties.update(key for key in source_data.keys()
                             if not key.startswith("Source") and
                                not key.startswith("Calculation") and
                                not key.startswith("Error") and
                                not key.startswith("Standardization") and
                                not key.startswith("Merge"))
    
    # Process each property
    for prop in all_properties:
        # Get priority list for this property
        priority_list = priority_config.get(prop, list(data_sources.keys()))
        
        # Try to get property value from sources in priority order
        value = None
        selected_source = None
        
        for source_name in priority_list:
            if source_name in data_sources and prop in data_sources[source_name]:
                source_value = data_sources[source_name][prop]
                
                # Skip None or empty values
                if source_value is None or (isinstance(source_value, str) and not source_value.strip()):
                    continue
                
                # Use this value
                value = source_value
                selected_source = source_name
                break
        
        # Add property to merged data if a value was found
        if value is not None:
            merged_data[prop] = value
            merged_data["Merge Details"][prop] = selected_source
    
    # Check for missing critical properties and try to fill them
    critical_properties = [
        "Molecular Formula", "Molecular Weight", "SMILES", "InChI", "InChIKey"
    ]
    
    for prop in critical_properties:
        if prop not in merged_data:
            # Try to find any non-None value for this property
            for source_name, source_data in data_sources.items():
                if prop in source_data and source_data[prop] is not None:
                    merged_data[prop] = source_data[prop]
                    merged_data["Merge Details"][prop] = f"{source_name} (fallback)"
                    break
    
    # Add validation information
    merged_data["Validation"] = {
        "Missing Properties": [prop for prop in critical_properties if prop not in merged_data],
        "Complete": all(prop in merged_data for prop in critical_properties)
    }
    
    return merged_data

def get_compound_properties(
    cid: Union[str, int],
    client: Optional[Any] = None,
    use_cache: bool = True,
    fallback_to_rdkit: bool = True,
    rdkit_fallback_identifier: str = None,
    include_additional_properties: bool = False,
    priority_config: Dict[str, List[str]] = None
) -> Dict[str, Any]:
    """
    Get compound properties with RDKit fallback.
    
    Args:
        cid: PubChem Compound ID
        client: PubChem client instance (optional)
        use_cache: Whether to use cached data
        fallback_to_rdkit: Whether to use RDKit as fallback
        rdkit_fallback_identifier: Identifier to use for RDKit (SMILES, InChI)
        include_additional_properties: Whether to include additional RDKit properties
        priority_config: Configuration for property source prioritization
        
    Returns:
        Compound properties
    """
    # Convert CID to string for consistency
    cid_str = str(cid)
    
    # Create cache key for this request
    cache_key = f"rdkit_fallback_{cid_str}"
    if include_additional_properties:
        cache_key += "_extended"
    
    # Check cache if enabled
    if use_cache:
        cached_data = get_compound(int(cid_str))
        if cached_data and cached_data.get("Source") == "merged":
            logger.debug(f"Cache hit for compound properties with RDKit fallback for CID {cid_str}")
            return cached_data
    
    # Initialize data sources dictionary
    data_sources = {}
    
    # Try to get data from PubChem API if client is provided
    pubchem_data = None
    if client is not None:
        try:
            pubchem_data = client.get_molecule_properties(
                cid=cid_str,
                use_cache=use_cache,
                fallback_to_cache=True
            )
            
            # Check if PubChem data has error
            if "Error" not in pubchem_data:
                data_sources["pubchem"] = pubchem_data
                
                # If we don't need to use RDKit fallback, return PubChem data directly
                if not fallback_to_rdkit:
                    return pubchem_data
                
                # Extract SMILES or InChI for RDKit if not provided
                if rdkit_fallback_identifier is None:
                    if "SMILES" in pubchem_data and pubchem_data["SMILES"]:
                        rdkit_fallback_identifier = pubchem_data["SMILES"]
                    elif "InChI" in pubchem_data and pubchem_data["InChI"]:
                        rdkit_fallback_identifier = pubchem_data["InChI"]
            else:
                logger.warning(f"Error in PubChem data for CID {cid_str}: {pubchem_data.get('Error')}")
        except Exception as e:
            logger.warning(f"Error fetching PubChem data for CID {cid_str}: {str(e)}")
    
    # Calculate properties with RDKit if fallback is enabled and we have an identifier
    rdkit_data = None
    if fallback_to_rdkit and rdkit_fallback_identifier:
        try:
            # Determine input format
            input_format = 'smiles'
            if rdkit_fallback_identifier.startswith('InChI='):
                input_format = 'inchi'
            
            # Calculate properties with RDKit
            rdkit_data = calculate_rdkit_properties(
                molecule_data=rdkit_fallback_identifier,
                input_format=input_format,
                include_additional_properties=include_additional_properties
            )
            
            # Check if RDKit calculation was successful
            if "Error" not in rdkit_data:
                data_sources["rdkit"] = rdkit_data
            else:
                logger.warning(f"Error calculating RDKit properties for CID {cid_str}: {rdkit_data.get('Error')}")
        except Exception as e:
            logger.warning(f"Error calculating RDKit properties for CID {cid_str}: {str(e)}")
    
    # If we have no data sources, return error
    if not data_sources:
        error_msg = "No data available from PubChem or RDKit"
        logger.error(f"{error_msg} for CID {cid_str}")
        return {
            "CID": cid_str,
            "Error": error_msg,
            "Source": "none"
        }
    
    # Merge data from all sources
    merged_data = merge_compound_data(
        data_sources=data_sources,
        priority_config=priority_config
    )
    
    # Add CID to merged data
    merged_data["CID"] = cid_str
    
    # Add PubChem link
    merged_data["PubChem Link"] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid_str}"
    
    # Cache the merged result
    if use_cache:
        try:
            store_compound(int(cid_str), merged_data, source="merged")
            logger.debug(f"Cached merged compound data for CID {cid_str}")
        except Exception as e:
            logger.warning(f"Error caching merged compound data for CID {cid_str}: {str(e)}")
    
    return merged_data

def get_compound_properties_batch(
    cids: List[Union[str, int]],
    client: Optional[Any] = None,
    use_cache: bool = True,
    fallback_to_rdkit: bool = True,
    include_additional_properties: bool = False,
    priority_config: Dict[str, List[str]] = None,
    max_workers: int = 4
) -> Dict[str, Dict[str, Any]]:
    """
    Get compound properties with RDKit fallback for multiple compounds.
    
    Args:
        cids: List of PubChem Compound IDs
        client: PubChem client instance (optional)
        use_cache: Whether to use cached data
        fallback_to_rdkit: Whether to use RDKit as fallback
        include_additional_properties: Whether to include additional RDKit properties
        priority_config: Configuration for property source prioritization
        max_workers: Maximum number of worker threads
        
    Returns:
        Dictionary mapping CIDs to compound properties
    """
    import concurrent.futures
    
    # Convert all CIDs to strings
    cid_strs = [str(cid) for cid in cids]
    
    # Initialize results dictionary
    results = {}
    
    # Define worker function
    def process_cid(cid_str):
        try:
            return cid_str, get_compound_properties(
                cid=cid_str,
                client=client,
                use_cache=use_cache,
                fallback_to_rdkit=fallback_to_rdkit,
                include_additional_properties=include_additional_properties,
                priority_config=priority_config
            )
        except Exception as e:
            logger.error(f"Error processing CID {cid_str}: {str(e)}")
            return cid_str, {
                "CID": cid_str,
                "Error": str(e),
                "Source": "error"
            }
    
    # Process CIDs in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        for cid_str, result in executor.map(process_cid, cid_strs):
            results[cid_str] = result
    
    return results