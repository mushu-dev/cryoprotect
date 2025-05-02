"""
Reference Compounds Module

This module defines a list of standard reference compounds used in the CryoProtect system.
These compounds serve as a baseline for comparison and include common cryoprotectants.
"""

from typing import List

def get_reference_compound_ids() -> List[str]:
    """
    Get list of reference cryoprotectant compound IDs.
    
    Returns:
        List of ChEMBL IDs for reference compounds
    """
    return [
        "CHEMBL388978",     # Glycerol
        "CHEMBL1098659",    # DMSO
        "CHEMBL66195",      # beta-Alanine
        "CHEMBL500033",     # tert-Butanol
        "CHEMBL1487",       # Urea
        "CHEMBL6196",       # Ethylene glycol
        "CHEMBL967",        # Propylene glycol
        "CHEMBL262548",     # Trehalose
        "CHEMBL6752"        # Glycine
    ]

# Extended reference compounds for more comprehensive coverage
def get_extended_reference_compound_ids():
    """
    Returns an extended list of ChEMBL IDs for reference compounds.
    
    This includes the standard reference compounds plus additional
    compounds of interest for cryoprotection research.
    
    Returns:
        list: A list of ChEMBL IDs as strings
    """
    # Start with standard references
    extended_ids = get_reference_compound_ids()
    
    # Add additional compounds
    additional_ids = [
        "CHEMBL71",     # Sucrose
        "CHEMBL1487",   # Mannitol
        "CHEMBL1229",   # Sorbitol
        "CHEMBL1200485",# Propylene glycol
        "CHEMBL1231",   # Urea
        "CHEMBL1308",   # Glycine
        "CHEMBL276624", # Betaine
        "CHEMBL1641",   # Proline
    ]
    
    extended_ids.extend(additional_ids)
    return extended_ids