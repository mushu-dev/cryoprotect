"""
Data standardization module for PubChem and RDKit property data.

This module provides functions to standardize compound data from different sources,
ensuring consistent units, formats, and naming conventions across the system.
"""

import logging
import re
from typing import Dict, Any, List, Optional, Union, Tuple

logger = logging.getLogger(__name__)

# Define expected property types and formats
PROPERTY_TYPES = {
    "Molecular Formula": str,
    "Molecular Weight": float,
    "LogP": float,
    "TPSA": float,
    "H-Bond Donors": int,
    "H-Bond Acceptors": int,
    "SMILES": str,
    "InChI": str,
    "InChIKey": str,
    "IUPACName": str,
    "Title": str,
    # Additional RDKit properties
    "Rotatable Bonds": int,
    "Ring Count": int,
    "Aromatic Ring Count": int,
    "Fraction CSP3": float,
    "Heavy Atom Count": int
}

# Define expected value ranges for numerical properties
PROPERTY_RANGES = {
    "Molecular Weight": (0, 5000),  # Da
    "LogP": (-20, 20),  # Typical range for drug-like compounds
    "TPSA": (0, 500),   # Å²
    "H-Bond Donors": (0, 50),
    "H-Bond Acceptors": (0, 50),
    "Rotatable Bonds": (0, 100),
    "Ring Count": (0, 50),
    "Aromatic Ring Count": (0, 30),
    "Fraction CSP3": (0, 1),
    "Heavy Atom Count": (0, 500)
}

# Define default values for properties when missing
DEFAULT_VALUES = {
    "Molecular Formula": "",
    "Molecular Weight": 0.0,
    "LogP": 0.0,
    "TPSA": 0.0,
    "H-Bond Donors": 0,
    "H-Bond Acceptors": 0,
    "SMILES": "",
    "InChI": "",
    "InChIKey": "",
    "IUPACName": "",
    "Title": "",
    "Rotatable Bonds": 0,
    "Ring Count": 0,
    "Aromatic Ring Count": 0,
    "Fraction CSP3": 0.0,
    "Heavy Atom Count": 0
}

# Regular expression patterns for validation
PATTERNS = {
    "Molecular Formula": r'^([A-Z][a-z]?\d*)+$',
    "InChI": r'^InChI=1S?/.*$',
    "InChIKey": r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'
}


def standardize_compound_data(
    compound_data: Dict[str, Any],
    source: str = 'pubchem'
) -> Dict[str, Any]:
    """
    Standardize compound data to ensure consistent format and values.
    
    Args:
        compound_data: Raw compound data
        source: Source of the data ('pubchem', 'rdkit', etc.)
        
    Returns:
        Standardized compound data
    """
    if not compound_data:
        logger.warning(f"Empty compound data from {source}")
        return {}
    
    # Initialize standardized data with source information
    standardized_data = {
        "Source": source,
        "Standardization": {
            "Warnings": [],
            "Modified": []
        }
    }
    
    # Map property names based on source
    mapped_data = _map_property_names(compound_data, source)
    
    # Process each property
    for prop_name, expected_type in PROPERTY_TYPES.items():
        # Skip properties not relevant to this source
        if prop_name not in mapped_data and prop_name not in DEFAULT_VALUES:
            continue
            
        # Get raw value or default
        raw_value = mapped_data.get(prop_name)
        
        # Standardize the property
        try:
            std_value = _standardize_property(prop_name, raw_value, expected_type)
            standardized_data[prop_name] = std_value
            
            # Check if value was modified
            if raw_value is not None and std_value != raw_value:
                standardized_data["Standardization"]["Modified"].append(prop_name)
            
            # Special handling for None values that were converted to defaults
            if raw_value is None and prop_name in DEFAULT_VALUES:
                standardized_data["Standardization"]["Modified"].append(prop_name)
                
        except ValueError as e:
            # Log warning and use default value
            warning = f"Invalid {prop_name}: {str(e)}"
            logger.warning(warning)
            standardized_data["Standardization"]["Warnings"].append(warning)
            standardized_data[prop_name] = DEFAULT_VALUES.get(prop_name)
    
    return standardized_data


def _map_property_names(
    compound_data: Dict[str, Any],
    source: str
) -> Dict[str, Any]:
    """
    Map property names from source-specific names to standard names.
    
    Args:
        compound_data: Raw compound data
        source: Source of the data ('pubchem', 'rdkit', etc.)
        
    Returns:
        Compound data with standardized property names
    """
    # Define mapping from source-specific names to standard names
    pubchem_to_standard = {
        "MolecularFormula": "Molecular Formula",
        "MolecularWeight": "Molecular Weight",
        "XLogP": "LogP",
        "TPSA": "TPSA",
        "HBondDonorCount": "H-Bond Donors",
        "HBondAcceptorCount": "H-Bond Acceptors",
        "IsomericSMILES": "SMILES",
        "InChI": "InChI",
        "InChIKey": "InChIKey",
        "IUPACName": "IUPACName",
        "Title": "Title"
    }
    
    # RDKit already uses standard names as defined in rdkit_fallback.py
    
    mapped_data = {}
    
    if source.lower() == 'pubchem':
        # Map PubChem property names to standard names
        for pubchem_name, std_name in pubchem_to_standard.items():
            if pubchem_name in compound_data:
                mapped_data[std_name] = compound_data[pubchem_name]
    else:
        # For other sources (including rdkit), copy properties directly
        for prop_name in PROPERTY_TYPES.keys():
            if prop_name in compound_data:
                mapped_data[prop_name] = compound_data[prop_name]
    
    return mapped_data


def _standardize_property(
    prop_name: str,
    value: Any,
    expected_type: type
) -> Any:
    """
    Standardize a single property value.
    
    Args:
        prop_name: Name of the property
        value: Raw property value
        expected_type: Expected type of the property
        
    Returns:
        Standardized property value
        
    Raises:
        ValueError: If the value cannot be standardized
    """
    # Handle None/empty values
    if value is None:
        return DEFAULT_VALUES.get(prop_name)
    
    # Convert to expected type
    try:
        if expected_type == int:
            # Convert to integer
            if isinstance(value, str) and value.strip() == '':
                return DEFAULT_VALUES.get(prop_name)
            std_value = int(float(value))  # Handle string or float inputs
        elif expected_type == float:
            # Convert to float and round to 2 decimal places
            if isinstance(value, str) and value.strip() == '':
                return DEFAULT_VALUES.get(prop_name)
            std_value = round(float(value), 2)
        elif expected_type == str:
            # Convert to string
            std_value = str(value).strip()
            if std_value == '':
                return DEFAULT_VALUES.get(prop_name)
        else:
            # Unsupported type
            raise ValueError(f"Unsupported type {expected_type} for {prop_name}")
    except (ValueError, TypeError):
        raise ValueError(f"Cannot convert {value} to {expected_type.__name__}")
    
    # Validate value range for numerical properties
    if prop_name in PROPERTY_RANGES:
        min_val, max_val = PROPERTY_RANGES[prop_name]
        if not min_val <= std_value <= max_val:
            raise ValueError(f"Value {std_value} outside expected range [{min_val}, {max_val}]")
    
    # Validate format for string properties with patterns
    if prop_name in PATTERNS and isinstance(std_value, str) and std_value:
        pattern = PATTERNS[prop_name]
        if not re.match(pattern, std_value):
            # For molecular formula, try to fix common formatting issues
            if prop_name == "Molecular Formula":
                std_value = _fix_molecular_formula(std_value)
                if not re.match(pattern, std_value):
                    raise ValueError(f"Invalid format: {std_value} does not match pattern {pattern}")
            else:
                raise ValueError(f"Invalid format: {std_value} does not match pattern {pattern}")
    
    # Special handling for specific properties
    if prop_name == "SMILES" and std_value:
        # Ensure SMILES is properly formatted (basic check)
        if '.' in std_value and not any(c in std_value for c in 'CNO'):
            raise ValueError(f"Invalid SMILES: {std_value}")
    
    return std_value


def _fix_molecular_formula(formula: str) -> str:
    """
    Fix common formatting issues in molecular formulas.
    
    Args:
        formula: Raw molecular formula
        
    Returns:
        Fixed molecular formula
    """
    # Remove spaces, parentheses, and brackets
    formula = re.sub(r'[\s\(\)\[\]]', '', formula)
    
    # Handle lowercase elements by capitalizing first letter of each element
    # First, identify element boundaries
    # This regex looks for patterns like 'c9h8o4' and converts to 'C9H8O4'
    result = ''
    i = 0
    while i < len(formula):
        # If we find a letter
        if formula[i].isalpha():
            # Capitalize the first letter of the element
            result += formula[i].upper()
            i += 1
            # Add any lowercase letters that follow (part of same element)
            while i < len(formula) and formula[i].isalpha() and formula[i].islower():
                result += formula[i]
                i += 1
            # Add any numbers that follow (count for this element)
            while i < len(formula) and formula[i].isdigit():
                result += formula[i]
                i += 1
        else:
            # Just add any other characters (like digits)
            result += formula[i]
            i += 1
    
    return result


def validate_compound_data(
    compound_data: Dict[str, Any],
    required_properties: List[str] = None
) -> Tuple[bool, List[str]]:
    """
    Validate compound data against expected types, ranges, and formats.
    
    Args:
        compound_data: Compound data to validate
        required_properties: List of properties that must be present
        
    Returns:
        Tuple of (is_valid, list_of_validation_errors)
    """
    if not compound_data:
        return False, ["Empty compound data"]
    
    validation_errors = []
    
    # Check required properties
    if required_properties:
        for prop in required_properties:
            if prop not in compound_data or compound_data[prop] in (None, ''):
                validation_errors.append(f"Missing required property: {prop}")
    
    # Validate property types and ranges
    for prop_name, value in compound_data.items():
        if prop_name in PROPERTY_TYPES and value is not None:
            expected_type = PROPERTY_TYPES[prop_name]
            
            # Check type
            if not isinstance(value, expected_type):
                validation_errors.append(
                    f"Invalid type for {prop_name}: expected {expected_type.__name__}, got {type(value).__name__}"
                )
                continue
            
            # Check range for numerical properties
            if prop_name in PROPERTY_RANGES:
                min_val, max_val = PROPERTY_RANGES[prop_name]
                if not min_val <= value <= max_val:
                    validation_errors.append(
                        f"Value for {prop_name} ({value}) outside expected range [{min_val}, {max_val}]"
                    )
            
            # Check format for string properties
            if prop_name in PATTERNS and isinstance(value, str) and value:
                pattern = PATTERNS[prop_name]
                if not re.match(pattern, value):
                    validation_errors.append(
                        f"Invalid format for {prop_name}: {value}"
                    )
    
    return len(validation_errors) == 0, validation_errors


def standardize_properties_batch(
    compounds_data: List[Dict[str, Any]],
    source: str = 'pubchem'
) -> List[Dict[str, Any]]:
    """
    Standardize multiple compound data entries in batch.
    
    Args:
        compounds_data: List of compound data dictionaries
        source: Source of the data ('pubchem', 'rdkit', etc.)
        
    Returns:
        List of standardized compound data dictionaries
    """
    standardized_compounds = []
    
    for compound_data in compounds_data:
        try:
            standardized_data = standardize_compound_data(compound_data, source)
            standardized_compounds.append(standardized_data)
        except Exception as e:
            logger.error(f"Error standardizing compound: {str(e)}")
            # Add error information to the compound data
            error_data = {
                "Error": f"Standardization failed: {str(e)}",
                "Source": source,
                "Raw Data": compound_data
            }
            standardized_compounds.append(error_data)
    
    return standardized_compounds