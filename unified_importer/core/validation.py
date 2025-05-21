"""
Validation utilities for molecular data imports.

This module provides validation functions for molecular structures,
properties, and other data types used in the import process.
"""

import re
import logging
from typing import Dict, Any, Optional, List, Tuple, Set, Union


class ValidationError(Exception):
    """Exception raised for data validation errors."""
    pass


class MoleculeValidator:
    """
    Validator for molecular data.
    
    This class provides methods to validate molecular structures,
    properties, and identifiers before database insertion.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize validator.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        
        # Compile regex patterns once
        self._inchikey_pattern = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')
        self._chembl_id_pattern = re.compile(r'^CHEMBL\d+$')
        self._pubchem_cid_pattern = re.compile(r'^\d+$')
    
    def validate_molecule(
        self,
        molecule_data: Dict[str, Any],
        required_fields: Optional[List[str]] = None
    ) -> Tuple[bool, Optional[str]]:
        """
        Validate a molecule record.
        
        Args:
            molecule_data: Dictionary containing molecule data
            required_fields: List of fields that must be present
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if required_fields is None:
            required_fields = ['name', 'smiles']
            
        # Check required fields
        for field in required_fields:
            if field not in molecule_data or not molecule_data[field]:
                return False, f"Missing required field: {field}"
        
        # Validate SMILES if present
        if 'smiles' in molecule_data and molecule_data['smiles']:
            is_valid, error = self.validate_smiles(molecule_data['smiles'])
            if not is_valid:
                return False, f"Invalid SMILES: {error}"
        
        # Validate InChI Key if present
        if 'inchikey' in molecule_data and molecule_data['inchikey']:
            is_valid, error = self.validate_inchikey(molecule_data['inchikey'])
            if not is_valid:
                return False, f"Invalid InChI Key: {error}"
        
        # Validate identifiers
        if 'pubchem_cid' in molecule_data and molecule_data['pubchem_cid']:
            is_valid, error = self.validate_pubchem_cid(molecule_data['pubchem_cid'])
            if not is_valid:
                return False, f"Invalid PubChem CID: {error}"
                
        if 'chembl_id' in molecule_data and molecule_data['chembl_id']:
            is_valid, error = self.validate_chembl_id(molecule_data['chembl_id'])
            if not is_valid:
                return False, f"Invalid ChEMBL ID: {error}"
        
        return True, None
    
    def validate_smiles(self, smiles: str) -> Tuple[bool, Optional[str]]:
        """
        Validate a SMILES string.
        
        Args:
            smiles: SMILES string to validate
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not smiles:
            return False, "Empty SMILES string"
            
        # Basic structural validation - this is minimal
        # For real validation, RDKit would be used to parse the SMILES
        if len(smiles) < 2:
            return False, "SMILES string too short"
            
        # Check for unbalanced parentheses
        if smiles.count('(') != smiles.count(')'):
            return False, "Unbalanced parentheses in SMILES"
            
        # Check for unbalanced brackets
        if smiles.count('[') != smiles.count(']'):
            return False, "Unbalanced brackets in SMILES"
        
        return True, None
    
    def validate_inchikey(self, inchikey: str) -> Tuple[bool, Optional[str]]:
        """
        Validate an InChI Key.
        
        Args:
            inchikey: InChI Key to validate
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not inchikey:
            return False, "Empty InChI Key"
            
        # Check format with regex
        if not self._inchikey_pattern.match(inchikey):
            return False, "InChI Key does not match expected format"
        
        return True, None
    
    def validate_chembl_id(self, chembl_id: Union[str, int]) -> Tuple[bool, Optional[str]]:
        """
        Validate a ChEMBL ID.
        
        Args:
            chembl_id: ChEMBL ID to validate
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not chembl_id:
            return False, "Empty ChEMBL ID"
            
        # Convert to string if needed
        chembl_id_str = str(chembl_id)
        
        # Check format with regex
        if not self._chembl_id_pattern.match(chembl_id_str):
            return False, "ChEMBL ID does not match expected format"
        
        return True, None
    
    def validate_pubchem_cid(self, pubchem_cid: Union[str, int]) -> Tuple[bool, Optional[str]]:
        """
        Validate a PubChem CID.
        
        Args:
            pubchem_cid: PubChem CID to validate
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not pubchem_cid:
            return False, "Empty PubChem CID"
            
        # Convert to string if needed
        cid_str = str(pubchem_cid)
        
        # Check format with regex
        if not self._pubchem_cid_pattern.match(cid_str):
            return False, "PubChem CID does not match expected format"
            
        # Additional numeric validation
        try:
            cid_int = int(cid_str)
            if cid_int <= 0:
                return False, "PubChem CID must be a positive integer"
        except ValueError:
            return False, "PubChem CID must be a valid integer"
        
        return True, None
    
    def validate_property(
        self,
        property_data: Dict[str, Any],
        required_fields: Optional[List[str]] = None
    ) -> Tuple[bool, Optional[str]]:
        """
        Validate a molecular property record.
        
        Args:
            property_data: Dictionary containing property data
            required_fields: List of fields that must be present
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if required_fields is None:
            required_fields = ['property_name', 'property_type']
            
        # Check required fields
        for field in required_fields:
            if field not in property_data or not property_data[field]:
                return False, f"Missing required field: {field}"
        
        # Check that at least one value field is set
        value_fields = ['numeric_value', 'text_value', 'boolean_value']
        if not any(field in property_data and property_data[field] is not None for field in value_fields):
            return False, "At least one value field must be set"
        
        # Validate value types
        if 'numeric_value' in property_data and property_data['numeric_value'] is not None:
            try:
                float(property_data['numeric_value'])
            except (ValueError, TypeError):
                return False, "numeric_value must be a valid number"
                
        if 'boolean_value' in property_data and property_data['boolean_value'] is not None:
            if not isinstance(property_data['boolean_value'], bool):
                return False, "boolean_value must be a boolean"
        
        return True, None
    
    def validate_synonyms(
        self,
        synonyms: List[str]
    ) -> Tuple[bool, Optional[str], List[str]]:
        """
        Validate and clean a list of synonyms.
        
        Args:
            synonyms: List of synonym strings
            
        Returns:
            Tuple of (is_valid, error_message, cleaned_synonyms)
        """
        if not synonyms:
            return True, None, []
            
        # Clean and filter synonyms
        cleaned = []
        seen = set()
        
        for synonym in synonyms:
            if not synonym or not synonym.strip():
                continue
                
            # Clean whitespace
            clean_synonym = synonym.strip()
            
            # Skip duplicates
            if clean_synonym.lower() in seen:
                continue
                
            seen.add(clean_synonym.lower())
            cleaned.append(clean_synonym)
        
        return True, None, cleaned
    
    def filter_valid_molecules(
        self,
        molecules: List[Dict[str, Any]],
        required_fields: Optional[List[str]] = None
    ) -> Tuple[List[Dict[str, Any]], List[Tuple[Dict[str, Any], str]]]:
        """
        Filter a list of molecules, keeping only valid ones.
        
        Args:
            molecules: List of molecule dictionaries
            required_fields: List of fields that must be present
            
        Returns:
            Tuple of (valid_molecules, invalid_molecules_with_reasons)
        """
        valid_molecules = []
        invalid_molecules = []
        
        for molecule in molecules:
            is_valid, error = self.validate_molecule(molecule, required_fields)
            
            if is_valid:
                valid_molecules.append(molecule)
            else:
                invalid_molecules.append((molecule, error or "Unknown validation error"))
        
        return valid_molecules, invalid_molecules