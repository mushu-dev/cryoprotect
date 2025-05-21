"""
Property transformation utilities.

This module provides functions for standardizing and transforming
property data across different chemical data sources.
"""

import logging
import re
from typing import Dict, List, Any, Optional, Tuple, Union, Set


class PropertyTransformer:
    """
    Transforms property data to standardized formats.

    This class provides methods for standardizing property names,
    units, and values across different data sources.
    """

    def __init__(self, logger: Optional[logging.Logger] = None, config: Optional[Dict[str, Any]] = None):
        """
        Initialize property transformer.

        Args:
            logger: Logger instance
            config: Configuration dictionary with options for transformation
        """
        self.logger = logger or logging.getLogger(__name__)
        self.config = config or {}

        # Configure transformation options
        self.standardize_units = self.config.get('standardize_units', True)
        self.add_missing_properties = self.config.get('add_missing_properties', True)
        self.preferred_units = self.config.get('preferred_units', {})

        # Initialize property name mappings
        self._init_property_mappings()

        # Initialize unit conversions
        self._init_unit_conversions()
    
    def _init_property_mappings(self) -> None:
        """Initialize mappings for standardizing property names."""
        # Map various property names to standard names
        self.property_name_map = {
            # Lipophilicity
            'alogp': 'LogP',
            'logp': 'LogP',
            'xlogp': 'LogP',
            'clogp': 'LogP',
            'partition_coefficient': 'LogP',
            'octanol_water_partition_coefficient': 'LogP',
            'logd': 'LogD',
            
            # Molecular weight
            'molecular_weight': 'MolecularWeight',
            'mol_weight': 'MolecularWeight',
            'mw': 'MolecularWeight',
            'monoisotopic_weight': 'MonoisotopicMass',
            'exact_mass': 'MonoisotopicMass',
            'average_mass': 'MolecularWeight',
            'weight': 'MolecularWeight',
            'full_mwt': 'MolecularWeight',
            
            # Topological properties
            'tpsa': 'TPSA',
            'polar_surface_area': 'TPSA',
            'topological_polar_surface_area': 'TPSA',
            
            # Hydrogen bonding
            'hba': 'HBondAcceptorCount',
            'hbd': 'HBondDonorCount',
            'h_bond_donor_count': 'HBondDonorCount',
            'h_bond_acceptor_count': 'HBondAcceptorCount',
            'hydrogen_bond_donors': 'HBondDonorCount',
            'hydrogen_bond_acceptors': 'HBondAcceptorCount',
            'num_h_acceptors': 'HBondAcceptorCount',
            'num_h_donors': 'HBondDonorCount',
            
            # Rotatable bonds
            'rtb': 'RotatableBondCount',
            'rotatable_bond_count': 'RotatableBondCount',
            'num_rotatable_bonds': 'RotatableBondCount',
            'nrot': 'RotatableBondCount',
            
            # Molecular formula
            'formula': 'MolecularFormula',
            'molecular_formula': 'MolecularFormula',
            'mol_formula': 'MolecularFormula',
            'full_molformula': 'MolecularFormula',
            
            # Complexity
            'complexity': 'Complexity',
            'molecular_complexity': 'Complexity',
            'bertz_complexity': 'Complexity',
            
            # Heavy atoms
            'heavy_atom_count': 'HeavyAtomCount',
            'num_heavy_atoms': 'HeavyAtomCount',
            'num_atoms': 'HeavyAtomCount',
            
            # Aromatic rings
            'aromatic_rings': 'AromaticRingsCount',
            'num_aromatic_rings': 'AromaticRingsCount',
            'ar': 'AromaticRingsCount',
            
            # Solubility
            'solubility': 'AqueousSolubility',
            'aqueous_solubility': 'AqueousSolubility',
            'water_solubility': 'AqueousSolubility',
            'logs': 'LogS',
            'solubility_classification': 'SolubilityClass',
            
            # Melting/boiling points
            'melting_point': 'MeltingPoint',
            'boiling_point': 'BoilingPoint',
            'mp': 'MeltingPoint',
            'bp': 'BoilingPoint',
            
            # pKa
            'pka': 'pKa',
            'pka_strongest_acidic': 'pKaStrongestAcidic',
            'pka_strongest_basic': 'pKaStrongestBasic',
            
            # Drug-likeness
            'lipinski_rule_of_5_violations': 'LipinskiViolations',
            'ro5_violations': 'LipinskiViolations',
            'num_lipinski_ro5_violations': 'LipinskiViolations',
            'rule_of_5_violations': 'LipinskiViolations',
            'drug_likeness': 'DrugLikeness',
            'qed': 'QED',
            'qed_weighted': 'QED',
        }
        
        # Map property names to types
        self.property_type_map = {
            'LogP': 'physicochemical',
            'LogD': 'physicochemical',
            'MolecularWeight': 'physicochemical',
            'MonoisotopicMass': 'physicochemical',
            'TPSA': 'physicochemical',
            'HBondAcceptorCount': 'structural',
            'HBondDonorCount': 'structural',
            'RotatableBondCount': 'structural',
            'MolecularFormula': 'structural',
            'Complexity': 'structural',
            'HeavyAtomCount': 'structural',
            'AromaticRingsCount': 'structural',
            'AqueousSolubility': 'physicochemical',
            'LogS': 'physicochemical',
            'SolubilityClass': 'physicochemical',
            'MeltingPoint': 'physicochemical',
            'BoilingPoint': 'physicochemical',
            'pKa': 'physicochemical',
            'pKaStrongestAcidic': 'physicochemical',
            'pKaStrongestBasic': 'physicochemical',
            'LipinskiViolations': 'druglikeness',
            'DrugLikeness': 'druglikeness',
            'QED': 'druglikeness',
        }
    
    def _init_unit_conversions(self) -> None:
        """Initialize unit conversion factors and names."""
        # Map property names to standard units
        self.standard_units = {
            'LogP': 'log units',
            'LogD': 'log units',
            'MolecularWeight': 'g/mol',
            'MonoisotopicMass': 'g/mol',
            'TPSA': 'Å²',
            'AqueousSolubility': 'mg/mL',
            'LogS': 'log(mol/L)',
            'MeltingPoint': '°C',
            'BoilingPoint': '°C',
            'pKa': 'pKa units',
            'pKaStrongestAcidic': 'pKa units',
            'pKaStrongestBasic': 'pKa units',
        }
        
        # Unit conversion factors
        self.unit_conversions = {
            # Temperature conversions
            ('°F', '°C'): lambda x: (x - 32) * 5/9,
            ('K', '°C'): lambda x: x - 273.15,
            ('°C', '°F'): lambda x: x * 9/5 + 32,
            ('°C', 'K'): lambda x: x + 273.15,

            # Solubility conversions
            ('µg/mL', 'mg/mL'): lambda x: x / 1000,
            ('g/L', 'mg/mL'): lambda x: x,
            ('mg/L', 'mg/mL'): lambda x: x / 1000,
            ('mol/L', 'mg/mL'): lambda x, mw: x * mw if mw else x,  # mw in g/mol -> mg/mL,

            # Area conversions
            ('Å²', 'nm²'): lambda x: x / 100,
            ('nm²', 'Å²'): lambda x: x * 100,

            # Molecular weight conversions
            ('g/mol', 'kDa'): lambda x: x / 1000,
            ('kDa', 'g/mol'): lambda x: x * 1000,
        }
        
        # Map source-specific unit names to standard names
        self.unit_name_map = {
            # Temperature
            'c': '°C',
            'celsius': '°C',
            'f': '°F',
            'fahrenheit': '°F',
            'k': 'K',
            'kelvin': 'K',
            
            # Molecular weight
            'g/mol': 'g/mol',
            'daltons': 'g/mol',
            'da': 'g/mol',
            'kda': 'kDa',
            'amu': 'g/mol',
            
            # Surface area
            'a^2': 'Å²',
            'a2': 'Å²',
            'angstrom^2': 'Å²',
            'angstrom2': 'Å²',
            'square angstroms': 'Å²',
            'nm^2': 'nm²',
            'nm2': 'nm²',
            'square nanometers': 'nm²',
            
            # Solubility
            'mg/ml': 'mg/mL',
            'g/l': 'g/L',
            'mg/l': 'mg/L',
            'µg/ml': 'µg/mL',
            'ug/ml': 'µg/mL',
            'mol/l': 'mol/L',
            'mmol/l': 'mmol/L',
            'µmol/l': 'µmol/L',
            'umol/l': 'µmol/L',
            
            # Other
            'log units': 'log units',
            'pka units': 'pKa units',
        }
    
    def standardize_property(
        self,
        property_data: Dict[str, Any],
        molecule_data: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Standardize a property data dictionary.

        Args:
            property_data: Property data dictionary
            molecule_data: Optional molecule data for context

        Returns:
            Standardized property data dictionary
        """
        # Make a copy to avoid modifying the original
        result = property_data.copy()

        # Standardize property name
        if 'property_name' in result:
            result['property_name'] = self.standardize_property_name(result['property_name'])

        # Set property type if missing
        if ('property_name' in result and
            ('property_type' not in result or not result['property_type'])):
            result['property_type'] = self.get_property_type(result['property_name'])

        # Standardize unit if present and standardization is enabled
        if self.standardize_units and 'unit' in result and result['unit']:
            result['unit'] = self.standardize_unit(result['unit'])

            # Get target unit - either from preferred_units, standard_units, or keep as is
            target_unit = None
            if 'property_name' in result:
                prop_name = result['property_name']
                if prop_name in self.preferred_units:
                    target_unit = self.preferred_units[prop_name]
                elif prop_name in self.standard_units:
                    target_unit = self.standard_units[prop_name]

            # Convert to target unit if needed and possible
            if (target_unit and result['unit'] != target_unit and
                'numeric_value' in result and result['numeric_value'] is not None):
                converted_value = self.convert_unit(
                    result['numeric_value'],
                    result['unit'],
                    target_unit,
                    molecule_data
                )
                if converted_value is not None:
                    result['numeric_value'] = converted_value
                    result['unit'] = target_unit
                    # Store original value and unit if configured
                    if self.config.get('preserve_original_values', False):
                        result['original_value'] = property_data.get('numeric_value')
                        result['original_unit'] = property_data.get('unit')

        # Add standard unit if missing and known
        if (self.standardize_units and
            ('unit' not in result or not result['unit']) and
            'property_name' in result):
            prop_name = result['property_name']
            if prop_name in self.preferred_units:
                result['unit'] = self.preferred_units[prop_name]
            elif prop_name in self.standard_units:
                result['unit'] = self.standard_units[prop_name]

        # Ensure data source is set
        if 'source' not in result or not result['source']:
            result['source'] = 'Unknown'

        # Add any additional metadata from config
        if 'additional_metadata' in self.config:
            for key, value in self.config['additional_metadata'].items():
                if key not in result:
                    result[key] = value

        return result
    
    def standardize_property_name(self, property_name: str) -> str:
        """
        Standardize a property name to a canonical form.
        
        Args:
            property_name: Property name to standardize
            
        Returns:
            Standardized property name
        """
        # Clean and lowercase the name
        clean_name = property_name.strip().lower()
        clean_name = re.sub(r'[^a-z0-9_]', '', clean_name.replace(' ', '_'))
        
        # Look up in mapping
        if clean_name in self.property_name_map:
            return self.property_name_map[clean_name]
        
        # Handle special cases with regex
        if re.match(r'^(xlogp|alogp|clogp)[0-9]?$', clean_name):
            return 'LogP'
        
        if re.match(r'^logd[0-9.]*$', clean_name):
            return 'LogD'
        
        # If no match, capitalize words
        words = property_name.split('_')
        return ''.join(word.capitalize() for word in words)
    
    def get_property_type(self, property_name: str) -> str:
        """
        Get the type of a property based on its standardized name.
        
        Args:
            property_name: Standardized property name
            
        Returns:
            Property type
        """
        if property_name in self.property_type_map:
            return self.property_type_map[property_name]
        
        # Default to 'other'
        return 'other'
    
    def standardize_unit(self, unit: str) -> str:
        """
        Standardize a unit to canonical form.
        
        Args:
            unit: Unit to standardize
            
        Returns:
            Standardized unit
        """
        # Clean and lowercase the unit
        clean_unit = unit.strip().lower()
        
        # Look up in mapping
        if clean_unit in self.unit_name_map:
            return self.unit_name_map[clean_unit]
        
        # If no match, return as is
        return unit
    
    def convert_unit(
        self,
        value: float,
        from_unit: str,
        to_unit: str,
        molecule_data: Optional[Dict[str, Any]] = None
    ) -> Optional[float]:
        """
        Convert a value from one unit to another.
        
        Args:
            value: Value to convert
            from_unit: Source unit
            to_unit: Target unit
            molecule_data: Optional molecule data for context (e.g., molecular weight)
            
        Returns:
            Converted value or None if conversion not possible
        """
        # If units are the same, no conversion needed
        if from_unit == to_unit:
            return value
        
        # Standardize units
        from_unit_std = self.standardize_unit(from_unit)
        to_unit_std = self.standardize_unit(to_unit)
        
        # Check if we have a conversion function
        if (from_unit_std, to_unit_std) in self.unit_conversions:
            converter = self.unit_conversions[(from_unit_std, to_unit_std)]
            
            # Some conversions need molecular properties (e.g., MW)
            if 'mol/L' in (from_unit_std, to_unit_std) and molecule_data:
                mw = molecule_data.get('molecular_weight')
                return converter(value, mw)
            else:
                try:
                    return converter(value)
                except TypeError:
                    # Converter needs additional context we don't have
                    return None
        
        # No conversion available
        return None
    
    def standardize_properties(
        self,
        properties: List[Dict[str, Any]],
        molecule_data: Optional[Dict[str, Any]] = None
    ) -> List[Dict[str, Any]]:
        """
        Standardize a list of property dictionaries.

        Args:
            properties: List of property dictionaries
            molecule_data: Optional molecule data for context

        Returns:
            List of standardized property dictionaries
        """
        standardized = [
            self.standardize_property(prop, molecule_data)
            for prop in properties
        ]

        # Add any missing standard properties if configured
        if self.add_missing_properties and molecule_data:
            standardized = self._add_missing_standard_properties(standardized, molecule_data)

        return standardized

    def _add_missing_standard_properties(
        self,
        properties: List[Dict[str, Any]],
        molecule_data: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """
        Add any missing standard properties using the molecule data.

        This uses any available information in the molecule data to add
        missing properties, particularly those that can be derived from
        molecular structure.

        Args:
            properties: List of property dictionaries
            molecule_data: Molecule data dictionary

        Returns:
            Augmented list of property dictionaries
        """
        result = properties.copy()

        # Create a set of existing property names
        existing_props = {
            prop['property_name']
            for prop in result
            if 'property_name' in prop
        }

        # Check for calculated properties in molecule data
        if 'properties' in molecule_data and isinstance(molecule_data['properties'], dict):
            mol_props = molecule_data['properties']

            # Add missing properties from molecule data
            for prop_name, value in mol_props.items():
                # Skip if property already exists
                if prop_name in existing_props:
                    continue

                # Create new property entry
                new_prop = {
                    'property_name': prop_name,
                    'property_type': self.get_property_type(prop_name),
                    'numeric_value': value,
                    'source': 'Derived'
                }

                # Add standard unit if known
                if prop_name in self.standard_units:
                    new_prop['unit'] = self.standard_units[prop_name]

                result.append(new_prop)
                existing_props.add(prop_name)

        return result
    
    def merge_properties(
        self,
        primary_properties: List[Dict[str, Any]],
        secondary_properties: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Merge two lists of property data, keeping the best values.
        
        If the same property appears in both lists, the primary list's
        value is kept unless it is missing or the secondary has a higher
        confidence score.
        
        Args:
            primary_properties: Primary list of property dictionaries
            secondary_properties: Secondary list of property dictionaries
            
        Returns:
            Merged list of property dictionaries
        """
        # Create a copy of the primary properties
        result = primary_properties.copy()
        
        # Index primary properties by name
        primary_by_name = {
            prop['property_name']: i
            for i, prop in enumerate(result)
            if 'property_name' in prop
        }
        
        # Process secondary properties
        for sec_prop in secondary_properties:
            if 'property_name' not in sec_prop:
                continue
                
            prop_name = sec_prop['property_name']
            
            # If property exists in primary, check if we should replace
            if prop_name in primary_by_name:
                primary_idx = primary_by_name[prop_name]
                primary_prop = result[primary_idx]
                
                # Replace if primary value is missing or has lower confidence
                primary_value_missing = (
                    'numeric_value' not in primary_prop or
                    primary_prop['numeric_value'] is None
                )
                
                primary_confidence = primary_prop.get('confidence_score', 0.0)
                secondary_confidence = sec_prop.get('confidence_score', 0.0)
                
                # Replace if primary is missing or secondary has higher confidence
                if (primary_value_missing or
                    (secondary_confidence > 0 and secondary_confidence > primary_confidence)):
                    result[primary_idx] = sec_prop
            else:
                # Property doesn't exist in primary, add it
                result.append(sec_prop)
                primary_by_name[prop_name] = len(result) - 1
        
        return result