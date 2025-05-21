"""
Molecule transformation utilities.

This module provides functions for transforming molecular data between
different formats and representations, including converting between
different chemical identifier types.
"""

import logging
import asyncio
import aiohttp
import time
import re
import json
from typing import Dict, List, Any, Optional, Tuple, Union, Set, Callable

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, DataStructs, MolStandardize
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


class MoleculeTransformer:
    """
    Transforms molecular data between different formats and representations.

    This class provides methods for converting between different chemical
    identifiers (SMILES, InChI, InChIKey), calculating molecular properties,
    resolving cross-references between databases, and handling molecular mixtures.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None, config: Optional[Dict[str, Any]] = None):
        """
        Initialize molecule transformer.

        Args:
            logger: Logger instance
            config: Configuration dictionary with options for transformation
        """
        self.logger = logger or logging.getLogger(__name__)
        self.config = config or {}

        # Configure rate limiting for API calls
        self.api_delay = self.config.get('api_delay', 0.5)
        self.last_request_time = time.time()
        self._session = None

        # Configure cross-reference resolution options
        self.resolve_cross_references = self.config.get('resolve_cross_references', True)
        self.pubchem_api_url = self.config.get('pubchem_api_url', 'https://pubchem.ncbi.nlm.nih.gov/rest/pug')
        self.chembl_api_url = self.config.get('chembl_api_url', 'https://www.ebi.ac.uk/chembl/api/data')

        # Configure mixture handling options
        self.handle_mixtures = self.config.get('handle_mixtures', True)
        self.separate_mixtures = self.config.get('separate_mixtures', False)

        # Set up batch processing options
        self.batch_size = self.config.get('batch_size', 50)
        self.max_retries = self.config.get('max_retries', 3)
        self.retry_delay = self.config.get('retry_delay', 2.0)

        # Initialize RDKit components if available
        if not RDKIT_AVAILABLE:
            self.logger.warning(
                "RDKit is not available. Some molecule transformation "
                "functionality will be limited."
            )
    
    async def _get_session(self) -> aiohttp.ClientSession:
        """
        Get or create the HTTP session.

        Returns:
            aiohttp.ClientSession instance
        """
        if self._session is None or self._session.closed:
            self._session = aiohttp.ClientSession()
        return self._session

    async def _close_session(self) -> None:
        """Close the HTTP session if open."""
        if self._session and not self._session.closed:
            await self._session.close()
            self._session = None

    async def _apply_rate_limiting(self) -> None:
        """Apply rate limiting to avoid overwhelming APIs."""
        now = time.time()
        time_since_last_request = now - self.last_request_time

        if time_since_last_request < self.api_delay:
            await asyncio.sleep(self.api_delay - time_since_last_request)

        self.last_request_time = time.time()

    async def standardize_molecule(
        self,
        molecule_data: Dict[str, Any],
        resolve_ids: bool = True
    ) -> Dict[str, Any]:
        """
        Standardize a molecule data dictionary.

        Ensures all required fields are present and correctly formatted.
        Converts between identifiers and adds calculated properties.
        Optionally resolves cross-database identifiers.

        Args:
            molecule_data: Molecule data dictionary
            resolve_ids: Whether to resolve cross-database IDs (ChEMBL <-> PubChem)

        Returns:
            Standardized molecule data dictionary
        """
        # Make a copy to avoid modifying the original
        result = molecule_data.copy()

        # Ensure name is present
        if 'name' not in result or not result['name']:
            if 'pubchem_cid' in result and result['pubchem_cid']:
                result['name'] = f"Compound {result['pubchem_cid']}"
            elif 'chembl_id' in result and result['chembl_id']:
                result['name'] = result['chembl_id']
            elif 'id' in result:
                result['name'] = f"Molecule {result['id']}"
            else:
                result['name'] = "Unknown Compound"

        # Set data source if not present
        if 'data_source' not in result:
            result['data_source'] = "Unknown"

        # Check if molecule might be a mixture
        is_mixture = False
        if self.handle_mixtures and RDKIT_AVAILABLE:
            is_mixture, components = await self._check_if_mixture(result)
            result['is_mixture'] = is_mixture
            if is_mixture:
                result['mixture_components'] = components

        # Add missing fields
        if RDKIT_AVAILABLE:
            result = await self._add_missing_identifiers(result)
            result = await self._add_calculated_properties(result)

            # Add molecular scaffold if not already present
            if ('smiles' in result and result['smiles'] and
                ('scaffold_smiles' not in result or not result['scaffold_smiles'])):
                result = await self._add_molecular_scaffold(result)

        # Resolve cross-database IDs if requested and enabled in config
        if resolve_ids and self.resolve_cross_references:
            try:
                result = await self._resolve_cross_database_ids(result)
            except Exception as e:
                self.logger.warning(f"Error resolving cross-database IDs: {str(e)}")

        return result
    
    async def _add_missing_identifiers(
        self,
        molecule_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Add missing chemical identifiers to molecule data.

        If SMILES is present but InChI or InChIKey is missing, calculates
        the missing identifiers. Similarly for other combinations.

        Args:
            molecule_data: Molecule data dictionary

        Returns:
            Molecule data with missing identifiers added
        """
        result = molecule_data.copy()

        try:
            # Try to get a valid RDKit molecule
            mol = None

            # Try from SMILES
            if 'smiles' in result and result['smiles']:
                mol = Chem.MolFromSmiles(result['smiles'])

                # If molecule conversion succeeded, standardize the SMILES
                if mol is not None:
                    result['smiles'] = normalize_smiles(result['smiles'])

            # Try from InChI if SMILES failed
            if mol is None and 'inchi' in result and result['inchi']:
                mol = Chem.MolFromInchi(result['inchi'])

            # If we have a valid molecule, calculate missing identifiers
            if mol is not None:
                # Add standardized SMILES if missing
                if 'smiles' not in result or not result['smiles']:
                    result['smiles'] = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)

                # Add canonical SMILES (non-isomeric) if missing
                if 'canonical_smiles' not in result or not result['canonical_smiles']:
                    result['canonical_smiles'] = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)

                # Add InChI if missing
                if 'inchi' not in result or not result['inchi']:
                    result['inchi'] = Chem.MolToInchi(mol)

                # Add InChIKey if missing
                if 'inchikey' not in result or not result['inchikey']:
                    if 'inchi' in result and result['inchi']:
                        result['inchikey'] = Chem.InchiToInchiKey(result['inchi'])
                    else:
                        result['inchikey'] = Chem.MolToInchiKey(mol)

                # Add molecular formula if missing
                if 'formula' not in result or not result['formula']:
                    result['formula'] = AllChem.CalcMolFormula(mol)

                # Add fingerprint for similarity searching if missing
                if 'fingerprint' not in result or not result['fingerprint']:
                    fp = get_molecule_fingerprint(result['smiles'])
                    if fp:
                        result['fingerprint'] = fp

                # Add standardized structure if requested in config
                if self.config.get('add_standardized_structure', False):
                    result = await self._add_standardized_structure(result, mol)

        except Exception as e:
            self.logger.warning(
                f"Error adding missing identifiers for {result.get('name', 'unknown')}: {str(e)}"
            )

        return result

    async def _add_standardized_structure(self, molecule_data: Dict[str, Any], mol=None) -> Dict[str, Any]:
        """
        Add standardized molecular structure to molecule data.

        Uses RDKit's structure standardization utilities to create a canonical
        representation of the molecular structure.

        Args:
            molecule_data: Molecule data dictionary
            mol: Optional RDKit molecule object (will be generated if not provided)

        Returns:
            Molecule data with standardized structure added
        """
        result = molecule_data.copy()

        try:
            if not mol and 'smiles' in result and result['smiles']:
                mol = Chem.MolFromSmiles(result['smiles'])

            if mol is not None:
                # Use RDKit's standardization utilities
                standardizer = MolStandardize.Standardizer()
                standard_mol = standardizer.standardize(mol)

                # Add standardized SMILES
                result['standardized_smiles'] = Chem.MolToSmiles(standard_mol, isomericSmiles=True, canonical=True)

                # Add standardized InChI and InChIKey
                std_inchi = Chem.MolToInchi(standard_mol)
                result['standardized_inchi'] = std_inchi
                result['standardized_inchikey'] = Chem.InchiToInchiKey(std_inchi)

        except Exception as e:
            self.logger.warning(
                f"Error adding standardized structure for {result.get('name', 'unknown')}: {str(e)}"
            )

        return result

    async def _add_molecular_scaffold(self, molecule_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Add molecular scaffold to molecule data.

        Calculates the Bemis-Murcko scaffold, which represents the core structure
        of the molecule without side chains.

        Args:
            molecule_data: Molecule data dictionary

        Returns:
            Molecule data with scaffold added
        """
        result = molecule_data.copy()

        try:
            if 'smiles' in result and result['smiles']:
                mol = Chem.MolFromSmiles(result['smiles'])

                if mol is not None:
                    # Get Bemis-Murcko scaffold
                    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                    scaffold_smiles = Chem.MolToSmiles(scaffold, isomericSmiles=True, canonical=True)
                    result['scaffold_smiles'] = scaffold_smiles

                    # Get Generic Scaffold (removing atom types)
                    generic_scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
                    generic_smiles = Chem.MolToSmiles(generic_scaffold, isomericSmiles=False, canonical=True)
                    result['generic_scaffold'] = generic_smiles
        except Exception as e:
            self.logger.warning(
                f"Error adding molecular scaffold for {result.get('name', 'unknown')}: {str(e)}"
            )

        return result
    
    async def _add_calculated_properties(
        self,
        molecule_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Add calculated molecular properties to molecule data.

        Calculates properties like molecular weight, LogP, etc. if missing.

        Args:
            molecule_data: Molecule data dictionary

        Returns:
            Molecule data with calculated properties added
        """
        result = molecule_data.copy()

        try:
            # Try to get a valid RDKit molecule
            mol = None

            # Try from SMILES
            if 'smiles' in result and result['smiles']:
                mol = Chem.MolFromSmiles(result['smiles'])

            # Try from InChI if SMILES failed
            if mol is None and 'inchi' in result and result['inchi']:
                mol = Chem.MolFromInchi(result['inchi'])

            # If we have a valid molecule, calculate properties
            if mol is not None:
                # Initialize properties list and properties dict for JSONB field
                properties = []
                if 'properties' not in result:
                    result['properties'] = {}

                # Add molecular weight if missing
                if 'molecular_weight' not in result or not result['molecular_weight']:
                    mol_weight = Descriptors.MolWt(mol)
                    result['molecular_weight'] = mol_weight
                    result['properties']['MolecularWeight'] = mol_weight

                # Calculate other properties
                # LogP
                logp = Descriptors.MolLogP(mol)
                properties.append({
                    'property_name': 'LogP',
                    'property_type': 'physicochemical',
                    'numeric_value': logp,
                    'source': 'Calculated'
                })
                result['properties']['LogP'] = logp

                # TPSA
                tpsa = Descriptors.TPSA(mol)
                properties.append({
                    'property_name': 'TPSA',
                    'property_type': 'physicochemical',
                    'numeric_value': tpsa,
                    'unit': 'Å²',
                    'source': 'Calculated'
                })
                result['properties']['TPSA'] = tpsa

                # HBond donors
                hbd = Descriptors.NumHDonors(mol)
                properties.append({
                    'property_name': 'HBondDonorCount',
                    'property_type': 'structural',
                    'numeric_value': hbd,
                    'source': 'Calculated'
                })
                result['properties']['HBondDonorCount'] = hbd

                # HBond acceptors
                hba = Descriptors.NumHAcceptors(mol)
                properties.append({
                    'property_name': 'HBondAcceptorCount',
                    'property_type': 'structural',
                    'numeric_value': hba,
                    'source': 'Calculated'
                })
                result['properties']['HBondAcceptorCount'] = hba

                # Rotatable bonds
                rtb = Descriptors.NumRotatableBonds(mol)
                properties.append({
                    'property_name': 'RotatableBondCount',
                    'property_type': 'structural',
                    'numeric_value': rtb,
                    'source': 'Calculated'
                })
                result['properties']['RotatableBondCount'] = rtb

                # Heavy atom count
                heavy_atoms = mol.GetNumHeavyAtoms()
                properties.append({
                    'property_name': 'HeavyAtomCount',
                    'property_type': 'structural',
                    'numeric_value': heavy_atoms,
                    'source': 'Calculated'
                })
                result['properties']['HeavyAtomCount'] = heavy_atoms

                # Ring count
                rings = Chem.GetSSSR(mol)
                ring_count = len(rings)
                properties.append({
                    'property_name': 'RingCount',
                    'property_type': 'structural',
                    'numeric_value': ring_count,
                    'source': 'Calculated'
                })
                result['properties']['RingCount'] = ring_count

                # Aromatic ring count
                aromatic_rings = sum(1 for ring in rings if Chem.MolToSmiles(mol, kekuleSmiles=False, rootedAtAtom=ring[0]).count('c') >= 6)
                properties.append({
                    'property_name': 'AromaticRingCount',
                    'property_type': 'structural',
                    'numeric_value': aromatic_rings,
                    'source': 'Calculated'
                })
                result['properties']['AromaticRingCount'] = aromatic_rings

                # Complexity
                complexity = Descriptors.BertzCT(mol)
                properties.append({
                    'property_name': 'Complexity',
                    'property_type': 'structural',
                    'numeric_value': complexity,
                    'source': 'Calculated'
                })
                result['properties']['Complexity'] = complexity

                # Fraction of sp3 hybridized carbons
                fsp3 = Descriptors.FractionCSP3(mol)
                properties.append({
                    'property_name': 'FractionCSP3',
                    'property_type': 'structural',
                    'numeric_value': fsp3,
                    'source': 'Calculated'
                })
                result['properties']['FractionCSP3'] = fsp3

                # QED (Quantitative Estimation of Drug-likeness)
                try:
                    from rdkit.Chem import QED
                    qed = QED.qed(mol)
                    properties.append({
                        'property_name': 'QED',
                        'property_type': 'druglikeness',
                        'numeric_value': qed,
                        'source': 'Calculated'
                    })
                    result['properties']['QED'] = qed
                except Exception as qed_error:
                    self.logger.debug(f"Error calculating QED: {str(qed_error)}")

                # Lipinski rule of 5 violations
                lipinski_violations = 0
                if mol_weight > 500: lipinski_violations += 1
                if logp > 5: lipinski_violations += 1
                if hbd > 5: lipinski_violations += 1
                if hba > 10: lipinski_violations += 1

                properties.append({
                    'property_name': 'LipinskiViolations',
                    'property_type': 'druglikeness',
                    'numeric_value': lipinski_violations,
                    'source': 'Calculated'
                })
                result['properties']['LipinskiViolations'] = lipinski_violations

                # Store calculated properties
                if 'calculated_properties' not in result:
                    result['calculated_properties'] = []

                result['calculated_properties'].extend(properties)
        except Exception as e:
            self.logger.warning(
                f"Error calculating properties for {result.get('name', 'unknown')}: {str(e)}"
            )

        return result
    
    async def merge_molecule_data(
        self,
        primary_data: Dict[str, Any],
        secondary_data: Dict[str, Any],
        prefer_source: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Merge data from two molecule records.

        Takes the primary molecule data and enhances it with any
        additional information from the secondary data. Can optionally
        prefer data from a specific source.

        Args:
            primary_data: Primary molecule data dictionary
            secondary_data: Secondary molecule data dictionary
            prefer_source: Data source to prefer (e.g., "ChEMBL", "PubChem")

        Returns:
            Merged molecule data dictionary
        """
        # Make a deep copy to avoid modifying the original
        result = primary_data.copy()

        # Fields to always keep from primary
        primary_only_fields = {'id', 'created_at', 'updated_at'}

        # Determine if we should prefer the secondary data
        prefer_secondary = False
        if prefer_source:
            primary_source = primary_data.get('data_source', '')
            secondary_source = secondary_data.get('data_source', '')
            if secondary_source == prefer_source and primary_source != prefer_source:
                prefer_secondary = True

        # Merge fields from secondary data
        for key, value in secondary_data.items():
            # Skip primary-only fields
            if key in primary_only_fields:
                continue

            # If field is missing in primary or we prefer secondary, add it
            if key not in result or result[key] is None or result[key] == '' or prefer_secondary:
                result[key] = value

            # Special handling for specific fields
            if key == 'synonyms' and 'synonyms' in result and isinstance(result['synonyms'], list):
                # Merge synonyms
                if isinstance(secondary_data.get('synonyms'), list):
                    primary_names = {name.lower() for name in result['synonyms'] if isinstance(name, str)}
                    for name in secondary_data['synonyms']:
                        if isinstance(name, str) and name.lower() not in primary_names:
                            result['synonyms'].append(name)
                            primary_names.add(name.lower())

            # Merge JSONB properties fields
            elif key == 'properties' and isinstance(result.get('properties'), dict):
                if isinstance(secondary_data.get('properties'), dict):
                    for prop_key, prop_value in secondary_data['properties'].items():
                        if prop_key not in result['properties'] or prefer_secondary:
                            result['properties'][prop_key] = prop_value

        # Merge calculated properties
        if ('calculated_properties' in secondary_data and
            secondary_data['calculated_properties']):

            if 'calculated_properties' not in result:
                result['calculated_properties'] = []

            # Get existing property names
            existing_props = {
                prop['property_name']
                for prop in result['calculated_properties']
                if isinstance(prop, dict) and 'property_name' in prop
            }

            # Add non-duplicate properties or replace if we prefer secondary
            for prop in secondary_data['calculated_properties']:
                if isinstance(prop, dict) and 'property_name' in prop:
                    if prop['property_name'] not in existing_props or prefer_secondary:
                        # If we prefer secondary and it already exists, find and replace it
                        if prefer_secondary and prop['property_name'] in existing_props:
                            for i, existing_prop in enumerate(result['calculated_properties']):
                                if ('property_name' in existing_prop and
                                    existing_prop['property_name'] == prop['property_name']):
                                    result['calculated_properties'][i] = prop
                                    break
                        else:
                            result['calculated_properties'].append(prop)
                            existing_props.add(prop['property_name'])

        return result


    async def _check_if_mixture(self, molecule_data: Dict[str, Any]) -> Tuple[bool, List[Dict[str, Any]]]:
        """
        Check if a molecule is a mixture and extract components if possible.

        This checks for special patterns in the SMILES string, like dots
        indicating disconnected structures (which often represent mixtures),
        and tries to extract and identify the individual components.

        Args:
            molecule_data: Molecule data dictionary

        Returns:
            Tuple of (is_mixture, components_list)
        """
        is_mixture = False
        components = []

        try:
            if 'smiles' in molecule_data and molecule_data['smiles']:
                smiles = molecule_data['smiles']

                # Check for dots in SMILES, which indicate disconnected components
                if '.' in smiles:
                    # Split the SMILES into components
                    component_smiles = smiles.split('.')

                    # If we have multiple components, it's a mixture
                    if len(component_smiles) > 1:
                        is_mixture = True

                        # Process each component
                        for i, comp_smiles in enumerate(component_smiles):
                            # Skip empty components
                            if not comp_smiles.strip():
                                continue

                            mol = Chem.MolFromSmiles(comp_smiles)
                            if mol is None:
                                continue

                            # Calculate component properties
                            component = {
                                'smiles': comp_smiles,
                                'canonical_smiles': Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True),
                                'molecular_weight': Descriptors.MolWt(mol),
                                'formula': AllChem.CalcMolFormula(mol),
                                'index': i
                            }

                            # Add additional identifiers
                            component['inchi'] = Chem.MolToInchi(mol)
                            component['inchikey'] = Chem.InchiToInchiKey(component['inchi'])

                            # Try to determine component name
                            component_name = None
                            # If the original molecule has synonyms, try to find one for this component
                            if 'synonyms' in molecule_data and molecule_data['synonyms']:
                                for synonym in molecule_data['synonyms']:
                                    # Look for patterns like "Component (50%)" or similar
                                    if '(' in synonym and any(c.isdigit() for c in synonym):
                                        component_name = synonym
                                        break

                            # If no specific name found, use generic name
                            if not component_name:
                                component_name = f"Component {i+1}"

                            component['name'] = component_name
                            components.append(component)
        except Exception as e:
            self.logger.warning(f"Error checking if molecule is a mixture: {str(e)}")

        return is_mixture, components

    async def standardize_molecules_batch(
        self,
        molecules_data: List[Dict[str, Any]],
        resolve_ids: bool = True
    ) -> List[Dict[str, Any]]:
        """
        Standardize multiple molecule data dictionaries in a batch.

        This is more efficient than calling standardize_molecule() repeatedly.

        Args:
            molecules_data: List of molecule data dictionaries
            resolve_ids: Whether to resolve cross-database IDs (ChEMBL <-> PubChem)

        Returns:
            List of standardized molecule data dictionaries
        """
        results = []

        # Process all molecules that don't need ID resolution first
        for mol_data in molecules_data:
            try:
                # Make a copy to avoid modifying the original
                mol_copy = mol_data.copy()

                # Ensure name is present
                if 'name' not in mol_copy or not mol_copy['name']:
                    if 'pubchem_cid' in mol_copy and mol_copy['pubchem_cid']:
                        mol_copy['name'] = f"Compound {mol_copy['pubchem_cid']}"
                    elif 'chembl_id' in mol_copy and mol_copy['chembl_id']:
                        mol_copy['name'] = mol_copy['chembl_id']
                    elif 'id' in mol_copy:
                        mol_copy['name'] = f"Molecule {mol_copy['id']}"
                    else:
                        mol_copy['name'] = "Unknown Compound"

                # Set data source if not present
                if 'data_source' not in mol_copy:
                    mol_copy['data_source'] = "Unknown"

                # Add missing identifiers and properties if RDKit is available
                if RDKIT_AVAILABLE:
                    mol_copy = await self._add_missing_identifiers(mol_copy)
                    mol_copy = await self._add_calculated_properties(mol_copy)

                results.append(mol_copy)
            except Exception as e:
                self.logger.error(f"Error standardizing molecule {mol_data.get('name', 'unknown')}: {str(e)}")
                results.append(mol_data)  # Keep the original in case of error

        # Now process ID resolution if requested
        if resolve_ids and self.resolve_cross_references:
            try:
                results = await self._resolve_batch_cross_references(results)
            except Exception as e:
                self.logger.warning(f"Error resolving batch cross-references: {str(e)}")

        return results

    async def _resolve_cross_database_ids(
        self,
        molecule_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Resolve cross-database identifiers (ChEMBL <-> PubChem).

        If the molecule has a ChEMBL ID but no PubChem CID, attempts to find
        the corresponding PubChem CID and vice versa.

        Args:
            molecule_data: Molecule data dictionary

        Returns:
            Molecule data with cross-references added
        """
        result = molecule_data.copy()

        try:
            # Check if we need to find PubChem CID
            if ('chembl_id' in result and result['chembl_id'] and
                ('pubchem_cid' not in result or not result['pubchem_cid'])):

                pubchem_cid = await self._resolve_chembl_to_pubchem(result['chembl_id'], result)
                if pubchem_cid:
                    result['pubchem_cid'] = pubchem_cid

            # Check if we need to find ChEMBL ID
            if ('pubchem_cid' in result and result['pubchem_cid'] and
                ('chembl_id' not in result or not result['chembl_id'])):

                chembl_id = await self._resolve_pubchem_to_chembl(result['pubchem_cid'], result)
                if chembl_id:
                    result['chembl_id'] = chembl_id
        except Exception as e:
            self.logger.warning(f"Error resolving cross-database IDs: {str(e)}")

        return result

    async def _resolve_batch_cross_references(
        self,
        molecules_data: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Resolve cross-database identifiers for a batch of molecules.

        This is more efficient than resolving each molecule individually.

        Args:
            molecules_data: List of molecule data dictionaries

        Returns:
            List of molecule data dictionaries with cross-references added
        """
        results = []
        chembl_tasks = []
        pubchem_tasks = []

        # Create tasks for resolving cross-references
        for mol_data in molecules_data:
            # Check if we need to find PubChem CID
            if ('chembl_id' in mol_data and mol_data['chembl_id'] and
                ('pubchem_cid' not in mol_data or not mol_data['pubchem_cid'])):

                chembl_tasks.append((mol_data, len(results)))

            # Check if we need to find ChEMBL ID
            elif ('pubchem_cid' in mol_data and mol_data['pubchem_cid'] and
                ('chembl_id' not in mol_data or not mol_data['chembl_id'])):

                pubchem_tasks.append((mol_data, len(results)))

            # No cross-reference needed
            else:
                pass

            # Add the molecule to results (will be updated later)
            results.append(mol_data.copy())

        # Process ChEMBL to PubChem tasks in batches
        for i in range(0, len(chembl_tasks), self.batch_size):
            batch = chembl_tasks[i:i + self.batch_size]
            try:
                # Collect ChEMBL IDs for the batch
                chembl_ids = [mol[0]['chembl_id'] for mol in batch]
                indices = [mol[1] for mol in batch]

                # Resolve the batch
                cid_mapping = await self._batch_resolve_chembl_to_pubchem(chembl_ids)

                # Update the results
                for chembl_id, idx in zip(chembl_ids, indices):
                    if chembl_id in cid_mapping and cid_mapping[chembl_id]:
                        results[idx]['pubchem_cid'] = cid_mapping[chembl_id]
            except Exception as e:
                self.logger.warning(f"Error in batch ChEMBL->PubChem resolution: {str(e)}")

            # Apply rate limiting
            await self._apply_rate_limiting()

        # Process PubChem to ChEMBL tasks in batches
        for i in range(0, len(pubchem_tasks), self.batch_size):
            batch = pubchem_tasks[i:i + self.batch_size]
            try:
                # Collect PubChem CIDs for the batch
                pubchem_cids = [mol[0]['pubchem_cid'] for mol in batch]
                indices = [mol[1] for mol in batch]

                # Resolve the batch
                chembl_mapping = await self._batch_resolve_pubchem_to_chembl(pubchem_cids)

                # Update the results
                for pubchem_cid, idx in zip(pubchem_cids, indices):
                    if pubchem_cid in chembl_mapping and chembl_mapping[pubchem_cid]:
                        results[idx]['chembl_id'] = chembl_mapping[pubchem_cid]
            except Exception as e:
                self.logger.warning(f"Error in batch PubChem->ChEMBL resolution: {str(e)}")

            # Apply rate limiting
            await self._apply_rate_limiting()

        return results

    async def _resolve_chembl_to_pubchem(
        self,
        chembl_id: str,
        molecule_data: Optional[Dict[str, Any]] = None
    ) -> Optional[str]:
        """
        Resolve a ChEMBL ID to a PubChem CID.

        Args:
            chembl_id: ChEMBL ID to resolve
            molecule_data: Optional molecule data for additional context

        Returns:
            PubChem CID or None if not found
        """
        try:
            # Try direct lookup first
            session = await self._get_session()
            await self._apply_rate_limiting()

            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceid/ChEMBL/{chembl_id}/cids/JSON"
            async with session.get(url) as response:
                if response.status == 200:
                    data = await response.json()
                    if 'InformationList' in data and 'Information' in data['InformationList']:
                        for item in data['InformationList']['Information']:
                            if 'CID' in item:
                                return str(item['CID'])

            # If direct lookup fails and we have molecule data, try looking up by InChIKey
            if molecule_data and 'inchikey' in molecule_data and molecule_data['inchikey']:
                inchikey = molecule_data['inchikey']
                await self._apply_rate_limiting()

                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                            return str(data['IdentifierList']['CID'][0])

            # As a last resort, try SMILES-based search
            if molecule_data and 'smiles' in molecule_data and molecule_data['smiles']:
                # This would require additional implementation for SMILES search in PubChem
                pass
        except Exception as e:
            self.logger.warning(f"Error resolving ChEMBL ID {chembl_id} to PubChem: {str(e)}")

        return None

    async def _batch_resolve_chembl_to_pubchem(
        self,
        chembl_ids: List[str]
    ) -> Dict[str, Optional[str]]:
        """
        Resolve multiple ChEMBL IDs to PubChem CIDs in a batch.

        Args:
            chembl_ids: List of ChEMBL IDs to resolve

        Returns:
            Dictionary mapping ChEMBL IDs to PubChem CIDs (or None if not found)
        """
        result = {chembl_id: None for chembl_id in chembl_ids}

        try:
            # Join IDs with commas for batch request
            batch_ids = ','.join(chembl_ids)

            session = await self._get_session()
            await self._apply_rate_limiting()

            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceid/ChEMBL/{batch_ids}/cids/JSON"
            async with session.get(url) as response:
                if response.status == 200:
                    data = await response.json()
                    if 'InformationList' in data and 'Information' in data['InformationList']:
                        for item in data['InformationList']['Information']:
                            if 'SourceID' in item and 'CID' in item:
                                source_id = item['SourceID']
                                # Extract the ChEMBL ID from the source ID
                                chembl_id = source_id.split('/')[-1] if '/' in source_id else source_id
                                if chembl_id in result:
                                    result[chembl_id] = str(item['CID'])
        except Exception as e:
            self.logger.warning(f"Error in batch resolution of ChEMBL IDs to PubChem: {str(e)}")

        return result

    async def _resolve_pubchem_to_chembl(
        self,
        pubchem_cid: str,
        molecule_data: Optional[Dict[str, Any]] = None
    ) -> Optional[str]:
        """
        Resolve a PubChem CID to a ChEMBL ID.

        Args:
            pubchem_cid: PubChem CID to resolve
            molecule_data: Optional molecule data for additional context

        Returns:
            ChEMBL ID or None if not found
        """
        try:
            # Try direct lookup first using PubChem's source links
            session = await self._get_session()
            await self._apply_rate_limiting()

            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/xrefs/SourceID/JSON"
            async with session.get(url) as response:
                if response.status == 200:
                    data = await response.json()
                    if 'InformationList' in data and 'Information' in data['InformationList']:
                        for item in data['InformationList']['Information']:
                            if 'SourceID' in item:
                                for source_id in item['SourceID']:
                                    if source_id.startswith('CHEMBL'):
                                        return source_id

            # If direct lookup fails, try ChEMBL API with InChIKey if available
            if molecule_data and 'inchikey' in molecule_data and molecule_data['inchikey']:
                inchikey = molecule_data['inchikey']
                await self._apply_rate_limiting()

                url = f"{self.chembl_api_url}/molecule?molecule_structures__standard_inchi_key={inchikey}"
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        if 'molecules' in data and data['molecules']:
                            return data['molecules'][0]['molecule_chembl_id']
        except Exception as e:
            self.logger.warning(f"Error resolving PubChem CID {pubchem_cid} to ChEMBL: {str(e)}")

        return None

    async def _batch_resolve_pubchem_to_chembl(
        self,
        pubchem_cids: List[str]
    ) -> Dict[str, Optional[str]]:
        """
        Resolve multiple PubChem CIDs to ChEMBL IDs in a batch.

        Args:
            pubchem_cids: List of PubChem CIDs to resolve

        Returns:
            Dictionary mapping PubChem CIDs to ChEMBL IDs (or None if not found)
        """
        result = {cid: None for cid in pubchem_cids}

        try:
            # PubChem doesn't have a good batch API for this, so we'll make individual requests
            # but process them concurrently
            session = await self._get_session()

            async def resolve_single_cid(cid):
                await self._apply_rate_limiting()
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/SourceID/JSON"
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        if 'InformationList' in data and 'Information' in data['InformationList']:
                            for item in data['InformationList']['Information']:
                                if 'SourceID' in item:
                                    for source_id in item['SourceID']:
                                        if source_id.startswith('CHEMBL'):
                                            return cid, source_id
                return cid, None

            # Create tasks for concurrent resolution (limited by semaphore to avoid overwhelming the API)
            tasks = [resolve_single_cid(cid) for cid in pubchem_cids]

            # Process in smaller concurrent batches
            for i in range(0, len(tasks), 5):  # Process 5 at a time
                batch_tasks = tasks[i:i+5]
                batch_results = await asyncio.gather(*batch_tasks)

                # Update the result dictionary
                for cid, chembl_id in batch_results:
                    if chembl_id:
                        result[cid] = chembl_id
        except Exception as e:
            self.logger.warning(f"Error in batch resolution of PubChem CIDs to ChEMBL: {str(e)}")

        return result

def normalize_smiles(smiles: str) -> str:
    """
    Normalize a SMILES string for consistent representation.

    Args:
        smiles: SMILES string

    Returns:
        Normalized SMILES string
    """
    if not RDKIT_AVAILABLE:
        return smiles

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception:
        return smiles


def get_molecule_fingerprint(smiles: str) -> Optional[str]:
    """
    Calculate a molecular fingerprint for similarity comparisons.

    Args:
        smiles: SMILES string

    Returns:
        Fingerprint string or None if calculation fails
    """
    if not RDKIT_AVAILABLE:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Calculate Morgan fingerprint
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)

        # Convert to hex string
        return fp.ToBitString()
    except Exception:
        return None


def calculate_similarity(
    smiles1: str,
    smiles2: str,
    fingerprint_type: str = 'morgan'
) -> float:
    """
    Calculate similarity between two molecules.

    Args:
        smiles1: First SMILES string
        smiles2: Second SMILES string
        fingerprint_type: Type of fingerprint to use ('morgan', 'maccs', 'topological')

    Returns:
        Similarity score (0-1) or 0 if calculation fails
    """
    if not RDKIT_AVAILABLE:
        return 0.0

    try:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)

        if mol1 is None or mol2 is None:
            return 0.0

        # Choose fingerprint type
        if fingerprint_type == 'maccs':
            from rdkit.Chem import MACCSkeys
            fp1 = MACCSkeys.GenMACCSKeys(mol1)
            fp2 = MACCSkeys.GenMACCSKeys(mol2)
        elif fingerprint_type == 'topological':
            fp1 = Chem.RDKFingerprint(mol1)
            fp2 = Chem.RDKFingerprint(mol2)
        else:  # Default to Morgan/ECFP
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, 2048)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, 2048)

        # Calculate Tanimoto similarity
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        return similarity
    except Exception:
        return 0.0