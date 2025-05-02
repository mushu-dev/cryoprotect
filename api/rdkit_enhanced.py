"""
CryoProtect Analyzer API - Enhanced RDKit Integration

This module extends the core RDKit utilities with additional functionality for the CryoProtect application.
It provides caching mechanisms, batch processing, and specialized functions for cryoprotectant analysis.

Features:
    - Molecular property calculation with caching
    - Substructure and similarity search optimized for cryoprotectants
    - Batch processing for multiple molecules
    - Specialized cryoprotectant property predictions
    - Molecular visualization with customization options
"""

import os
import json
import hashlib
import logging
import time
from typing import Dict, List, Optional, Tuple, Union, Any, Set
from functools import lru_cache
from datetime import datetime

# Import RDKit modules
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, MolSurf, AllChem, Fragments, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdFingerprintGenerator

# Import existing RDKit utilities
from api.rdkit_utils import (
    parse_molecule, calculate_hydrogen_bonding, calculate_logp, calculate_tpsa,
    calculate_molecular_properties, identify_functional_groups, estimate_permeability,
    calculate_all_properties, generate_molecule_svg, perform_substructure_search,
    calculate_similarity
)

# Set up logging
logger = logging.getLogger(__name__)

# Configure cache directory
CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'cache', 'rdkit')
os.makedirs(CACHE_DIR, exist_ok=True)

class PropertyCache:
    """
    Cache for molecular properties to avoid redundant calculations.
    
    This class provides both memory-based caching (using LRU cache) and disk-based
    caching for persistence between application restarts.
    """
    
    def __init__(self, cache_size=1000, cache_dir=CACHE_DIR):
        """
        Initialize the property cache.
        
        Args:
            cache_size: Maximum number of molecules to cache in memory
            cache_dir: Directory to store persistent cache files
        """
        self.cache_dir = cache_dir
        self.cache_size = cache_size
        os.makedirs(self.cache_dir, exist_ok=True)
        
        # Initialize in-memory cache
        self._init_memory_cache()
    
    def _init_memory_cache(self):
        """Initialize the in-memory LRU cache for molecular properties."""
        @lru_cache(maxsize=self.cache_size)
        def cached_calculate_properties(smiles_hash):
            # Check disk cache first
            cache_file = os.path.join(self.cache_dir, f"{smiles_hash}.json")
            if os.path.exists(cache_file):
                try:
                    with open(cache_file, 'r') as f:
                        return json.load(f)
                except (json.JSONDecodeError, IOError) as e:
                    logger.warning(f"Error reading cache file {cache_file}: {str(e)}")
            
            return None
        
        self.cached_calculate_properties = cached_calculate_properties
    
    def _get_smiles_hash(self, smiles: str) -> str:
        """
        Generate a hash for a SMILES string to use as a cache key.
        
        Args:
            smiles: SMILES string to hash
            
        Returns:
            Hash string for the SMILES
        """
        # Canonicalize the SMILES first to ensure consistent hashing
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            canonical_smiles = Chem.MolToSmiles(mol)
            return hashlib.md5(canonical_smiles.encode()).hexdigest()
        return hashlib.md5(smiles.encode()).hexdigest()
    
    def get(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Get cached properties for a molecule if available.
        
        Args:
            smiles: SMILES string of the molecule
            
        Returns:
            Dictionary of cached properties or None if not in cache
        """
        smiles_hash = self._get_smiles_hash(smiles)
        return self.cached_calculate_properties(smiles_hash)
    
    def set(self, smiles: str, properties: Dict[str, Any]) -> None:
        """
        Cache properties for a molecule.
        
        Args:
            smiles: SMILES string of the molecule
            properties: Dictionary of properties to cache
        """
        smiles_hash = self._get_smiles_hash(smiles)
        
        # Update disk cache
        cache_file = os.path.join(self.cache_dir, f"{smiles_hash}.json")
        try:
            with open(cache_file, 'w') as f:
                json.dump(properties, f)
        except IOError as e:
            logger.warning(f"Error writing to cache file {cache_file}: {str(e)}")
        
        # Update memory cache (just store the hash to trigger the cache)
        self.cached_calculate_properties(smiles_hash)
    
    def clear(self) -> None:
        """Clear both memory and disk caches."""
        # Clear memory cache
        self.cached_calculate_properties.cache_clear()
        
        # Clear disk cache
        for file in os.listdir(self.cache_dir):
            if file.endswith('.json'):
                try:
                    os.remove(os.path.join(self.cache_dir, file))
                except OSError as e:
                    logger.warning(f"Error removing cache file {file}: {str(e)}")

# Initialize global property cache
property_cache = PropertyCache()

def calculate_properties_with_cache(mol_data: str, input_format: str = 'smiles') -> Dict[str, Any]:
    """
    Calculate molecular properties with caching to avoid redundant calculations.
    
    Args:
        mol_data: Molecular data as a string (SMILES, MOL, SDF)
        input_format: Format of the input data ('smiles', 'mol', 'sdf')
        
    Returns:
        Dictionary of calculated properties
    """
    # For non-SMILES formats, convert to SMILES first
    if input_format != 'smiles':
        mol = parse_molecule(mol_data, input_format)
        if mol is None:
            return {"error": f"Failed to parse molecule from {input_format} data"}
        smiles = Chem.MolToSmiles(mol)
    else:
        smiles = mol_data
    
    # Check cache first
    cached_properties = property_cache.get(smiles)
    if cached_properties:
        logger.debug(f"Cache hit for SMILES: {smiles}")
        return cached_properties
    
    # Calculate properties if not in cache
    logger.debug(f"Cache miss for SMILES: {smiles}, calculating properties")
    properties = calculate_all_properties(mol_data, input_format)
    
    # Cache the results if calculation was successful
    if "error" not in properties:
        property_cache.set(smiles, properties)
    
    return properties

def calculate_cryoprotectant_properties(mol_data: str, input_format: str = 'smiles') -> Dict[str, Any]:
    """
    Calculate properties specifically relevant for cryoprotectants.
    
    This function extends the standard property calculations with additional
    properties that are particularly relevant for cryoprotectant applications.
    
    Args:
        mol_data: Molecular data as a string (SMILES, MOL, SDF)
        input_format: Format of the input data ('smiles', 'mol', 'sdf')
        
    Returns:
        Dictionary of calculated properties with cryoprotectant-specific additions
    """
    # Get base properties
    properties = calculate_properties_with_cache(mol_data, input_format)
    
    if "error" in properties:
        return properties
    
    # Parse molecule
    mol = parse_molecule(mol_data, input_format)
    if mol is None:
        return {"error": f"Failed to parse molecule from {input_format} data"}
    
    # Add cryoprotectant-specific properties
    cryo_properties = {}
    
    # 1. Glass transition temperature estimate (Tg)
    # This is a simplified model based on molecular weight and hydrogen bonding
    try:
        mw = properties["molecular_properties"]["molecular_weight"]
        hbond_total = properties["hydrogen_bonding"]["total"]
        logp = properties["logp"]
        
        # Simplified Tg estimation model (for demonstration)
        # Actual model would be based on experimental data and more sophisticated algorithms
        tg_estimate = -50 + (0.1 * mw) + (5 * hbond_total) - (10 * logp)
        cryo_properties["glass_transition_temp_estimate"] = tg_estimate
    except (KeyError, TypeError) as e:
        logger.warning(f"Error calculating glass transition temperature: {str(e)}")
        cryo_properties["glass_transition_temp_estimate"] = None
    
    # 2. Vitrification tendency score (0-100)
    try:
        # Factors that promote vitrification:
        # - High molecular weight
        # - Many hydrogen bond donors/acceptors
        # - Low logP (more hydrophilic)
        # - High TPSA
        
        mw_factor = min(100, properties["molecular_properties"]["molecular_weight"] / 10)
        hbond_factor = min(100, properties["hydrogen_bonding"]["total"] * 10)
        logp_factor = max(0, 100 - (properties["logp"] * 20))
        tpsa_factor = min(100, properties["tpsa"] / 2)
        
        vitrification_score = (mw_factor + hbond_factor + logp_factor + tpsa_factor) / 4
        cryo_properties["vitrification_tendency"] = vitrification_score
    except (KeyError, TypeError) as e:
        logger.warning(f"Error calculating vitrification tendency: {str(e)}")
        cryo_properties["vitrification_tendency"] = None
    
    # 3. Cell membrane permeability score for cryoprotection (0-100)
    try:
        # Ideal cryoprotectants balance permeability with stability
        # - Moderate logP (not too hydrophobic or hydrophilic)
        # - Moderate molecular weight (not too large)
        # - Moderate TPSA
        
        logp_opt = max(0, 100 - abs(properties["logp"] - 1.5) * 25)
        mw_opt = max(0, 100 - abs(properties["molecular_properties"]["molecular_weight"] - 150) / 3)
        tpsa_opt = max(0, 100 - abs(properties["tpsa"] - 60) / 1.5)
        
        permeability_score = (logp_opt + mw_opt + tpsa_opt) / 3
        cryo_properties["cryo_permeability_score"] = permeability_score
    except (KeyError, TypeError) as e:
        logger.warning(f"Error calculating cryo permeability score: {str(e)}")
        cryo_properties["cryo_permeability_score"] = None
    
    # 4. Toxicity risk assessment
    try:
        # Simple toxicity risk model based on structural features
        toxicity_risk = 0
        
        # Check for potentially toxic functional groups
        if properties.get("functional_groups", {}).get("nitro", 0) > 0:
            toxicity_risk += 20
        if properties.get("functional_groups", {}).get("alkyl_halide", 0) > 0:
            toxicity_risk += 15
        
        # High logP can indicate bioaccumulation potential
        if properties["logp"] > 5:
            toxicity_risk += 25
        
        # Lipinski violations can indicate poor drug-likeness
        if properties.get("permeability", {}).get("rule_of_5_violations", 0) > 2:
            toxicity_risk += 15
        
        cryo_properties["toxicity_risk"] = min(100, toxicity_risk)
    except (KeyError, TypeError) as e:
        logger.warning(f"Error calculating toxicity risk: {str(e)}")
        cryo_properties["toxicity_risk"] = None
    
    # 5. Overall cryoprotectant score (0-100)
    try:
        # Weight the different factors according to importance for cryoprotection
        vitrification_weight = 0.4
        permeability_weight = 0.3
        toxicity_weight = 0.3
        
        overall_score = (
            cryo_properties["vitrification_tendency"] * vitrification_weight +
            cryo_properties["cryo_permeability_score"] * permeability_weight +
            (100 - cryo_properties["toxicity_risk"]) * toxicity_weight
        )
        
        cryo_properties["overall_cryoprotectant_score"] = overall_score
    except (KeyError, TypeError) as e:
        logger.warning(f"Error calculating overall cryoprotectant score: {str(e)}")
        cryo_properties["overall_cryoprotectant_score"] = None
    
    # Add the cryoprotectant properties to the result
    properties["cryoprotectant_properties"] = cryo_properties
    
    return properties

def batch_calculate_properties(molecules: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """
    Calculate properties for multiple molecules in batch.
    
    Args:
        molecules: List of dictionaries, each containing 'data' and 'format' keys
        
    Returns:
        List of property dictionaries for each molecule
    """
    results = []
    
    for i, mol_info in enumerate(molecules):
        try:
            mol_data = mol_info.get('data', '')
            input_format = mol_info.get('format', 'smiles')
            
            # Calculate properties with caching
            properties = calculate_properties_with_cache(mol_data, input_format)
            
            # Add molecule index and original data for reference
            properties['molecule_index'] = i
            properties['original_data'] = mol_data
            
            results.append(properties)
        except Exception as e:
            logger.error(f"Error calculating properties for molecule {i}: {str(e)}")
            results.append({
                'molecule_index': i,
                'original_data': mol_info.get('data', ''),
                'error': str(e)
            })
    
    return results

def find_similar_molecules(query_mol_data: str, target_molecules: List[str], 
                          query_format: str = 'smiles', 
                          fingerprint_type: str = 'morgan',
                          similarity_threshold: float = 0.7) -> List[Dict[str, Any]]:
    """
    Find molecules similar to a query molecule from a list of target molecules.
    
    Args:
        query_mol_data: Query molecule data
        target_molecules: List of SMILES strings to search
        query_format: Format of the query data ('smiles', 'mol', 'sdf')
        fingerprint_type: Type of fingerprint to use ('morgan', 'maccs', 'topological')
        similarity_threshold: Minimum similarity score to include in results
        
    Returns:
        List of dictionaries with similarity results, sorted by similarity score
    """
    # Parse query molecule
    query_mol = parse_molecule(query_mol_data, query_format)
    if query_mol is None:
        return [{"error": f"Failed to parse query molecule from {query_format} data"}]
    
    # Get canonical SMILES for the query molecule
    query_smiles = Chem.MolToSmiles(query_mol)
    
    # Generate query fingerprint
    if fingerprint_type == 'morgan':
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
    elif fingerprint_type == 'maccs':
        query_fp = AllChem.GetMACCSKeysFingerprint(query_mol)
    elif fingerprint_type == 'topological':
        query_fp = FingerprintMols.FingerprintMol(query_mol)
    else:
        return [{"error": f"Unsupported fingerprint type: {fingerprint_type}"}]
    
    # Calculate similarity for each target molecule
    results = []
    
    for i, target_smiles in enumerate(target_molecules):
        try:
            target_mol = Chem.MolFromSmiles(target_smiles)
            if target_mol is None:
                results.append({
                    "molecule_index": i,
                    "smiles": target_smiles,
                    "error": "Failed to parse target molecule"
                })
                continue
            
            # Get canonical SMILES for the target molecule
            canonical_target = Chem.MolToSmiles(target_mol)
            
            # If the molecules are identical (same canonical SMILES), set similarity to 1.0
            if canonical_target == query_smiles:
                tanimoto = 1.0
                dice = 1.0
            else:
                # Generate target fingerprint
                if fingerprint_type == 'morgan':
                    target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, 2, nBits=2048)
                elif fingerprint_type == 'maccs':
                    target_fp = AllChem.GetMACCSKeysFingerprint(target_mol)
                elif fingerprint_type == 'topological':
                    target_fp = FingerprintMols.FingerprintMol(target_mol)
                
                # Calculate similarity
                tanimoto = DataStructs.TanimotoSimilarity(query_fp, target_fp)
                dice = DataStructs.DiceSimilarity(query_fp, target_fp)
            tanimoto = DataStructs.TanimotoSimilarity(query_fp, target_fp)
            dice = DataStructs.DiceSimilarity(query_fp, target_fp)
            
            # Only include results above threshold
            if tanimoto >= similarity_threshold:
                results.append({
                    "molecule_index": i,
                    "smiles": target_smiles,
                    "tanimoto": tanimoto,
                    "dice": dice,
                    "fingerprint_type": fingerprint_type
                })
        except Exception as e:
            logger.error(f"Error calculating similarity for molecule {i}: {str(e)}")
            results.append({
                "molecule_index": i,
                "smiles": target_smiles,
                "error": str(e)
            })
    
    # Sort results by similarity score (descending)
    results = sorted(results, key=lambda x: x.get("tanimoto", 0), reverse=True)
    
    return results

def batch_substructure_search(query_mol_data: str, target_molecules: List[str],
                             query_format: str = 'smarts') -> List[Dict[str, Any]]:
    """
    Perform substructure search on multiple target molecules.
    
    Args:
        query_mol_data: Query molecule data (usually SMARTS pattern)
        target_molecules: List of SMILES strings to search
        query_format: Format of the query data ('smarts', 'smiles')
        
    Returns:
        List of dictionaries with search results for each target molecule
    """
    # Parse query molecule
    if query_format.lower() == 'smarts':
        query_mol = Chem.MolFromSmarts(query_mol_data)
    else:
        query_mol = parse_molecule(query_mol_data, query_format)
    
    if query_mol is None:
        return [{"error": f"Failed to parse query molecule from {query_format} data"}]
    
    # Perform search on each target molecule
    results = []
    
    for i, target_smiles in enumerate(target_molecules):
        try:
            target_mol = Chem.MolFromSmiles(target_smiles)
            if target_mol is None:
                results.append({
                    "molecule_index": i,
                    "smiles": target_smiles,
                    "match": False,
                    "error": "Failed to parse target molecule"
                })
                continue
            
            # Perform substructure search
            matches = target_mol.GetSubstructMatches(query_mol)
            
            result = {
                "molecule_index": i,
                "smiles": target_smiles,
                "match": len(matches) > 0,
                "match_count": len(matches),
                "matches": [list(match) for match in matches]
            }
            
            # Generate visualization if there are matches
            if matches:
                # Generate SVG with highlighted matches
                all_match_atoms = set()
                for match in matches:
                    all_match_atoms.update(match)
                
                result["visualization"] = generate_molecule_svg(
                    target_smiles, 
                    'smiles', 
                    highlight_atoms=list(all_match_atoms)
                )
            
            results.append(result)
        except Exception as e:
            logger.error(f"Error performing substructure search for molecule {i}: {str(e)}")
            results.append({
                "molecule_index": i,
                "smiles": target_smiles,
                "match": False,
                "error": str(e)
            })
    
    return results

def generate_molecule_grid(mol_data_list: List[str], labels: List[str] = None,
                          mol_per_row: int = 3, sub_img_size: Tuple[int, int] = (200, 200),
                          input_format: str = 'smiles') -> str:
    """
    Generate a grid of molecule visualizations.
    
    Args:
        mol_data_list: List of molecular data strings
        labels: Optional list of labels for each molecule
        mol_per_row: Number of molecules per row in the grid
        sub_img_size: Size of each molecule image (width, height)
        input_format: Format of the input data ('smiles', 'mol', 'sdf')
        
    Returns:
        SVG string of the molecule grid
    """
    # Parse molecules
    mols = []
    for mol_data in mol_data_list:
        mol = parse_molecule(mol_data, input_format)
        if mol is not None:
            # Remove hydrogens for cleaner visualization
            mol = Chem.RemoveHs(mol)
            mols.append(mol)
        else:
            # Add a placeholder for failed molecules
            mols.append(None)
    
    # Generate grid image
    if not labels:
        labels = [f"Molecule {i+1}" for i in range(len(mols))]
    
    try:
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=mol_per_row,
            subImgSize=sub_img_size,
            legends=labels,
            useSVG=True
        )
        return img
    except Exception as e:
        logger.error(f"Error generating molecule grid: {str(e)}")
        return ""

def analyze_scaffold(mol_data: str, input_format: str = 'smiles') -> Dict[str, Any]:
    """
    Analyze the molecular scaffold and return scaffold information.
    
    Args:
        mol_data: Molecular data as a string
        input_format: Format of the input data ('smiles', 'mol', 'sdf')
        
    Returns:
        Dictionary with scaffold information
    """
    mol = parse_molecule(mol_data, input_format)
    if mol is None:
        return {"error": f"Failed to parse molecule from {input_format} data"}
    
    try:
        # Get Murcko scaffold
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smiles = Chem.MolToSmiles(scaffold) if scaffold else ""
        
        # Get framework (scaffold including ring substituents)
        framework = MurckoScaffold.MurckoScaffoldSmiles(
            Chem.MolToSmiles(mol),
            includeChirality=False
        )
        
        # Calculate fraction of atoms in scaffold
        if mol.GetNumHeavyAtoms() > 0:
            scaffold_fraction = scaffold.GetNumHeavyAtoms() / mol.GetNumHeavyAtoms()
        else:
            scaffold_fraction = 0.0
        
        # Generate visualizations
        mol_svg = generate_molecule_svg(mol_data, input_format)
        scaffold_svg = generate_molecule_svg(scaffold_smiles, 'smiles') if scaffold_smiles else ""
        
        return {
            "scaffold_smiles": scaffold_smiles,
            "framework_smiles": framework,
            "scaffold_fraction": scaffold_fraction,
            "molecule_svg": mol_svg,
            "scaffold_svg": scaffold_svg
        }
    except Exception as e:
        logger.error(f"Error analyzing scaffold: {str(e)}")
        return {"error": f"Error analyzing scaffold: {str(e)}"}

def clear_property_cache() -> Dict[str, Any]:
    """
    Clear the property cache.
    
    Returns:
        Status message
    """
    try:
        property_cache.clear()
        return {"status": "success", "message": "Property cache cleared"}
    except Exception as e:
        logger.error(f"Error clearing property cache: {str(e)}")
        return {"status": "error", "message": f"Error clearing property cache: {str(e)}"}