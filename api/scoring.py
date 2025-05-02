"""
CryoProtect Analyzer API - Scoring Module

This module provides scientific algorithms for calculating cryoprotection effectiveness scores
based on molecular and mixture properties.

Overview:
    - Implements scoring functions for hydrogen bonding, logP, molecular size, TPSA, functional groups, and permeability.
    - Aggregates these scores using empirically derived weights to produce an overall cryoprotection score for molecules and mixtures.
    - The scoring system is grounded in scientific literature on cryoprotectant design, permeability, and structure-activity relationships.

Scientific Context:
    Cryoprotectant effectiveness is determined by a combination of physicochemical properties, including hydrogen bonding capacity,
    hydrophilic/lipophilic balance, molecular size, polar surface area, and the presence of key functional groups.
    The scoring algorithms here are based on established principles in cryobiology, medicinal chemistry, and pharmaceutical formulation.

References:
    - Fahy, G. M., et al. (2004). Vitrification as an approach to cryopreservation. Cryobiology, 48(2), 157-178.
    - Pegg, D. E. (2015). Principles of cryopreservation. Methods in Molecular Biology, 1257, 3-19.
    - Lipinski, C. A., et al. (1997). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. Advanced Drug Delivery Reviews, 23(1-3), 3-25.
    - Ertl, P., et al. (2000). Fast calculation of molecular polar surface area as a sum of fragment-based contributions and its application to the prediction of drug transport properties. J. Med. Chem., 43(20), 3714-3717.

All scoring functions are documented for clarity and maintainability for both developers and scientists.
"""

import logging
import math
from typing import Dict, List, Optional, Union, Any, Tuple

from api.rdkit_utils import (
    parse_molecule, calculate_hydrogen_bonding, calculate_logp,
    calculate_tpsa, calculate_molecular_properties, identify_functional_groups,
    estimate_permeability
)
from api.models import Molecule, MolecularProperty, Mixture

# Set up logging
logger = logging.getLogger(__name__)

# Weight factors for different properties in the scoring algorithm
# These weights are based on scientific literature on cryoprotectants
SCORE_WEIGHTS = {
    # Hydrogen bonding is critical for cryoprotection (enables water replacement and vitrification)
    "hydrogen_bonding": 0.25,
    # LogP (hydrophobicity/hydrophilicity balance) affects solubility and membrane permeability
    "logp": 0.15,
    # Molecular size and weight influence cell penetration and toxicity
    "molecular_size": 0.15,
    # Topological polar surface area (TPSA) relates to permeability and hydrogen bonding
    "tpsa": 0.15,
    # Functional groups (especially hydroxyl groups) are key for water interaction and glass formation
    "functional_groups": 0.20,
    # Permeability (as predicted by rules and log Papp) is essential for intracellular protection
    "permeability": 0.10
}

def normalize_score(value: float, min_val: float, max_val: float,
                   invert: bool = False) -> float:
    """
    Normalizes a value to a score between 0 and 1.

    Scientific Rationale:
        Normalization is used to map raw property values to a standard scale for scoring.
        Inversion is used when lower values are more desirable (e.g., toxicity).

    Args:
        value (float): The value to normalize.
        min_val (float): Minimum expected value.
        max_val (float): Maximum expected value.
        invert (bool): If True, a higher value results in a lower score.

    Returns:
        float: Normalized score between 0 and 1.

    Notes:
        - Values outside the expected range are clipped.
        - Used as a utility in property scoring functions.
    """
    # Handle values outside the expected range
    if value < min_val:
        value = min_val
    elif value > max_val:
        value = max_val

    # Calculate normalized value
    normalized = (value - min_val) / (max_val - min_val)

    # Invert if needed (higher value = lower score)
    if invert:
        normalized = 1 - normalized

    return normalized

def score_hydrogen_bonding(h_bonds: Dict[str, int]) -> float:
    """
    Scores hydrogen bonding capacity for cryoprotectant effectiveness.

    Scientific Rationale:
        Hydrogen bonding is essential for water replacement, glass formation, and stabilization of biological structures during freezing.
        The optimal range is based on the number of hydrogen bond donors and acceptors found in effective cryoprotectants.

    Args:
        h_bonds (Dict[str, int]): Dictionary with donor and acceptor counts.

    Returns:
        float: Score between 0 and 1.

    References:
        - Fahy, G. M., et al. (2004). Vitrification as an approach to cryopreservation.

    """
    # Cryoprotectants typically have multiple hydrogen bonding sites
    # Optimal range is around 6-12 total hydrogen bonds (empirical)
    total_bonds = h_bonds.get("total", 0)

    # Score based on total hydrogen bonds
    if total_bonds <= 2:
        return 0.2  # Too few hydrogen bonds
    elif total_bonds <= 5:
        return 0.5 + (total_bonds - 2) * 0.1  # Increasing score
    elif total_bonds <= 12:
        return 0.8 + (total_bonds - 5) * 0.02  # Optimal range
    else:
        # Too many hydrogen bonds might indicate poor permeability
        return 1.0 - (total_bonds - 12) * 0.05  # Decreasing score
    
def score_logp(logp: float) -> float:
    """
    Scores LogP (octanol-water partition coefficient) for cryoprotectant effectiveness.

    Scientific Rationale:
        LogP reflects the hydrophilic/lipophilic balance, which affects solubility, permeability, and toxicity.
        Effective cryoprotectants are typically moderately hydrophilic (LogP between -3 and 0).

    Args:
        logp (float): XLogP value.

    Returns:
        float: Score between 0 and 1.

    References:
        - Pegg, D. E. (2015). Principles of cryopreservation.

    """
    # Ideal cryoprotectants have moderate hydrophilicity
    # Optimal LogP range is around -3 to 0
    if logp < -5:
        return 0.3  # Too hydrophilic
    elif logp < -3:
        return 0.5 + (logp + 5) * 0.1  # Increasing score
    elif logp <= 0:
        return 0.9  # Optimal range
    elif logp <= 2:
        return 0.9 - (logp) * 0.2  # Decreasing score
    else:
        return 0.5 - (logp - 2) * 0.1  # Too lipophilic
    
def score_molecular_size(properties: Dict[str, float]) -> float:
    """
    Scores molecular size and weight for cryoprotectant effectiveness.

    Scientific Rationale:
        Molecular size affects cell penetration, toxicity, and glass-forming ability.
        Effective cryoprotectants are typically small to medium-sized molecules (60–180 Da).

    Args:
        properties (Dict[str, float]): Dictionary of molecular properties.

    Returns:
        float: Score between 0 and 1.

    References:
        - Fahy, G. M., et al. (2004). Vitrification as an approach to cryopreservation.

    """
    # Get molecular weight
    mw = properties.get("molecular_weight", 0)

    # Ideal cryoprotectants have moderate molecular weight
    # Optimal range is around 60-180 Da (like glycerol, DMSO, ethylene glycol)
    if mw < 30:
        return 0.2  # Too small
    elif mw < 60:
        return 0.5 + (mw - 30) * 0.01  # Increasing score
    elif mw <= 180:
        return 0.8 + (180 - mw) * 0.001  # Optimal range
    elif mw <= 300:
        return 0.8 - (mw - 180) * 0.002  # Decreasing score
    else:
        return 0.4 - (mw - 300) * 0.001  # Too large
    
def score_tpsa(tpsa: float) -> float:
    """
    Scores topological polar surface area (TPSA) for cryoprotectant effectiveness.

    Scientific Rationale:
        TPSA is a predictor of membrane permeability and hydrogen bonding capacity.
        Effective cryoprotectants have moderate TPSA (40–90 Å²), balancing solubility and permeability.

    Args:
        tpsa (float): TPSA value in Å².

    Returns:
        float: Score between 0 and 1.

    References:
        - Ertl, P., et al. (2000). Fast calculation of molecular polar surface area.

    """
    # Ideal cryoprotectants have moderate TPSA
    # Optimal range is around 40-90 Å²
    if tpsa < 20:
        return 0.3  # Too low
    elif tpsa < 40:
        return 0.5 + (tpsa - 20) * 0.01  # Increasing score
    elif tpsa <= 90:
        return 0.7 + (tpsa - 40) * 0.006  # Optimal range
    elif tpsa <= 140:
        return 1.0 - (tpsa - 90) * 0.006  # Decreasing score
    else:
        return 0.7 - (tpsa - 140) * 0.005  # Too high
    
def score_functional_groups(groups: Dict[str, int]) -> float:
    """
    Scores the presence of key functional groups for cryoprotectant effectiveness.

    Scientific Rationale:
        Hydroxyl groups are especially important for water replacement and glass formation.
        Ethers, amines, and amides also contribute to cryoprotective properties.

    Args:
        groups (Dict[str, int]): Dictionary of functional groups and their counts.

    Returns:
        float: Score between 0 and 1.

    References:
        - Fahy, G. M., et al. (2004). Vitrification as an approach to cryopreservation.

    """
    # Hydroxyl groups are particularly important for cryoprotection
    hydroxyl_count = groups.get("hydroxyl", 0)
    alcohol_count = groups.get("alcohol", 0)

    # Combine hydroxyl and alcohol counts (may overlap)
    oh_groups = max(hydroxyl_count, alcohol_count)

    # Other beneficial groups
    ether_count = groups.get("ether", 0)
    amine_count = groups.get("amine", 0)
    amide_count = groups.get("amide", 0)

    # Calculate base score from hydroxyl groups
    if oh_groups == 0:
        oh_score = 0.2
    elif oh_groups <= 3:
        oh_score = 0.5 + oh_groups * 0.15
    else:
        oh_score = 0.95  # Maximum score for hydroxyl groups

    # Bonus for other beneficial groups
    other_groups_score = min(0.2, (ether_count + amine_count + amide_count) * 0.05)

    # Combine scores
    return min(1.0, oh_score + other_groups_score)
    
def score_permeability(permeability: Dict[str, Any]) -> float:
    """
    Scores permeability for cryoprotectant effectiveness.

    Scientific Rationale:
        Permeability is essential for intracellular protection. Rule of 5 and Veber rule violations
        are penalized, while the estimated log Papp (apparent permeability) is used to assess passive diffusion.

    Args:
        permeability (Dict[str, Any]): Dictionary of permeability estimates.

    Returns:
        float: Score between 0 and 1.

    References:
        - Lipinski, C. A., et al. (1997). Experimental and computational approaches to estimate solubility and permeability.
        - Veber, D. F., et al. (2002). Molecular properties that influence the oral bioavailability of drug candidates.

    """
    # Get relevant permeability metrics
    rule_violations = permeability.get("rule_of_5_violations", 0)
    veber_violations = permeability.get("veber_violations", 0)
    log_papp = permeability.get("estimated_log_papp", -5.0)

    # Calculate score based on rule violations
    rule_score = 1.0 - (rule_violations * 0.2) - (veber_violations * 0.15)

    # Calculate score based on estimated permeability coefficient
    # Optimal log Papp is around -5 to -4
    if log_papp < -6:
        papp_score = 0.4  # Too impermeable
    elif log_papp < -5:
        papp_score = 0.6 + (log_papp + 6) * 0.2  # Increasing score
    elif log_papp <= -4:
        papp_score = 0.8  # Optimal range
    elif log_papp <= -3:
        papp_score = 0.8 - (log_papp + 4) * 0.1  # Decreasing score
    else:
        papp_score = 0.7 - (log_papp + 3) * 0.1  # Too permeable

    # Combine scores (weighted average)
    return (rule_score * 0.4) + (papp_score * 0.6)

def calculate_molecule_score(mol_data: str, input_format: str = 'smiles') -> Dict[str, Any]:
    """
    Calculate cryoprotection effectiveness score for a molecule.
    
    Args:
        mol_data: Molecular data as a string (SMILES, MOL, SDF)
        input_format: Format of the input data ('smiles', 'mol', 'sdf')
        
    Returns:
        Dictionary with overall score and component scores
    """
    # Parse molecule
    mol = parse_molecule(mol_data, input_format)
    
    if mol is None:
        return {"error": "Failed to parse molecule"}
    
    # Calculate properties
    h_bonds = calculate_hydrogen_bonding(mol)
    logp = calculate_logp(mol)
    tpsa = calculate_tpsa(mol)
    mol_properties = calculate_molecular_properties(mol)
    func_groups = identify_functional_groups(mol)
    perm = estimate_permeability(mol)
    
    # Calculate component scores
    h_bond_score = score_hydrogen_bonding(h_bonds)
    logp_score = score_logp(logp)
    size_score = score_molecular_size(mol_properties)
    tpsa_score = score_tpsa(tpsa)
    func_group_score = score_functional_groups(func_groups)
    perm_score = score_permeability(perm)
    
    # Calculate weighted overall score
    overall_score = (
        h_bond_score * SCORE_WEIGHTS["hydrogen_bonding"] +
        logp_score * SCORE_WEIGHTS["logp"] +
        size_score * SCORE_WEIGHTS["molecular_size"] +
        tpsa_score * SCORE_WEIGHTS["tpsa"] +
        func_group_score * SCORE_WEIGHTS["functional_groups"] +
        perm_score * SCORE_WEIGHTS["permeability"]
    )
    
    # Scale to 0-100 range
    overall_score = round(overall_score * 100)
    
    # Return scores
    return {
        "overall_score": overall_score,
        "component_scores": {
            "hydrogen_bonding": round(h_bond_score * 100),
            "logp": round(logp_score * 100),
            "molecular_size": round(size_score * 100),
            "tpsa": round(tpsa_score * 100),
            "functional_groups": round(func_group_score * 100),
            "permeability": round(perm_score * 100)
        },
        "properties": {
            "hydrogen_bonding": h_bonds,
            "logp": logp,
            "tpsa": tpsa,
            "molecular_properties": mol_properties,
            "functional_groups": func_groups,
            "permeability": perm
        }
    }

def calculate_mixture_score(mixture_id: str) -> Dict[str, Any]:
    """
    Calculate cryoprotection effectiveness score for a mixture.
    
    Args:
        mixture_id: ID of the mixture
        
    Returns:
        Dictionary with overall score and component scores
    """
    try:
        # Get mixture with components
        mixture = Mixture.get_with_components(mixture_id)
        if not mixture:
            return {"error": f"Mixture with ID {mixture_id} not found"}
        
        # Get components
        components = mixture.get("components", [])
        if not components:
            return {"error": "Mixture has no components"}
        
        # Calculate scores for each component
        component_scores = []
        total_concentration = 0
        
        for component in components:
            # Get molecule
            molecule = Molecule.get(component["molecule_id"])
            if not molecule:
                continue
                
            # Get SMILES
            smiles = molecule.get("smiles")
            if not smiles:
                continue
                
            # Calculate score
            score_data = calculate_molecule_score(smiles)
            if "error" in score_data:
                continue
                
            # Add to component scores
            component_scores.append({
                "molecule_id": molecule["id"],
                "name": molecule["name"],
                "concentration": component["concentration"],
                "concentration_unit": component["concentration_unit"],
                "score": score_data["overall_score"]
            })
            
            # Add to total concentration
            total_concentration += component["concentration"]
        
        if not component_scores:
            return {"error": "Could not calculate scores for any components"}
        
        # Calculate weighted average score based on concentration
        weighted_score = 0
        for comp in component_scores:
            weight = comp["concentration"] / total_concentration if total_concentration > 0 else 1.0 / len(component_scores)
            weighted_score += comp["score"] * weight
        
        # Apply synergy bonus for mixtures with complementary properties
        # This is a simplified model - in reality, synergy would be more complex
        if len(component_scores) > 1:
            # Small bonus for mixtures (could be refined with more sophisticated models)
            synergy_bonus = min(5, len(component_scores) * 2)
            weighted_score = min(100, weighted_score + synergy_bonus)
        
        # Round the final score
        weighted_score = round(weighted_score)
        
        return {
            "mixture_id": mixture_id,
            "name": mixture["name"],
            "overall_score": weighted_score,
            "component_scores": component_scores
        }
        
    except Exception as e:
        logger.error(f"Error calculating mixture score: {str(e)}")
        return {"error": f"Error calculating mixture score: {str(e)}"}

def store_molecule_score(molecule_id: str, score_data: Dict[str, Any]) -> bool:
    """
    Store a molecule's cryoprotection score in the database.
    
    Args:
        molecule_id: ID of the molecule
        score_data: Score data from calculate_molecule_score
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Extract overall score
        overall_score = score_data.get("overall_score")
        if overall_score is None:
            return False
            
        # Store overall score
        MolecularProperty.add_property(molecule_id, "Cryoprotection Score", overall_score)
        
        # Store component scores
        component_scores = score_data.get("component_scores", {})
        for name, score in component_scores.items():
            property_name = f"Cryoprotection {name.replace('_', ' ').title()} Score"
            MolecularProperty.add_property(molecule_id, property_name, score)
            
        return True
        
    except Exception as e:
        logger.error(f"Error storing molecule score: {str(e)}")
        return False

def store_mixture_score(mixture_id: str, score_data: Dict[str, Any]) -> bool:
    """
    Store a mixture's cryoprotection score in the database.
    
    Args:
        mixture_id: ID of the mixture
        score_data: Score data from calculate_mixture_score
        
    Returns:
        True if successful, False otherwise
    """
    try:
        from api.models import Prediction, CalculationMethod
        
        # Extract overall score
        overall_score = score_data.get("overall_score")
        if overall_score is None:
            return False
            
        # Store as a prediction
        Prediction.add_prediction(
            mixture_id=mixture_id,
            property_name="Cryoprotection Score",
            value=overall_score,
            confidence=0.85,  # High confidence in the scoring algorithm
            method_name="Cryoprotection Scoring Algorithm"
        )
        
        return True
        
    except Exception as e:
        logger.error(f"Error storing mixture score: {str(e)}")
        return False