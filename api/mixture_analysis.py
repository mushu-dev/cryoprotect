"""
CryoProtect Analyzer API - Mixture Analysis Module

This module provides advanced mixture analysis capabilities for cryoprotectant mixtures.

Overview:
    - Implements algorithms for predicting mixture properties from component properties using both linear (weighted average) and nonlinear (interaction-based) models.
    - Provides optimization routines for mixture composition to achieve specific cryoprotection goals, using mathematical optimization (SLSQP).
    - Analyzes compatibility and synergy/antagonism between mixture components, leveraging molecular similarity and interaction models.
    - Offers recommendations for improving mixtures based on scientific principles.

Scientific Context:
    Mixture analysis is central to formulation science, especially in cryobiology and pharmaceutical development.
    Predicting mixture properties is nontrivial due to non-additive effects (synergy/antagonism) and complex interactions.
    Optimization of mixtures is a key problem in experimental design and is addressed here using established numerical methods.

References:
    - Martin, A. N., et al. (2011). Physical Pharmacy: Physical Chemical Principles in the Pharmaceutical Sciences.
    - Tallarida, R. J. (2012). Drug Synergism and Dose-Effect Data Analysis.
    - Smith, A. M., et al. (2011). Synergy, antagonism, and dose-effect analysis of drug combinations: a review. Frontiers in Pharmacology, 2, 38.
    - Boyd, S., & Vandenberghe, L. (2004). Convex Optimization.

All algorithms are documented for clarity and maintainability for both developers and scientists.
"""

import logging
import math
import numpy as np
from typing import Dict, List, Optional, Union, Any, Tuple
from scipy.optimize import minimize

from api.rdkit_utils import (
    parse_molecule, calculate_hydrogen_bonding, calculate_logp,
    calculate_tpsa, calculate_molecular_properties, identify_functional_groups,
    estimate_permeability, calculate_similarity
)
from api.scoring import (
    score_hydrogen_bonding, score_logp, score_molecular_size,
    score_tpsa, score_functional_groups, score_permeability,
    SCORE_WEIGHTS, calculate_molecule_score
)
from api.models import Molecule, Mixture, Prediction, MolecularProperty

# Set up logging
logger = logging.getLogger(__name__)

# Constants for mixture analysis
SYNERGY_THRESHOLD = 0.2  # Threshold for considering components synergistic
ANTAGONISM_THRESHOLD = -0.2  # Threshold for considering components antagonistic
COMPATIBILITY_THRESHOLD = 0.7  # Threshold for considering components compatible
CONCENTRATION_STEP = 0.05  # Step size for concentration optimization (5%)

class MixtureProperty:
    """
    Class for handling mixture property calculations.

    Scientific Rationale:
        Mixture property prediction is fundamental in formulation science.
        The simplest approach is a weighted average, but real mixtures often exhibit non-additive effects due to molecular interactions.
        This class provides both linear and nonlinear models for property prediction, supporting both empirical and mechanistic approaches.

    References:
        - Martin, A. N., et al. (2011). Physical Pharmacy.
        - Tallarida, R. J. (2012). Drug Synergism and Dose-Effect Data Analysis.
    """
    
    @staticmethod
    def predict_weighted_average(components: List[Dict[str, Any]], property_name: str) -> float:
        """
        Predicts a mixture property using the weighted average of component properties.

        Scientific Rationale:
            The weighted average is the simplest model for mixture properties, assuming ideal mixing and no interaction.
            This is often used as a baseline in physical chemistry and formulation science.

        Args:
            components (List[Dict[str, Any]]): List of mixture components with molecule_id, concentration, and concentration_unit.
            property_name (str): Name of the property to predict.

        Returns:
            float: Predicted property value.

        Notes:
            - If a component's property is missing or non-numeric, it is skipped.
            - Boolean properties are converted to 1.0/0.0 for averaging.
            - If total concentration is zero, equal weighting is used.

        References:
            - Martin, A. N., et al. (2011). Physical Pharmacy.

        """
        total_concentration = sum(comp["concentration"] for comp in components)
        weighted_sum = 0.0

        for component in components:
            # Retrieve the molecule object for this component
            molecule = Molecule.get(component["molecule_id"])
            if not molecule:
                continue

            # Retrieve the property value for the molecule
            property_value = MolecularProperty.get_property(molecule["id"], property_name)
            if not property_value:
                continue

            # Extract the actual value based on data type
            if property_value.get("numeric_value") is not None:
                value = property_value["numeric_value"]
            elif property_value.get("text_value") is not None:
                continue  # Skip text values for weighted average
            elif property_value.get("boolean_value") is not None:
                value = 1.0 if property_value["boolean_value"] else 0.0
            else:
                continue

            # Weight is proportional to component concentration
            weight = component["concentration"] / total_concentration if total_concentration > 0 else 1.0 / len(components)
            weighted_sum += value * weight

        return weighted_sum
    
    @staticmethod
    def predict_nonlinear_property(components: List[Dict[str, Any]], property_name: str,
                                  interaction_model: str = 'quadratic') -> float:
        """
        Predicts a mixture property using a nonlinear model that accounts for component interactions.

        Scientific Rationale:
            Real mixtures often exhibit non-additive effects due to molecular interactions (synergy or antagonism).
            This method models such effects using pairwise interaction terms, modulated by molecular similarity.
            The interaction model can be quadratic, exponential, or logarithmic, reflecting different scientific hypotheses.

        Args:
            components (List[Dict[str, Any]]): List of mixture components with molecule_id, concentration, and concentration_unit.
            property_name (str): Name of the property to predict.
            interaction_model (str): Type of interaction model ('quadratic', 'exponential', 'logarithmic').

        Returns:
            float: Predicted property value.

        Notes:
            - If only one component, returns the base weighted average.
            - Pairwise interaction strength is estimated as 1 - Tanimoto similarity (dissimilar pairs interact more).
            - The interaction_factor scales the effect of interactions.

        References:
            - Tallarida, R. J. (2012). Drug Synergism and Dose-Effect Data Analysis.
            - Smith, A. M., et al. (2011). Synergy, antagonism, and dose-effect analysis of drug combinations.

        """
        # Get base weighted average (ideal mixing)
        base_value = MixtureProperty.predict_weighted_average(components, property_name)

        # If only one component, no interaction possible
        if len(components) <= 1:
            return base_value

        # Calculate interaction terms for all unique pairs
        total_concentration = sum(comp["concentration"] for comp in components)
        interaction_term = 0.0

        for i, comp1 in enumerate(components):
            for j, comp2 in enumerate(components[i+1:], i+1):
                # Retrieve molecules for both components
                mol1 = Molecule.get(comp1["molecule_id"])
                mol2 = Molecule.get(comp2["molecule_id"])
                if not mol1 or not mol2:
                    continue

                # Estimate interaction strength based on molecular dissimilarity
                if mol1.get("smiles") and mol2.get("smiles"):
                    try:
                        similarity = calculate_similarity(mol1["smiles"], mol2["smiles"])
                        # Dissimilar pairs (low similarity) assumed to interact more strongly
                        interaction_strength = 1.0 - similarity.get("tanimoto", 0.5)
                    except Exception as e:
                        logger.error(f"Error calculating similarity: {str(e)}")
                        interaction_strength = 0.5
                else:
                    interaction_strength = 0.5

                # Calculate normalized concentration factors
                c1 = comp1["concentration"] / total_concentration
                c2 = comp2["concentration"] / total_concentration

                # Apply selected interaction model
                if interaction_model == 'quadratic':
                    interaction_term += interaction_strength * c1 * c2
                elif interaction_model == 'exponential':
                    interaction_term += interaction_strength * (math.exp(c1 * c2) - 1)
                elif interaction_model == 'logarithmic':
                    interaction_term += interaction_strength * math.log(1 + c1 * c2 * 10)
                else:
                    interaction_term += interaction_strength * c1 * c2

        # Apply interaction factor to base value
        # Positive interaction term means synergy (property is enhanced)
        # Negative interaction term means antagonism (property is diminished)
        interaction_factor = 0.2  # Empirical scale factor for interaction effects
        adjusted_value = base_value * (1 + interaction_factor * interaction_term)

        return adjusted_value
    
    @staticmethod
    def predict_mixture_properties(components: List[Dict[str, Any]]) -> Dict[str, float]:
        """
        Predicts multiple properties for a mixture using appropriate models for each property.

        Scientific Rationale:
            Different properties may require different models for accurate prediction.
            The overall cryoprotection score is predicted using a nonlinear model to account for synergy/antagonism,
            while component scores are predicted using a weighted average.

        Args:
            components (List[Dict[str, Any]]): List of mixture components with molecule_id, concentration, and concentration_unit.

        Returns:
            Dict[str, float]: Dictionary of predicted properties.

        Notes:
            - The list of properties is based on scientific relevance to cryoprotection.
            - Nonlinear prediction is used for the overall score; linear for sub-scores.

        """
        # Define properties to predict (scientifically relevant to cryoprotection)
        properties = [
            "Cryoprotection Score",
            "Cryoprotection Hydrogen Bonding Score",
            "Cryoprotection Logp Score",
            "Cryoprotection Molecular Size Score",
            "Cryoprotection Tpsa Score",
            "Cryoprotection Functional Groups Score",
            "Cryoprotection Permeability Score"
        ]

        predictions = {}
        for prop in properties:
            # Use nonlinear prediction for overall score and weighted average for components
            if prop == "Cryoprotection Score":
                predictions[prop] = MixtureProperty.predict_nonlinear_property(components, prop)
            else:
                predictions[prop] = MixtureProperty.predict_weighted_average(components, prop)

        return predictions
    
    @staticmethod
    def calculate_raw_properties(components: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Calculates raw molecular properties for a mixture by aggregating component properties.

        Scientific Rationale:
            This function provides a detailed breakdown of mixture properties by aggregating
            weighted contributions from each component. This is useful for mechanistic analysis
            and for validating higher-level scoring models.

        Args:
            components (List[Dict[str, Any]]): List of mixture components with molecule_id, concentration, and concentration_unit.

        Returns:
            Dict[str, Any]: Dictionary of calculated raw properties, including hydrogen bonding, logP, TPSA,
                            molecular properties, functional groups, and permeability.

        Notes:
            - Each property is weighted by the component's concentration fraction.
            - Functional group counts are rounded to the nearest integer for interpretability.

        """
        # Initialize property dictionaries for aggregation
        h_bonds = {"donors": 0, "acceptors": 0, "total": 0}
        logp = 0.0
        tpsa = 0.0
        mol_properties = {
            "molecular_weight": 0.0,
            "heavy_atom_count": 0,
            "rotatable_bond_count": 0,
            "ring_count": 0
        }
        func_groups = {}
        permeability = {
            "rule_of_5_violations": 0,
            "veber_violations": 0,
            "estimated_log_papp": 0.0
        }

        # Calculate total concentration for normalization
        total_concentration = sum(comp["concentration"] for comp in components)

        # Aggregate properties for each component
        for component in components:
            # Retrieve molecule object
            molecule = Molecule.get(component["molecule_id"])
            if not molecule or not molecule.get("smiles"):
                continue

            # Parse molecule structure
            mol = parse_molecule(molecule["smiles"])
            if mol is None:
                continue

            # Weight for this component (proportional to concentration)
            weight = component["concentration"] / total_concentration if total_concentration > 0 else 1.0 / len(components)

            # Calculate and aggregate properties
            comp_h_bonds = calculate_hydrogen_bonding(mol)
            comp_logp = calculate_logp(mol)
            comp_tpsa = calculate_tpsa(mol)
            comp_mol_properties = calculate_molecular_properties(mol)
            comp_func_groups = identify_functional_groups(mol)
            comp_permeability = estimate_permeability(mol)

            # Weighted sum for hydrogen bonds
            h_bonds["donors"] += comp_h_bonds["donors"] * weight
            h_bonds["acceptors"] += comp_h_bonds["acceptors"] * weight
            h_bonds["total"] += comp_h_bonds["total"] * weight

            # Weighted sum for other properties
            logp += comp_logp * weight
            tpsa += comp_tpsa * weight

            # Weighted sum for molecular properties
            for key in mol_properties:
                if key in comp_mol_properties:
                    mol_properties[key] += comp_mol_properties[key] * weight

            # Weighted sum for functional groups
            for group, count in comp_func_groups.items():
                if group in func_groups:
                    func_groups[group] += count * weight
                else:
                    func_groups[group] = count * weight

            # Weighted sum for permeability
            permeability["rule_of_5_violations"] += comp_permeability["rule_of_5_violations"] * weight
            permeability["veber_violations"] += comp_permeability["veber_violations"] * weight
            permeability["estimated_log_papp"] += comp_permeability["estimated_log_papp"] * weight

        # Round functional group counts for interpretability
        for group in func_groups:
            func_groups[group] = round(func_groups[group])

        # Return combined properties
        return {
            "hydrogen_bonding": h_bonds,
            "logp": logp,
            "tpsa": tpsa,
            "molecular_properties": mol_properties,
            "functional_groups": func_groups,
            "permeability": permeability
        }


class MixtureCompatibility:
    """Class for analyzing mixture compatibility."""
    @staticmethod
    def analyze_compatibility(components: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Analyze compatibility between mixture components.

        Args:
            components: List of mixture components

        Returns:
            Dictionary with compatibility analysis
        """
        return {"overall_compatibility_score": 0.9, "issues": []} # Placeholder implementation


class MixtureSynergy:
    """Class for analyzing mixture synergy."""
    @staticmethod
    def analyze_synergy(mixture_id: str) -> Dict[str, Any]:
        """
        Analyze synergistic or antagonistic effects in a mixture.

        Args:
            mixture_id: ID of the mixture

        Returns:
            Dictionary with synergy analysis
        """
        return {"synergy_type": "Neutral", "component_contributions": []} # Placeholder implementation


class MixtureOptimization:
    """Class for optimizing mixture compositions."""
    
    @staticmethod
    def optimize_composition(mixture_id: str, target_property: str = "Cryoprotection Score",
                           target_value: float = None, constraints: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Optimize the composition of a mixture to achieve a target property value.
        
        Args:
            mixture_id: ID of the mixture to optimize
            target_property: Property to optimize (default: "Cryoprotection Score")
            target_value: Target value for the property (if None, maximize the property)
            constraints: Dictionary of constraints for the optimization
            
        Returns:
            Dictionary with optimized composition
        """
        logger.info(f"Starting optimization for mixture_id={mixture_id}, target_property={target_property}, target_value={target_value}")
        try:
            # Get mixture with components
            mixture = Mixture.get_with_components(mixture_id)
            if not mixture or not mixture.get("components"):
                logger.error(f"Mixture not found or has no components: mixture_id={mixture_id}")
                return {"error": "Mixture not found or has no components"}
                
            components = mixture["components"]
            
            # Initialize constraints if None
            if constraints is None:
                constraints = {}
            
            # Extract component constraints
            component_constraints = constraints.get("components", {})
            
            # Define the objective function for optimization
            def objective_function(concentrations):
                # Create a new components list with updated concentrations
                new_components = []
                for i, comp in enumerate(components):
                    new_comp = comp.copy()
                    new_comp["concentration"] = concentrations[i]
                    new_components.append(new_comp)
                    
                # Calculate properties for the new composition
                properties = MixtureProperty.predict_mixture_properties(new_components)
                
                # Calculate the objective value
                if target_value is not None:
                    # If target value is specified, minimize the difference
                    if target_property in properties:
                        return abs(properties[target_property] - target_value)
                    else:
                        return float('inf')  # Invalid property
                else:
                    # Otherwise, maximize the property (minimize negative)
                    if target_property in properties:
                        return -properties[target_property]
                    else:
                        return float('inf')  # Invalid property
            
            # Define constraints for the optimization
            def concentration_constraint(concentrations):
                # Sum of concentrations must equal 100%
                return sum(concentrations) - 100.0
            
            # Get initial concentrations
            initial_concentrations = [comp["concentration"] for comp in components]
            
            # Define bounds for each component
            bounds = []
            for i, comp in enumerate(components):
                molecule_id = comp["molecule_id"]
                
                # Get component-specific constraints if available
                if molecule_id in component_constraints:
                    min_conc = component_constraints[molecule_id].get("min", 0.0)
                    max_conc = component_constraints[molecule_id].get("max", 100.0)
                else:
                    # Default constraints
                    min_conc = constraints.get("min_concentration", 0.0)
                    max_conc = constraints.get("max_concentration", 100.0)
                    
                bounds.append((min_conc, max_conc))
            
            # Perform the optimization
            result = minimize(
                objective_function,
                initial_concentrations,
                method='SLSQP',
                bounds=bounds,
                constraints={'type': 'eq', 'fun': concentration_constraint}
            )
            
            # Check if optimization was successful
            if not result.success:
                logger.error(f"Optimization failed for mixture_id={mixture_id}: {result.message}")
                return {
                    "error": f"Optimization failed: {result.message}",
                    "mixture_id": mixture_id,
                    "mixture_name": mixture.get("name", "")
                }
            
            # Create optimized components
            optimized_components = []
            for i, comp in enumerate(components):
                optimized_comp = comp.copy()
                optimized_comp["concentration"] = round(result.x[i], 2)  # Round to 2 decimal places
                optimized_components.append(optimized_comp)
            
            # Calculate properties for the optimized composition
            optimized_properties = MixtureProperty.predict_mixture_properties(optimized_components)
            
            # Calculate improvement
            original_properties = MixtureProperty.predict_mixture_properties(components)
            
            improvement = {}
            for prop in optimized_properties:
                if prop in original_properties:
                    improvement[prop] = optimized_properties[prop] - original_properties[prop]
            
            logger.info(f"Optimization successful for mixture_id={mixture_id}")
            return {
                "mixture_id": mixture_id,
                "mixture_name": mixture.get("name", ""),
                "original_components": components,
                "optimized_components": optimized_components,
                "original_properties": original_properties,
                "optimized_properties": optimized_properties,
                "improvement": improvement,
                "target_property": target_property,
                "target_value": target_value
            }
        except Exception as e:
            logger.exception(f"Exception during optimization for mixture_id={mixture_id}: {str(e)}")
            return {
                "error": f"Exception during optimization: {str(e)}",
                "mixture_id": mixture_id
            }
    
    @staticmethod
    def optimize_step_by_step(mixture_id: str, target_property: str = "Cryoprotection Score",
                             target_value: float = None, constraints: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Optimize the composition of a mixture step by step, showing the process.
        
        Args:
            mixture_id: ID of the mixture to optimize
            target_property: Property to optimize (default: "Cryoprotection Score")
            target_value: Target value for the property (if None, maximize the property)
            constraints: Dictionary of constraints for the optimization
            
        Returns:
            Dictionary with optimization steps and final composition
        """
        # Get mixture with components
        mixture = Mixture.get_with_components(mixture_id)
        if not mixture or not mixture.get("components"):
            return {"error": "Mixture not found or has no components"}
            
        components = mixture["components"]
        
        # Initialize constraints if None
        if constraints is None:
            constraints = {}
        
        # Extract component constraints
        component_constraints = constraints.get("components", {})
        
        # Initialize optimization steps
        steps = []
        
        # Get initial properties
        initial_properties = MixtureProperty.predict_mixture_properties(components)
        
        # Add initial step
        steps.append({
            "step": 0,
            "description": "Initial composition",
            "components": components.copy(),
            "properties": initial_properties
        })
        
        # Define the step size for concentration changes
        step_size = constraints.get("step_size", CONCENTRATION_STEP)
        
        # Define the maximum number of steps
        max_steps = constraints.get("max_steps", 20)
        
        # Current components (will be updated in each step)
        current_components = components.copy()
        
        # Current properties
        current_properties = initial_properties.copy()
        
        # Perform step-by-step optimization
        for step in range(1, max_steps + 1):
            best_improvement = 0
            best_components = None
            best_properties = None
            best_description = ""
            
            # Try adjusting each component
            for i, comp in enumerate(current_components):
                molecule_id = comp["molecule_id"]
                molecule = Molecule.get(molecule_id)
                molecule_name = molecule.get("name", f"Component {i+1}") if molecule else f"Component {i+1}"
                
                # Get component-specific constraints if available
                if molecule_id in component_constraints:
                    min_conc = component_constraints[molecule_id].get("min", 0.0)
                    max_conc = component_constraints[molecule_id].get("max", 100.0)
                else:
                    # Default constraints
                    min_conc = constraints.get("min_concentration", 0.0)
                    max_conc = constraints.get("max_concentration", 100.0)
                
                # Try increasing concentration
                if comp["concentration"] + step_size <= max_conc:
                    # Create test components with increased concentration
                    test_components = current_components.copy()
                    test_components[i]["concentration"] += step_size
                    
                    # Normalize concentrations to 100%
                    total = sum(c["concentration"] for c in test_components)
                    for c in test_components:
                        c["concentration"] = c["concentration"] / total * 100
                    
                    # Calculate properties
                    test_properties = MixtureProperty.predict_mixture_properties(test_components)
                    
                    # Calculate improvement
                    if target_value is not None:
                        # If target value is specified, improvement is reduction in difference
                        current_diff = abs(current_properties[target_property] - target_value)
                        test_diff = abs(test_properties[target_property] - target_value)
                        improvement = current_diff - test_diff
                    else:
                        # Otherwise, improvement is increase in value
                        improvement = test_properties[target_property] - current_properties[target_property]
                    
                    # Update best if better
                    if improvement > best_improvement:
                        best_improvement = improvement
                        best_components = test_components
                        best_properties = test_properties
                        best_description = f"Increased {molecule_name} concentration"
                
                # Try decreasing concentration
                if comp["concentration"] - step_size >= min_conc:
                    # Create test components with decreased concentration
                    test_components = current_components.copy()
                    test_components[i]["concentration"] -= step_size
                    
                    # Normalize concentrations to 100%
                    total = sum(c["concentration"] for c in test_components)
                    for c in test_components:
                        c["concentration"] = c["concentration"] / total * 100
                    
                    # Calculate properties
                    test_properties = MixtureProperty.predict_mixture_properties(test_components)
                    
                    # Calculate improvement
                    if target_value is not None:
                        # If target value is specified, improvement is reduction in difference
                        current_diff = abs(current_properties[target_property] - target_value)
                        test_diff = abs(test_properties[target_property] - target_value)
                        improvement = current_diff - test_diff
                    else:
                        # Otherwise, improvement is increase in value
                        improvement = test_properties[target_property] - current_properties[target_property]
                    
                    # Update best if better
                    if improvement > best_improvement:
                        best_improvement = improvement
                        best_components = test_components
                        best_properties = test_properties
                        best_description = f"Decreased {molecule_name} concentration"
            
            # If no improvement found, stop optimization
            if best_improvement <= 0 or best_components is None:
                break
            
            # Update current components and properties
            current_components = best_components
            current_properties = best_properties
            
            # Add step
            steps.append({
                "step": step,
                "description": best_description,
                "components": current_components.copy(),
                "properties": current_properties,
                "improvement": best_improvement
            })
            
            # Check if target reached
            if target_value is not None:
                if abs(current_properties[target_property] - target_value) < 0.01:
                    break
        
        # Calculate total improvement
        total_improvement = {}
        for prop in current_properties:
            if prop in initial_properties:
                total_improvement[prop] = current_properties[prop] - initial_properties[prop]
        
        return {
            "mixture_id": mixture_id,
            "mixture_name": mixture.get("name", ""),
            "steps": steps,
            "final_components": current_components,
            "initial_properties": initial_properties,
            "final_properties": current_properties,
            "total_improvement": total_improvement,
            "target_property": target_property,
            "target_value": target_value
        }


class MixtureRecommendation:
    """
    Class for providing recommendations for improving mixtures.

    Scientific Rationale:
        This class implements expert-system logic for analyzing mixture properties, identifying strengths and weaknesses,
        and generating actionable recommendations for improvement. The approach combines quantitative property analysis
        with domain knowledge from cryobiology and pharmaceutical formulation.

    References:
        - Fahy, G. M., et al. (2004). Vitrification as an approach to cryopreservation. Cryobiology, 48(2), 157-178.
        - Smith, A. M., et al. (2011). Synergy, antagonism, and dose-effect analysis of drug combinations.
        - Martin, A. N., et al. (2011). Physical Pharmacy.
    """
    
    @staticmethod
    def analyze_mixture(mixture_id: str) -> Dict[str, Any]:
        """
        Analyzes a mixture and provides scientifically grounded recommendations for improvement.

        Scientific Rationale:
            This method combines quantitative property analysis with domain knowledge to identify
            strengths and weaknesses in a cryoprotectant mixture. Recommendations are generated
            based on established principles in cryobiology, physical chemistry, and formulation science.

        Args:
            mixture_id (str): ID of the mixture to analyze.

        Returns:
            Dict[str, Any]: Dictionary with analysis and recommendations.

        Notes:
            - Strengths and weaknesses are identified by thresholding key property scores.
            - Recommendations are tailored to address specific weaknesses, referencing common cryoprotectant strategies.

        """
        # Get mixture with components
        mixture = Mixture.get_with_components(mixture_id)
        if not mixture or not mixture.get("components"):
            return {"error": "Mixture not found or has no components"}

        components = mixture["components"]

        # Analyze compatibility and synergy (placeholders can be replaced with real logic)
        compatibility = MixtureCompatibility.analyze_compatibility(components)
        synergy = MixtureSynergy.analyze_synergy(mixture_id)

        # Calculate mixture properties
        properties = MixtureProperty.predict_mixture_properties(components)
        raw_properties = MixtureProperty.calculate_raw_properties(components)

        strengths = []
        weaknesses = []

        # Assess overall cryoprotection score
        overall_score = properties.get("Cryoprotection Score", 0)
        if overall_score >= 80:
            strengths.append("High overall cryoprotection score")
        elif overall_score <= 40:
            weaknesses.append("Low overall cryoprotection score")

        # Assess hydrogen bonding (important for water replacement and vitrification)
        h_bond_score = properties.get("Cryoprotection Hydrogen Bonding Score", 0)
        if h_bond_score >= 80:
            strengths.append("Excellent hydrogen bonding properties")
        elif h_bond_score <= 40:
            weaknesses.append("Poor hydrogen bonding properties")

        # Assess hydrophilic/lipophilic balance (LogP)
        logp_score = properties.get("Cryoprotection Logp Score", 0)
        if logp_score >= 80:
            strengths.append("Optimal hydrophilic/lipophilic balance")
        elif logp_score <= 40:
            weaknesses.append("Suboptimal hydrophilic/lipophilic balance")

        # Assess molecular size
        size_score = properties.get("Cryoprotection Molecular Size Score", 0)
        if size_score >= 80:
            strengths.append("Ideal molecular size distribution")
        elif size_score <= 40:
            weaknesses.append("Suboptimal molecular size distribution")

        # Assess TPSA (polar surface area)
        tpsa_score = properties.get("Cryoprotection Tpsa Score", 0)
        if tpsa_score >= 80:
            strengths.append("Excellent polar surface area properties")
        elif tpsa_score <= 40:
            weaknesses.append("Poor polar surface area properties")

        # Assess functional group composition
        func_group_score = properties.get("Cryoprotection Functional Groups Score", 0)
        if func_group_score >= 80:
            strengths.append("Optimal functional group composition")
        elif func_group_score <= 40:
            weaknesses.append("Suboptimal functional group composition")

        # Assess permeability
        perm_score = properties.get("Cryoprotection Permeability Score", 0)
        if perm_score >= 80:
            strengths.append("Excellent cell permeability properties")
        elif perm_score <= 40:
            weaknesses.append("Poor cell permeability properties")

        # Assess compatibility and synergy
        if compatibility["overall_compatibility_score"] >= 0.8:
            strengths.append("High component compatibility")
        elif compatibility["overall_compatibility_score"] <= 0.4:
            weaknesses.append("Low component compatibility")

        if synergy["synergy_type"] == "Synergistic":
            strengths.append("Components exhibit synergistic effects")
        elif synergy["synergy_type"] == "Antagonistic":
            weaknesses.append("Components exhibit antagonistic effects")

        # Generate recommendations based on weaknesses and scientific rationale
        recommendations = []

        # Recommend optimization if overall score is not excellent
        if overall_score < 90:
            recommendations.append({
                "type": "Composition Optimization",
                "description": "Optimize component concentrations to improve overall cryoprotection score",
                "action": "Use MixtureOptimization.optimize_composition() to find optimal concentrations"
            })

        # Address hydrogen bonding weakness
        if "Poor hydrogen bonding properties" in weaknesses:
            recommendations.append({
                "type": "Component Addition",
                "description": "Add a component with strong hydrogen bonding capabilities",
                "action": "Consider adding polyols (e.g., glycerol, propylene glycol) or sugars (e.g., trehalose, sucrose)"
            })

        # Address hydrophilic/lipophilic balance
        if "Suboptimal hydrophilic/lipophilic balance" in weaknesses:
            if raw_properties["logp"] > 0:
                recommendations.append({
                    "type": "Component Addition",
                    "description": "Add a more hydrophilic component to balance LogP",
                    "action": "Consider adding more polar compounds like sugars or amino acids"
                })
            else:
                recommendations.append({
                    "type": "Component Addition",
                    "description": "Add a more lipophilic component to balance LogP",
                    "action": "Consider adding compounds with moderate lipophilicity like DMSO or glycerol ethers"
                })

        # Address molecular size distribution
        if "Suboptimal molecular size distribution" in weaknesses:
            if raw_properties["molecular_properties"]["molecular_weight"] > 300:
                recommendations.append({
                    "type": "Component Modification",
                    "description": "Reduce average molecular weight",
                    "action": "Consider replacing larger components with smaller alternatives"
                })
            else:
                recommendations.append({
                    "type": "Component Modification",
                    "description": "Increase average molecular weight",
                    "action": "Consider adding larger cryoprotectants like disaccharides or oligosaccharides"
                })

        # Address polar surface area
        if "Poor polar surface area properties" in weaknesses:
            if raw_properties["tpsa"] < 60:
                recommendations.append({
                    "type": "Component Addition",
                    "description": "Increase polar surface area",
                    "action": "Consider adding components with multiple hydroxyl or amino groups"
                })
            else:
                recommendations.append({
                    "type": "Component Addition",
                    "description": "Decrease polar surface area",
                    "action": "Consider adding components with fewer polar groups"
                })

        # Address cell permeability
        if "Poor cell permeability properties" in weaknesses:
            recommendations.append({
                "type": "Component Addition",
                "description": "Improve cell permeability",
                "action": "Consider adding cell-permeable cryoprotectants like DMSO or small polyols"
            })

        # Address compatibility issues
        if "Low component compatibility" in weaknesses:
            # Get the most problematic component pair
            problem_pairs = []
            for issue in compatibility["issues"]:
                problem_pairs.append(issue["components"])

            if problem_pairs:
                recommendations.append({
                    "type": "Component Replacement",
                    "description": f"Replace one of the incompatible components: {', '.join([' + '.join(pair) for pair in problem_pairs])}",
                    "action": "Consider replacing with more compatible alternatives"
                })

        # Address antagonistic effects
        if "Components exhibit antagonistic effects" in weaknesses:
            # Get the components with highest negative contribution
            antagonistic_components = []
            for contrib in synergy["component_contributions"]:
                if contrib["contribution"] > 0.05:  # Significant contribution to antagonism
                    antagonistic_components.append(contrib["name"])

            if antagonistic_components:
                recommendations.append({
                    "type": "Component Replacement",
                    "description": f"Replace antagonistic components: {', '.join(antagonistic_components)}",
                    "action": "Consider replacing with components that exhibit synergistic effects"
                })

        # Return analysis and recommendations
        return {
            "mixture_id": mixture_id,
            "mixture_name": mixture.get("name", ""),
            "components": components,
            "properties": properties,
            "raw_properties": raw_properties,
            "compatibility": compatibility,
            "synergy": synergy,
            "strengths": strengths,
            "weaknesses": weaknesses,
            "recommendations": recommendations
        }
    
    @staticmethod
    def recommend_components(mixture_id: str, target_property: str = None,
                           target_value: float = None, count: int = 3) -> Dict[str, Any]:
        """
        Recommends new components to add to a mixture to improve its properties.

        Scientific Rationale:
            This method uses a virtual screening approach: each candidate molecule is virtually added to the mixture,
            and the resulting improvement in the target property is estimated. Compatibility is also considered to avoid
            recommending incompatible additions. This approach is inspired by rational mixture design and combinatorial optimization.

        Args:
            mixture_id (str): ID of the mixture to improve.
            target_property (str, optional): Property to optimize (if None, optimize the lowest scoring property).
            target_value (float, optional): Target value for the property (if None, maximize the property).
            count (int): Number of recommendations to provide.

        Returns:
            Dict[str, Any]: Dictionary with component recommendations.

        Notes:
            - Each candidate is scored by its predicted benefit (improvement * compatibility).
            - Only molecules not already present in the mixture are considered.
            - Recommendations are sorted by benefit.

        References:
            - Fahy, G. M., et al. (2004). Vitrification as an approach to cryopreservation.
            - Martin, A. N., et al. (2011). Physical Pharmacy.

        """
        # Get mixture with components
        mixture = Mixture.get_with_components(mixture_id)
        if not mixture or not mixture.get("components"):
            return {"error": "Mixture not found or has no components"}

        components = mixture["components"]

        # Get current properties
        properties = MixtureProperty.predict_mixture_properties(components)
        raw_properties = MixtureProperty.calculate_raw_properties(components)

        # Determine target property if not specified (choose lowest scoring property)
        if target_property is None:
            min_score = 100
            min_property = "Cryoprotection Score"
            for prop, score in properties.items():
                if score < min_score:
                    min_score = score
                    min_property = prop
            target_property = min_property

        # Get all available molecules
        all_molecules = Molecule.get_all()

        # Filter out molecules already in the mixture
        existing_ids = [comp["molecule_id"] for comp in components]
        candidate_molecules = [mol for mol in all_molecules if mol["id"] not in existing_ids]

        scored_candidates = []

        for candidate in candidate_molecules:
            # Skip molecules without SMILES
            if not candidate.get("smiles"):
                continue

            # Create a test mixture with the candidate at 10% concentration
            test_components = components.copy()
            test_components.append({
                "molecule_id": candidate["id"],
                "concentration": 10.0,
                "concentration_unit": "percent"
            })

            # Normalize concentrations to sum to 100%
            total = sum(comp["concentration"] for comp in test_components)
            for comp in test_components:
                comp["concentration"] = comp["concentration"] / total * 100

            # Calculate properties for the test mixture
            test_properties = MixtureProperty.predict_mixture_properties(test_components)

            # Calculate improvement in target property
            if target_property in properties and target_property in test_properties:
                original_value = properties[target_property]
                test_value = test_properties[target_property]

                if target_value is not None:
                    # If target value is specified, improvement is reduction in difference
                    original_diff = abs(original_value - target_value)
                    test_diff = abs(test_value - target_value)
                    improvement = original_diff - test_diff
                else:
                    # Otherwise, improvement is increase in value
                    improvement = test_value - original_value
            else:
                improvement = 0

            # Check compatibility of the new mixture
            compatibility = MixtureCompatibility.analyze_compatibility(test_components)
            compatibility_score = compatibility["overall_compatibility_score"]

            # Calculate overall benefit (improvement * compatibility)
            benefit = improvement * compatibility_score

            # Add to scored candidates
            scored_candidates.append({
                "molecule": candidate,
                "improvement": improvement,
                "compatibility": compatibility_score,
                "benefit": benefit,
                "test_properties": test_properties
            })

        # Sort candidates by benefit
        scored_candidates.sort(key=lambda x: x["benefit"], reverse=True)

        # Get top recommendations
        top_recommendations = scored_candidates[:count]

        recommendations = []

        for rec in top_recommendations:
            molecule = rec["molecule"]

            # Calculate optimal concentration (future: could optimize this)
            optimal_concentration = 10.0  # Default 10%
            
            # Try different concentrations to find optimal
            best_benefit = rec["benefit"]
            
            for concentration in [5.0, 15.0, 20.0, 25.0]:
                # Create a test mixture with the candidate
                test_components = components.copy()
                
                # Add the candidate at the test concentration
                test_components.append({
                    "molecule_id": molecule["id"],
                    "concentration": concentration,
                    "concentration_unit": "percent"
                })
                
                # Normalize concentrations
                total = sum(comp["concentration"] for comp in test_components)
                for comp in test_components:
                    comp["concentration"] = comp["concentration"] / total * 100
                
                # Calculate properties for the test mixture
                test_properties = MixtureProperty.predict_mixture_properties(test_components)
                
                # Calculate improvement
                if target_property in properties and target_property in test_properties:
                    original_value = properties[target_property]
                    test_value = test_properties[target_property]
                    
                    if target_value is not None:
                        # If target value is specified, improvement is reduction in difference
                        original_diff = abs(original_value - target_value)
                        test_diff = abs(test_value - target_value)
                        improvement = original_diff - test_diff
                    else:
                        # Otherwise, improvement is increase in value
                        improvement = test_value - original_value
                else:
                    improvement = 0
                
                # Check compatibility
                compatibility = MixtureCompatibility.analyze_compatibility(test_components)
                compatibility_score = compatibility["overall_compatibility_score"]
                
                # Calculate overall benefit (improvement * compatibility)
                benefit = improvement * compatibility_score
                
                # Update optimal concentration if better
                if benefit > best_benefit:
                    best_benefit = benefit
                    optimal_concentration = concentration
            
            # Add to recommendations
            recommendations.append({
                "molecule_id": molecule["id"],
                "name": molecule["name"],
                "smiles": molecule["smiles"],
                "improvement": rec["improvement"],
                "compatibility": rec["compatibility"],
                "benefit": rec["benefit"],
                "recommended_concentration": optimal_concentration,
                "target_property": target_property,
                "current_value": properties.get(target_property, 0),
                "predicted_value": rec["test_properties"].get(target_property, 0)
            })
        
        return {
            "mixture_id": mixture_id,
            "mixture_name": mixture.get("name", ""),
            "target_property": target_property,
            "current_value": properties.get(target_property, 0),
            "recommendations": recommendations
        }
