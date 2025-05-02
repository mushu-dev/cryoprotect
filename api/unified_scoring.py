"""
CryoProtect Analyzer API - Unified Scoring Module

This module provides a comprehensive scoring system that combines:
1. Cryoprotection efficacy (from scoring.py)
2. Toxicity data (from Tox21)
3. Glass transition temperature (Tg) predictions

The unified scoring system provides a balanced assessment of cryoprotectant candidates,
considering both effectiveness and safety profiles.

Scientific Context:
    The ideal cryoprotectant should balance several key properties:
    - High cryoprotective efficacy (ability to prevent ice formation and protect cells)
    - Low toxicity (minimal adverse effects on biological systems)
    - Appropriate glass transition temperature (Tg) for the specific application
    - Good permeability and biocompatibility

    Different applications may require different balancing of these properties:
    - Cell preservation may prioritize efficacy and permeability
    - Organ preservation may prioritize low toxicity and optimal Tg
    - Long-term storage may prioritize glass-forming ability and stability

References:
    - Elliott, G. D., et al. (2017). Cryoprotectants: A review of the actions and applications
      of cryoprotective solutes that modulate cell recovery from ultra-low temperatures.
      Cryobiology, 76, 74-91.
    - Best, B. P. (2015). Cryoprotectant toxicity: facts, issues, and questions.
      Rejuvenation research, 18(5), 422-436.
    - Wowk, B. (2010). Thermodynamic aspects of vitrification. Cryobiology, 60(1), 11-22.
"""

import logging
import math
import numpy as np
from typing import Dict, List, Optional, Union, Any, Tuple
from datetime import datetime

from api.scoring import (
    calculate_molecule_score, calculate_mixture_score,
    normalize_score, SCORE_WEIGHTS
)
from api.predictive_models import (
    predict_vitrification_tendency, predict_cryoprotection_effectiveness,
    predict_cell_membrane_protection, model_manager
)
from api.models import Molecule, MolecularProperty, Mixture, Prediction
from chemical_data.toxicity.toxicity_scorer import ToxicityScorer

# Set up logging
logger = logging.getLogger(__name__)

# Application contexts with different weighting profiles
APPLICATION_CONTEXTS = {
    "general": {
        "efficacy": 0.40,
        "toxicity": 0.30,
        "tg": 0.30,
        "description": "Balanced profile suitable for general cryopreservation applications"
    },
    "cell_preservation": {
        "efficacy": 0.50,
        "toxicity": 0.25,
        "tg": 0.25,
        "description": "Prioritizes efficacy and permeability for cell preservation"
    },
    "organ_preservation": {
        "efficacy": 0.30,
        "toxicity": 0.40,
        "tg": 0.30,
        "description": "Prioritizes low toxicity and optimal glass transition for organ preservation"
    },
    "long_term_storage": {
        "efficacy": 0.30,
        "toxicity": 0.20,
        "tg": 0.50,
        "description": "Prioritizes glass-forming ability and stability for long-term storage"
    },
    "sensitive_tissues": {
        "efficacy": 0.30,
        "toxicity": 0.50,
        "tg": 0.20,
        "description": "Prioritizes safety for sensitive tissues like neural or reproductive cells"
    }
}

# Glass transition temperature (Tg) optimal ranges for different applications
TG_OPTIMAL_RANGES = {
    "general": (-120, -80),
    "cell_preservation": (-130, -100),
    "organ_preservation": (-110, -90),
    "long_term_storage": (-100, -70),
    "sensitive_tissues": (-120, -100)
}

class UnifiedScorer:
    """
    Class for calculating unified scores that combine efficacy, toxicity, and Tg.
    
    Scientific Rationale:
        This class implements a comprehensive scoring system that balances the three
        key aspects of cryoprotectant performance: efficacy, safety, and physical properties.
        
        The scoring approach is based on the concept that different cryopreservation
        applications have different requirements, and the ideal balance between efficacy,
        toxicity, and glass transition properties varies by context.
    """
    
    def __init__(self):
        """Initialize the UnifiedScorer."""
        self.toxicity_scorer = ToxicityScorer()
    
    def _score_tg_value(self, tg_value: float, application_context: str = "general") -> float:
        """
        Score a glass transition temperature (Tg) value based on the optimal range for an application.
        
        Scientific Rationale:
            Different cryopreservation applications have different optimal Tg ranges.
            This method scores a Tg value based on how close it is to the optimal range
            for the specified application context.
            
        Args:
            tg_value: Glass transition temperature in °C
            application_context: Application context for optimal Tg range
            
        Returns:
            Score between 0 and 100
        """
        # Get optimal range for the application
        optimal_range = TG_OPTIMAL_RANGES.get(application_context, TG_OPTIMAL_RANGES["general"])
        min_tg, max_tg = optimal_range
        
        # Score based on distance from optimal range
        if min_tg <= tg_value <= max_tg:
            # Within optimal range - high score
            # Higher score for values in the middle of the range
            normalized_pos = (tg_value - min_tg) / (max_tg - min_tg)
            # Parabolic function with maximum at the center of the range
            score = 90 + 10 * (1 - 4 * (normalized_pos - 0.5) ** 2)
        elif tg_value < min_tg:
            # Below optimal range
            # Score decreases as distance from range increases
            distance = min_tg - tg_value
            score = max(0, 90 - distance * 0.5)
        else:  # tg_value > max_tg
            # Above optimal range
            # Score decreases as distance from range increases
            distance = tg_value - max_tg
            score = max(0, 90 - distance * 0.5)
        
        return round(score)
    
    def _store_unified_score(
        self, 
        molecule_id: Optional[str], 
        mixture_id: Optional[str],
        score: float,
        application_context: str,
        details: Dict[str, Any]
    ) -> bool:
        """
        Store a unified score in the database.
        
        Args:
            molecule_id: ID of the molecule (None for mixtures)
            mixture_id: ID of the mixture (None for molecules)
            score: Unified score
            application_context: Application context used for scoring
            details: Score details
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Store as a prediction
            Prediction.add_prediction(
                molecule_id=molecule_id,
                mixture_id=mixture_id,
                property_name=f"Unified Score ({application_context})",
                value=score,
                confidence=0.8,
                method_name="Unified Scoring System",
                metadata={
                    "application_context": application_context,
                    "component_scores": {
                        "efficacy": details["component_scores"]["efficacy"]["score"],
                        "toxicity": details["component_scores"]["toxicity"]["score"],
                        "glass_transition": details["component_scores"]["glass_transition"]["score"]
                    },
                    "weights": {
                        "efficacy": details["component_scores"]["efficacy"]["weight"],
                        "toxicity": details["component_scores"]["toxicity"]["weight"],
                        "glass_transition": details["component_scores"]["glass_transition"]["weight"]
                    }
                }
            )
            
            return True
            
        except Exception as e:
            logger.error(f"Error storing unified score: {str(e)}")
            return False


# Create a global instance of the UnifiedScorer
unified_scorer = UnifiedScorer()


def calculate_unified_molecule_score(
    molecule_id: str, 
    application_context: str = "general",
    algorithm: str = "random_forest",
    recalculate: bool = False
) -> Dict[str, Any]:
    """
    Calculate a unified score for a molecule combining efficacy, toxicity, and Tg.
    
    Scientific Rationale:
        This function provides a comprehensive assessment of a cryoprotectant molecule
        by combining three key aspects:
        1. Efficacy (ability to prevent ice formation and protect cells)
        2. Safety (low toxicity profile based on Tox21 data)
        3. Physical properties (appropriate glass transition temperature)
        
        The weighting of these aspects is adjusted based on the application context,
        as different cryopreservation applications have different requirements.
    
    Args:
        molecule_id: ID of the molecule
        application_context: Application context for weighting ("general", "cell_preservation", etc.)
        algorithm: Algorithm to use for predictive models
        recalculate: Whether to recalculate scores even if they already exist
        
    Returns:
        Dictionary with unified score and component scores
        
    References:
        - Elliott, G. D., et al. (2017). Cryoprotectants: A review of the actions and applications
          of cryoprotective solutes that modulate cell recovery from ultra-low temperatures.
          Cryobiology, 76, 74-91.
        - Best, B. P. (2015). Cryoprotectant toxicity: facts, issues, and questions.
          Rejuvenation research, 18(5), 422-436.
    """
    return unified_scorer.calculate_unified_molecule_score(
        molecule_id, application_context, algorithm, recalculate
    )


def calculate_unified_mixture_score(
    mixture_id: str, 
    application_context: str = "general",
    algorithm: str = "random_forest",
    recalculate: bool = False
) -> Dict[str, Any]:
    """
    Calculate a unified score for a mixture combining efficacy, toxicity, and Tg.
    
    Scientific Rationale:
        This function provides a comprehensive assessment of a cryoprotectant mixture
        by combining three key aspects:
        1. Efficacy (ability to prevent ice formation and protect cells)
        2. Safety (low toxicity profile based on Tox21 data)
        3. Physical properties (appropriate glass transition temperature)
        
        The weighting of these aspects is adjusted based on the application context,
        as different cryopreservation applications have different requirements.
        
        Mixtures are evaluated based on their component properties, with consideration
        for synergistic effects between components.
    
    Args:
        mixture_id: ID of the mixture
        application_context: Application context for weighting ("general", "cell_preservation", etc.)
        algorithm: Algorithm to use for predictive models
        recalculate: Whether to recalculate scores even if they already exist
        
    Returns:
        Dictionary with unified score and component scores
        
    References:
        - Elliott, G. D., et al. (2017). Cryoprotectants: A review of the actions and applications
          of cryoprotective solutes that modulate cell recovery from ultra-low temperatures.
          Cryobiology, 76, 74-91.
        - Wowk, B. (2010). Thermodynamic aspects of vitrification. Cryobiology, 60(1), 11-22.
    """
    return unified_scorer.calculate_unified_mixture_score(
        mixture_id, application_context, algorithm, recalculate
    )


def get_available_application_contexts() -> Dict[str, Dict[str, Any]]:
    """
    Get available application contexts with their descriptions and weight profiles.
    
    Returns:
        Dictionary of application contexts with descriptions and weights
    """
    return {
        context: {
            "description": data["description"],
            "weights": {
                "efficacy": data["efficacy"],
                "toxicity": data["toxicity"],
                "glass_transition": data["tg"]
            },
            "tg_optimal_range": TG_OPTIMAL_RANGES.get(context, TG_OPTIMAL_RANGES["general"])
        }
        for context, data in APPLICATION_CONTEXTS.items()
    }


def batch_calculate_unified_scores(
    entity_ids: List[str],
    entity_type: str = "molecule",
    application_context: str = "general",
    algorithm: str = "random_forest",
    recalculate: bool = False
) -> Dict[str, Any]:
    """
    Calculate unified scores for multiple molecules or mixtures.
    
    Args:
        entity_ids: List of molecule or mixture IDs
        entity_type: Type of entities ("molecule" or "mixture")
        application_context: Application context for weighting
        algorithm: Algorithm to use for predictive models
        recalculate: Whether to recalculate scores even if they already exist
        
    Returns:
        Dictionary with results for each entity
    """
    results = {}
    errors = []
    
    for entity_id in entity_ids:
        try:
            if entity_type.lower() == "molecule":
                result = calculate_unified_molecule_score(
                    entity_id, application_context, algorithm, recalculate
                )
            elif entity_type.lower() == "mixture":
                result = calculate_unified_mixture_score(
                    entity_id, application_context, algorithm, recalculate
                )
            else:
                errors.append({
                    "entity_id": entity_id,
                    "error": f"Invalid entity type: {entity_type}"
                })
                continue
                
            if "error" in result:
                errors.append({
                    "entity_id": entity_id,
                    "error": result["error"]
                })
            else:
                results[entity_id] = result
                
        except Exception as e:
            logger.error(f"Error calculating unified score for {entity_type} {entity_id}: {str(e)}")
            errors.append({
                "entity_id": entity_id,
                "error": str(e)
            })
    
    return {
        "results": results,
        "errors": errors,
        "summary": {
            "total": len(entity_ids),
            "successful": len(results),
            "failed": len(errors),
            "application_context": application_context
        }
    }
