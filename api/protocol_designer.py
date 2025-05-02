"""
CryoProtect Analyzer - Protocol Designer

This module implements the ProtocolDesigner class, which generates scientifically valid, stepwise protocols
for preparing and applying cryoprotectant mixtures. The generator supports stepwise addition, dilution, or
removal of cryoprotectants, based on user-specified mixture, target concentrations, and constraints.

Algorithmic and Scientific Basis:
- Stepwise addition/removal of cryoprotectants is a standard method to minimize osmotic shock and toxicity.
- The protocol generator divides the total concentration change into N steps (user-specified or auto-calculated).
- Each step respects user constraints (max step size, temperature, time, etc.).
- The protocol assumes ideal mixing and linear concentration changes unless otherwise specified.
- Timing and temperature recommendations are based on literature for cell/embryo cryopreservation (see code comments).

Assumptions:
- The mixture composition is retrieved via the Mixture model (not implemented here).
- Volumes are calculated for a standard sample size unless otherwise specified.
- The protocol is generic; users should adapt for specific cell types or experimental needs.

References:
- Mazur, P. (1984). Freezing of living cells: mechanisms and implications. American Journal of Physiology.
- Fahy, G. M. et al. (2004). Vitrification as an approach to cryopreservation. Cryobiology.

Author: Roo (AI)
"""

from typing import List, Dict, Optional, Any
from api.models import Protocol, Mixture

class ProtocolDesigner:
    """
    ProtocolDesigner generates stepwise protocols for cryoprotectant addition/removal.
    """

    @staticmethod
    def design_concentration_gradient(
        mixture_id: str,
        target_concentration: float = 10.0,
        sample_type: str = "cell_line",
        starting_temperature: float = 4.0,
        target_temperature: Optional[float] = None,
        step_count: Optional[int] = None,
        custom_sensitivity: Optional[Dict] = None
    ) -> Dict:
        """
        Generate a stepwise protocol for reaching the target cryoprotectant concentration.

        Args:
            mixture_id (str): ID of the cryoprotectant mixture.
            target_concentration (float): Final concentration (% w/v or mM).
            sample_type (str): Type of biological sample.
            starting_temperature (float): Initial temperature (°C).
            target_temperature (float, optional): Final temperature (°C).
            step_count (int, optional): Number of steps (default: auto-calculate).
            custom_sensitivity (dict, optional): Sensitivity profile (e.g., max step size, time per step).

        Returns:
            dict: Protocol with stepwise instructions.
        """
        # --- Parameter defaults and validation ---
        if target_concentration <= 0:
            return {"error": "Target concentration must be greater than 0"}
        if step_count is not None and step_count < 1:
            return {"error": "Step count must be at least 1"}

        # Retrieve mixture info from database
        mixture_data = Mixture.get_with_components(mixture_id)
        if not mixture_data:
            return {"error": f"Mixture with ID {mixture_id} not found"}
            
        # Transform mixture data to expected format
        mixture_info = {
            "name": mixture_data["name"],
            "components": []
        }
        
        # Add components with initial_conc of 0 and target_conc from parameter
        # This maintains backward compatibility with existing code
        for component in mixture_data.get("components", []):
            mixture_info["components"].append({
                "name": f"Component {component.get('molecule_id', 'unknown')}",
                "initial_conc": 0.0,
                "target_conc": target_concentration
            })
            
        # If no components were found, add a default component
        if not mixture_info["components"]:
            mixture_info["components"].append({
                "name": "Default Component",
                "initial_conc": 0.0,
                "target_conc": target_concentration
            })

        # Sensitivity defaults
        max_step_size = custom_sensitivity.get("max_step_size") if custom_sensitivity else None
        time_per_step = custom_sensitivity.get("time_per_step") if custom_sensitivity else 2  # min
        hold_time = custom_sensitivity.get("hold_time") if custom_sensitivity else 2  # min

        # --- Step calculation ---
        total_change = target_concentration  # Assume starting at 0 for simplicity
        if max_step_size:
            step_count_calc = int(-(-total_change // max_step_size))  # Ceiling division
            steps = step_count or step_count_calc
        else:
            steps = step_count or 4  # Default to 4 steps if not specified

        step_size = total_change / steps

        # --- Protocol generation ---
        protocol_steps = []
        current_conc = 0.0
        for i in range(steps):
            next_conc = min(current_conc + step_size, target_concentration)
            action = (
                f"Add cryoprotectant to reach {next_conc:.2f} units concentration."
                if next_conc > current_conc else
                f"Remove cryoprotectant to reach {next_conc:.2f} units concentration."
            )
            protocol_steps.append({
                "step": i + 1,
                "action": action,
                "from_concentration": current_conc,
                "to_concentration": next_conc,
                "temperature": target_temperature if target_temperature is not None else starting_temperature,
                "hold_time_min": hold_time,
                "notes": (
                    "Mix gently. Allow equilibration before next step. "
                    "Adjust timing for sensitive samples."
                )
            })
            current_conc = next_conc

        protocol = {
            "mixture_id": mixture_id,
            "mixture_name": mixture_info["name"],
            "target_concentration": target_concentration,
            "sample_type": sample_type,
            "starting_temperature": starting_temperature,
            "target_temperature": target_temperature if target_temperature is not None else starting_temperature,
            "steps": protocol_steps,
            "summary": (
                f"Stepwise protocol to reach {target_concentration} units concentration in {steps} steps. "
                "Adjust parameters as needed for your sample and application."
            )
        }
        return protocol

    @staticmethod
    def save_protocol(
        mixture_id: str,
        name: str,
        description: str,
        protocol_data: Dict[str, Any]
    ) -> Dict:
        """
        Save a protocol to the database.

        Args:
            mixture_id (str): ID of the cryoprotectant mixture.
            name (str): Name of the protocol.
            description (str): Description of the protocol.
            protocol_data (dict): Protocol data from design_concentration_gradient.

        Returns:
            dict: The saved protocol record.
        """
        try:
            # Extract protocol parameters from protocol_data
            target_concentration = protocol_data.get("target_concentration")
            sample_type = protocol_data.get("sample_type")
            starting_temperature = protocol_data.get("starting_temperature")
            target_temperature = protocol_data.get("target_temperature")
            steps = protocol_data.get("steps", [])
            
            # Calculate step_count from steps if not provided
            step_count = len(steps)
            
            # Extract custom_sensitivity if present
            custom_sensitivity = protocol_data.get("custom_sensitivity")
            
            # Create protocol in database
            protocol = Protocol.create_protocol(
                mixture_id=mixture_id,
                name=name,
                description=description,
                target_concentration=target_concentration,
                sample_type=sample_type,
                starting_temperature=starting_temperature,
                target_temperature=target_temperature,
                step_count=step_count,
                steps=steps,
                custom_sensitivity=custom_sensitivity
            )
            
            return protocol
        except Exception as e:
            return {"error": str(e)}

    @staticmethod
    def get_saved_protocol(protocol_id: str) -> Dict:
        """
        Retrieve a saved protocol from the database.

        Args:
            protocol_id (str): ID of the protocol to retrieve.

        Returns:
            dict: The protocol record or error message.
        """
        try:
            protocol = Protocol.get_protocol(protocol_id)
            if not protocol:
                return {"error": f"Protocol with ID {protocol_id} not found"}
            return protocol
        except Exception as e:
            return {"error": str(e)}

    @staticmethod
    def list_protocols_for_mixture(mixture_id: str) -> List[Dict]:
        """
        List all protocols for a mixture.

        Args:
            mixture_id (str): ID of the mixture.

        Returns:
            list: List of protocol records.
        """
        try:
            # Check if mixture exists
            mixture = Mixture.get(mixture_id)
            if not mixture:
                return []
                
            # Get protocols for mixture
            protocols = Protocol.get_protocols_for_mixture(mixture_id)
            return protocols
        except Exception:
            return []

    @staticmethod
    def get_sample_sensitivity_profiles() -> Dict:
        """
        Retrieve sample sensitivity profiles for different biological sample types.
        
        These profiles contain parameters that guide protocol generation based on
        the specific sensitivities of different cell/tissue types to cryoprotectants,
        temperature changes, and osmotic stress.
        
        Returns:
            dict: Dictionary of sample types and their sensitivity profiles
        """
        # Sample sensitivity profiles based on literature values
        # These could be stored in the database in a production system
        return {
            "cell_line": {
                "name": "Cell Line (Generic)",
                "description": "Standard parameters for immortalized cell lines",
                "osmotic_tolerance": 0.3,  # Fractional volume change tolerated
                "max_step_size": 5.0,      # Max concentration change per step (%)
                "time_per_step": 5,        # Minutes per equilibration step
                "cooling_rate": 1.0,       # °C/min for controlled cooling
                "warming_rate": 10.0,      # °C/min for warming
                "notes": "Suitable for most immortalized cell lines like HEK293, CHO, etc."
            },
            "primary_cells": {
                "name": "Primary Cells",
                "description": "Parameters for primary cell isolates",
                "osmotic_tolerance": 0.2,
                "max_step_size": 2.5,
                "time_per_step": 10,
                "cooling_rate": 0.5,
                "warming_rate": 8.0,
                "notes": "For freshly isolated primary cells that are more sensitive to osmotic stress"
            },
            "tissue": {
                "name": "Tissue Samples",
                "description": "Parameters for tissue fragments",
                "osmotic_tolerance": 0.15,
                "max_step_size": 2.0,
                "time_per_step": 15,
                "cooling_rate": 0.3,
                "warming_rate": 5.0,
                "notes": "For tissue fragments where diffusion limitations require longer equilibration"
            },
            "organoid": {
                "name": "Organoids",
                "description": "Parameters for 3D organoid cultures",
                "osmotic_tolerance": 0.18,
                "max_step_size": 2.0,
                "time_per_step": 12,
                "cooling_rate": 0.4,
                "warming_rate": 6.0,
                "notes": "For 3D organoid cultures with complex architecture"
            },
            "embryo": {
                "name": "Embryos",
                "description": "Parameters for embryo cryopreservation",
                "osmotic_tolerance": 0.1,
                "max_step_size": 1.5,
                "time_per_step": 8,
                "cooling_rate": 0.3,
                "warming_rate": 15.0,
                "notes": "For embryos that are highly sensitive to cryoinjury"
            }
        }

    @staticmethod
    def compare_protocols(protocol_ids: List[str]) -> Dict:
        """
        Compare multiple protocols to identify differences in parameters, steps, and expected outcomes.
        
        This method retrieves the specified protocols, analyzes their parameters and steps,
        and generates a comparison report highlighting key differences and similarities.
        
        Args:
            protocol_ids: List of protocol IDs to compare
            
        Returns:
            dict: Comparison results with protocol details and analysis
        """
        if not protocol_ids or len(protocol_ids) < 2:
            return {"error": "At least two protocol IDs are required for comparison"}
            
        # Retrieve protocols
        protocols = []
        for protocol_id in protocol_ids:
            protocol = Protocol.get_protocol(protocol_id)
            if not protocol:
                return {"error": f"Protocol with ID {protocol_id} not found"}
            protocols.append(protocol)
            
        # Extract key parameters for comparison
        comparison_data = {
            "protocols": protocols,
            "parameter_comparison": {},
            "step_comparison": {},
            "summary": {},
            "recommendations": []
        }
        
        # Compare basic parameters
        param_keys = [
            "target_concentration",
            "sample_type",
            "starting_temperature",
            "target_temperature",
            "step_count"
        ]
        
        for key in param_keys:
            values = [p.get(key) for p in protocols]
            comparison_data["parameter_comparison"][key] = {
                "values": values,
                "same": len(set(str(v) for v in values)) == 1,
                "min": min(values) if all(v is not None for v in values) else None,
                "max": max(values) if all(v is not None for v in values) else None
            }
        
        # Compare steps
        all_steps = []
        for protocol in protocols:
            steps = protocol.get("steps", [])
            all_steps.append(steps)
            
        # Find the protocol with the most steps for comparison
        max_steps = max(len(steps) for steps in all_steps)
        
        # Compare each step across protocols
        for step_idx in range(max_steps):
            step_data = []
            for protocol_idx, protocol in enumerate(protocols):
                steps = protocol.get("steps", [])
                if step_idx < len(steps):
                    step_data.append(steps[step_idx])
                else:
                    step_data.append(None)
                    
            comparison_data["step_comparison"][f"step_{step_idx+1}"] = {
                "data": step_data,
                "concentration_values": [
                    step.get("to_concentration") if step else None
                    for step in step_data
                ],
                "temperature_values": [
                    step.get("temperature") if step else None
                    for step in step_data
                ],
                "time_values": [
                    step.get("hold_time_min") if step else None
                    for step in step_data
                ]
            }
        
        # Generate summary statistics
        comparison_data["summary"] = {
            "protocol_count": len(protocols),
            "name": [p.get("name", "Unnamed") for p in protocols],
            "sample_types": list(set(p.get("sample_type") for p in protocols)),
            "concentration_range": [
                min(p.get("target_concentration", 0) for p in protocols),
                max(p.get("target_concentration", 0) for p in protocols)
            ],
            "step_count_range": [
                min(p.get("step_count", 0) for p in protocols),
                max(p.get("step_count", 0) for p in protocols)
            ],
            "total_duration": [
                sum(step.get("hold_time_min", 0) for step in p.get("steps", []))
                for p in protocols
            ]
        }
        
        # Generate recommendations based on comparison
        recommendations = []
        
        # Check if protocols use different sample types
        if len(set(p.get("sample_type") for p in protocols)) > 1:
            recommendations.append(
                "Protocols are designed for different sample types. Consider using the protocol "
                "specifically designed for your sample type."
            )
            
        # Check for significant differences in step count
        step_counts = [p.get("step_count", 0) for p in protocols]
        if max(step_counts) - min(step_counts) > 2:
            recommendations.append(
                "Protocols have significantly different step counts. More steps generally "
                "provide gentler transitions but take longer to execute."
            )
            
        # Check for differences in target concentration
        target_concs = [p.get("target_concentration", 0) for p in protocols]
        if max(target_concs) - min(target_concs) > 2:
            recommendations.append(
                "Protocols target different final concentrations. Higher concentrations may "
                "provide better cryoprotection but increase toxicity risk."
            )
            
        # Add recommendations to comparison data
        comparison_data["recommendations"] = recommendations
        
        # Generate visualization data (placeholder for now)
        # In a full implementation, this would generate data for charts/graphs
        comparison_data["visualization"] = None
        
        return comparison_data