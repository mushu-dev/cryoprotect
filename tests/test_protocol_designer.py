"""
Tests for the Protocol Designer module.
"""

import json
import sys
import os
from datetime import datetime
import unittest
from unittest.mock import patch, MagicMock

# Add parent directory to path to import app modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tests.base_test_case import BaseTestCase, MockSupabaseBaseTestCase

from api.protocol_designer import ProtocolDesigner
from api.models import Mixture, Protocol

class MockMixture:
    """Mock Mixture class for testing."""
    
    @staticmethod
    def get(mixture_id):
        """Mock get method."""
        return {
            "id": mixture_id,
            "name": "Test Mixture",
            "description": "A test mixture for unit tests"
        }
    
    @staticmethod
    def get_with_components(mixture_id):
        """Mock get_with_components method."""
        return {
            "id": mixture_id,
            "name": "Test Mixture",
            "description": "A test mixture for unit tests",
            "components": [
                {
                    "molecule_id": "mol1",
                    "concentration": 70,
                    "concentration_unit": "%"
                },
                {
                    "molecule_id": "mol2",
                    "concentration": 30,
                    "concentration_unit": "%"
                }
            ]
        }

class MockProtocol:
    """Mock Protocol class for testing."""
    
    @staticmethod
    def create_protocol(mixture_id, name, description, target_concentration, sample_type,
                       starting_temperature, target_temperature=None, step_count=None,
                       steps=None, custom_sensitivity=None):
        """Mock create_protocol method."""
        return {
            "id": "protocol1",
            "mixture_id": mixture_id,
            "name": name,
            "description": description,
            "target_concentration": target_concentration,
            "sample_type": sample_type,
            "starting_temperature": starting_temperature,
            "target_temperature": target_temperature,
            "step_count": step_count or 4,
            "steps": steps or [],
            "custom_sensitivity": custom_sensitivity,
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat()
        }
    
    @staticmethod
    def get_protocol(protocol_id):
        """Mock get_protocol method."""
        if protocol_id == "protocol1":
            return {
                "id": protocol_id,
                "mixture_id": "mix1",
                "mixture_name": "Test Mixture",
                "name": "Test Protocol",
                "description": "A test protocol",
                "target_concentration": 10.0,
                "sample_type": "cell_line",
                "starting_temperature": 4.0,
                "target_temperature": 0.0,
                "step_count": 4,
                "steps": [
                    {
                        "step": 1,
                        "action": "Add cryoprotectant to reach 2.50 units concentration.",
                        "from_concentration": 0.0,
                        "to_concentration": 2.5,
                        "temperature": 4.0,
                        "hold_time_min": 2,
                        "notes": "Mix gently. Allow equilibration before next step."
                    },
                    {
                        "step": 2,
                        "action": "Add cryoprotectant to reach 5.00 units concentration.",
                        "from_concentration": 2.5,
                        "to_concentration": 5.0,
                        "temperature": 4.0,
                        "hold_time_min": 2,
                        "notes": "Mix gently. Allow equilibration before next step."
                    },
                    {
                        "step": 3,
                        "action": "Add cryoprotectant to reach 7.50 units concentration.",
                        "from_concentration": 5.0,
                        "to_concentration": 7.5,
                        "temperature": 4.0,
                        "hold_time_min": 2,
                        "notes": "Mix gently. Allow equilibration before next step."
                    },
                    {
                        "step": 4,
                        "action": "Add cryoprotectant to reach 10.00 units concentration.",
                        "from_concentration": 7.5,
                        "to_concentration": 10.0,
                        "temperature": 4.0,
                        "hold_time_min": 2,
                        "notes": "Mix gently. Allow equilibration before next step."
                    }
                ],
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat()
            }
        elif protocol_id == "protocol2":
            return {
                "id": protocol_id,
                "mixture_id": "mix1",
                "mixture_name": "Test Mixture",
                "name": "Test Protocol 2",
                "description": "A second test protocol",
                "target_concentration": 15.0,
                "sample_type": "primary_cells",
                "starting_temperature": 4.0,
                "target_temperature": -1.0,
                "step_count": 6,
                "steps": [
                    {
                        "step": 1,
                        "action": "Add cryoprotectant to reach 2.50 units concentration.",
                        "from_concentration": 0.0,
                        "to_concentration": 2.5,
                        "temperature": 4.0,
                        "hold_time_min": 5,
                        "notes": "Mix gently. Allow equilibration before next step."
                    },
                    {
                        "step": 2,
                        "action": "Add cryoprotectant to reach 5.00 units concentration.",
                        "from_concentration": 2.5,
                        "to_concentration": 5.0,
                        "temperature": 4.0,
                        "hold_time_min": 5,
                        "notes": "Mix gently. Allow equilibration before next step."
                    },
                    {
                        "step": 3,
                        "action": "Add cryoprotectant to reach 7.50 units concentration.",
                        "from_concentration": 5.0,
                        "to_concentration": 7.5,
                        "temperature": 4.0,
                        "hold_time_min": 5,
                        "notes": "Mix gently. Allow equilibration before next step."
                    },
                    {
                        "step": 4,
                        "action": "Add cryoprotectant to reach 10.00 units concentration.",
                        "from_concentration": 7.5,
                        "to_concentration": 10.0,
                        "temperature": 4.0,
                        "hold_time_min": 5,
                        "notes": "Mix gently. Allow equilibration before next step."
                    },
                    {
                        "step": 5,
                        "action": "Add cryoprotectant to reach 12.50 units concentration.",
                        "from_concentration": 10.0,
                        "to_concentration": 12.5,
                        "temperature": 2.0,
                        "hold_time_min": 5,
                        "notes": "Mix gently. Allow equilibration before next step."
                    },
                    {
                        "step": 6,
                        "action": "Add cryoprotectant to reach 15.00 units concentration.",
                        "from_concentration": 12.5,
                        "to_concentration": 15.0,
                        "temperature": 0.0,
                        "hold_time_min": 5,
                        "notes": "Mix gently. Allow equilibration before next step."
                    }
                ],
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat()
            }
        return None
    
    @staticmethod
    def get_protocols_for_mixture(mixture_id):
        """Mock get_protocols_for_mixture method."""
        if mixture_id == "mix1":
            return [
                {
                    "id": "protocol1",
                    "mixture_id": mixture_id,
                    "name": "Test Protocol",
                    "description": "A test protocol",
                    "target_concentration": 10.0,
                    "sample_type": "cell_line",
                    "starting_temperature": 4.0,
                    "target_temperature": 0.0,
                    "step_count": 4,
                    "created_at": datetime.now().isoformat(),
                    "updated_at": datetime.now().isoformat()
                },
                {
                    "id": "protocol2",
                    "mixture_id": mixture_id,
                    "name": "Test Protocol 2",
                    "description": "A second test protocol",
                    "target_concentration": 15.0,
                    "sample_type": "primary_cells",
                    "starting_temperature": 4.0,
                    "target_temperature": -1.0,
                    "step_count": 6,
                    "created_at": datetime.now().isoformat(),
                    "updated_at": datetime.now().isoformat()
                }
            ]
        return []

class TestProtocolDesigner(MockSupabaseBaseTestCase):
    """Test cases for the Protocol Designer module."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Mock the Mixture and Protocol classes
        self.original_mixture_get = Mixture.get
        self.original_mixture_get_with_components = Mixture.get_with_components
        self.original_protocol_create_protocol = Protocol.create_protocol
        self.original_protocol_get_protocol = Protocol.get_protocol
        self.original_protocol_get_protocols_for_mixture = Protocol.get_protocols_for_mixture
        
        Mixture.get = MockMixture.get
        Mixture.get_with_components = MockMixture.get_with_components
        Protocol.create_protocol = MockProtocol.create_protocol
        Protocol.get_protocol = MockProtocol.get_protocol
        Protocol.get_protocols_for_mixture = MockProtocol.get_protocols_for_mixture
    
    def tearDown(self):
        """Tear down test fixtures."""
        # Restore the original methods
        Mixture.get = self.original_mixture_get
        Mixture.get_with_components = self.original_mixture_get_with_components
        Protocol.create_protocol = self.original_protocol_create_protocol
        Protocol.get_protocol = self.original_protocol_get_protocol
        Protocol.get_protocols_for_mixture = self.original_protocol_get_protocols_for_mixture
    
    def test_protocol_designer_basic(self):
        """Test the basic functionality of the ProtocolDesigner class."""
        # Test design_concentration_gradient method
        protocol = ProtocolDesigner.design_concentration_gradient(
            mixture_id="mix1",
            target_concentration=10.0,
            sample_type="cell_line",
            starting_temperature=4.0,
            target_temperature=0.0,
            step_count=4
        )
        
        self.assertEqual(protocol["mixture_id"], "mix1")
        self.assertEqual(protocol["target_concentration"], 10.0)
        self.assertEqual(protocol["sample_type"], "cell_line")
        self.assertEqual(protocol["starting_temperature"], 4.0)
        self.assertEqual(protocol["target_temperature"], 0.0)
        self.assertEqual(len(protocol["steps"]), 4)
        
        # Test save_protocol method
        saved_protocol = ProtocolDesigner.save_protocol(
            mixture_id="mix1",
            name="Test Protocol",
            description="A test protocol",
            protocol_data=protocol
        )
        
        self.assertEqual(saved_protocol["id"], "protocol1")
        self.assertEqual(saved_protocol["mixture_id"], "mix1")
        self.assertEqual(saved_protocol["name"], "Test Protocol")
        self.assertEqual(saved_protocol["description"], "A test protocol")
        
        # Test get_saved_protocol method
        retrieved_protocol = ProtocolDesigner.get_saved_protocol("protocol1")
        self.assertEqual(retrieved_protocol["id"], "protocol1")
        self.assertEqual(retrieved_protocol["mixture_id"], "mix1")
        self.assertEqual(retrieved_protocol["mixture_name"], "Test Mixture")
        self.assertEqual(retrieved_protocol["name"], "Test Protocol")
        
        # Test list_protocols_for_mixture method
        protocols = ProtocolDesigner.list_protocols_for_mixture("mix1")
        self.assertEqual(len(protocols), 2)
        self.assertEqual(protocols[0]["id"], "protocol1")
        self.assertEqual(protocols[0]["mixture_id"], "mix1")
        self.assertEqual(protocols[0]["name"], "Test Protocol")
        self.assertEqual(protocols[1]["id"], "protocol2")
        self.assertEqual(protocols[1]["name"], "Test Protocol 2")

    def test_design_concentration_gradient_with_custom_sensitivity(self):
        """Test the design_concentration_gradient method with custom sensitivity profile."""
        custom_sensitivity = {
            "max_step_size": 2.0,
            "time_per_step": 5,
            "hold_time": 3
        }
        
        protocol = ProtocolDesigner.design_concentration_gradient(
            mixture_id="mix1",
            target_concentration=10.0,
            sample_type="primary_cells",
            starting_temperature=4.0,
            target_temperature=0.0,
            custom_sensitivity=custom_sensitivity
        )
        
        # With max_step_size of 2.0 and target of 10.0, we should get 5 steps
        self.assertEqual(len(protocol["steps"]), 5)
        
        # Check that hold_time is applied from custom sensitivity
        self.assertEqual(protocol["steps"][0]["hold_time_min"], 3)
        
        # Check that steps respect max_step_size
        for step in protocol["steps"]:
            concentration_change = step["to_concentration"] - step["from_concentration"]
            self.assertLessEqual(concentration_change, 2.0)

    def test_design_concentration_gradient_edge_cases(self):
        """Test edge cases for the design_concentration_gradient method."""
        # Test with invalid target concentration
        protocol = ProtocolDesigner.design_concentration_gradient(
            mixture_id="mix1",
            target_concentration=0.0,
            sample_type="cell_line",
            starting_temperature=4.0
        )
        self.assertIn("error", protocol)
        self.assertEqual(protocol["error"], "Target concentration must be greater than 0")
        
        # Test with invalid step count
        protocol = ProtocolDesigner.design_concentration_gradient(
            mixture_id="mix1",
            target_concentration=10.0,
            sample_type="cell_line",
            starting_temperature=4.0,
            step_count=0
        )
        self.assertIn("error", protocol)
        self.assertEqual(protocol["error"], "Step count must be at least 1")
        
        # Test with non-existent mixture
        with patch.object(Mixture, 'get_with_components', return_value=None):
            protocol = ProtocolDesigner.design_concentration_gradient(
                mixture_id="nonexistent",
                target_concentration=10.0,
                sample_type="cell_line",
                starting_temperature=4.0
            )
            self.assertIn("error", protocol)
            self.assertEqual(protocol["error"], "Mixture with ID nonexistent not found")

    def test_get_sample_sensitivity_profiles(self):
        """Test the get_sample_sensitivity_profiles method."""
        profiles = ProtocolDesigner.get_sample_sensitivity_profiles()
        
        # Check that all expected profiles are present
        expected_profiles = ["cell_line", "primary_cells", "tissue", "organoid", "embryo"]
        for profile_name in expected_profiles:
            self.assertIn(profile_name, profiles)
            
        # Check that profiles have the expected structure
        for profile_name, profile in profiles.items():
            self.assertIn("name", profile)
            self.assertIn("description", profile)
            self.assertIn("osmotic_tolerance", profile)
            self.assertIn("max_step_size", profile)
            self.assertIn("time_per_step", profile)
            self.assertIn("cooling_rate", profile)
            self.assertIn("warming_rate", profile)
            self.assertIn("notes", profile)
            
        # Check specific values for a profile
        cell_line_profile = profiles["cell_line"]
        self.assertEqual(cell_line_profile["name"], "Cell Line (Generic)")
        self.assertEqual(cell_line_profile["osmotic_tolerance"], 0.3)
        self.assertEqual(cell_line_profile["max_step_size"], 5.0)
        self.assertEqual(cell_line_profile["time_per_step"], 5)

    def test_compare_protocols(self):
        """Test the compare_protocols method."""
        # Test with valid protocol IDs
        comparison = ProtocolDesigner.compare_protocols(["protocol1", "protocol2"])
        
        # Check that the comparison has the expected structure
        self.assertIn("protocols", comparison)
        self.assertIn("parameter_comparison", comparison)
        self.assertIn("step_comparison", comparison)
        self.assertIn("summary", comparison)
        self.assertIn("recommendations", comparison)
        
        # Check that parameter comparison includes all expected parameters
        param_keys = [
            "target_concentration",
            "sample_type",
            "starting_temperature",
            "target_temperature",
            "step_count"
        ]
        for key in param_keys:
            self.assertIn(key, comparison["parameter_comparison"])
            
        # Check that step comparison includes all steps
        self.assertIn("step_1", comparison["step_comparison"])
        self.assertIn("step_2", comparison["step_comparison"])
        self.assertIn("step_3", comparison["step_comparison"])
        self.assertIn("step_4", comparison["step_comparison"])
        
        # Check that summary includes expected fields
        self.assertIn("protocol_count", comparison["summary"])
        self.assertIn("sample_types", comparison["summary"])
        self.assertIn("concentration_range", comparison["summary"])
        self.assertIn("step_count_range", comparison["summary"])
        self.assertIn("total_duration", comparison["summary"])
        
        # Check that recommendations are generated
        self.assertTrue(len(comparison["recommendations"]) > 0)
        
        # Test with invalid protocol IDs
        with self.assertRaises(Exception):
            ProtocolDesigner.compare_protocols(["nonexistent"])
            
        # Test with insufficient protocol IDs
        result = ProtocolDesigner.compare_protocols(["protocol1"])
        self.assertIn("error", result)
        self.assertEqual(result["error"], "At least two protocol IDs are required for comparison")

class TestProtocolStep(unittest.TestCase):
    """Test cases for the Protocol.ProtocolStep class."""
    
    def test_validate_step(self):
        """Test the validate_step method."""
        # Valid step data
        valid_step = {
            "step": 1,
            "action": "Add cryoprotectant to reach 2.50 units concentration.",
            "from_concentration": 0.0,
            "to_concentration": 2.5,
            "temperature": 4.0,
            "hold_time_min": 2,
            "notes": "Mix gently. Allow equilibration before next step."
        }
        
        # Validate should return the same data for valid step
        result = Protocol.ProtocolStep.validate_step(valid_step)
        self.assertEqual(result, valid_step)
        
        # Test with missing required field
        invalid_step = valid_step.copy()
        del invalid_step["action"]
        with self.assertRaises(ValueError) as context:
            Protocol.ProtocolStep.validate_step(invalid_step)
        self.assertIn("Missing required field: action", str(context.exception))
        
        # Test with negative from_concentration
        invalid_step = valid_step.copy()
        invalid_step["from_concentration"] = -1.0
        with self.assertRaises(ValueError) as context:
            Protocol.ProtocolStep.validate_step(invalid_step)
        self.assertIn("From concentration must be non-negative", str(context.exception))
        
        # Test with negative to_concentration
        invalid_step = valid_step.copy()
        invalid_step["to_concentration"] = -1.0
        with self.assertRaises(ValueError) as context:
            Protocol.ProtocolStep.validate_step(invalid_step)
        self.assertIn("To concentration must be non-negative", str(context.exception))
        
        # Test with non-positive hold_time
        invalid_step = valid_step.copy()
        invalid_step["hold_time_min"] = 0
        with self.assertRaises(ValueError) as context:
            Protocol.ProtocolStep.validate_step(invalid_step)
        self.assertIn("Hold time must be positive", str(context.exception))

    def test_calculate_step_parameters(self):
        """Test the calculate_step_parameters method."""
        # Test with default sample type
        params = Protocol.ProtocolStep.calculate_step_parameters(
            from_concentration=0.0,
            to_concentration=5.0,
            temperature=4.0
        )
        
        self.assertIn("hold_time_min", params)
        self.assertIn("notes", params)
        
        # Test with different sample types
        sample_types = ["cell_line", "primary_cells", "tissue", "organoid", "embryo"]
        for sample_type in sample_types:
            params = Protocol.ProtocolStep.calculate_step_parameters(
                from_concentration=0.0,
                to_concentration=5.0,
                temperature=4.0,
                sample_type=sample_type
            )
            self.assertIn("hold_time_min", params)
            self.assertIn("notes", params)
        
        # Test with custom sensitivity
        custom_sensitivity = {
            "time_per_step": 10
        }
        params = Protocol.ProtocolStep.calculate_step_parameters(
            from_concentration=0.0,
            to_concentration=5.0,
            temperature=4.0,
            custom_sensitivity=custom_sensitivity
        )
        
        # Hold time should be adjusted based on concentration change and custom sensitivity
        self.assertGreater(params["hold_time_min"], 10)
        
        # Test with small concentration change
        params = Protocol.ProtocolStep.calculate_step_parameters(
            from_concentration=0.0,
            to_concentration=1.0,
            temperature=4.0
        )
        
        # Hold time should be smaller for smaller concentration changes
        self.assertLessEqual(params["hold_time_min"],
                            Protocol.ProtocolStep.calculate_step_parameters(
                                from_concentration=0.0,
                                to_concentration=10.0,
                                temperature=4.0
                            )["hold_time_min"])

class TestProtocolComparison(unittest.TestCase):
    """Test cases for the Protocol.ProtocolComparison class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Mock Protocol.get_protocol
        self.original_get_protocol = Protocol.get_protocol
        Protocol.get_protocol = MockProtocol.get_protocol
    
    def tearDown(self):
        """Tear down test fixtures."""
        # Restore original method
        Protocol.get_protocol = self.original_get_protocol
    
    def test_create_comparison(self):
        """Test the create_comparison method."""
        # Test with valid protocol IDs
        comparison = Protocol.ProtocolComparison.create_comparison(["protocol1", "protocol2"])
        
        # Check that the comparison has the expected structure
        self.assertIn("id", comparison)
        self.assertIn("protocol_ids", comparison)
        self.assertIn("protocols", comparison)
        self.assertIn("created_at", comparison)
        self.assertIn("parameter_comparison", comparison)
        self.assertIn("step_comparison", comparison)
        self.assertIn("summary", comparison)
        self.assertIn("recommendations", comparison)
        
        # Check that protocols are included
        self.assertEqual(len(comparison["protocols"]), 2)
        
        # Check that parameter comparison includes all expected parameters
        param_keys = [
            "target_concentration",
            "sample_type",
            "starting_temperature",
            "target_temperature",
            "step_count"
        ]
        for key in param_keys:
            self.assertIn(key, comparison["parameter_comparison"])
            
        # Check that step comparison includes all steps
        max_steps = max(len(p.get("steps", [])) for p in comparison["protocols"])
        for i in range(max_steps):
            self.assertIn(f"step_{i+1}", comparison["step_comparison"])
            
        # Check that summary includes expected fields
        self.assertEqual(comparison["summary"]["protocol_count"], 2)
        self.assertEqual(len(comparison["summary"]["sample_types"]), 2)
        self.assertEqual(comparison["summary"]["concentration_range"], [10.0, 15.0])
        self.assertEqual(comparison["summary"]["step_count_range"], [4, 6])
        
        # Test with insufficient protocol IDs
        with self.assertRaises(ValueError) as context:
            Protocol.ProtocolComparison.create_comparison(["protocol1"])
        self.assertIn("At least two protocol IDs are required for comparison", str(context.exception))
        
        # Test with non-existent protocol ID
        with self.assertRaises(ValueError) as context:
            Protocol.ProtocolComparison.create_comparison(["protocol1", "nonexistent"])
        self.assertIn("Protocol with ID nonexistent not found", str(context.exception))

    def test_save_comparison(self):
        """Test the save_comparison method."""
        # Create a comparison
        comparison = Protocol.ProtocolComparison.create_comparison(["protocol1", "protocol2"])
        
        # Save the comparison
        saved_comparison = Protocol.ProtocolComparison.save_comparison(comparison)
        
        # Check that the saved comparison has the expected structure
        self.assertIn("id", saved_comparison)
        self.assertIn("created_at", saved_comparison)
        
        # ID and created_at should be different from the original
        self.assertNotEqual(saved_comparison["id"], comparison["id"])
        self.assertNotEqual(saved_comparison["created_at"], comparison["created_at"])

class TestProtocolDesignerIntegration(MockSupabaseBaseTestCase):
    """Integration tests for the Protocol Designer module."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Mock the Mixture and Protocol classes
        self.original_mixture_get = Mixture.get
        self.original_mixture_get_with_components = Mixture.get_with_components
        self.original_protocol_create_protocol = Protocol.create_protocol
        self.original_protocol_get_protocol = Protocol.get_protocol
        self.original_protocol_get_protocols_for_mixture = Protocol.get_protocols_for_mixture
        
        Mixture.get = MockMixture.get
        Mixture.get_with_components = MockMixture.get_with_components
        Protocol.create_protocol = MockProtocol.create_protocol
        Protocol.get_protocol = MockProtocol.get_protocol
        Protocol.get_protocols_for_mixture = MockProtocol.get_protocols_for_mixture
    
    def tearDown(self):
        """Tear down test fixtures."""
        # Restore the original methods
        Mixture.get = self.original_mixture_get
        Mixture.get_with_components = self.original_mixture_get_with_components
        Protocol.create_protocol = self.original_protocol_create_protocol
        Protocol.get_protocol = self.original_protocol_get_protocol
        Protocol.get_protocols_for_mixture = self.original_protocol_get_protocols_for_mixture
    
    def test_end_to_end_protocol_workflow(self):
        """Test the end-to-end workflow of creating, saving, and comparing protocols."""
        # 1. Get sensitivity profiles
        profiles = ProtocolDesigner.get_sample_sensitivity_profiles()
        cell_line_profile = profiles["cell_line"]
        primary_cells_profile = profiles["primary_cells"]
        
        # 2. Create two protocols with different sensitivity profiles
        protocol1 = ProtocolDesigner.design_concentration_gradient(
            mixture_id="mix1",
            target_concentration=10.0,
            sample_type="cell_line",
            starting_temperature=4.0,
            target_temperature=0.0,
            custom_sensitivity={
                "max_step_size": cell_line_profile["max_step_size"],
                "time_per_step": cell_line_profile["time_per_step"]
            }
        )
        
        protocol2 = ProtocolDesigner.design_concentration_gradient(
            mixture_id="mix1",
            target_concentration=15.0,
            sample_type="primary_cells",
            starting_temperature=4.0,
            target_temperature=-1.0,
            custom_sensitivity={
                "max_step_size": primary_cells_profile["max_step_size"],
                "time_per_step": primary_cells_profile["time_per_step"]
            }
        )
        
        # 3. Save the protocols
        saved_protocol1 = ProtocolDesigner.save_protocol(
            mixture_id="mix1",
            name="Cell Line Protocol",
            description="Protocol for cell lines",
            protocol_data=protocol1
        )
        
        saved_protocol2 = ProtocolDesigner.save_protocol(
            mixture_id="mix1",
            name="Primary Cells Protocol",
            description="Protocol for primary cells",
            protocol_data=protocol2
        )
        
        # 4. Compare the protocols
        comparison = ProtocolDesigner.compare_protocols([saved_protocol1["id"], saved_protocol2["id"]])
        
        # 5. Verify the comparison results
        self.assertEqual(len(comparison["protocols"]), 2)
        self.assertIn("recommendations", comparison)
        self.assertTrue(len(comparison["recommendations"]) > 0)
        
        # 6. Check that the comparison identifies key differences
        sample_type_comparison = comparison["parameter_comparison"]["sample_type"]
        self.assertFalse(sample_type_comparison["same"])
        
        target_conc_comparison = comparison["parameter_comparison"]["target_concentration"]
        self.assertFalse(target_conc_comparison["same"])
        
        # 7. Verify that step comparison works correctly
        step1_comparison = comparison["step_comparison"]["step_1"]
        self.assertEqual(len(step1_comparison["concentration_values"]), 2)
        
        # Both protocols should have the same first step concentration
        self.assertEqual(step1_comparison["concentration_values"][0],
                         step1_comparison["concentration_values"][1])

if __name__ == '__main__':
    unittest.main()
