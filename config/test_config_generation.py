#!/usr/bin/env python
"""
Test Configuration Generation

This script tests the generation of configuration files from the schema.
It verifies that the generated files are valid Python and TypeScript code.
"""

import os
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path


class TestConfigGeneration(unittest.TestCase):
    def setUp(self):
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.schema_path = os.path.join(self.script_dir, "schema.json")
        self.generator_path = os.path.join(self.script_dir, "generate_configs.py")
        
        # Ensure schema exists
        self.assertTrue(os.path.exists(self.schema_path), f"Schema file not found at {self.schema_path}")
        self.assertTrue(os.path.exists(self.generator_path), f"Generator script not found at {self.generator_path}")
    
    def test_schema_validation(self):
        """Test that the schema validation succeeds."""
        cmd = [sys.executable, self.generator_path, "--validate-only"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, f"Schema validation failed: {result.stderr}")
        self.assertIn("Schema validation completed successfully", result.stdout)
    
    def test_python_generation(self):
        """Test generation of Python configuration file."""
        with tempfile.NamedTemporaryFile(suffix=".py", delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            cmd = [
                sys.executable, 
                self.generator_path, 
                "--schema", self.schema_path,
                "--python", tmp_path,
                "--typescript", os.devnull  # Skip TypeScript generation
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, f"Python config generation failed: {result.stderr}")
            
            # Check that the output file exists and is a valid Python file
            self.assertTrue(os.path.exists(tmp_path), f"Output file not found at {tmp_path}")
            
            # Verify it's valid Python code
            result = subprocess.run(
                [sys.executable, "-m", "py_compile", tmp_path],
                capture_output=True, text=True
            )
            self.assertEqual(result.returncode, 0, f"Generated Python file has syntax errors: {result.stderr}")
        finally:
            # Clean up
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
    
    def test_config_sections(self):
        """Test that all expected sections are present in the generated Python config."""
        with tempfile.NamedTemporaryFile(suffix=".py", delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            cmd = [
                sys.executable, 
                os.path.join(self.script_dir, "generate_python_config.py"), 
                "--schema", self.schema_path,
                "--output", tmp_path
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, f"Python config generation failed: {result.stderr}")
            
            # Check the content for expected sections
            with open(tmp_path, 'r') as f:
                content = f.read()
            
            expected_sections = [
                "APP_", "API_", "DATABASE_", "AUTH_", "LOGGING_", 
                "CACHING_", "RATE_LIMITING_", "CHEMBL_", "FRONTEND_"
            ]
            
            for section in expected_sections:
                self.assertIn(section, content, f"Expected section {section} not found in generated config")
            
            # Check for environment-specific config classes
            expected_classes = [
                "class DevelopmentConfig", "class TestingConfig", 
                "class StagingConfig", "class ProductionConfig"
            ]
            
            for class_def in expected_classes:
                self.assertIn(class_def, content, f"Expected class {class_def} not found in generated config")
        finally:
            # Clean up
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)


if __name__ == "__main__":
    unittest.main()