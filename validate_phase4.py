#!/usr/bin/env python3
"""
Phase 4 Validation Script - Validates the implementation of Phase 4 components.

This script checks that all components of Phase 4 have been properly implemented
and are functioning correctly. It runs a series of validation checks on:
1. Scientific models
2. Enhanced RDKit integration
3. API extensions
4. Database schema
"""
import sys
import os
import json
import logging
import importlib
import unittest
from pathlib import Path
import inspect
import subprocess
import importlib.util

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("phase4_validation")

# Define validation results class
class ValidationResult:
    def __init__(self, component):
        self.component = component
        self.checks = []
        self.passed = True
    
    def add_check(self, name, passed, message=None):
        self.checks.append({
            "name": name,
            "passed": passed,
            "message": message
        })
        if not passed:
            self.passed = False
    
    def get_summary(self):
        total = len(self.checks)
        passed = sum(1 for check in self.checks if check["passed"])
        return {
            "component": self.component,
            "status": "PASSED" if self.passed else "FAILED",
            "checks_passed": passed,
            "checks_total": total,
            "percentage": round((passed / total) * 100) if total > 0 else 0,
            "checks": self.checks
        }

# Validation functions
def validate_file_exists(file_path):
    """Check if a file exists."""
    return os.path.exists(file_path)

def validate_module_imports(module_name):
    """Check if a Python module can be imported."""
    try:
        importlib.import_module(module_name)
        return True, None
    except ImportError as e:
        return False, str(e)

def validate_class_exists(module_name, class_name):
    """Check if a class exists in a module."""
    try:
        module = importlib.import_module(module_name)
        return hasattr(module, class_name), None
    except ImportError as e:
        return False, str(e)

def validate_method_exists(module_name, class_name, method_name):
    """Check if a method exists in a class."""
    try:
        module = importlib.import_module(module_name)
        if class_name is None:
            # We're checking for a function in the module instead of a method in a class
            return hasattr(module, method_name), None
        else:
            # We're checking for a method in a class
            cls = getattr(module, class_name)
            return hasattr(cls, method_name), None
    except (ImportError, AttributeError) as e:
        return False, str(e)

def run_unit_tests(test_module):
    """Run unit tests for a module."""
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromName(test_module)
    runner = unittest.TextTestRunner(verbosity=0)
    result = runner.run(suite)
    return result.wasSuccessful(), f"Failures: {len(result.failures)}, Errors: {len(result.errors)}"

def validate_scientific_models():
    """Validate scientific models implementation."""
    result = ValidationResult("Scientific Models")
    
    # Check base files
    result.add_check(
        "scientific_models package exists",
        validate_file_exists("scientific_models/__init__.py")
    )
    
    result.add_check(
        "scientific_models.base module exists",
        validate_file_exists("scientific_models/base.py")
    )
    
    # Check model implementations
    model_modules = [
        "concentration", "temperature", "mixtures", "rdkit_enhanced"
    ]
    for module in model_modules:
        exists = validate_file_exists(f"scientific_models/{module}.py")
        result.add_check(
            f"scientific_models.{module} module exists",
            exists
        )
    
    # Check base classes
    base_classes = [
        "ScientificModel", "MoleculePropertyModel", "MixtureModel"
    ]
    for cls in base_classes:
        passed, error = validate_class_exists("scientific_models.base", cls)
        result.add_check(
            f"scientific_models.base.{cls} class exists",
            passed,
            error
        )
    
    # Check if core methods exist in ScientificModel
    core_methods = ["validate_inputs", "calculate"]
    for method in core_methods:
        passed, error = validate_method_exists(
            "scientific_models.base", "ScientificModel", method
        )
        result.add_check(
            f"scientific_models.base.ScientificModel.{method} method exists",
            passed,
            error
        )
    
    # Check concentration models
    concentration_models = [
        "LinearConcentrationModel", "ExponentialConcentrationModel", 
        "SigmoidConcentrationModel"
    ]
    for model in concentration_models:
        try:
            passed, error = validate_class_exists(
                "scientific_models.concentration", model
            )
            result.add_check(
                f"scientific_models.concentration.{model} class exists",
                passed,
                error
            )
        except:
            result.add_check(
                f"scientific_models.concentration.{model} class exists",
                False,
                "Error checking class"
            )
    
    # Check temperature models
    temperature_models = [
        "LinearTemperatureModel", "ArrheniusTemperatureModel",
        "GlassTransitionModel"
    ]
    for model in temperature_models:
        try:
            passed, error = validate_class_exists(
                "scientific_models.temperature", model
            )
            result.add_check(
                f"scientific_models.temperature.{model} class exists",
                passed,
                error
            )
        except:
            result.add_check(
                f"scientific_models.temperature.{model} class exists",
                False,
                "Error checking class"
            )
    
    # Check mixture models
    mixture_models = [
        "MixtureOptimizationModel", "GeneticOptimizationModel",
        "GridSearchModel", "SynergyPredictionModel"
    ]
    for model in mixture_models:
        try:
            passed, error = validate_class_exists(
                "scientific_models.mixtures", model
            )
            result.add_check(
                f"scientific_models.mixtures.{model} class exists",
                passed,
                error
            )
        except:
            result.add_check(
                f"scientific_models.mixtures.{model} class exists",
                False,
                "Error checking class"
            )
    
    return result

def validate_enhanced_rdkit():
    """Validate enhanced RDKit implementation."""
    result = ValidationResult("Enhanced RDKit")
    
    # Check files
    result.add_check(
        "Dockerfile.rdkit-enhanced exists",
        validate_file_exists("Dockerfile.rdkit-enhanced")
    )
    
    result.add_check(
        "requirements_rdkit_enhanced.txt exists",
        validate_file_exists("requirements_rdkit_enhanced.txt")
    )
    
    result.add_check(
        "rdkit_enhanced_service.py exists",
        validate_file_exists("rdkit_enhanced_service.py")
    )
    
    result.add_check(
        "run_rdkit_enhanced.sh exists",
        validate_file_exists("run_rdkit_enhanced.sh")
    )
    
    # Check if run_rdkit_enhanced.sh is executable
    if validate_file_exists("run_rdkit_enhanced.sh"):
        try:
            is_executable = os.access("run_rdkit_enhanced.sh", os.X_OK)
            result.add_check(
                "run_rdkit_enhanced.sh is executable",
                is_executable
            )
        except:
            result.add_check(
                "run_rdkit_enhanced.sh is executable",
                False,
                "Error checking execution permissions"
            )
    
    # Check RDKit Python interface
    result.add_check(
        "scientific_models/rdkit_enhanced.py exists",
        validate_file_exists("scientific_models/rdkit_enhanced.py")
    )
    
    rdkit_classes = [
        "EnhancedRDKitCalculator", "EnhancedRDKitError",
        "DescriptorBasedModel", "MolecularSimilarityModel",
        "ConformerAnalysisModel", "PharmacophoreModel"
    ]
    for cls in rdkit_classes:
        try:
            passed, error = validate_class_exists(
                "scientific_models.rdkit_enhanced", cls
            )
            result.add_check(
                f"scientific_models.rdkit_enhanced.{cls} class exists",
                passed,
                error
            )
        except:
            result.add_check(
                f"scientific_models.rdkit_enhanced.{cls} class exists",
                False,
                "Error checking class"
            )
    
    # Check RDKit API resources
    result.add_check(
        "api/rdkit_enhanced_resources.py exists",
        validate_file_exists("api/rdkit_enhanced_resources.py")
    )
    
    # Direct function check instead of using validate_method_exists
    if validate_file_exists("api/rdkit_enhanced_resources.py"):
        try:
            import api.rdkit_enhanced_resources
            result.add_check(
                "api.rdkit_enhanced_resources.register_resources function exists",
                hasattr(api.rdkit_enhanced_resources, "register_resources"),
                None
            )
        except Exception as e:
            result.add_check(
                "api.rdkit_enhanced_resources.register_resources function exists",
                False,
                str(e)
            )
    else:
        result.add_check(
            "api.rdkit_enhanced_resources.register_resources function exists",
            False,
            "RDKit enhanced resources file does not exist"
        )
    
    # Check API resource classes
    rdkit_resources = [
        "RDKitDescriptorsResource", "RDKitConformersResource",
        "RDKitSimilarityResource", "RDKitPharmacophoreResource",
        "RDKitConformerAnalysisResource", "RDKitPropertyPredictionResource",
        "RDKitServiceInfoResource"
    ]
    for resource in rdkit_resources:
        try:
            passed, error = validate_class_exists(
                "api.rdkit_enhanced_resources", resource
            )
            result.add_check(
                f"api.rdkit_enhanced_resources.{resource} class exists",
                passed,
                error
            )
        except:
            result.add_check(
                f"api.rdkit_enhanced_resources.{resource} class exists",
                False,
                "Error checking class"
            )
    
    # Check test files
    result.add_check(
        "tests/test_rdkit_enhanced.py exists",
        validate_file_exists("tests/test_rdkit_enhanced.py")
    )
    
    result.add_check(
        "tests/test_rdkit_enhanced_integration.py exists",
        validate_file_exists("tests/test_rdkit_enhanced_integration.py")
    )
    
    # Try to run unit tests if they exist
    if validate_file_exists("tests/test_rdkit_enhanced.py"):
        try:
            passed, message = run_unit_tests("tests.test_rdkit_enhanced")
            result.add_check(
                "Unit tests for rdkit_enhanced pass",
                passed,
                message
            )
        except Exception as e:
            result.add_check(
                "Unit tests for rdkit_enhanced pass",
                False,
                str(e)
            )
    
    # Check documentation
    result.add_check(
        "RDKIT_ENHANCED_ENVIRONMENT.md exists",
        validate_file_exists("RDKIT_ENHANCED_ENVIRONMENT.md")
    )
    
    result.add_check(
        "ENHANCED_RDKIT_GUIDE.md exists",
        validate_file_exists("ENHANCED_RDKIT_GUIDE.md")
    )
    
    return result

def validate_api_extensions():
    """Validate API extensions implementation."""
    result = ValidationResult("API Extensions")
    
    # Check scientific resources
    result.add_check(
        "api/scientific_resources.py exists",
        validate_file_exists("api/scientific_resources.py")
    )
    
    # Direct function check instead of using validate_method_exists
    if validate_file_exists("api/scientific_resources.py"):
        try:
            import api.scientific_resources
            result.add_check(
                "api.scientific_resources.register_resources function exists",
                hasattr(api.scientific_resources, "register_resources"),
                None
            )
            result.add_check(
                "api.scientific_resources.register_scientific_docs function exists",
                hasattr(api.scientific_resources, "register_scientific_docs"),
                None
            )
        except Exception as e:
            result.add_check(
                "api.scientific_resources.register_resources function exists",
                False,
                str(e)
            )
            result.add_check(
                "api.scientific_resources.register_scientific_docs function exists",
                False,
                str(e)
            )
    else:
        result.add_check(
            "api.scientific_resources.register_resources function exists",
            False,
            "Scientific resources file does not exist"
        )
        result.add_check(
            "api.scientific_resources.register_scientific_docs function exists",
            False,
            "Scientific resources file does not exist"
        )
    
    # Check if API extensions are registered in __init__.py
    if validate_file_exists("api/__init__.py"):
        try:
            with open("api/__init__.py", "r") as f:
                init_content = f.read()
            
            result.add_check(
                "Scientific resources imported in api/__init__.py",
                "from api.scientific_resources import" in init_content
            )
            
            result.add_check(
                "RDKit Enhanced resources imported in api/__init__.py",
                "from api.rdkit_enhanced_resources import" in init_content
            )
            
            result.add_check(
                "Scientific resources registered in api/__init__.py",
                "register_scientific_resources(api)" in init_content
            )
            
            result.add_check(
                "Enhanced RDKit resources registered in api/__init__.py",
                "register_rdkit_enhanced_resources(api)" in init_content
            )
        except:
            result.add_check(
                "API extensions registered in api/__init__.py",
                False,
                "Error checking API registration"
            )
    
    return result

def validate_database_schema():
    """Validate database schema extensions implementation."""
    result = ValidationResult("Database Schema")
    
    # Check migration files
    migration_files = [
        "migrations/030_scientific_models_schema.sql"
    ]
    for file in migration_files:
        result.add_check(
            f"{file} exists",
            validate_file_exists(file)
        )
    
    # Check if database schema migration has appropriate content
    if validate_file_exists("migrations/030_scientific_models_schema.sql"):
        try:
            with open("migrations/030_scientific_models_schema.sql", "r") as f:
                migration_content = f.read().lower()
            
            schema_elements = [
                "create table", "concentration_models", 
                "temperature_models", "mixture_optimizations", 
                "create index", "create materialized view"
            ]
            
            for element in schema_elements:
                result.add_check(
                    f"Schema migration contains '{element}'",
                    element in migration_content
                )
        except:
            result.add_check(
                "Schema migration contains required elements",
                False,
                "Error checking migration file content"
            )
    
    return result

def run_validation():
    """Run all validation checks and generate a report."""
    logger.info("Starting Phase 4 validation...")
    
    # Run validations
    results = [
        validate_scientific_models(),
        validate_enhanced_rdkit(),
        validate_api_extensions(),
        validate_database_schema()
    ]
    
    # Generate summary
    summaries = [result.get_summary() for result in results]
    
    total_checks = sum(summary["checks_total"] for summary in summaries)
    total_passed = sum(summary["checks_passed"] for summary in summaries)
    overall_percentage = round((total_passed / total_checks) * 100) if total_checks > 0 else 0
    
    overall_status = "PASSED" if all(summary["status"] == "PASSED" for summary in summaries) else "FAILED"
    
    report = {
        "overall_status": overall_status,
        "checks_passed": total_passed,
        "checks_total": total_checks,
        "percentage": overall_percentage,
        "components": summaries
    }
    
    # Print report
    logger.info(f"Validation completed. Overall status: {overall_status} ({overall_percentage}%)")
    for summary in summaries:
        logger.info(f"  {summary['component']}: {summary['status']} ({summary['percentage']}%)")
    
    # Save report
    with open("phase4_validation_report.json", "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Validation report saved to 'phase4_validation_report.json'")
    
    return overall_status == "PASSED"

if __name__ == "__main__":
    success = run_validation()
    sys.exit(0 if success else 1)