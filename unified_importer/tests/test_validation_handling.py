"""
Unit tests for the validation error handling system.
"""

import unittest
import tempfile
import os
import time
import json
from unittest.mock import Mock, patch
import logging

from unified_importer.core.error_handling.error_classification import (
    ErrorCategory, ErrorSeverity, RecoveryStrategy, ErrorContext
)
from unified_importer.core.error_handling.validation_handling import (
    ValidationResult, ValidationError, ValidationReport,
    ValidationErrorHandler, ChemicalValidationRules
)

class TestValidationHandling(unittest.TestCase):
    """Tests for the validation error handling system."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.logger = logging.getLogger("test")
        self.logger.setLevel(logging.DEBUG)
    
    def test_validation_result(self):
        """Test ValidationResult class."""
        # Valid result
        valid_result = ValidationResult(
            is_valid=True,
            data="test data"
        )
        self.assertTrue(valid_result.is_valid)
        self.assertTrue(valid_result)  # Bool context
        
        # Invalid result
        invalid_result = ValidationResult(
            is_valid=False,
            data="invalid data",
            errors=["Error 1", "Error 2"]
        )
        self.assertFalse(invalid_result.is_valid)
        self.assertFalse(invalid_result)  # Bool context
        self.assertEqual(len(invalid_result.errors), 2)
    
    def test_validation_error(self):
        """Test ValidationError class."""
        error = ValidationError(
            message="Invalid SMILES",
            data="C1CC1C",
            field="smiles",
            code="INVALID_SMILES",
            severity=ErrorSeverity.MEDIUM,
            context={"extra": "info"}
        )
        
        self.assertEqual(error.message, "Invalid SMILES")
        self.assertEqual(error.data, "C1CC1C")
        self.assertEqual(error.field, "smiles")
        self.assertEqual(error.code, "INVALID_SMILES")
        self.assertEqual(error.severity, ErrorSeverity.MEDIUM)
        self.assertEqual(error.context["extra"], "info")
        
        # Test to_dict method
        error_dict = error.to_dict()
        self.assertEqual(error_dict["message"], "Invalid SMILES")
        self.assertEqual(error_dict["data"], "C1CC1C")
        self.assertEqual(error_dict["field"], "smiles")
        self.assertEqual(error_dict["code"], "INVALID_SMILES")
        self.assertEqual(error_dict["severity"], "MEDIUM")
    
    def test_validation_report(self):
        """Test ValidationReport class."""
        report = ValidationReport()
        
        # Initially empty
        self.assertEqual(report.valid_count, 0)
        self.assertEqual(report.invalid_count, 0)
        self.assertEqual(report.error_count, 0)
        self.assertEqual(report.warning_count, 0)
        self.assertTrue(report.is_empty())
        
        # Add errors
        error1 = ValidationError(
            message="Invalid SMILES",
            data="C1CC1C",
            field="smiles",
            code="INVALID_SMILES",
            severity=ErrorSeverity.MEDIUM
        )
        
        error2 = ValidationError(
            message="Name too short",
            data="X",
            field="name",
            code="INVALID_NAME",
            severity=ErrorSeverity.LOW
        )
        
        report.add_error(error1)
        report.add_error(error2)
        
        # Add warnings
        warning1 = ValidationError(
            message="Unusual atom",
            data="C1CCU1C",
            field="smiles",
            code="UNUSUAL_ATOM",
            severity=ErrorSeverity.LOW
        )
        
        report.add_warning(warning1)
        
        # Mark items
        report.mark_valid()
        report.mark_valid()
        report.mark_invalid()
        report.mark_invalid()
        report.mark_invalid()
        
        # Check counts
        self.assertEqual(report.valid_count, 2)
        self.assertEqual(report.invalid_count, 3)
        self.assertEqual(report.error_count, 2)
        self.assertEqual(report.warning_count, 1)
        self.assertFalse(report.is_empty())
        
        # Check indexes
        self.assertEqual(len(report.by_field["smiles"]), 2)  # 1 error, 1 warning
        self.assertEqual(len(report.by_field["name"]), 1)
        self.assertEqual(len(report.by_code["INVALID_SMILES"]), 1)
        self.assertEqual(len(report.by_code["INVALID_NAME"]), 1)
        self.assertEqual(len(report.by_code["UNUSUAL_ATOM"]), 1)
        self.assertEqual(len(report.by_severity["MEDIUM"]), 1)
        self.assertEqual(len(report.by_severity["LOW"]), 2)
        
        # Test get_summary
        summary = report.get_summary()
        self.assertEqual(summary["total_items"], 5)
        self.assertEqual(summary["valid_items"], 2)
        self.assertEqual(summary["invalid_items"], 3)
        self.assertEqual(summary["error_count"], 2)
        self.assertEqual(summary["warning_count"], 1)
        self.assertIn("smiles", summary["fields_with_errors"])
        self.assertIn("name", summary["fields_with_errors"])
        self.assertIn("INVALID_SMILES", summary["error_codes"])
        self.assertEqual(summary["error_counts_by_severity"]["MEDIUM"], 1)
        
        # Test to_dict and file operations
        report_dict = report.to_dict()
        self.assertEqual(report_dict["valid_count"], 2)
        self.assertEqual(report_dict["invalid_count"], 3)
        self.assertEqual(len(report_dict["errors"]), 2)
        self.assertEqual(len(report_dict["warnings"]), 1)
        
        # Save to file and reload
        with tempfile.NamedTemporaryFile(delete=False, suffix=".json") as temp_file:
            report_path = temp_file.name
        
        try:
            report.save_to_file(report_path)
            self.assertTrue(os.path.exists(report_path))
            
            loaded_report = ValidationReport.load_from_file(report_path)
            self.assertEqual(loaded_report.valid_count, 2)
            self.assertEqual(loaded_report.invalid_count, 3)
            self.assertEqual(loaded_report.error_count, 2)
            self.assertEqual(loaded_report.warning_count, 1)
            self.assertEqual(len(loaded_report.errors), 2)
            self.assertEqual(len(loaded_report.warnings), 1)
            self.assertEqual(len(loaded_report.by_field["smiles"]), 2)
        finally:
            if os.path.exists(report_path):
                os.remove(report_path)
    
    def test_validation_report_merge(self):
        """Test merging validation reports."""
        report1 = ValidationReport()
        report1.mark_valid()
        report1.mark_invalid()
        report1.add_error(ValidationError(
            message="Error 1",
            data="data1",
            field="field1",
            code="CODE1"
        ))
        
        report2 = ValidationReport()
        report2.mark_valid()
        report2.mark_valid()
        report2.mark_invalid()
        report2.add_error(ValidationError(
            message="Error 2",
            data="data2",
            field="field2",
            code="CODE2"
        ))
        report2.add_warning(ValidationError(
            message="Warning 1",
            data="data3",
            field="field1",
            code="WARN1"
        ))
        
        # Merge reports
        report1.merge(report2)
        
        # Check merged results
        self.assertEqual(report1.valid_count, 3)
        self.assertEqual(report1.invalid_count, 2)
        self.assertEqual(report1.error_count, 2)
        self.assertEqual(report1.warning_count, 1)
        self.assertEqual(len(report1.by_field["field1"]), 2)
        self.assertEqual(len(report1.by_field["field2"]), 1)
        self.assertEqual(len(report1.by_code["CODE1"]), 1)
        self.assertEqual(len(report1.by_code["CODE2"]), 1)
        self.assertEqual(len(report1.by_code["WARN1"]), 1)
    
    def test_validation_error_handler_basic(self):
        """Test basic ValidationErrorHandler functionality."""
        handler = ValidationErrorHandler(logger=self.logger)
        
        # Test validation function
        def validate_positive(value):
            if value <= 0:
                raise ValueError("Value must be positive")
            return value
        
        # Test with valid input
        result = handler.handle_validation(
            validate_positive, 42,
            component="Number Validator",
            operation="validate_positive"
        )
        self.assertEqual(result, 42)
        
        # Check report
        report = handler.get_report()
        self.assertEqual(report.valid_count, 1)
        self.assertEqual(report.invalid_count, 0)
        self.assertEqual(report.error_count, 0)
        
        # Test with invalid input and SKIP strategy
        result = handler.handle_validation(
            validate_positive, -1,
            component="Number Validator",
            operation="validate_positive",
            error_strategy=RecoveryStrategy.SKIP,
            default_value=0
        )
        self.assertEqual(result, 0)
        
        # Check updated report
        report = handler.get_report()
        self.assertEqual(report.valid_count, 1)
        self.assertEqual(report.invalid_count, 1)
        self.assertEqual(report.error_count, 1)
        
        # Reset report
        handler.reset_report()
        report = handler.get_report()
        self.assertEqual(report.valid_count, 0)
        self.assertEqual(report.error_count, 0)
    
    def test_validation_error_handler_strategies(self):
        """Test ValidationErrorHandler with different strategies."""
        handler = ValidationErrorHandler(logger=self.logger)
        
        def validate_in_range(value):
            if not (0 <= value <= 100):
                raise ValueError(f"Value {value} must be between 0 and 100")
            return value
        
        # Test FALLBACK strategy
        result = handler.handle_validation(
            validate_in_range, 150,
            component="Range Validator",
            operation="validate_in_range",
            error_strategy=RecoveryStrategy.FALLBACK,
            fallback_func=lambda x: min(max(0, x), 100)  # Clamp value to range
        )
        self.assertEqual(result, 100)
        
        # Test ABORT strategy
        with self.assertRaises(ValueError):
            handler.handle_validation(
                validate_in_range, 150,
                component="Range Validator",
                operation="validate_in_range",
                error_strategy=RecoveryStrategy.ABORT
            )
        
        # Test LOG_ONLY strategy
        result = handler.handle_validation(
            validate_in_range, 150,
            component="Range Validator",
            operation="validate_in_range",
            error_strategy=RecoveryStrategy.LOG_ONLY
        )
        self.assertEqual(result, 150)  # Original value returned despite being invalid
    
    def test_validate_batch(self):
        """Test batch validation."""
        handler = ValidationErrorHandler(logger=self.logger)
        
        # Validation function
        def validate_positive(value):
            if not isinstance(value, (int, float)) or value <= 0:
                raise ValueError(f"Value must be a positive number: {value}")
            return value
        
        # Create test data with some invalid items
        items = [5, 10, -3, "not a number", 20, 0]
        
        # Test batch validation with SKIP strategy
        results, report = handler.validate_batch(
            validate_positive, items,
            component="Number Validator",
            operation="validate_batch",
            error_strategy=RecoveryStrategy.SKIP,
            default_value=1
        )
        
        # Check results
        self.assertEqual(results[0], 5)
        self.assertEqual(results[1], 10)
        self.assertEqual(results[2], 1)  # Default for invalid
        self.assertEqual(results[3], 1)  # Default for invalid
        self.assertEqual(results[4], 20)
        self.assertEqual(results[5], 1)  # Default for invalid
        
        # Check report
        self.assertEqual(report.valid_count, 3)
        self.assertEqual(report.invalid_count, 3)
        self.assertEqual(report.error_count, 3)
        
        # Test batch validation with FALLBACK strategy
        results, report = handler.validate_batch(
            validate_positive, items,
            component="Number Validator",
            operation="validate_batch",
            error_strategy=RecoveryStrategy.FALLBACK,
            fallback_func=lambda x: 1 if not isinstance(x, (int, float)) or x <= 0 else x
        )
        
        # Check results
        self.assertEqual(results[0], 5)
        self.assertEqual(results[1], 10)
        self.assertEqual(results[2], 1)  # Fallback for invalid
        self.assertEqual(results[3], 1)  # Fallback for invalid
        self.assertEqual(results[4], 20)
        self.assertEqual(results[5], 1)  # Fallback for invalid
    
    def test_validate_batch_parallel(self):
        """Test parallel batch validation."""
        handler = ValidationErrorHandler(logger=self.logger)
        
        # Validation function
        def validate_positive(value):
            if not isinstance(value, (int, float)) or value <= 0:
                raise ValueError(f"Value must be a positive number: {value}")
            return value
        
        # Create test data with some invalid items
        items = [5, 10, -3, "not a number", 20, 0]
        
        # Test parallel batch validation
        results, report = handler.validate_batch(
            validate_positive, items,
            component="Number Validator",
            operation="validate_batch_parallel",
            error_strategy=RecoveryStrategy.SKIP,
            default_value=1,
            parallel=True,
            max_workers=3
        )
        
        # Check results (order may vary due to parallel execution)
        self.assertEqual(len(results), 6)
        valid_results = [r for r in results if r in (5, 10, 20)]
        default_results = [r for r in results if r == 1]
        self.assertEqual(len(valid_results), 3)
        self.assertEqual(len(default_results), 3)
        
        # Check report
        self.assertEqual(report.valid_count, 3)
        self.assertEqual(report.invalid_count, 3)
        self.assertEqual(report.error_count, 3)
    
    def test_validation_abort_on_error(self):
        """Test batch validation with abort on error."""
        handler = ValidationErrorHandler(logger=self.logger)
        
        # Validation function
        def validate_positive(value):
            if not isinstance(value, (int, float)) or value <= 0:
                raise ValueError(f"Value must be a positive number: {value}")
            return value
        
        # Create test data with some invalid items
        items = [5, 10, -3, "not a number", 20, 0]
        
        # Test batch validation with continue_on_error=False
        with self.assertRaises(ValueError):
            results, report = handler.validate_batch(
                validate_positive, items,
                component="Number Validator",
                operation="validate_batch_abort",
                error_strategy=RecoveryStrategy.SKIP,
                default_value=1,
                continue_on_error=False
            )
    
    def test_chemical_validation_rules_smiles(self):
        """Test SMILES validation rules."""
        # Valid SMILES
        self.assertEqual(ChemicalValidationRules.validate_smiles("CC(=O)OC1=CC=CC=C1C(=O)O"), "CC(=O)OC1=CC=CC=C1C(=O)O")
        self.assertEqual(ChemicalValidationRules.validate_smiles("C"), "C")
        self.assertEqual(ChemicalValidationRules.validate_smiles("c1ccccc1"), "c1ccccc1")
        
        # Invalid SMILES
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_smiles("")
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_smiles(None)
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_smiles("XYZ123")  # No element symbols
    
    def test_chemical_validation_rules_inchi(self):
        """Test InChI validation rules."""
        # Valid InChI
        self.assertEqual(
            ChemicalValidationRules.validate_inchi("InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"),
            "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
        )
        
        # Invalid InChI
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_inchi("")
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_inchi(None)
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_inchi("1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)")
    
    def test_chemical_validation_rules_molecule_name(self):
        """Test molecule name validation rules."""
        # Valid names
        self.assertEqual(ChemicalValidationRules.validate_molecule_name("Aspirin"), "Aspirin")
        self.assertEqual(ChemicalValidationRules.validate_molecule_name("2-Acetoxybenzoic acid"), "2-Acetoxybenzoic acid")
        
        # Invalid names
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecule_name("")
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecule_name(None)
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecule_name("A")  # Too short
    
    def test_chemical_validation_rules_cid(self):
        """Test CID validation rules."""
        # Valid CIDs
        self.assertEqual(ChemicalValidationRules.validate_cid("2244"), "2244")
        self.assertEqual(ChemicalValidationRules.validate_cid(2244), "2244")
        
        # Invalid CIDs
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_cid("")
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_cid(None)
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_cid("CID2244")
    
    def test_chemical_validation_rules_chembl_id(self):
        """Test ChEMBL ID validation rules."""
        # Valid ChEMBL IDs
        self.assertEqual(ChemicalValidationRules.validate_chembl_id("CHEMBL25"), "CHEMBL25")
        self.assertEqual(ChemicalValidationRules.validate_chembl_id("CHEMBL1201"), "CHEMBL1201")
        
        # Invalid ChEMBL IDs
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_chembl_id("")
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_chembl_id(None)
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_chembl_id("25")
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_chembl_id("CHEMBL-25")
    
    def test_chemical_validation_rules_molecular_formula(self):
        """Test molecular formula validation rules."""
        # Valid formulas
        self.assertEqual(ChemicalValidationRules.validate_molecular_formula("C9H8O4"), "C9H8O4")
        self.assertEqual(ChemicalValidationRules.validate_molecular_formula("H2O"), "H2O")
        
        # Invalid formulas
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecular_formula("")
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecular_formula(None)
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecular_formula("123456")  # No element symbols
    
    def test_chemical_validation_rules_molecular_weight(self):
        """Test molecular weight validation rules."""
        # Valid weights
        self.assertEqual(ChemicalValidationRules.validate_molecular_weight(180.16), 180.16)
        self.assertEqual(ChemicalValidationRules.validate_molecular_weight("180.16"), 180.16)
        self.assertEqual(ChemicalValidationRules.validate_molecular_weight(18), 18.0)
        
        # Invalid weights
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecular_weight("")
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecular_weight(None)
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecular_weight("abc")
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecular_weight(0)
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecular_weight(-10)
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_molecular_weight(6000)
    
    def test_chemical_validation_rules_logp(self):
        """Test logP validation rules."""
        # Valid logP values
        self.assertEqual(ChemicalValidationRules.validate_logp(1.43), 1.43)
        self.assertEqual(ChemicalValidationRules.validate_logp("1.43"), 1.43)
        self.assertEqual(ChemicalValidationRules.validate_logp(-2), -2.0)
        self.assertEqual(ChemicalValidationRules.validate_logp(0), 0.0)
        
        # Invalid logP values
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_logp("")
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_logp(None)
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_logp("abc")
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_logp(-15)
        
        with self.assertRaises(ValueError):
            ChemicalValidationRules.validate_logp(25)


if __name__ == '__main__':
    unittest.main()