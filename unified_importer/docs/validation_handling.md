# Validation Error Handling Documentation

## Overview

The Validation Error Handling system provides a comprehensive framework for validating data, handling validation errors with configurable strategies, and generating detailed validation reports. It is particularly designed for chemical data validation in the CryoProtect Unified Importer but can be used for any type of data validation.

## Key Components

- **ValidationResult**: Represents the result of a validation operation
- **ValidationError**: Captures detailed information about validation errors
- **ValidationReport**: Aggregates validation errors and provides analysis
- **ValidationErrorHandler**: Implements validation logic with different strategies
- **ChemicalValidationRules**: Provides domain-specific validation for chemical data

## Core Classes

### ValidationResult

Represents the outcome of a validation operation:

```python
@dataclass
class ValidationResult:
    is_valid: bool           # Whether the data is valid
    data: Any                # The validated (or invalid) data
    errors: List[str]        # List of error messages
    warnings: List[str]      # List of warning messages
    context: Dict[str, Any]  # Additional context information
```

Usage example:

```python
# Create a validation result
result = ValidationResult(
    is_valid=True,
    data="CC(=O)OC1=CC=CC=C1C(=O)O",
    warnings=["Contains aromatic rings"]
)

# Use in boolean context
if result:
    process_valid_data(result.data)
else:
    handle_errors(result.errors)
```

### ValidationError

Captures detailed information about a validation error:

```python
@dataclass
class ValidationError:
    message: str                     # Human-readable error message
    data: Any                        # The invalid data
    field: Optional[str] = None      # Field name with the error
    code: Optional[str] = None       # Error code for programmatic handling
    severity: ErrorSeverity = ...    # Severity level from ErrorSeverity enum
    context: Dict[str, Any] = ...    # Additional context information
    timestamp: float = ...           # Time when the error occurred
```

Usage example:

```python
# Create a validation error
error = ValidationError(
    message="Invalid SMILES string",
    data="C1CC1CC",
    field="smiles",
    code="INVALID_SMILES",
    severity=ErrorSeverity.MEDIUM,
    context={"source": "ChEMBL"}
)

# Access error information
print(f"Error in field {error.field}: {error.message}")
```

### ValidationReport

Aggregates validation errors and provides analysis:

```python
@dataclass
class ValidationReport:
    valid_count: int = 0                     # Count of valid items
    invalid_count: int = 0                   # Count of invalid items
    error_count: int = 0                     # Count of errors
    warning_count: int = 0                   # Count of warnings
    errors: List[ValidationError] = ...      # List of errors
    warnings: List[ValidationError] = ...    # List of warnings
    by_field: Dict[str, List[ValidationError]] = ...    # Errors indexed by field
    by_code: Dict[str, List[ValidationError]] = ...     # Errors indexed by code
    by_severity: Dict[str, List[ValidationError]] = ... # Errors indexed by severity
```

Key methods:

```python
# Add an error to the report
report.add_error(error)

# Add a warning to the report
report.add_warning(warning)

# Mark an item as valid
report.mark_valid()

# Mark an item as invalid
report.mark_invalid()

# Get a summary of the validation
summary = report.get_summary()

# Convert to dictionary for serialization
report_dict = report.to_dict()

# Save to a file
report.save_to_file("validation_report.json")

# Load from a file
report = ValidationReport.load_from_file("validation_report.json")

# Merge with another report
report.merge(other_report)
```

### ValidationErrorHandler

Implements validation logic with different strategies:

```python
class ValidationErrorHandler:
    # Initialize with optional logger
    def __init__(self, logger=None)
    
    # Handle validation for a single item
    def handle_validation(
        self, 
        validation_func,         # Function to validate data
        data,                    # Data to validate
        component,               # Component name for context
        operation,               # Operation name for context
        error_strategy=RecoveryStrategy.SKIP,  # Strategy for handling errors
        default_value=None,      # Default value for SKIP strategy
        fallback_func=None,      # Function for FALLBACK strategy
        record_errors=True       # Whether to record errors in report
    )
    
    # Validate a batch of items
    def validate_batch(
        self, 
        validation_func,         # Function to validate each item
        items,                   # List of items to validate
        component,               # Component name for context
        operation,               # Operation name for context
        error_strategy=RecoveryStrategy.SKIP,  # Strategy for handling errors
        default_value=None,      # Default value for SKIP strategy
        fallback_func=None,      # Function for FALLBACK strategy
        parallel=False,          # Whether to validate in parallel
        max_workers=None,        # Max workers for parallel validation
        continue_on_error=True   # Whether to continue after errors
    )
    
    # Reset the current validation report
    def reset_report()
    
    # Get the current validation report
    def get_report()
    
    # Save the current report to a file
    def save_report(file_path)
```

### ChemicalValidationRules

Provides domain-specific validation for chemical data:

```python
class ChemicalValidationRules:
    # Validate SMILES string
    @staticmethod
    def validate_smiles(smiles)
    
    # Validate InChI string
    @staticmethod
    def validate_inchi(inchi)
    
    # Validate molecule name
    @staticmethod
    def validate_molecule_name(name)
    
    # Validate PubChem CID
    @staticmethod
    def validate_cid(cid)
    
    # Validate ChEMBL ID
    @staticmethod
    def validate_chembl_id(chembl_id)
    
    # Validate molecular formula
    @staticmethod
    def validate_molecular_formula(formula)
    
    # Validate molecular weight
    @staticmethod
    def validate_molecular_weight(weight)
    
    # Validate logP value
    @staticmethod
    def validate_logp(logp)
```

## Recovery Strategies

The validation error handling system supports multiple recovery strategies:

- **SKIP**: Skip the invalid data and return a default value
- **FALLBACK**: Apply a fallback function to transform invalid data
- **ABORT**: Raise the validation error and stop processing
- **LOG_ONLY**: Log the error but return the original data

## Usage Examples

### Basic Validation

```python
from unified_importer.core.error_handling import (
    ValidationErrorHandler, RecoveryStrategy, ChemicalValidationRules
)

# Create a validation handler
handler = ValidationErrorHandler()

# Validate a SMILES string
try:
    smiles = handler.handle_validation(
        ChemicalValidationRules.validate_smiles,
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        component="Molecule Importer",
        operation="validate_smiles"
    )
    print(f"Valid SMILES: {smiles}")
except Exception as e:
    print(f"Invalid SMILES: {e}")

# Validate with SKIP strategy for invalid data
smiles = handler.handle_validation(
    ChemicalValidationRules.validate_smiles,
    "invalid_smiles",
    component="Molecule Importer",
    operation="validate_smiles",
    error_strategy=RecoveryStrategy.SKIP,
    default_value="C"  # Default to methane
)
print(f"SMILES (with default if invalid): {smiles}")

# Validate with FALLBACK strategy
def normalize_smiles(smiles):
    """Simple fallback to normalize SMILES"""
    if not isinstance(smiles, str):
        return "C"  # Default to methane
    return smiles.strip().upper()

smiles = handler.handle_validation(
    ChemicalValidationRules.validate_smiles,
    " c1ccccc1 ",  # Benzene with whitespace
    component="Molecule Importer",
    operation="validate_smiles",
    error_strategy=RecoveryStrategy.FALLBACK,
    fallback_func=normalize_smiles
)
print(f"Normalized SMILES: {smiles}")
```

### Batch Validation

```python
from unified_importer.core.error_handling import (
    ValidationErrorHandler, RecoveryStrategy, ChemicalValidationRules
)

# Create a validation handler
handler = ValidationErrorHandler()

# Batch of SMILES to validate
smiles_batch = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin (valid)
    "C1CCCCC1",                   # Cyclohexane (valid)
    "invalid_smiles",             # Invalid
    "C1=CC=CC=C1",                # Benzene (valid)
    "C1CC2CC1CC2",                # Norbornane (valid)
    "XYZ"                         # Invalid
]

# Validate batch with SKIP strategy
results, report = handler.validate_batch(
    ChemicalValidationRules.validate_smiles,
    smiles_batch,
    component="Molecule Importer",
    operation="validate_smiles_batch",
    error_strategy=RecoveryStrategy.SKIP,
    default_value="C"  # Default to methane
)

# Print results
print(f"Valid: {report.valid_count}, Invalid: {report.invalid_count}")
for i, result in enumerate(results):
    print(f"Result {i+1}: {result}")

# Get a summary of validation errors
summary = report.get_summary()
print(f"Total errors: {summary['error_count']}")

# Save the report for later analysis
handler.save_report("validation_report.json")
```

### Parallel Batch Validation

```python
from unified_importer.core.error_handling import (
    ValidationErrorHandler, RecoveryStrategy, ChemicalValidationRules
)

# Create a validation handler
handler = ValidationErrorHandler()

# Large batch of data to validate
large_batch = [f"C{'C' * i}" for i in range(1000)]

# Validate in parallel
results, report = handler.validate_batch(
    ChemicalValidationRules.validate_smiles,
    large_batch,
    component="Molecule Importer",
    operation="validate_large_batch",
    error_strategy=RecoveryStrategy.SKIP,
    parallel=True,
    max_workers=8  # Use 8 worker threads
)

print(f"Validated {len(results)} items")
print(f"Valid: {report.valid_count}, Invalid: {report.invalid_count}")
```

### Creating Custom Validation Rules

```python
from unified_importer.core.error_handling import ValidationErrorHandler, RecoveryStrategy

# Create a custom validation function
def validate_molecule_descriptor(descriptor):
    """Validate a molecule descriptor dictionary."""
    if not isinstance(descriptor, dict):
        raise ValueError("Descriptor must be a dictionary")
    
    required_fields = ["smiles", "name", "molecular_weight"]
    for field in required_fields:
        if field not in descriptor:
            raise ValueError(f"Missing required field: {field}")
    
    # Validate molecular weight
    weight = descriptor["molecular_weight"]
    if not isinstance(weight, (int, float)) or weight <= 0:
        raise ValueError(f"Invalid molecular weight: {weight}")
    
    return descriptor

# Create a validation handler
handler = ValidationErrorHandler()

# Validate a descriptor
valid_descriptor = {
    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "name": "Aspirin",
    "molecular_weight": 180.16
}

invalid_descriptor = {
    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "name": "Aspirin"
    # Missing molecular_weight
}

# Validate with different strategies
result1 = handler.handle_validation(
    validate_molecule_descriptor,
    valid_descriptor,
    component="Descriptor Validator",
    operation="validate_descriptor"
)

result2 = handler.handle_validation(
    validate_molecule_descriptor,
    invalid_descriptor,
    component="Descriptor Validator",
    operation="validate_descriptor",
    error_strategy=RecoveryStrategy.SKIP,
    default_value={"smiles": "C", "name": "Unknown", "molecular_weight": 0}
)

print(f"Valid result: {result1['name']}")
print(f"Default result: {result2['name']}")
```

## Integration with Error Manager

The validation error handling system integrates seamlessly with the unified `ErrorManager`:

```python
from unified_importer.core.error_handling import (
    ErrorManager, RecoveryStrategy, ChemicalValidationRules
)

# Create an error manager
error_manager = ErrorManager()

# Validate a SMILES string
smiles = error_manager.validate(
    ChemicalValidationRules.validate_smiles,
    "CC(=O)OC1=CC=CC=C1C(=O)O",
    component="Molecule Importer",
    operation="validate_smiles"
)

# Validate a batch of SMILES
smiles_batch = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
    "C1CCCCC1",                   # Cyclohexane
    "invalid_smiles",             # Invalid
    "C1=CC=CC=C1"                 # Benzene
]

results, report = error_manager.validate_batch(
    ChemicalValidationRules.validate_smiles,
    smiles_batch,
    component="Molecule Importer",
    operation="validate_smiles_batch",
    error_strategy=RecoveryStrategy.SKIP,
    default_value="C"
)

# Get the validation report
report = error_manager.get_validation_report()
print(f"Valid: {report.valid_count}, Invalid: {report.invalid_count}")

# Reset the validation report
error_manager.reset_validation_report()
```

## Best Practices

### 1. Choose Appropriate Error Strategies

Select the error strategy based on the importance of the data:

- Use **ABORT** for critical data where invalid values are unacceptable
- Use **SKIP** for non-critical data where defaults are acceptable
- Use **FALLBACK** when data can be recovered or transformed
- Use **LOG_ONLY** when you want to track issues but continue with the original data

### 2. Provide Meaningful Default Values

When using the SKIP strategy, provide meaningful default values:

```python
# Good: Use reasonable defaults for missing properties
result = handler.handle_validation(
    validate_property,
    value,
    component="Property Validator",
    operation="validate_logP",
    error_strategy=RecoveryStrategy.SKIP,
    default_value=0.0  # Reasonable default for logP
)

# Bad: Don't use None unless appropriate
result = handler.handle_validation(
    validate_property,
    value,
    component="Property Validator",
    operation="validate_logP",
    error_strategy=RecoveryStrategy.SKIP,
    default_value=None  # May cause problems downstream
)
```

### 3. Use Fallback Functions for Complex Transformations

For complex data recovery, implement custom fallback functions:

```python
def sanitize_smiles(smiles):
    """Attempt to sanitize invalid SMILES."""
    if not isinstance(smiles, str):
        return "C"
    
    # Remove unsupported characters
    smiles = re.sub(r'[^A-Za-z0-9\(\)\[\]\=\#\-\+]', '', smiles)
    
    # Add proper valence if needed
    if smiles and not any(c in smiles for c in "()[]"):
        smiles = "C"
    
    return smiles

result = handler.handle_validation(
    validate_smiles,
    "C@1CC!!1",  # Invalid SMILES
    component="SMILES Validator",
    operation="validate_smiles",
    error_strategy=RecoveryStrategy.FALLBACK,
    fallback_func=sanitize_smiles
)
```

### 4. Use Batch Validation for Performance

For validating multiple items, use the batch validation with parallel processing:

```python
# Efficient: Use batch validation with parallel processing
results, report = handler.validate_batch(
    validation_func,
    large_dataset,
    component="Batch Validator",
    operation="validate_batch",
    error_strategy=RecoveryStrategy.SKIP,
    parallel=True,
    max_workers=8
)

# Inefficient: Don't loop manually for large datasets
results = []
for item in large_dataset:
    result = handler.handle_validation(
        validation_func,
        item,
        component="Single Validator",
        operation="validate_item",
        error_strategy=RecoveryStrategy.SKIP
    )
    results.append(result)
```

### 5. Analyze Validation Reports

Use the validation reports to understand and improve data quality:

```python
# Get summary statistics
summary = report.get_summary()
print(f"Valid: {summary['valid_items']}, Invalid: {summary['invalid_items']}")
print(f"Error types: {summary['error_codes']}")
print(f"Fields with errors: {summary['fields_with_errors']}")

# Analyze errors by severity
critical_errors = report.by_severity.get("CRITICAL", [])
high_errors = report.by_severity.get("HIGH", [])
medium_errors = report.by_severity.get("MEDIUM", [])

print(f"Critical errors: {len(critical_errors)}")
print(f"High severity errors: {len(high_errors)}")
print(f"Medium severity errors: {len(medium_errors)}")

# Analyze errors by field
smiles_errors = report.by_field.get("smiles", [])
name_errors = report.by_field.get("name", [])

print(f"SMILES errors: {len(smiles_errors)}")
print(f"Name errors: {len(name_errors)}")

# Save for further analysis
report.save_to_file("validation_analysis.json")
```

## Conclusion

The Validation Error Handling system provides a comprehensive framework for validating data, handling validation errors with configurable strategies, and generating detailed validation reports. It integrates seamlessly with the broader error handling framework and supports both individual and batch validation operations.

By using this system consistently throughout the application, we ensure robust handling of invalid data, maintain data integrity, and capture detailed information for debugging and analysis.