# Database Verification Module

This module provides tools for verifying database integrity, schema consistency, and data quality for the CryoProtect v2 project.

## Usage

### Command Line Interface

The verification module provides a command line interface for running verification checks:

```bash
# Run standard verification
python -m database.verification

# Run comprehensive verification
python -m database.verification --level comprehensive

# Run specific verification modules
python -m database.verification --include schema data_quality

# Generate HTML report
python -m database.verification --format html --output report.html
```

### Programmatic Usage

The verification module can also be used programmatically:

```python
from database.verification import verify_database, verify_schema, verify_constraints, verify_data_quality

# Run full database verification
results = verify_database(level='standard')

# Run specific verification
schema_results = verify_schema()
constraint_results = verify_constraints()
data_quality_results = verify_data_quality()

# Check success
if results['success']:
    print("Verification passed!")
else:
    print("Verification failed!")
    for module, module_results in results['results'].items():
        for issue in module_results.get('issues', []):
            print(f"{module}: {issue['message']}")
```

## Verification Levels

The verification module supports three verification levels:

- **basic**: Fast, essential checks only
- **standard**: Default level, balanced thoroughness and performance
- **comprehensive**: Thorough checks, may be slower

## Report Formats

Verification results can be output in several formats:

- **console**: Print results to the console
- **json**: Generate a JSON report
- **html**: Generate an HTML report

## Adding Custom Verifications

To add a custom verification, create a new function in the appropriate module and add it to the verification runner.

```python
# In database/verification/data.py
def _check_custom_requirement(conn):
    issues = []
    # Implement custom check
    return issues

# In database/verification/data.py
def run_data_quality_verification(conn, level):
    issues = []
    # Add custom check to appropriate level
    if level == 'comprehensive':
        issues.extend(_check_custom_requirement(conn))
    # ...