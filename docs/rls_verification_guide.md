# RLS Verification and Remediation Guide

This guide explains how to use the RLS (Row Level Security) verification and remediation utilities provided in the `rls_utils.py` module. These utilities ensure that RLS is always enabled and all required policies are present for the `molecule` and `molecular_property` tables.

## Overview

The RLS verification and remediation utilities provide the following functionality:

1. **Verification**: Check if RLS is enabled and all required policies are present for the `molecule` and `molecular_property` tables.
2. **Remediation**: Automatically fix any issues found during verification.
3. **Decorator**: A Python decorator that automatically verifies and remediates RLS before and after a function call.
4. **Audit Logging**: All verification and remediation actions are logged to a dedicated audit log.

## Usage

### Verifying RLS Settings

To verify RLS settings for the `molecule` and `molecular_property` tables:

```python
from rls_utils import verify_rls

# Verify RLS settings
issues = verify_rls()

if not issues:
    print("No RLS issues found.")
else:
    print(f"RLS issues found: {issues}")
```

The `verify_rls()` function returns a dictionary of issues found for each table. If no issues are found, an empty dictionary is returned.

### Remediating RLS Issues

To remediate RLS issues:

```python
from rls_utils import verify_rls, remediate_rls

# Verify RLS settings
issues = verify_rls()

if issues:
    # Remediate issues
    success = remediate_rls(issues)
    if success:
        print("All RLS issues remediated successfully.")
    else:
        print("Some RLS issues could not be remediated.")
```

The `remediate_rls()` function takes a dictionary of issues (as returned by `verify_rls()`) and attempts to fix them. It returns `True` if all issues were fixed, `False` otherwise.

### Using the Decorator

The `@ensure_rls_restored` decorator can be applied to any function that might affect RLS settings. It automatically verifies and remediates RLS before and after the function call, even if an exception occurs:

```python
from rls_utils import ensure_rls_restored

@ensure_rls_restored
def import_data_to_database():
    # This function might affect RLS settings
    # ...
```

The decorator will:

1. Verify and remediate RLS before the function call.
2. Call the function.
3. Verify and remediate RLS after the function call.
4. If an exception occurs, verify and remediate RLS before re-raising the exception.

### Audit Logging

All verification and remediation actions are logged to a dedicated audit log at `logs/rls_audit.jsonl`. Each log entry includes:

- Timestamp
- Action (verify or remediate)
- Table
- Details of the action
- User
- Script

To view the audit log:

```bash
cat logs/rls_audit.jsonl | jq
```

## RLS Policies

The following RLS policies are verified and maintained for the `molecule` and `molecular_property` tables:

### Molecule Table

1. **Select molecules for project members**: Allows project members to view molecules in their projects.
2. **Insert molecules for project members**: Allows project members to add molecules to their projects.
3. **Update molecules for project members**: Allows project members to update molecules in their projects.
4. **Delete molecules for project members**: Allows project members to delete molecules in their projects.
5. **Allow service role inserts on molecule**: Allows service role to insert molecules regardless of project membership.

### Molecular Property Table

1. **Select molecular_properties for project members**: Allows project members to view molecular properties in their projects.
2. **Insert molecular_properties for project members**: Allows project members to add molecular properties to their projects.
3. **Update molecular_properties for project members**: Allows project members to update molecular properties in their projects.
4. **Delete molecular_properties for project members**: Allows project members to delete molecular properties in their projects.
5. **Allow service role inserts on molecular_property**: Allows service role to insert molecular properties regardless of project membership.

## Testing

A comprehensive test suite is provided in `test_rls_utils.py`. To run the tests:

```bash
python test_rls_utils.py
```

The tests verify that:

1. RLS verification correctly identifies issues.
2. RLS remediation correctly fixes issues.
3. The decorator correctly verifies and remediates RLS before and after function calls.
4. Audit logging works correctly.

## Integration with ChEMBL Import

The `@ensure_rls_restored` decorator has been applied to the `import_compounds_to_database()` function in `ChEMBL_Integrated_Import.py`. This ensures that RLS is always verified and restored before and after importing compounds from ChEMBL.

## Troubleshooting

If you encounter issues with RLS verification or remediation:

1. Check the audit log at `logs/rls_audit.jsonl` for details of the verification and remediation actions.
2. Run the tests in `test_rls_utils.py` to verify that the utilities are working correctly.
3. Manually verify RLS settings using the `verify_rls()` function.
4. If necessary, manually remediate RLS issues using the `remediate_rls()` function.

## References

- [PostgreSQL Row Level Security](https://www.postgresql.org/docs/current/ddl-rowsecurity.html)
- [Supabase Row Level Security](https://supabase.com/docs/guides/auth/row-level-security)
- [RLS Restoration Protocol](../.specs/rls_restoration_protocol.md)