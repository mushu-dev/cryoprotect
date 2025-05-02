# Running Fix Scripts

This guide provides detailed instructions for running the CryoProtect fix scripts to address various issues with the database schema, security implementation, relationship design, and API integration.

## Overview of Fix Scripts

The CryoProtect project includes several fix scripts that address specific issues:

1. **Database Backup** (`create_database_backup.py`): Creates a backup of the database before making changes
2. **Schema Standardization** (`standardize_schema.py`): Standardizes table names and ensures proper primary/foreign keys
3. **Security Implementation** (`implement_security.py`): Implements Row Level Security and app-specific roles
4. **Relationship Fixes** (`fix_relationships.py`): Fixes relationship design issues and adds proper junction tables
5. **API Integration Fixes** (`fix_api_integration.py`): Fixes API endpoint registration and JSON serialization issues
6. **Testing** (`tests/run_tests.py`): Runs tests to verify that all fixes have been applied correctly

These scripts can be run individually or together using the master integration script.

## Master Integration Script

The master integration script (`run_cryoprotect_fixes.py`) orchestrates the execution of all fix scripts in the correct order. This is the recommended approach for applying all fixes at once.

### Running the Master Script

To run the master integration script:

```bash
# Windows
python run_cryoprotect_fixes.py [options]

# Unix/Linux/macOS
python3 run_cryoprotect_fixes.py [options]
```

### Command-line Options

The master script supports several command-line options:

| Option | Description |
|--------|-------------|
| `--all` | Run all fixes in sequence (default) |
| `--backup-only` | Run only database backup |
| `--schema-only` | Run only schema standardization |
| `--security-only` | Run only security implementation |
| `--relationships-only` | Run only relationship fixes |
| `--api-only` | Run only API integration fixes |
| `--test-only` | Run only tests |
| `--verify` | Verify the state after each fix |
| `--rollback` | Rollback changes if a fix fails |
| `--dry-run` | Show what would be done without making changes |
| `--report` | Generate a comprehensive final report |
| `--help` | Show help message and exit |

### Examples

Run all fixes with verification:
```bash
python run_cryoprotect_fixes.py --verify
```

Run only schema standardization:
```bash
python run_cryoprotect_fixes.py --schema-only
```

Run in dry-run mode to see what would be done:
```bash
python run_cryoprotect_fixes.py --dry-run
```

## Running Individual Fix Scripts

You can also run each fix script individually if you need to apply specific fixes.

### 1. Database Backup

```bash
python create_database_backup.py [--dry-run]
```

This script creates a backup of the database before making any changes. The backup is stored in the `backups` directory with a timestamp.

### 2. Schema Standardization

```bash
python standardize_schema.py [--dry-run] [--verify] [--rollback]
```

This script standardizes table names (singular to plural), ensures proper UUID primary keys, and adds appropriate foreign key constraints.

### 3. Security Implementation

```bash
python implement_security.py [--dry-run] [--verify-only] [--rollback]
```

This script implements Row Level Security (RLS) on all public tables, creates RLS policies, and sets up app-specific database roles.

### 4. Relationship Fixes

```bash
python fix_relationships.py [--dry-run] [--verify-only] [--rollback]
```

This script fixes relationship design issues by replacing fan traps with proper junction tables, applying 3NF normalization, and adding missing foreign key constraints.

### 5. API Integration Fixes

```bash
python fix_api_integration.py [--dry-run]
```

This script fixes API endpoint registration issues, database configuration issues, and JSON serialization problems.

### 6. Testing

```bash
python tests/run_tests.py [--report] [--verbose]
```

This script runs tests to verify that all fixes have been applied correctly.

## Batch Files and Shell Scripts

For convenience, batch files and shell scripts are provided for common operations:

### Windows Batch Files

- `run_cryoprotect_fixes.bat`: Runs the master integration script
- `run_database_backup.bat`: Runs the database backup script
- `run_schema_standardization.bat`: Runs the schema standardization script
- `run_api_fixes.bat`: Runs the API integration fixes script
- `run_performance_tests.bat`: Runs performance tests

### Unix/Linux/macOS Shell Scripts

- `run_cryoprotect_fixes.sh`: Runs the master integration script
- `run_schema_standardization.sh`: Runs the schema standardization script
- `run_api_fixes.sh`: Runs the API integration fixes script
- `run_performance_tests.sh`: Runs performance tests

## Environment Setup

Before running the fix scripts, ensure that your environment is properly set up:

1. Create a `.env` file with your Supabase credentials:
   ```
   SUPABASE_URL=your-supabase-url
   SUPABASE_KEY=your-supabase-key
   SUPABASE_USER=your-supabase-user
   SUPABASE_PASSWORD=your-supabase-password
   ```

2. Install required Python packages:
   ```bash
   pip install -r requirements.txt
   ```

3. Ensure you have appropriate permissions to modify the Supabase database.

## Logging

All fix scripts generate detailed logs that can be used to troubleshoot issues:

- Master integration script: `cryoprotect_fixes_YYYYMMDD_HHMMSS.log`
- Schema standardization: `schema_migration_YYYYMMDD_HHMMSS.log`
- Security implementation: `security_implementation.log`
- Relationship fixes: `fix_relationships.log`
- API integration fixes: `api_fixes.log`

These logs are stored in the project root directory and contain information about the operations performed, any errors encountered, and the overall status of the fixes.

## Best Practices

When running the fix scripts, follow these best practices:

1. **Always create a backup** before making changes
2. **Run in dry-run mode first** to see what would be done
3. **Run in a test environment** before applying to production
4. **Verify changes** after running each script
5. **Check logs** for any errors or warnings
6. **Have a rollback plan** in case something goes wrong

## Troubleshooting

If you encounter issues when running the fix scripts:

1. Check the log files for error messages
2. Verify that your Supabase credentials are correct
3. Ensure you have appropriate permissions
4. Try running in dry-run mode to identify potential issues
5. Use the `--verify` option to check the state after each fix
6. If needed, use the `--rollback` option to revert changes

For more detailed troubleshooting information, see the [Troubleshooting Guide](../appendix/troubleshooting-guide.md).