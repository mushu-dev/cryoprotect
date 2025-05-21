# Property Data Standardization Guide

This guide explains the property data standardization implementation for the CryoProtect application, including the schema changes, data migration strategy, and usage examples.

## Overview

The property standardization system ensures consistent formatting and units for molecular properties in the database. It provides:

1. Standardized property types and units
2. Automatic value validation and data quality enforcement
3. Unit conversion capabilities
4. Text-to-numeric value extraction with unit detection
5. A standardized view for property access

## Schema Changes

The standardization implementation adds the following tables and columns:

### New Tables

1. **property_types**: Reference table defining standard property types
   - Includes value type, allowed ranges, validation rules, and default units
   - Ensures consistent property naming and validation

2. **units**: Reference table defining standard units of measurement
   - Includes categories, symbols, and conversion factors
   - Enables unit conversion and standardization

### Column Additions

The following columns were added to the `molecular_properties` table:

- `property_type_id`: Links to a standardized property type
- `unit_id`: Links to a standardized unit of measurement

## Key Functions

### Data Quality Functions

- `validate_property_value()`: Validates property values against their defined constraints
- `standardize_property_name()`: Normalizes property names to a consistent format
- `convert_unit()`: Converts values between compatible units
- `extract_unit_from_text()`: Extracts numeric values and units from text fields

### Data Access Functions

- `add_molecular_property()`: Adds a property with validation and automatic type creation
- `standardized_properties` view: Provides a unified view of all properties with their types and units

## Migration Process

The property standardization process migrates existing data to the new schema:

1. **Analysis Phase**:
   - Identifies existing property names
   - Matches properties to standardized types
   - Detects units in text values
   - Creates a migration plan

2. **Migration Phase**:
   - Creates missing property types
   - Links properties to their standardized types
   - Converts text values with units to numeric values
   - Updates unit references

## Usage Examples

### Adding a New Property

```sql
-- Add a new molecular property with automatic validation
SELECT public.add_molecular_property(
    molecule_id := '123e4567-e89b-12d3-a456-426614174000',
    property_name := 'melting_point',
    numeric_value := 78.5,
    unit_name := 'celsius'
);
```

### Accessing Standardized Properties

```sql
-- Query all properties with standardized names and units
SELECT 
    molecule_id,
    property_name,
    numeric_value,
    unit_symbol
FROM 
    standardized_properties
WHERE 
    value_type = 'numeric'
ORDER BY
    property_name;
```

### Converting Units

```sql
-- Convert temperature from Celsius to Fahrenheit
SELECT
    public.convert_unit(
        value := 25.0,
        from_unit_id := (SELECT id FROM units WHERE name = 'celsius'),
        to_unit_id := (SELECT id FROM units WHERE name = 'fahrenheit')
    ) AS temperature_fahrenheit;
```

## Running the Standardization Process

The standardization process can be run using the provided script:

```bash
# Run the standardization with dry-run (no changes)
./run_property_standardization.sh --dry-run

# Analyze properties without migration
./run_property_standardization.sh --analyze-only

# Run the full standardization
./run_property_standardization.sh
```

### Offline Mode

If the database is not accessible, you can use the offline mode script for testing and development:

```bash
# Run a simulation in dry-run mode
./run_property_standardization_offline.sh --dry-run

# Run analysis-only simulation
./run_property_standardization_offline.sh --analyze-only

# Run full standardization simulation
./run_property_standardization_offline.sh
```

The offline script generates sample reports and simulates the standardization process without requiring a database connection. This is useful for:

- Development without database access
- Testing the migration process structure
- Reviewing the expected output format
- Training new team members

### Prerequisites

Before running the standardization process, ensure:

1. Database connection is properly configured in `.env` or environment variables:
   ```
   DB_HOST=<your_db_host>
   DB_PORT=<your_db_port>
   DB_NAME=<your_db_name>
   DB_USER=<your_db_user>
   DB_PASSWORD=<your_db_password>
   ```

2. The user has sufficient privileges to:
   - Create new tables
   - Modify existing tables
   - Create functions and triggers
   - Update data

3. The `molecular_properties` table exists and contains data to standardize

## Testing the Migration

It's recommended to test the migration in a development environment before applying to production:

1. Create a database backup:
   ```bash
   ./create_production_backup.sh
   ```

2. Restore the backup to a test environment:
   ```bash
   ./restore_backup.py --source <backup_file> --target <test_db>
   ```

3. Run the migration in dry-run mode:
   ```bash
   ./run_property_standardization.sh --dry-run
   ```

4. Review the generated report and fix any issues

5. Run the full migration on the test database and verify results

## Monitoring and Verification

After running the standardization process, verify the results:

1. Check the generated report file for details on the migration
2. Query the `property_types` and `units` tables to see the standardized references
3. Use the `standardized_properties` view to verify data quality
4. Run data quality checks to ensure all properties are linked to property types

```sql
-- Check for properties without property_type_id
SELECT property_name, COUNT(*) 
FROM molecular_properties 
WHERE property_type_id IS NULL 
GROUP BY property_name;

-- Verify text-to-numeric conversions
SELECT property_name, COUNT(*) 
FROM molecular_properties 
WHERE unit_id IS NOT NULL 
GROUP BY property_name;
```

## Troubleshooting

### Common Issues

If issues occur during standardization:

1. Review the error log in the generated report
2. Check for incompatible unit conversions
3. Verify property type mappings for unmapped properties
4. Use the `--analyze-only` flag to diagnose issues without making changes
5. Check database permissions and connection settings
6. Examine the SQL migration file for syntax errors

### Database Connection Problems

If you encounter database connection issues:

1. **Connection Timeout**:
   - Check if the database server is accessible from your network
   - Verify firewall settings aren't blocking PostgreSQL (port 5432)
   - Try connecting with `psql` to isolate connection issues

2. **Authentication Failures**:
   - Verify your database credentials in the `.env` file
   - Check if the database user has the necessary permissions
   - For "Ident authentication failed" errors, check PostgreSQL's `pg_hba.conf` authentication method

3. **SSL/TLS Issues**:
   - Ensure `USE_SSL=true` when connecting to cloud databases like Supabase
   - Check if your PostgreSQL client has SSL support enabled
   - Verify your SSL certificates if using client certificate authentication

4. **Offline Work Option**:
   - If database connectivity can't be resolved immediately, use the offline mode:
     ```bash
     ./run_property_standardization_offline.sh
     ```
   - Continue development and apply changes when connectivity is restored

### Testing the Database Connection

Use the provided test script to diagnose connection issues:

```bash
python3 test_db_connection.py
```

This will show detailed configuration information and attempt to connect to the database.