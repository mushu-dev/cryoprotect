# CryoProtect Configuration Management System

This directory contains the unified configuration management system for the CryoProtect application. It provides a single source of truth for all configuration parameters across both backend and frontend components.

## Overview

The configuration system is designed around a central schema definition (`schema.json`) that describes all configuration parameters, their types, default values, and environment-specific overrides. From this schema, language-specific configuration files are generated for both the Python backend and TypeScript frontend.

### Key Features

- **Single source of truth**: All configuration parameters are defined in one place
- **Type safety**: All parameters have defined types that are enforced in both Python and TypeScript
- **Environment-specific overrides**: Different environments (development, testing, staging, production) can have different default values
- **Hierarchical structure**: Configuration is organized in logical sections
- **Validation**: Both runtime and build-time validation to catch configuration errors early
- **Secure handling of secrets**: Sensitive values are properly marked and handled securely

## Directory Structure

- `schema.json` - The central configuration schema definition
- `generate_configs.py` - Main script to generate all language-specific configuration files
- `generate_python_config.py` - Script to generate Python configuration
- `generate_typescript_config.js` - Script to generate TypeScript configuration
- `consolidate_backend_config.py` - Script to merge backend configuration files
- `README.md` - This documentation file

## Getting Started

### Prerequisites

- Python 3.7+
- Node.js 14+ (for TypeScript generation)

### Generating Configuration Files

To generate all configuration files from the schema:

```bash
./generate_configs.py
```

This will:
1. Validate the schema for correctness
2. Generate a new Python configuration file at `../config.py.new`
3. Generate a new TypeScript configuration file at `../frontend/src/config/config.ts`

### Customizing the Generation

You can specify custom paths for the input schema and output files:

```bash
./generate_configs.py --schema custom_schema.json --python ../my_config.py --typescript ../frontend/src/my_config.ts
```

To only validate the schema without generating files:

```bash
./generate_configs.py --validate-only
```

### Consolidating Backend Configuration Files

If you have multiple backend configuration files (e.g., `config.py`, `db_config.py`), you can consolidate them:

```bash
./consolidate_backend_config.py --main ../config.py --secondary ../frontend/db_config.py ../another_config.py
```

This will analyze the configurations and generate a consolidated file at `../config.consolidated.py`.

## Schema Format

The schema follows a hierarchical structure:

```json
{
  "properties": {
    "section_name": {
      "type": "object",
      "description": "Section description",
      "properties": {
        "property_name": {
          "type": "string",
          "description": "Property description",
          "default": "default value"
        }
      }
    }
  },
  "environmentMapping": {
    "development": {
      "section_name.property_name": "development value"
    },
    "production": {
      "section_name.property_name": "production value"
    }
  }
}
```

### Top-level Sections

- `properties`: Contains all configuration sections
- `environmentMapping`: Defines environment-specific overrides using dot notation

### Property Attributes

Each property can have the following attributes:

- `type`: Data type (string, integer, number, boolean, object, array)
- `description`: Human-readable description of the property
- `default`: Default value
- `enum`: List of allowed values (for string types)
- `format`: Format hint (e.g., "uri", "password")
- `examples`: Example values

## Using the Configuration

### In Python (Backend)

```python
from config import active_config

# Access configuration values
database_url = active_config.DATABASE_URL
debug_mode = active_config.DEBUG

# Or use the convenience exports
from config import DATABASE_URL, DEBUG
```

### In TypeScript (Frontend)

```typescript
import config from 'src/config/config';

// Access configuration values
const apiUrl = config.api.baseUrl;
const useMockData = config.app.featureFlags.useMockData;
```

## Best Practices

1. **Always update the schema first**: When adding new configuration parameters, add them to `schema.json` first
2. **Regenerate configuration files**: After updating the schema, regenerate all configuration files
3. **Commit the schema and generated files**: Keep the schema and generated files in version control
4. **Use environment variables for sensitive values**: Never hardcode sensitive values in the schema
5. **Document all parameters**: Always add a description for each parameter in the schema
6. **Validate configuration at startup**: Use the validation functions to catch configuration errors early

## Extending the Schema

When adding new configuration parameters:

1. Add the parameter to the appropriate section in `schema.json`
2. Add environment-specific overrides if needed
3. Regenerate configuration files with `./generate_configs.py`
4. Update any code that depends on the new parameter

## Troubleshooting

### Validation Errors

If validation fails, check the error message for details. Common issues include:

- Missing required properties
- Type mismatches
- Invalid default values

### Generation Failures

If generation fails:

- Ensure the schema is valid
- Check that Python and Node.js are installed and available
- Verify that all paths exist and are writable
- Check for syntax errors in the generators

## Contributing

When contributing to the configuration system:

1. Keep the schema clean and well-organized
2. Document all changes to the schema
3. Test generation and validation after making changes
4. Update this documentation if you add new features or change behavior