# CryoProtect v2 Developer Guide

This guide provides an overview of the code structure, logging, error handling, testing, and best practices for maintaining and extending the CryoProtect Analyzer project.

---

## Project Structure

- **app.py**: Main Flask application entry point
- **api/**: Core backend modules (REST API endpoints, business logic)
  - `scoring.py`, `predictive_models.py`, `rdkit_utils.py`, etc.
- **static/**: Frontend static files (JavaScript, CSS)
- **templates/**: HTML templates for the web interface
- **migrations/**: Database migration scripts and schema documentation
- **tests/**: Comprehensive test suite (unit, integration, schema)
- **memory-bank/**: Persistent storage for RooFlow memory and cached data
- **config.py**, **config_staging.py**, **config_production.py**: Environment-specific configuration
- **logging_config.py**: Logging setup and configuration

---

## Logging

- Logging is configured in `logging_config.py` using Python's `logging` module.
- Log messages are written to both the console and log files.
- Use appropriate log levels (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`) throughout the codebase.
- Example usage:
  ```python
  import logging
  logger = logging.getLogger(__name__)
  logger.info("This is an info message")
  logger.error("This is an error message")
  ```

---

## Error Handling

- Use try/except blocks to handle exceptions gracefully.
- Return informative error messages in API responses.
- For Flask routes and API endpoints, use Flask's error handlers to provide consistent error responses.
- Log all exceptions using the logging system for traceability.

---

## Testing

- All core modules and endpoints are covered by tests in the `tests/` directory.
- See `tests/README.md` for details on running and contributing to tests.
- Use the test runner script (`python tests/run_tests.py`) to execute all tests.
- Add new tests for any new features, bug fixes, or schema changes.

---

## Best Practices

- Follow PEP 8 for Python code style.
- Document all functions and classes with clear docstrings.
- Use type hints where appropriate for clarity.
- Keep modules focused and maintainable; avoid large monolithic files.
- Use environment variables and configuration files for secrets and deployment-specific settings.
- Regularly run tests and validation scripts to ensure data and code integrity.
- When making database changes, update migration scripts and schema documentation.

---

## Contributing

- When adding new features:
  - Update or add relevant tests.
  - Update documentation as needed.
  - Follow the established code structure and naming conventions.
- For bug fixes:
  - Add regression tests to prevent recurrence.
  - Document the fix in the memory-bank/issue_fixes.json file if relevant.

---

## Additional Resources

- [API Documentation](README_API.md)
- [Database Schema](README_Database_Schema.md)
- [Scoring System](README_Scoring.md)
- [Predictive Models](README_Predictive_Models.md)
- [RDKit Integration](README_RDKit.md)
- [Testing Guide](tests/README.md)