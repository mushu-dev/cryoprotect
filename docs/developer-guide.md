# Developer Guide

This Developer Guide provides information for developers working with the CryoProtect system, including how to extend the system, best practices for future development, and API usage examples.

## Overview

The CryoProtect system has been significantly improved with fixes to the database schema, security implementation, relationship design, and API integration. This guide will help developers understand how to:

1. Extend the system with new features
2. Follow best practices for future development
3. Use the API effectively

## Table of Contents

1. [Extending the System](./developer/extending-the-system.md)
   - Adding new tables
   - Extending existing functionality
   - Creating new API endpoints
   - Adding new molecular properties
   - Integrating with external systems

2. [Best Practices](./developer/best-practices.md)
   - Database design principles
   - Security considerations
   - API design guidelines
   - Testing strategies
   - Code organization
   - Documentation standards

3. [API Usage Examples](./developer/api-usage-examples.md)
   - Authentication
   - Working with molecules
   - Working with mixtures
   - Running predictions
   - Managing experiments
   - Error handling

## System Architecture

The CryoProtect system follows a modular architecture with the following components:

### Backend Components

- **Flask Application**: The main web application that serves the UI and API
- **Database Layer**: Supabase PostgreSQL database with standardized schema
- **Authentication**: Supabase Auth for user authentication and authorization
- **Molecular Analysis**: RDKit for molecular property calculations and visualization
- **API Layer**: RESTful API for accessing and manipulating data

### Frontend Components

- **Web Interface**: HTML, CSS, and JavaScript with Bootstrap 5
- **Visualization**: Molecular visualization using RDKit and JavaScript libraries
- **Forms**: User input forms for data entry and search
- **Results Display**: Tables and charts for displaying results

### Integration Points

- **Supabase**: Database and authentication
- **RDKit**: Molecular analysis
- **PubChem**: Optional integration for retrieving molecular data

## Development Environment Setup

To set up a development environment for CryoProtect:

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/cryoprotect-analyzer.git
   cd cryoprotect-analyzer
   ```

2. Set up the environment:
   ```bash
   # Windows
   setup_environment.bat
   
   # Unix/Linux/macOS
   ./setup_environment.sh
   ```

3. Configure Supabase:
   - Create a `.env` file with your Supabase credentials
   - Apply the database migration

4. Run the application in development mode:
   ```bash
   # Windows
   run_app.bat
   
   # Unix/Linux/macOS
   ./run_app.sh
   ```

## Development Workflow

The recommended development workflow for CryoProtect is:

1. **Create a Feature Branch**: Create a branch for your feature or fix
2. **Implement Changes**: Make your changes to the code
3. **Run Tests**: Ensure all tests pass
4. **Verify Changes**: Verify that your changes work as expected
5. **Document Changes**: Update documentation to reflect your changes
6. **Submit Pull Request**: Submit a pull request for review

## Key Files and Directories

- `app.py`: Main Flask application
- `api/`: API endpoints and resources
- `static/`: Static files (CSS, JavaScript, images)
- `templates/`: HTML templates
- `models/`: Data models and database interaction
- `migrations/`: Database migration scripts
- `tests/`: Test scripts
- `docs/`: Documentation

For detailed information on extending the system, best practices, and API usage, please refer to the specific sections linked in the Table of Contents.