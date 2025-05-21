# CryoProtect Minimal Test Guide

This guide explains how to use the minimal test scripts to verify basic functionality of your CryoProtect application on Fedora Linux.

## Overview

Since the full application has many dependencies and containerization introduces SELinux challenges, these minimal scripts help you verify the core functionality:

1. Flask server running
2. Supabase connectivity
3. Environment configuration
4. Python dependencies

## Quick Start

### Run the Minimal App

```bash
./run_minimal.sh
```

Or use the simplified app:

```bash
./run_simplified_app.sh
```

These scripts:
- Create a dedicated Python virtual environment
- Install only the essential dependencies (Flask, requests, python-dotenv)
- Run a simplified Flask application with diagnostic endpoints

### Verify Functionality

Once running, open these URLs in your browser:

1. `http://localhost:5000/` - Basic verification that Flask is running
2. `http://localhost:5000/health` - Simple health check
3. `http://localhost:5000/env` - View environment variables
4. `http://localhost:5000/test-supabase` - Test Supabase connectivity
5. `http://localhost:5000/dependencies` - Check installed Python dependencies
6. `http://localhost:5000/api/sample/cryoprotectants` - View sample data

The `/test-supabase` endpoint is especially important as it confirms that your application can connect to your Supabase instance.

## Features

The simplified app provides:
- Basic Flask web server
- Supabase connectivity testing
- Comprehensive dependency checking
- Sample cryoprotectant data endpoint

## Database Structure

The CryoProtect Supabase database contains several important tables:

- `molecules` - 733 molecules with their chemical identifiers
- `molecular_properties` - 741 property values for molecules
- `property_types` - Definitions of property types
- `mixtures` - 8 cryoprotectant mixtures 
- `mixture_components` - Components of mixtures with concentrations

See `DATABASE_SCHEMA_SUMMARY.md` for more details on the database structure.

## Example VS55 Vitrification Solution

The database includes information about the VS55 Vitrification Solution, which consists of:
- Propylene glycol (concentration: 2.2)
- Formamide (concentration: 1.4)
- Dimethyl sulfoxide (concentration: 8.4)

## Troubleshooting

### Supabase Connection Failures

If the Supabase connection test fails:

1. Check your `.env` file:
   ```
   SUPABASE_URL=https://your-project.supabase.co
   SUPABASE_KEY=your-anon-key
   ```

2. Try the connection test from outside containers:
   ```
   python scripts/test_db_connection_simple.py
   ```

3. Check DNS resolution:
   ```
   host your-project.supabase.co
   ```

### SELinux Issues

If you encounter SELinux issues when trying to run containers:

1. Try disabling SELinux for testing purposes (not recommended for production):
   ```
   sudo setenforce 0  # Temporarily disable
   ```

2. Use the SELinux context setup script:
   ```
   sudo ./setup_selinux_sudo.sh
   ```

### Dependency Issues

If you encounter missing dependencies:

1. Check the `/dependencies` endpoint to see what's missing
2. Install the missing dependencies:
   ```bash
   source quick_env/bin/activate
   pip install missing_package_name
   ```

3. For RDKit-specific issues, use:
   ```bash
   conda install -c conda-forge rdkit
   ```

## Scripts

- `run_minimal.sh` - Runs the original minimal app
- `run_simplified_app.sh` - Runs the enhanced simplified app
- `run_simplified_app_background.sh` - Runs the app in the background
- `setup_selinux_sudo.sh` - Sets up SELinux contexts
- `fedora_setup.sh` - Complete Fedora environment setup

## Full Application Dependencies

To analyze all the dependencies required by the full app:

```bash
python scripts/debug_dependencies.py
```

This will scan app.py and report all the imported modules.

## Going Forward

Once the minimal app confirms your Supabase connectivity, you can proceed by:

1. Gradually adding more dependencies to the minimal version
2. Running the application directly with Python without containers 
3. Using Podman with SELinux configurations once basic functionality is confirmed

See `FEDORA_SETUP_SUMMARY.md` for more information on the next steps.