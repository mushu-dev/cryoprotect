# Heroku Deployment Guide for CryoProtect

This guide documents the deployment of the CryoProtect application to Heroku. The deployment uses a simplified version of the application due to the limitations of Heroku's environment (specifically the lack of RDKit support).

## Deployed Application

The application is deployed at: https://cryoprotect-8030e4025428.herokuapp.com/

## Architecture Overview

The Heroku deployment consists of:

1. **Database**: PostgreSQL database provisioned by Heroku
2. **API Server**: A simplified Flask application that provides basic API endpoints
3. **Release Phase**: Scripts for database setup and sample data population

## Deployment Components

### Database Setup

The database is initialized during the Heroku release phase using the `setup_database.py` script. This script:

1. Connects to the PostgreSQL database provided by Heroku
2. Creates the necessary tables if they don't exist
3. Sets up any required indexes

### Sample Data Population

Sample data is loaded into the database during the release phase using the `populate_sample_data.py` script. This includes:

1. Cryoprotectant molecules (Glycerol, DMSO, etc.)
2. Mixtures of cryoprotectants
3. Molecular properties
4. Property types

### Simplified API

Since RDKit is not available on Heroku by default, we created a simplified API version `simple_app.py` that provides the following endpoints:

- `/` - Health check and basic information
- `/health` - Detailed health check including database connection status
- `/api/molecules` - List molecules with pagination
- `/api/molecules/<id>` - Get details of a specific molecule

## Deployment Steps

1. Created database setup scripts compatible with Heroku PostgreSQL
2. Created a simplified application without RDKit dependencies
3. Updated the Procfile to specify both release phase commands and web server
4. Added sample data population script for initial data

## Mock RDKit Implementation

As RDKit is not available on Heroku, we created a mock implementation (`mock_rdkit.py`) that:

1. Provides stub implementations of RDKit classes and functions
2. Returns reasonable default values for molecular properties
3. Allows the code to run without the actual RDKit library

This approach allows the deployment of a minimal version of the application without requiring the complex setup needed for RDKit.

## Limitations

The Heroku deployment has the following limitations:

1. No RDKit-based molecular property calculations
2. No molecular visualization features
3. Limited to basic CRUD operations on molecules and mixtures
4. No complex analysis or prediction features

## Future Improvements

For a full-featured deployment, consider:

1. Using Heroku container deployment with a custom container that includes RDKit
2. Setting up a separate RDKit service that the Heroku app can call
3. Migrating to a platform that supports more complex scientific computation (AWS, GCP, etc.)

## API Documentation

### GET /api/molecules

Returns a list of molecules with pagination.

**Query Parameters:**
- `limit` (optional): Maximum number of results to return (default: 10)
- `offset` (optional): Number of results to skip (default: 0)

**Example Response:**
```json
{
  "status": "success",
  "data": [
    {
      "id": 1,
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "pubchem_cid": "753",
      "molecular_formula": "C3H8O3",
      "molecular_weight": 92.09
    },
    // More molecules...
  ],
  "pagination": {
    "total": 10,
    "limit": 5,
    "offset": 0
  }
}
```

### GET /api/molecules/{id}

Returns details of a specific molecule.

**Parameters:**
- `id`: The ID of the molecule to retrieve

**Example Response:**
```json
{
  "status": "success",
  "data": {
    "id": 1,
    "name": "Glycerol",
    "smiles": "C(C(CO)O)O",
    "pubchem_cid": "753",
    "molecular_formula": "C3H8O3",
    "molecular_weight": 92.09
  }
}
```