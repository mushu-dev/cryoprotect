# Phase 3.1: Database Population Guide

This guide details how to populate the Supabase database with scientifically accurate cryoprotectant data as part of the deployment infrastructure implementation.

## Database Population Overview

The CryoProtect v2 application requires a populated database with scientific data to function properly. The `populate_database_supabase.py` script provides a comprehensive solution for populating all required tables with scientifically accurate data:

- Molecules (cryoprotectants)
- Molecular properties
- Mixtures
- Mixture components
- Calculation methods
- Predictions
- Experiments
- Experiment properties
- Property types

## Key Files

| File | Purpose | Status |
|------|---------|--------|
| `populate_database_supabase.py` | Main script for database population | ✅ Complete |
| `service_role_helper.py` | Helper for service role authentication | ⚠️ Needs Implementation |
| `.env` | Environment variables for Supabase | ⚠️ Needs Configuration |
| `.env.template` | Template for environment configuration | ⚠️ Needs Creation |

## Implementation Steps

### 1. Create Service Role Helper

Create or update `service_role_helper.py` to provide secure authentication for database operations:

```python
#!/usr/bin/env python3
"""
Service Role Helper for Supabase Authentication

This module provides helper functions for authenticating with Supabase using the service role.
"""

import os
from supabase import create_client, Client
from dotenv import load_dotenv
import uuid

# Load environment variables
load_dotenv()

# Supabase connection
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_SERVICE_KEY = os.getenv("SUPABASE_SERVICE_KEY")
DEFAULT_USER_ID = os.getenv("DEFAULT_USER_ID", "00000000-0000-0000-0000-000000000000")

def get_supabase_client() -> Client:
    """
    Get a Supabase client authenticated with the service role key.
    
    Returns:
        supabase.Client: Authenticated Supabase client
    
    Raises:
        ValueError: If SUPABASE_URL or SUPABASE_SERVICE_KEY is not set
    """
    if not SUPABASE_URL or not SUPABASE_SERVICE_KEY:
        raise ValueError("SUPABASE_URL and SUPABASE_SERVICE_KEY must be set in .env file")
    
    return create_client(SUPABASE_URL, SUPABASE_SERVICE_KEY)

def get_user_id() -> str:
    """
    Get the default user ID for service role operations.
    
    Returns:
        str: Default user ID
    """
    return DEFAULT_USER_ID

def ensure_user_profile(supabase: Client) -> str:
    """
    Ensure a user profile exists for the default user ID.
    
    Args:
        supabase (supabase.Client): Authenticated Supabase client
    
    Returns:
        str: User profile ID
    """
    # Check if profile already exists
    response = supabase.table("user_profile").select("*").eq("auth_user_id", DEFAULT_USER_ID).execute()
    if hasattr(response, 'data') and response.data:
        return response.data[0]["id"]
    
    # Create new profile if it doesn't exist
    profile_id = str(uuid.uuid4())
    profile_data = {
        "id": profile_id,
        "auth_user_id": DEFAULT_USER_ID,
        "display_name": "System",
        "email": "system@cryoprotect.example",
        "affiliation": "CryoProtect System",
        "created_at": datetime.now().isoformat(),
        "updated_at": datetime.now().isoformat()
    }
    
    supabase.table("user_profile").insert(profile_data).execute()
    return profile_id
```

### 2. Create Environment Template

Create `.env.template` for Supabase configuration:

```
# Supabase Configuration
SUPABASE_URL=https://xxxxxxxxxxxxxxxxxxxx.supabase.co
SUPABASE_KEY=your-anon-key
SUPABASE_SERVICE_KEY=your-service-role-key
SUPABASE_USER=admin@example.com
SUPABASE_PASSWORD=secure-password
DEFAULT_USER_ID=00000000-0000-0000-0000-000000000000
```

### 3. Integrate Database Population into Deployment Process

Add database population to your deployment script (`scripts/deploy_blue_green.sh`):

```bash
#!/bin/bash
# Deployment script with database population

# Database Schema and Population
echo "Checking database schema and populating data if needed..."
python populate_database_supabase.py

# Continue with deployment steps
# ...
```

## Data Included

The population script includes scientifically accurate data:

1. **Cryoprotectants**: DMSO, glycerol, ethylene glycol, trehalose, etc.
2. **Properties**: Glass transition temperature, cell membrane permeability, etc.
3. **Mixtures**: VS55, M22, and other standard vitrification solutions
4. **Experiments**: Literature-based experimental protocols and results
5. **Predictions**: AI-predicted properties using different methods

## Security Considerations

- The service role key should be used only for deployment and database migration
- Create separate roles for application access
- Implement proper RLS policies for data access
- Store all secrets in secure environment variables, never in code

## Usage in CI/CD Pipeline

In your GitHub Actions workflow, include database population step after schema migration:

```yaml
- name: Apply database schema
  run: python apply_schema.py
  env:
    SUPABASE_URL: ${{ secrets.SUPABASE_URL }}
    SUPABASE_SERVICE_KEY: ${{ secrets.SUPABASE_SERVICE_KEY }}

- name: Populate database
  run: python populate_database_supabase.py
  env:
    SUPABASE_URL: ${{ secrets.SUPABASE_URL }}
    SUPABASE_SERVICE_KEY: ${{ secrets.SUPABASE_SERVICE_KEY }}
    DEFAULT_USER_ID: ${{ secrets.DEFAULT_USER_ID }}
```

## Troubleshooting

- **Authentication Issues**: Verify service role key has proper permissions
- **Data Already Exists**: The script checks for existing data to avoid duplicates
- **Missing Tables**: Ensure schema migrations run successfully before population
- **Permission Denied**: Check RLS policies and roles

---

This guide provides all necessary information for implementing database population as part of the deployment infrastructure. Include database population in your CI/CD pipeline to ensure a properly configured environment for CryoProtect v2.