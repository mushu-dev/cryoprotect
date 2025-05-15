# RDKit Service Integration Guide

This guide explains how to set up and integrate the RDKit microservice with the main CryoProtect application.

## Architecture Overview

The solution uses a microservice architecture:

1. **RDKit Microservice**: Deployed on Fly.io, provides RDKit molecular calculations via REST API
2. **Main Application**: Deployed on Heroku, communicates with the RDKit microservice

This architecture allows us to:
- Use conda and RDKit on Fly.io, bypassing Heroku's limitations
- Keep the main application on Heroku, utilizing student credits
- Scale each component independently based on demand

## Deployment Instructions

### 1. Deploy RDKit Microservice to Fly.io

```bash
# Install Fly CLI
curl -L https://fly.io/install.sh | sh

# Navigate to the RDKit microservice directory
cd rdkit-service

# Authenticate with Fly.io
fly auth login

# Create a new Fly.io app
fly launch --name cryoprotect-rdkit

# Deploy the app
fly deploy
```

### 2. Set Up Environment Variables on Heroku

```bash
# Set the RDKit service URL in Heroku
heroku config:set RDKIT_SERVICE_URL=https://cryoprotect-rdkit.fly.dev

# Optional: Set the API key for authentication
heroku config:set RDKIT_API_KEY=your-secure-api-key
```

### 3. Integrate the Client in Heroku Application

1. Copy the `heroku-client.py` file to your Heroku application
2. Import and use the client in your application:

```python
from heroku_client import RDKitServiceClient

# Initialize the client
rdkit_client = RDKitServiceClient()

# Use the client in your application
@app.route('/api/molecules/<int:molecule_id>/properties')
def get_molecule_properties(molecule_id):
    # Get the molecule from the database
    molecule = get_molecule_by_id(molecule_id)
    
    # Calculate properties using the RDKit service
    properties = rdkit_client.calculate_properties(molecule.smiles)
    
    return jsonify(properties)
```

## API Reference

The RDKit microservice provides the following endpoints:

### Calculate Molecular Properties

**Endpoint:** `/api/calculate-properties`  
**Method:** POST  
**Description:** Calculate molecular properties using RDKit  

**Request Body:**
```json
{
  "molecule_data": "CCO",
  "input_format": "smiles"
}
```

### Generate Molecular Visualization

**Endpoint:** `/api/visualization`  
**Method:** POST  
**Description:** Generate a visualization of a molecule as SVG  

**Request Body:**
```json
{
  "molecule_data": "CCO",
  "input_format": "smiles",
  "width": 400,
  "height": 300,
  "highlight_atoms": [0, 1]
}
```

### Perform Substructure Search

**Endpoint:** `/api/substructure-search`  
**Method:** POST  
**Description:** Search for a substructure within a molecule  

**Request Body:**
```json
{
  "query_mol_data": "[OH]",
  "target_mol_data": "CCO",
  "query_format": "smarts",
  "target_format": "smiles"
}
```

### Calculate Molecular Similarity

**Endpoint:** `/api/similarity`  
**Method:** POST  
**Description:** Calculate similarity between two molecules  

**Request Body:**
```json
{
  "mol1_data": "CCO",
  "mol2_data": "CC(=O)O",
  "mol1_format": "smiles",
  "mol2_format": "smiles",
  "fingerprint_type": "morgan"
}
```

## Client Configuration Options

The RDKit client can be configured with the following options:

- `service_url`: URL of the RDKit service. If not provided, uses the `RDKIT_SERVICE_URL` environment variable.
- `api_key`: API key for authentication. If not provided, uses the `RDKIT_API_KEY` environment variable.
- `timeout`: Request timeout in seconds. Defaults to 30 seconds.
- `max_retries`: Maximum number of retries for failed requests. Defaults to 3.

## Fallback Behavior

If the RDKit service is unavailable or not configured, the client falls back to a mock implementation. This ensures that your application can continue to function even if the RDKit service is down.

## Monitoring

- Monitor the RDKit service using Fly.io's built-in metrics dashboard:
  ```bash
  fly dashboard
  ```

- Monitor the client in your Heroku application by checking the logs:
  ```bash
  heroku logs --tail
  ```

## Scaling

The RDKit service can be scaled on Fly.io as needed:

```bash
# Scale the number of instances
fly scale count 3

# Scale the memory allocation
fly scale memory 1024
```

## Security Considerations

- The service includes optional API key authentication
- All traffic is encrypted via HTTPS
- The service can be restricted to only accept requests from your Heroku application's IP address

## Troubleshooting

If you encounter issues with the RDKit service:

1. Check the Fly.io logs:
   ```bash
   fly logs
   ```

2. Verify the service is running:
   ```bash
   fly status
   ```

3. Check the health endpoint:
   ```bash
   curl https://cryoprotect-rdkit.fly.dev/health
   ```

4. Test the client with debug logging:
   ```python
   import logging
   logging.getLogger('rdkit-client').setLevel(logging.DEBUG)
   
   client = RDKitServiceClient()
   result = client.calculate_properties("CCO")
   print(result)
   ```