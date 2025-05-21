# CryoProtect Container Setup Guide

This guide explains how to set up and use the containerized environment for CryoProtect.

## Container Architecture

The CryoProtect application is deployed using two containers:

1. **cryoprotect-app**: The main Flask application that handles API requests and serves the web interface.
2. **cryoprotect-rdkit**: A service container that provides RDKit functionality for molecular property calculations.

These containers communicate via a dedicated network (`cryoprotect-net`).

## Setup

### Prerequisites

- Podman (version 3.0+) installed on your system
- Git repository cloned locally

### Quick Setup

Run the setup script to create and start the containers:

```bash
./setup_containers.sh
```

This script will:
1. Create the `cryoprotect-net` network if it doesn't exist
2. Stop and remove existing containers if they exist
3. Start the RDKit service container
4. Start the main application container
5. Verify that both containers are running

### Verification

To verify that the container setup is working correctly:

```bash
./check_container_setup.sh
```

This script checks:
- If both containers are running
- If they are on the same network
- If the RDKit service is healthy
- If the app can communicate with the RDKit service
- If property calculations work correctly

## Accessing the Services

- Main application: http://localhost:5001
- RDKit service: http://localhost:5002

## Cleanup

To stop and remove the containers:

```bash
./cleanup_containers.sh
```

Add the `--all` flag to also remove the network:

```bash
./cleanup_containers.sh --all
```

## Using the Property Search Client

A command-line client is provided to interact with the CryoProtect API:

```bash
# Check RDKit service status
./property_search_client.py --check

# Get properties for a specific molecule
./property_search_client.py --smiles CCO

# Search for molecules by property criteria
./property_search_client.py --min-mw 100 --max-mw 500 --max-logp 5
```

## Troubleshooting

### Containers cannot communicate

1. Verify that both containers are on the same network:
   ```bash
   podman network inspect cryoprotect-net
   ```

2. Verify that container names are resolved correctly:
   ```bash
   podman exec cryoprotect-app ping cryoprotect-rdkit
   ```

### Services not accessible

1. Check if containers are running:
   ```bash
   podman ps
   ```

2. Check the container logs:
   ```bash
   podman logs cryoprotect-app
   podman logs cryoprotect-rdkit
   ```

### RDKit service errors

If the RDKit service fails to start or returns errors, try using the mock RDKit service:

1. Make sure `mock_rdkit_service.py` is present in the project root
2. Restart the RDKit service container:
   ```bash
   podman stop cryoprotect-rdkit
   podman rm cryoprotect-rdkit
   podman run -d --name=cryoprotect-rdkit --network=cryoprotect-net -p 5002:5000 -v ./mock_rdkit_service.py:/app/mock_rdkit_service.py:z python:3.10-slim sh -c "cd /app && pip install flask && python mock_rdkit_service.py"
   ```