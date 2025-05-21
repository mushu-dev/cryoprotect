# CryoProtect Container Setup Summary

This document summarizes the current status of the Podman container setup for CryoProtect and provides guidance on using the available options.

## Container Setup Options

We've created several approaches to containerize the CryoProtect application, each with different trade-offs:

1. **Comprehensive Conda Setup** (`create_cryoprotect_conda_container.sh`)
   - Full-featured setup with conda environment
   - Includes RDKit and all scientific libraries
   - Handles SELinux permissions properly
   - Provides systemd integration
   - **Status**: Created but installation takes longer (30+ minutes)

2. **Simplified Container** (`simplified_container_setup.sh`)
   - Lightweight Python container with basic dependencies
   - Uses mock RDKit for testing
   - Quick to set up (under 5 minutes)
   - **Status**: Working, accessible at http://localhost:5001

3. **Production Dockerfile** (`Dockerfile`)
   - Optimized multi-stage build
   - Security hardened with non-root user
   - Includes healthchecks
   - **Status**: Available but requires build time

## Using the Simplified Container

The simplified container is the quickest way to get started:

```bash
# Set up the container
./simplified_container_setup.sh

# Run commands in the container
./run_slim_container.sh "python -c 'import sys; print(sys.version)'"

# Start an interactive shell
./run_slim_container.sh

# Check container logs
podman logs CryoProtect-Slim
```

This container uses a mock RDKit implementation for testing purposes. It has a minimal Flask application running with these endpoints:

- `GET /`: Basic status endpoint
- `GET /health`: Health check endpoint

## Using the Full Conda Container

The full conda container provides a complete environment for development:

```bash
# Set up the container (takes time)
./create_cryoprotect_conda_container.sh

# Run commands in the container
./run_in_cryoprotect.sh "python -c 'import rdkit; print(rdkit.__version__)'"

# Start an interactive shell
./run_in_cryoprotect.sh

# Mount project files securely
./mount_project_in_container.sh

# Run the application
./run_app_in_container.sh
```

## SELinux Considerations

Both container setups handle SELinux permissions:

1. Using the `:Z` suffix for volumes to properly label mounted volumes
2. Setting `container_file_t` context on directories when possible
3. Using temporary directories with proper contexts for file mounting
4. Using `--security-opt label=type:container_t` to set the proper container type

## Next Steps

### For Using the Simplified Container
1. Test basic API functionality
2. Use mock RDKit for simpler tests
3. Implement minimal required features for testing

### For Using the Full Conda Container
1. Complete the conda environment installation (may take 30+ minutes)
2. Mount the full project for development
3. Implement all features with real RDKit

### For Production Deployment
1. Build the production Dockerfile
2. Set up proper volumes for data persistence
3. Configure systemd services for automatic startup

## Troubleshooting

### SELinux Issues
If you encounter "Permission denied" errors when mounting volumes:
```bash
sudo chcon -Rt container_file_t /path/to/volume
```

### Container Access Issues
If you can't access the container's API:
```bash
# Check if the container is running
podman ps

# Check the container logs
podman logs CryoProtect-Slim

# Check firewall settings
sudo firewall-cmd --list-all
```

### Port Conflicts
If you see "Address already in use" errors:
```bash
# Find what's using the port
sudo lsof -i :5000

# Modify the port in the setup script
# Change: podman run -d -p 5000:5000
# To: podman run -d -p 5001:5000
```