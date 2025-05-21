# SELinux Configuration for CryoProtect Containers

This guide provides a step-by-step approach to configuring SELinux for CryoProtect containers, ensuring proper security while allowing necessary container operations.

## Understanding SELinux and Containers

SELinux provides Mandatory Access Control (MAC) for Linux systems. When working with containers, SELinux enforces isolation between the host system and containers, as well as between different containers.

Key concepts:
- **SELinux Context**: Security labels applied to files, processes, and network ports
- **Type Enforcement**: Restricting access based on the type of the resource and the domain of the process
- **MCS (Multi-Category Security)**: Additional separation for processes with the same type

## Common SELinux Issues with Containers

1. **Volume Mounting**: Inability to read/write files in mounted volumes
2. **Network Access**: Restricted network communication between containers
3. **Socket Access**: Issues with accessing Unix sockets
4. **Port Binding**: Problems binding to network ports

## Solution 1: SELinux Volume Labels

For proper volume mounts, use volume labeling options:

- `:Z` - Relabel with a private unshared label
- `:z` - Relabel with a shared label

Example:
```bash
podman run -v /path/on/host:/path/in/container:Z ...
```

## Solution 2: Create a Custom SELinux Policy Module

For CryoProtect's specific needs, create a custom policy:

1. First, create a policy file:

```bash
cat > cryoprotect.te << EOF
module cryoprotect 1.0;

require {
    type container_t;
    type container_file_t;
    type container_runtime_t;
    class file { read write getattr open };
    class dir { read search open };
}

#============= Rules for CryoProtect containers ===============
allow container_t container_file_t:file { read write getattr open };
allow container_t container_file_t:dir { read search open };
EOF
```

2. Compile the policy:

```bash
checkmodule -M -m -o cryoprotect.mod cryoprotect.te
semodule_package -o cryoprotect.pp -m cryoprotect.mod
```

3. Install the policy:

```bash
sudo semodule -i cryoprotect.pp
```

## Solution 3: Configure Container File Contexts

For CryoProtect's application directories:

```bash
# Create context for CryoProtect directories
sudo semanage fcontext -a -t container_file_t "/home/mushu/Projects/CryoProtect(/.*)?"

# Apply the context
sudo restorecon -Rv /home/mushu/Projects/CryoProtect
```

## Solution 4: Container Runtime Options

Use specific SELinux options with Podman:

```bash
# For development environments only
podman run --security-opt label=disable ...

# For production, with custom process label
podman run --security-opt label=type:container_runtime_t ...
```

## Recommended CryoProtect Container Setup

For the most reliable configuration:

### Development Environment

```bash
#!/bin/bash
# Development environment setup with SELinux considerations

# Create the network if it doesn't exist
podman network create cryoprotect-net 2>/dev/null || true

# Set up data directories with proper context
mkdir -p ./data/rdkit ./data/app
chcon -Rt container_file_t ./data

# RDKit service container
podman run -d --name=cryoprotect-rdkit --replace \
    --network=cryoprotect-net \
    -p 5002:5000 \
    -v "./data/rdkit:/data:Z" \
    python:3.10-slim \
    sh -c "cd /app && pip install flask && python mock_rdkit_service.py"

# App container
podman run -d --name=cryoprotect-app --replace \
    --network=cryoprotect-net \
    -p 5001:5000 \
    -e RDKIT_SERVICE_URL=http://cryoprotect-rdkit:5000 \
    -v "./data/app:/data:Z" \
    -v "$(pwd)/app.py:/app/app.py:z" \
    python:3.10-slim \
    sh -c "cd /app && pip install -r requirements_essential.txt && python app.py"
```

### Production Environment

For production, use more restrictive SELinux configurations:

```bash
#!/bin/bash
# Production environment setup with SELinux

# Create the network if it doesn't exist
podman network create cryoprotect-net 2>/dev/null || true

# Set up data directories with proper context
mkdir -p /var/cryoprotect/data/rdkit /var/cryoprotect/data/app
chcon -Rt container_file_t /var/cryoprotect/data

# RDKit service container with resource limits
podman run -d --name=cryoprotect-rdkit --replace \
    --network=cryoprotect-net \
    --security-opt label=type:container_runtime_t \
    --cpus=1 --memory=1g \
    -p 5002:5000 \
    -v "/var/cryoprotect/data/rdkit:/data:Z" \
    cryoprotect/rdkit:latest

# App container with resource limits
podman run -d --name=cryoprotect-app --replace \
    --network=cryoprotect-net \
    --security-opt label=type:container_runtime_t \
    --cpus=2 --memory=2g \
    -p 5001:5000 \
    -e RDKIT_SERVICE_URL=http://cryoprotect-rdkit:5000 \
    -v "/var/cryoprotect/data/app:/data:Z" \
    cryoprotect/app:latest
```

## Troubleshooting SELinux Container Issues

### 1. Check SELinux Status

```bash
getenforce
sestatus
```

### 2. View SELinux Denials

```bash
sudo ausearch -m avc -ts recent
```

### 3. Generate Policy from Denials

```bash
sudo ausearch -m avc -ts recent | audit2allow -M cryoprotect
```

### 4. Verify File Contexts

```bash
ls -Z /path/to/directory
```

### 5. Temporarily Set SELinux to Permissive

For testing only, not recommended for production:

```bash
sudo setenforce 0  # Set to permissive mode
# Run tests
sudo setenforce 1  # Return to enforcing mode
```

## Best Practices

1. **Use Volume Labeling**: Always use `:Z` or `:z` for volume mounts
2. **Minimal Privileges**: Follow the principle of least privilege
3. **Separate Data Directories**: Use different directories for different containers
4. **Audit Regularly**: Regularly check SELinux denials
5. **Document Context Changes**: Keep track of all SELinux context changes

By following this guide, CryoProtect containers should function properly while maintaining SELinux protection.