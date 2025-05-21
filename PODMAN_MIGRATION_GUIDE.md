# Docker to Podman Migration Guide for CryoProtect

This guide provides a comprehensive overview of migrating CryoProtect from Docker to Podman on Fedora Linux. Podman is a daemonless container engine that provides several advantages over Docker, especially on SELinux-enabled systems like Fedora.

## Table of Contents

1. [Why Migrate to Podman?](#why-migrate-to-podman)
2. [Migration Process Overview](#migration-process-overview)
3. [Prerequisites](#prerequisites)
4. [Automated Migration](#automated-migration)
5. [Manual Migration Steps](#manual-migration-steps)
6. [Troubleshooting](#troubleshooting)
7. [Fedora-Specific Considerations](#fedora-specific-considerations)
8. [Podman vs Docker Command Reference](#podman-vs-docker-command-reference)
9. [Best Practices](#best-practices)

## Why Migrate to Podman?

Podman offers several advantages over Docker, particularly on Fedora:

- **Daemonless Architecture**: No need for a persistent daemon process
- **Rootless Containers**: Run containers as a non-privileged user
- **Enhanced Security**: Better SELinux integration and OCI compliance
- **Drop-in Replacement**: Compatible API with Docker for smooth transition
- **Kubernetes Integration**: Native support for Kubernetes YAML format
- **Fedora Official Support**: Recommended container engine for Fedora

## Migration Process Overview

The migration from Docker to Podman involves these high-level steps:

1. Install Podman and related tools
2. Test system readiness for migration
3. Convert Docker configuration files to Podman format
4. Update deployment scripts
5. Test application on Podman
6. Document changes and update operational procedures

## Prerequisites

Before migrating to Podman, ensure your system meets these requirements:

- Fedora Linux (version 35 or higher recommended)
- Sufficient disk space for container images
- Basic understanding of container concepts
- Administrator (sudo) privileges for installation

## Automated Migration

We provide two scripts to assist with migration:

### 1. System Readiness Test

```bash
# Run the readiness test
./test_podman_readiness.sh
```

This script performs comprehensive checks to determine if your system is ready for Podman migration, including:
- Fedora version verification
- Podman and dependencies installation status
- SELinux configuration
- Resource availability
- Basic Podman functionality

### 2. Migration Script

```bash
# Run the migration script
./migrate_to_podman.sh
```

The migration script automates the following:
- Installing Podman if not already installed
- Converting Docker compose files to Podman format
- Creating Podman-specific configuration files
- Generating a quick-start script
- Creating a comprehensive deployment guide

## Manual Migration Steps

If you prefer to migrate manually, follow these steps:

### 1. Install Podman and Related Tools

```bash
# Install Podman, podman-compose, and SELinux integration
sudo dnf install -y podman podman-compose container-selinux

# Verify installation
podman --version
podman-compose --version
```

### 2. Convert Docker Compose Files

Podman requires modifications to Docker Compose files, especially for SELinux:

```bash
# Create a backup
cp docker-compose.yml docker-compose.yml.bak

# Modify volume mounts for SELinux
# Add :Z suffix to dedicated volumes (one container)
# Add :z suffix to shared volumes (multiple containers)
```

Example modification:
```yaml
# Original Docker entry
volumes:
  - ./data:/app/data

# Modified Podman entry
volumes:
  - ./data:/app/data:Z
```

### 3. Update Container Management Commands

Update any scripts or documentation to use Podman commands:

```bash
# Replace Docker commands with Podman equivalents
# docker run -> podman run
# docker-compose up -> podman-compose up
```

### 4. Address SELinux Contexts

```bash
# Set appropriate SELinux contexts for volume mount directories
chcon -Rt container_file_t ./logs
chcon -Rt container_file_t ./data
```

### 5. Test the Application

```bash
# Start the application with Podman
podman-compose up -d

# Check logs
podman logs -f cryoprotect
```

## Troubleshooting

### Common Issues and Solutions

#### Permission Denied on Volume Mounts

```bash
# Check SELinux contexts
ls -lZ ./data

# Reset context if needed
chcon -Rt container_file_t ./data
```

#### Network Connectivity Issues

```bash
# Check firewalld settings
sudo firewall-cmd --list-all

# Allow required ports
sudo firewall-cmd --permanent --add-port=5000/tcp
sudo firewall-cmd --reload
```

#### Rootless Podman Issues

```bash
# Check user namespace support
cat /proc/sys/kernel/unprivileged_userns_clone

# If result is 0, enable user namespaces
sudo sysctl -w kernel.unprivileged_userns_clone=1
```

## Fedora-Specific Considerations

### SELinux Integration

Podman on Fedora enforces SELinux policies by default, which enhances security but requires specific configuration:

- Volume mounts need `:Z` or `:z` suffix
- Host directories need appropriate SELinux context
- Container processes run with specific SELinux types

### Firewalld Configuration

Fedora's default firewall may require configuration for container networking:

```bash
# Allow pod network traffic
sudo firewall-cmd --permanent --zone=trusted --add-interface=cni-podman0

# Allow application ports
sudo firewall-cmd --permanent --add-port=5000/tcp
sudo firewall-cmd --reload
```

### Rootless Mode Considerations

Rootless Podman works well on Fedora but has some limitations:

- Port binding below 1024 requires additional configuration
- Some system-level operations may not be available
- Resource limits may be more restrictive

## Podman vs Docker Command Reference

| Task | Docker Command | Podman Command |
|------|---------------|---------------|
| Pull image | `docker pull image:tag` | `podman pull image:tag` |
| List images | `docker images` | `podman images` |
| List containers | `docker ps` | `podman ps` |
| Build image | `docker build -t name .` | `podman build -t name .` |
| Run container | `docker run image` | `podman run image` |
| Execute in container | `docker exec -it cont cmd` | `podman exec -it cont cmd` |
| View logs | `docker logs container` | `podman logs container` |
| Compose up | `docker-compose up -d` | `podman-compose up -d` |
| Compose down | `docker-compose down` | `podman-compose down` |
| System pruning | `docker system prune` | `podman system prune` |

## Best Practices

### Security Considerations

- Run containers as non-root users both inside and outside containers
- Use SELinux type enforcement for enhanced container isolation
- Implement container resource limits
- Regularly scan container images for vulnerabilities

### Performance Optimization

- Use volume mounts carefully, especially with large directories
- Consider using podman generate kube for Kubernetes compatibility
- Use --security-opt label=disable only when absolutely necessary
- Optimize image size using multi-stage builds

### Operations Management

- Create shell aliases for common commands
- Use systemd units for managing persistent containers
- Document podman-specific configurations
- Update CI/CD pipelines to work with Podman

---

This guide should help you successfully migrate CryoProtect from Docker to Podman on Fedora. For additional assistance, refer to the included scripts or consult the Podman documentation.