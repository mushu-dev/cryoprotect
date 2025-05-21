# Fedora Podman Migration Guide for CryoProtect

This document provides a comprehensive guide for migrating the CryoProtect application from Docker to Podman on Fedora Linux systems. The migration process addresses Fedora-specific considerations, SELinux integration, IPv4/IPv6 connectivity issues, and ensuring a smooth transition with minimal downtime.

## Table of Contents

1. [Migration Overview](#migration-overview)
2. [Included Resources](#included-resources)
3. [Migration Steps](#migration-steps)
4. [Fedora-Specific Considerations](#fedora-specific-considerations)
5. [Troubleshooting Common Issues](#troubleshooting-common-issues)
6. [Post-Migration Verification](#post-migration-verification)
7. [Reverting if Necessary](#reverting-if-necessary)

## Migration Overview

This migration package provides a comprehensive set of tools to transition CryoProtect from Docker to Podman on Fedora Linux. The approach focuses on:

- **Minimal disruption**: Testing system readiness before proceeding
- **Gradual transition**: Providing both Docker and Podman configurations during transition
- **SELinux integration**: Properly configuring SELinux contexts for security
- **IPv4 compatibility**: Ensuring Supabase connections work properly
- **Comprehensive documentation**: Clear guidance for all migration steps

## Included Resources

The migration package includes the following resources:

1. **test_podman_readiness.sh**: Pre-migration system assessment
2. **migrate_to_podman.sh**: Main migration script
3. **podman-compose.yml**: Podman-compatible compose configuration
4. **podman-compose.minimal.yml**: Simplified configuration for easy testing
5. **quickstart_podman.sh**: Simplified startup script
6. **podman_post_installation_check.sh**: Post-migration verification
7. **PODMAN_DEPLOYMENT_GUIDE.md**: Comprehensive deployment documentation
8. **PODMAN_MIGRATION_GUIDE.md**: General migration documentation

## Migration Steps

### Step 1: Test System Readiness

Begin by assessing if your system is ready for Podman migration:

```bash
./test_podman_readiness.sh
```

This script checks:
- Fedora version compatibility
- Podman and podman-compose installation
- SELinux configuration
- System resources
- Container prerequisites

Address any issues identified by the readiness test before proceeding.

### Step 2: Run Migration Script

Once your system passes the readiness test, run the migration script:

```bash
./migrate_to_podman.sh
```

This script:
- Installs Podman if needed
- Converts Docker Compose files to Podman format
- Adds SELinux context labels to volume mounts
- Creates minimal configurations for testing
- Generates comprehensive documentation

### Step 3: Test with Minimal Configuration

Test the application with the minimal Podman configuration:

```bash
./quickstart_podman.sh
```

This provides a simplified startup experience to verify basic functionality.

### Step 4: Perform Post-Installation Verification

After testing basic functionality, run the post-installation verification:

```bash
./podman_post_installation_check.sh
```

This script performs comprehensive checks including:
- Container health verification
- SELinux context validation
- Network configuration testing
- IPv4 Supabase connectivity
- Resource usage assessment
- Application logs analysis

### Step 5: Transition to Full Deployment

Once verification passes, transition to the full deployment:

```bash
podman-compose -f podman-compose.yml up -d
```

## Fedora-Specific Considerations

### SELinux Integration

Fedora uses SELinux by default, which requires special handling for containers:

1. **Volume Labels**: All volume mounts need `:Z` (dedicated) or `:z` (shared) suffixes
2. **Directory Contexts**: Host directories need appropriate SELinux contexts:
   ```bash
   chcon -Rt container_file_t ./logs
   chcon -Rt container_file_t ./backup/data
   ```
3. **Container Security Options**: The `security_opt` setting in compose files:
   ```yaml
   security_opt:
     - label=type:container_file_t
   ```

### IPv4 Compatibility for Supabase

To ensure Supabase connections work properly on Fedora (which may use IPv6 by default):

1. **Force IPv4 Connection**: Use the connection test script in `scripts/test_connectivity.sh`
2. **Network Configuration**: If IPv6 causes issues, force IPv4 when needed:
   ```bash
   podman run --network=ip4only cryoprotect
   ```
3. **DNS Resolution**: Configure proper DNS resolution in `/etc/hosts`

### Firewall Configuration

Fedora's firewalld may require adjustments:

```bash
sudo firewall-cmd --permanent --add-port=5000/tcp
sudo firewall-cmd --reload
```

## Troubleshooting Common Issues

### Permission Denied on Volume Mounts

If you encounter permission issues with volume mounts:

```bash
# Check SELinux contexts
ls -lZ ./logs ./backup/data

# Reset context if needed
chcon -Rt container_file_t ./logs
chcon -Rt container_file_t ./backup/data
```

### IPv6/IPv4 Connectivity Issues

If you experience connectivity issues with Supabase:

```bash
# Test IPv4 connectivity
ping -4 your-project.supabase.co

# Add IPv4 address to /etc/hosts
echo "123.45.67.89 your-project.supabase.co" | sudo tee -a /etc/hosts
```

### Container Failed to Start

If containers fail to start:

```bash
# Check container logs
podman logs cryoprotect

# Check SELinux denials
sudo ausearch -m avc -ts recent

# Run with SELinux debugging
podman run --log-level=debug cryoprotect
```

## Post-Migration Verification

After completing migration, verify functionality:

1. **Application Functionality**: Test all application features
2. **Performance Testing**: Compare performance with Docker deployment
3. **Security Assessment**: Verify SELinux policies are working as expected
4. **Resource Utilization**: Monitor resource usage and compare to Docker

## Reverting if Necessary

If issues are encountered that cannot be resolved, you can revert to Docker:

1. Stop Podman containers:
   ```bash
   podman-compose -f podman-compose.yml down
   ```

2. Restart Docker containers:
   ```bash
   docker-compose up -d
   ```

3. Reset SELinux contexts if needed:
   ```bash
   restorecon -Rv ./logs ./backup/data
   ```

## Further Resources

- [Podman Documentation](https://podman.io/docs/)
- [Fedora Container Documentation](https://docs.fedoraproject.org/en-US/fedora-silverblue/docker-podman/)
- [SELinux User's Guide](https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/8/html/using_selinux/index)
- [Podman-Compose GitHub](https://github.com/containers/podman-compose)