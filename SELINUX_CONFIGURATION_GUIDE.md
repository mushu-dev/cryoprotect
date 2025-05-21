# SELinux Configuration Guide for CryoProtect

This comprehensive guide covers SELinux configuration for the CryoProtect application, providing security hardening for Fedora environments. SELinux adds an additional layer of security beyond standard Linux permissions, preventing unauthorized access and minimizing potential damage from security breaches.

## Table of Contents

1. [Introduction to SELinux](#introduction-to-selinux)
2. [SELinux Configuration Tools](#selinux-configuration-tools)
3. [Quick Start](#quick-start)
4. [Container Security](#container-security)
5. [Directory and File Contexts](#directory-and-file-contexts)
6. [Network Port Configuration](#network-port-configuration)
7. [Custom SELinux Policies](#custom-selinux-policies)
8. [Troubleshooting SELinux Issues](#troubleshooting-selinux-issues)
9. [Best Practices](#best-practices)
10. [Additional Resources](#additional-resources)

## Introduction to SELinux

Security-Enhanced Linux (SELinux) is a mandatory access control (MAC) security mechanism implemented in the Linux kernel. It provides fine-grained control over which processes and users can access specific files, directories, ports, and other resources.

### SELinux Modes

- **Enforcing**: SELinux policy is enforced. Access denied events are logged and prohibited.
- **Permissive**: SELinux policy is not enforced, but access denied events are logged.
- **Disabled**: SELinux is fully disabled.

For production environments, SELinux should be set to **Enforcing** mode.

### Key SELinux Concepts

- **Contexts**: Labels assigned to processes, files, and resources
- **Types**: The part of the context that defines what the resource is (e.g., `container_file_t`)
- **Domains**: Process types that define what a process can do
- **Booleans**: On/off switches that modify policy behavior
- **Policies**: Rules defining allowed interactions between types

## SELinux Configuration Tools

The CryoProtect project includes several tools for configuring and verifying SELinux:

1. **`setup_selinux.sh`**: Main configuration script for SELinux
   - Sets container booleans
   - Configures directory contexts
   - Sets port contexts
   - Creates custom policy module

2. **`setup_podman_selinux.sh`**: Podman-specific SELinux configuration
   - Sets Podman-specific booleans
   - Configures volume labeling
   - Creates helper scripts for container management

3. **`verify_selinux_config.sh`**: Verification script to test SELinux configuration
   - Checks SELinux status
   - Verifies booleans are set correctly
   - Checks directory contexts
   - Checks port contexts
   - Generates verification report

4. **`monitor_selinux_denials.sh`**: Tool to monitor and report SELinux denials
   - Displays recent denials affecting CryoProtect
   - Generates policy recommendations

## Quick Start

To quickly set up SELinux for CryoProtect:

```bash
# 1. Configure SELinux (requires root)
sudo ./setup_selinux.sh

# 2. Configure Podman SELinux settings (requires root)
sudo ./setup_podman_selinux.sh --apply

# 3. Verify the configuration
./verify_selinux_config.sh --generate-report

# 4. Monitor for SELinux denials (requires root)
sudo ./monitor_selinux_denials.sh
```

## Container Security

### Key Container SELinux Booleans

| Boolean | Description | Recommended |
|---------|-------------|-------------|
| `container_manage_cgroup` | Allow containers to manage cgroup files | ON |
| `container_use_devices` | Allow containers to use devices | ON |
| `container_connect_any` | Allow containers to connect to any port | OFF |
| `container_user_exec_content` | Allow container to execute content | ON |

### Container File Contexts

Container files should use the `container_file_t` context:

```bash
# Set context for a directory
sudo semanage fcontext -a -t container_file_t "/path/to/dir(/.*)?"
sudo restorecon -Rv /path/to/dir

# View context
ls -ldZ /path/to/dir
```

### Podman Volume Mounts

When mounting volumes in Podman, use these suffixes:

- `:Z` - For dedicated volume (used by a single container)
- `:z` - For shared volume (used by multiple containers)

Example:
```yaml
volumes:
  - ./database:/var/lib/postgresql/data:Z
```

## Directory and File Contexts

### Critical Directories for CryoProtect

| Directory | Recommended Context | Purpose |
|-----------|---------------------|---------|
| `database/` | `container_file_t` | Database files |
| `backups/` | `container_file_t` | Backup files |
| `postgres_data/` | `container_file_t` | PostgreSQL data |
| `/etc/systemd/system/cryoprotect*.service` | `systemd_unit_file_t` | Systemd service files |

### Setting Directory Contexts

```bash
# Set context for database directory
sudo semanage fcontext -a -t container_file_t "/home/mushu/Projects/CryoProtect/database(/.*)?"
sudo restorecon -Rv /home/mushu/Projects/CryoProtect/database

# Set context for backup directory
sudo semanage fcontext -a -t container_file_t "/home/mushu/Projects/CryoProtect/backups(/.*)?"
sudo restorecon -Rv /home/mushu/Projects/CryoProtect/backups
```

### Temporary Context Changes

For temporary changes (not persisting across restorecon):

```bash
# Temporarily change context
sudo chcon -Rt container_file_t /path/to/directory

# Change context of a single file
sudo chcon -t container_file_t /path/to/file
```

## Network Port Configuration

### Important Ports for CryoProtect

| Port | Service | Recommended Context |
|------|---------|---------------------|
| 5432 | PostgreSQL | `postgresql_port_t` |
| 5433 | PostgreSQL Replica | `postgresql_port_t` |
| 8000 | Flask Application | `http_port_t` |

### Managing Port Contexts

```bash
# Add PostgreSQL port context
sudo semanage port -a -t postgresql_port_t -p tcp 5432

# Modify existing port context
sudo semanage port -m -t postgresql_port_t -p tcp 5432

# List port contexts
sudo semanage port -l | grep -E "5432|http_port"
```

## Custom SELinux Policies

Custom policies allow precise control over SELinux behavior for the application.

### CryoProtect Custom Policy

The CryoProtect custom policy (`cryoprotect_custom.te`) allows:

- Containers to connect to PostgreSQL
- Containers to read/write container_file_t directories
- Containers to read user home directories for configuration

### Creating and Applying Custom Policies

```bash
# Generate policy from denials
sudo ausearch -m AVC -ts recent | audit2allow -M my_policy

# Apply the policy
sudo semodule -i my_policy.pp
```

### Managing SELinux Modules

```bash
# List all installed modules
sudo semodule -l

# Remove a module
sudo semodule -r module_name
```

## Troubleshooting SELinux Issues

### Common SELinux Issues

1. **"Permission denied" for container operations**
   - Check SELinux booleans: `getsebool -a | grep container`
   - Enable container management: `sudo setsebool -P container_manage_cgroup on`

2. **Container can't connect to database**
   - Check port context: `sudo semanage port -l | grep 5432`
   - Set port context: `sudo semanage port -a -t postgresql_port_t -p tcp 5432`

3. **"Permission denied" accessing volume**
   - Check file context: `ls -ldZ /path/to/volume`
   - Set correct context: `sudo chcon -Rt container_file_t /path/to/volume`

### Finding SELinux Denials

```bash
# View recent denials
sudo ausearch -m AVC -ts recent

# Filter for specific application
sudo ausearch -m AVC -ts recent | grep -i cryoprotect

# Get policy recommendations for denials
sudo ausearch -m AVC -ts recent | audit2allow -a
```

### Using the Monitor Script

The included monitoring script provides targeted information about CryoProtect-related SELinux denials:

```bash
sudo ./monitor_selinux_denials.sh
```

## Best Practices

1. **Always run in enforcing mode** in production environments
2. **Use volume mounts with `:Z` or `:z`** suffixes with Podman
3. **Run regular verification** with `verify_selinux_config.sh`
4. **Monitor for denials** with `monitor_selinux_denials.sh`
5. **Automate SELinux configuration** using the provided scripts
6. **Document policies** in a central location
7. **Never disable SELinux** as a solution to problems
8. **Test in permissive mode** before applying changes in production
9. **Label new directories** when created
10. **Update policies** when application behavior changes

## Additional Resources

- [Fedora SELinux Guide](https://docs.fedoraproject.org/en-US/quick-docs/selinux/)
- [Red Hat SELinux Users and Administrators Guide](https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/7/html/selinux_users_and_administrators_guide/)
- [Container SELinux Context Types](https://www.redhat.com/sysadmin/container-selinux-types)
- [Udica: Custom SELinux policies for containers](https://github.com/containers/udica)

---

This guide was created for the CryoProtect project to ensure proper SELinux configuration in Fedora environments. See also the Podman-specific guide (`PODMAN_SELINUX_GUIDE.md`) for additional container-specific information.