# Podman Conda Setup for CryoProtect

This guide explains how to set up a Podman container with conda for CryoProtect development and testing. The setup addresses SELinux permissions issues and provides a reliable environment for development.

## Overview

The CryoProtect project requires several specialized dependencies, including RDKit, which is best installed through conda. This setup provides two approaches:

1. **Comprehensive Setup**: Full-featured setup with SELinux configuration and systemd integration
2. **Quick Test Setup**: Simplified setup for quick testing without all the SELinux complexities

## Prerequisites

- Fedora Linux with Podman installed
- SELinux in enforcing mode
- Sufficient disk space (conda environments can be large)

## Comprehensive Setup

The comprehensive setup creates a fully configured Podman container with conda and all required dependencies, properly configured for SELinux. It includes several helper scripts and systemd integration.

### Setup Steps

1. **Create the Container**:
   ```bash
   ./create_cryoprotect_conda_container.sh
   ```

2. **Run Commands in the Container**:
   ```bash
   ./run_in_cryoprotect.sh "python -c 'import rdkit; print(rdkit.__version__)'"
   ```

3. **Start Interactive Shell**:
   ```bash
   ./run_in_cryoprotect.sh
   ```

4. **Run the Application**:
   ```bash
   ./run_app_in_container.sh
   ```

5. **Mount Project Files**:
   ```bash
   ./mount_project_in_container.sh
   ```

6. **Install as a Systemd Service**:
   ```bash
   systemctl --user daemon-reload
   systemctl --user enable cryoprotect-container.service
   systemctl --user start cryoprotect-container.service
   ```

### Helper Scripts

The comprehensive setup creates several helper scripts:

- **run_in_cryoprotect.sh**: Run commands in the container with conda environment
- **run_app_in_container.sh**: Run the CryoProtect application in the container
- **mount_project_in_container.sh**: Safely mount project files into the container
- **setup_mock_rdkit.sh**: Create a mock RDKit module for testing

## Quick Test Setup

The quick test setup is a simplified approach for testing without all the SELinux complexities. It's useful for quickly verifying that the container and dependencies are working correctly.

### Quick Test Steps

1. **Create the Test Container**:
   ```bash
   ./quick_conda_container.sh
   ```

2. **Run Interactive Shell**:
   ```bash
   podman exec -it CryoProtect-Test bash -c 'conda activate cryoprotect && bash'
   ```

3. **Run Tests**:
   ```bash
   podman exec -it CryoProtect-Test bash -c 'cd /app && conda activate cryoprotect && python test_import.py'
   ```

4. **Clean Up**:
   ```bash
   podman stop CryoProtect-Test
   podman rm CryoProtect-Test
   ```

## SELinux Considerations

SELinux can cause permission issues when working with container volume mounts. The comprehensive setup handles these issues by:

1. Using temporary directories with proper SELinux contexts
2. Using volume mounts with `:Z` suffix for proper labeling
3. Running containers with `--security-opt label=type:container_t`
4. Applying custom SELinux policies for container operations

For more detailed SELinux information, refer to:
- SELINUX_CONFIGURATION_GUIDE.md
- PODMAN_SELINUX_GUIDE.md

## Troubleshooting

### Common Issues

1. **"Permission denied" for container operations**
   - Check SELinux booleans: `getsebool -a | grep container`
   - Enable container management: `sudo setsebool -P container_manage_cgroup on`

2. **"Permission denied" accessing volume**
   - Check file context: `ls -ldZ /path/to/volume`
   - Set correct context: `sudo chcon -Rt container_file_t /path/to/volume`

3. **Conda environment not found**
   - Make sure to activate it: `conda activate cryoprotect`
   - Check if it was created: `conda env list`

4. **RDKit import errors**
   - Try using the mock RDKit: `./setup_mock_rdkit.sh`

### Checking SELinux Denials

To see if SELinux is blocking operations:
```bash
sudo ausearch -m AVC -ts recent
```

### Get Recommendations for Denials

```bash
sudo ausearch -m AVC -ts recent | audit2allow -a
```

## Best Practices

1. Use the comprehensive setup for development, the quick setup for testing
2. Keep the container running in the background to avoid setup time
3. Use the helper scripts rather than direct podman commands
4. Regularly check for SELinux denials
5. Use the mock RDKit for basic testing without full environment setup
6. For production, consider building a custom image with Dockerfile

## Additional Resources

- [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/)
- [Podman Documentation](https://podman.io/docs)
- [SELinux User's and Administrator's Guide](https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/7/html/selinux_users_and_administrators_guide/)
- [RDKit Documentation](https://www.rdkit.org/docs/)