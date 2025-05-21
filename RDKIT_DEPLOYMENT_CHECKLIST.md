# RDKit Deployment Checklist

This checklist helps ensure a smooth transition to the new RDKit integration approach in the CryoProtect project.

## Pre-Deployment Tasks

- [ ] Run comprehensive tests with `rdkit_integration_test.py` in host environment
- [ ] Build RDKit container with `build_rdkit_container.sh`
- [ ] Test with RDKit container using `run_with_rdkit_container.sh`
- [ ] Update existing code with `update_rdkit_integration.sh`
- [ ] Review update summary in `rdkit_integration_update_summary.md`
- [ ] Manually review files marked as requiring manual attention
- [ ] Run unit tests for all updated files

## Deployment Procedure

1. **Staging Environment Deployment**
   - [ ] Deploy updated files to staging environment
   - [ ] Deploy RDKit container to staging environment
   - [ ] Run integration tests in staging environment
   - [ ] Monitor for errors in logs

2. **Production Deployment**
   - [ ] Schedule maintenance window if needed
   - [ ] Back up production environment
   - [ ] Deploy updated files to production
   - [ ] Deploy RDKit container to production
   - [ ] Run quick tests to verify functionality
   - [ ] Monitor for errors in logs

## Container Configuration

Ensure the following environment variables are properly set for the RDKit container:

```
CONTAINER_NAME=CryoProtect-RDKit
IMAGE_NAME=cryoprotect-rdkit
HOST_PORT=5001
CONTAINER_PORT=5000
APP_VOLUME=/path/to/app/dir
```

## RDKit Container Commands

```bash
# Build the container
./build_rdkit_container.sh

# Run a script in the container
./run_with_rdkit_container.sh your_script.py

# Restart the container
podman restart $CONTAINER_NAME

# View container logs
podman logs $CONTAINER_NAME

# Get a shell in the container
podman exec -it $CONTAINER_NAME bash
```

## Post-Deployment Verification

- [ ] Run integration tests in production environment
- [ ] Verify all critical functionality works as expected
- [ ] Monitor resource usage (memory, CPU)
- [ ] Check logs for warnings or errors
- [ ] Verify container auto-restarts if it crashes

## Rollback Procedure

If issues are encountered:

1. Stop and remove the RDKit container:
   ```bash
   podman stop $CONTAINER_NAME
   podman rm $CONTAINER_NAME
   ```

2. Restore backup files:
   ```bash
   find . -name "*.bak" | while read f; do cp "$f" "${f%.bak}"; done
   ```

3. Restart services with original configuration

## Long-term Maintenance

- Schedule regular updates to the RDKit container (quarterly)
- Monitor RDKit version releases for new features or bug fixes
- Enhance the wrapper as needed to support new functionality
- Test with new RDKit versions before upgrading

## Notes on Environment-Specific Considerations

### Fedora Environment

- SELinux context adjustments may be needed:
  ```bash
  sudo chcon -Rt container_file_t /path/to/app/dir
  ```

### Container Orchestration

- For Kubernetes deployments, use the Kubernetes manifests in `k8s/rdkit-deployment.yaml`
- For Docker Swarm, use the stack configuration in `docker/rdkit-stack.yml`

## Support and Resources

- RDKit documentation: https://www.rdkit.org/docs/index.html
- RDKit GitHub repository: https://github.com/rdkit/rdkit
- Project-specific documentation: `RDKIT_INTEGRATION_GUIDE.md`