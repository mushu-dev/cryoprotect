# Podman Deployment Guide for CryoProtect on Fedora

This guide provides instructions for deploying CryoProtect using Podman on Fedora Linux systems. Podman is a daemonless container engine that provides a Docker-compatible command line interface.

## Prerequisites

- Fedora Linux (35 or newer)
- Podman installed (`sudo dnf install -y podman`)
- Podman Compose (will be installed by the deployment script if not present)
- Python 3 and pip (for podman-compose if not installed)

## Quick Start

1. Clone the repository:
   ```bash
   git clone https://github.com/yourorg/cryoprotect.git
   cd cryoprotect
   ```

2. Run the setup script:
   ```bash
   ./deploy_with_podman.sh setup
   ```
   This will create necessary directories, set proper SELinux contexts, and check dependencies.

3. Edit the generated `.env` file with your Supabase credentials and other configuration options.

4. Start CryoProtect in development mode:
   ```bash
   ./deploy_with_podman.sh dev
   ```

5. For production deployment:
   ```bash
   ./deploy_with_podman.sh prod
   ```

## Configuration

The deployment is configured through the `.env` file, which is created by the setup script if it doesn't exist. Key configuration options include:

- `FLASK_ENV`: Set to `development` or `production`
- `HOST_PORT`: The port to expose the application on (default: 5000)
- `SUPABASE_URL`: Your Supabase project URL
- `SUPABASE_KEY`: Your Supabase anon key
- `SECRET_KEY`: Secret key for Flask session encryption
- Volume paths for persistent data

## Directory Structure

The deployment script creates and configures the following directories:

- `./logs`: Application logs with proper SELinux context
- `./backup/data`: Backup data storage
- `./redis-data`: Redis data persistence

## SELinux Integration

This deployment is SELinux-aware and automatically sets the appropriate context for volume mount points:

```bash
chcon -Rt container_file_t ./logs
chcon -Rt container_file_t ./backup/data
chcon -Rt container_file_t ./redis-data
```

If you encounter permission issues with volume mounts, verify that SELinux contexts are properly set.

## Automatic Startup

To configure CryoProtect to start automatically with your user session:

1. Create a systemd service file:
   ```bash
   ./deploy_with_podman.sh systemd
   ```

2. Enable the service:
   ```bash
   systemctl --user enable cryoprotect.service
   ```

3. Start the service:
   ```bash
   systemctl --user start cryoprotect.service
   ```

## Deployment Commands

The deployment script supports the following commands:

- `setup`: Perform initial setup (directories, SELinux contexts, dependencies)
- `dev`: Start in development mode (interactive, with logs to console)
- `prod`: Start in production mode (detached, as a background service)
- `stop`: Stop all CryoProtect services
- `systemd`: Create systemd user service for automatic startup
- `help`: Show help message

## Troubleshooting

### Volume Mount Permission Issues

If you encounter errors related to volume mounts:

1. Check SELinux status: `getenforce`
2. If SELinux is enforcing, make sure contexts are set: `ls -laZ ./logs ./backup/data ./redis-data`
3. Manually set contexts if needed: `sudo chcon -Rt container_file_t ./logs`

### Connection Refused Errors

If you see "connection refused" errors:

1. Check if Podman is running: `systemctl --user status podman.socket`
2. Start the socket if needed: `systemctl --user start podman.socket`
3. Verify the container is running: `podman ps`

### Networking Issues

If containers can't communicate:

1. Check network status: `podman network ls`
2. Inspect the app network: `podman network inspect app-network`
3. Restart the network if needed: `podman network rm app-network && podman network create app-network`

## Differences from Docker Deployment

Key differences between Podman and Docker deployments:

1. **Rootless by default**: Podman runs as a non-root user, enhancing security
2. **SELinux integration**: Volume mounts require proper SELinux contexts
3. **User systemd services**: Services run at user level rather than system level
4. **No daemon**: Podman doesn't use a background daemon like Docker

## Monitoring

Access application metrics:
- CryoProtect metrics: http://localhost:5000/metrics
- Container metrics: `podman stats`

## Updating

To update the application:

1. Pull latest code: `git pull`
2. Stop the service: `./deploy_with_podman.sh stop`
3. Rebuild containers: `podman-compose -f podman-compose.yml build`
4. Start the service: `./deploy_with_podman.sh prod`