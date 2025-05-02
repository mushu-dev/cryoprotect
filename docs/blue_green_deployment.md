# Blue/Green Deployment Guide for CryoProtect v2

This document provides comprehensive instructions for implementing and managing blue/green deployments for the CryoProtect v2 application.

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Prerequisites](#prerequisites)
4. [Deployment Process](#deployment-process)
5. [Cutover Process](#cutover-process)
6. [Rollback Process](#rollback-process)
7. [Health Checks](#health-checks)
8. [Monitoring](#monitoring)
9. [Troubleshooting](#troubleshooting)

## Overview

Blue/green deployment is a release strategy that reduces downtime and risk by running two identical production environments called "blue" and "green". At any time, only one of the environments is live and serving production traffic. This approach allows for seamless updates and quick rollbacks if issues are detected.

## Architecture

The CryoProtect v2 blue/green deployment architecture consists of:

- **NGINX Load Balancer**: Routes traffic to either blue or green environment
- **Blue Environment**: One version of the application
- **Green Environment**: Another version of the application
- **Shared Resources**: Database, Redis, and other stateful services
- **Monitoring**: Prometheus and Grafana for metrics collection and visualization

### Key Files

- `docker-compose.blue-green.yml`: Defines the blue/green deployment services
- `nginx/nginx.conf`: Main NGINX configuration
- `nginx/conf.d/active.conf`: Points to the currently active environment
- `nginx/conf.d/blue.conf`: Configuration for blue environment
- `nginx/conf.d/green.conf`: Configuration for green environment
- `nginx/scripts/cutover.sh`: Script for managing cutover between environments

## Prerequisites

- Docker and Docker Compose installed
- Access to the container registry (GitHub Container Registry)
- Proper environment variables set (see below)

### Required Environment Variables

```
# Version tags
BLUE_VERSION=<version-tag>
GREEN_VERSION=<version-tag>

# Resource allocation
BLUE_MEM_LIMIT=512M
BLUE_CPU_LIMIT=0.5
GREEN_MEM_LIMIT=512M
GREEN_CPU_LIMIT=0.5

# Secret management
USE_EXTERNAL_SECRETS=true
SUPABASE_URL_SECRET=cryoprotect_supabase_url
SUPABASE_KEY_SECRET=cryoprotect_supabase_key
SECRET_KEY_SECRET=cryoprotect_secret_key
REDIS_URL_SECRET=cryoprotect_redis_url
```

## Deployment Process

### Initial Deployment

1. Set up the initial environment (typically blue):

```bash
# Set environment variables
export BLUE_VERSION=v1.0.0
export GREEN_VERSION=v1.0.0

# Start the services
docker-compose -f docker-compose.blue-green.yml up -d
```

2. Verify the deployment:

```bash
# Check container status
docker-compose -f docker-compose.blue-green.yml ps

# Check health endpoints
curl http://localhost:8080/health/blue
curl http://localhost:8080/health/green
```

### Deploying a New Version

1. Update the version tag for the inactive environment:

```bash
# If blue is active, update green
export GREEN_VERSION=v1.1.0

# If green is active, update blue
export BLUE_VERSION=v1.1.0
```

2. Update the inactive environment:

```bash
docker-compose -f docker-compose.blue-green.yml up -d
```

3. Wait for the new environment to be healthy:

```bash
# If deploying to green
curl http://localhost:8080/health/green

# If deploying to blue
curl http://localhost:8080/health/blue
```

## Cutover Process

The cutover process switches traffic from the currently active environment to the newly deployed environment.

### Using the Cutover Script

```bash
# Switch to green environment
docker exec -it nginx /usr/local/bin/cutover.sh --target green

# Switch to blue environment
docker exec -it nginx /usr/local/bin/cutover.sh --target blue
```

### Manual Cutover

If you need to perform a manual cutover:

1. Copy the target configuration to active.conf:

```bash
# For green
docker exec -it nginx cp /etc/nginx/conf.d/green.conf /etc/nginx/conf.d/active.conf

# For blue
docker exec -it nginx cp /etc/nginx/conf.d/blue.conf /etc/nginx/conf.d/active.conf
```

2. Reload NGINX configuration:

```bash
docker exec -it nginx nginx -s reload
```

3. Verify the cutover:

```bash
curl http://localhost
```

## Rollback Process

If issues are detected after a cutover, you can quickly roll back to the previous environment.

### Using the Cutover Script

```bash
# Rollback to the previous environment
docker exec -it nginx /usr/local/bin/cutover.sh --rollback
```

### Manual Rollback

1. Determine the previous environment (if current is green, previous is blue, and vice versa)
2. Copy the previous configuration to active.conf:

```bash
# If current is green, rollback to blue
docker exec -it nginx cp /etc/nginx/conf.d/blue.conf /etc/nginx/conf.d/active.conf

# If current is blue, rollback to green
docker exec -it nginx cp /etc/nginx/conf.d/green.conf /etc/nginx/conf.d/active.conf
```

3. Reload NGINX configuration:

```bash
docker exec -it nginx nginx -s reload
```

## Health Checks

Health checks are critical for ensuring the reliability of blue/green deployments.

### Application Health Checks

The CryoProtect application exposes a `/health` endpoint that returns:
- HTTP 200: Application is healthy
- HTTP 503: Application is unhealthy

### NGINX Health Checks

NGINX provides health check endpoints:
- `/health`: NGINX status
- `/health/blue`: Blue environment status
- `/health/green`: Green environment status

### Docker Health Checks

Docker Compose health checks are configured for all services to ensure they are running properly.

## Monitoring

Monitoring is essential for detecting issues during and after deployment.

### Prometheus Metrics

Access Prometheus metrics at:
```
http://localhost:9090
```

### Grafana Dashboards

Access Grafana dashboards at:
```
http://localhost:3000
```

Default credentials:
- Username: admin
- Password: admin (change on first login)

## Troubleshooting

### Common Issues

1. **Health Check Failures**

   If health checks are failing:
   - Check application logs: `docker-compose -f docker-compose.blue-green.yml logs cryoprotect-blue`
   - Verify the application is running: `docker exec -it cryoprotect-blue ps aux`
   - Check the health endpoint directly: `docker exec -it cryoprotect-blue curl localhost:5000/health`

2. **NGINX Configuration Issues**

   If NGINX is not routing traffic correctly:
   - Check NGINX logs: `docker-compose -f docker-compose.blue-green.yml logs nginx`
   - Verify active.conf: `docker exec -it nginx cat /etc/nginx/conf.d/active.conf`
   - Test NGINX configuration: `docker exec -it nginx nginx -t`

3. **Cutover Script Failures**

   If the cutover script fails:
   - Check script permissions: `docker exec -it nginx ls -la /usr/local/bin/cutover.sh`
   - Run with verbose output: `docker exec -it nginx /usr/local/bin/cutover.sh --target blue --force`
   - Perform manual cutover as described above

### Recovery Steps

If both environments are unhealthy:

1. Check application logs for both environments
2. Verify database connectivity
3. Ensure secrets are properly configured
4. Restart the unhealthy services:
   ```bash
   docker-compose -f docker-compose.blue-green.yml restart cryoprotect-blue cryoprotect-green
   ```
5. If necessary, roll back to a known good version:
   ```bash
   export BLUE_VERSION=v1.0.0
   export GREEN_VERSION=v1.0.0
   docker-compose -f docker-compose.blue-green.yml up -d
   ```

## Conclusion

The blue/green deployment strategy provides a robust mechanism for deploying new versions of the CryoProtect application with minimal downtime and risk. By following the procedures outlined in this document, you can ensure smooth deployments and quick recovery from any issues that may arise.