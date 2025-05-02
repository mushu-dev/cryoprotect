# Blue/Green Deployment for CryoProtect v2

This document provides comprehensive information on the blue/green deployment system implemented for CryoProtect v2, including setup, usage, and troubleshooting.

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Setup](#setup)
4. [Deployment Process](#deployment-process)
5. [Database Migrations](#database-migrations)
6. [Health Checking](#health-checking)
7. [Rollback Procedures](#rollback-procedures)
8. [Troubleshooting](#troubleshooting)
9. [Best Practices](#best-practices)

## Overview

Blue/green deployment is a technique that reduces downtime and risk by running two identical production environments called Blue and Green. At any time, only one of the environments is live, serving all production traffic. The other environment remains idle.

For CryoProtect v2, we've implemented a blue/green deployment system that allows for:

- Zero-downtime deployments
- Automated health checking
- Easy rollbacks
- Safe database migrations
- Traffic switching between environments

## Architecture

The blue/green deployment architecture consists of:

1. **NGINX Load Balancer**: Routes traffic to either the blue or green environment
2. **Blue Environment**: One of the two identical production environments
3. **Green Environment**: The other identical production environment
4. **Database**: Shared between both environments
5. **Deployment Scripts**: Automate the deployment process

```
                   ┌─────────────┐
                   │    NGINX    │
                   │Load Balancer│
                   └──────┬──────┘
                          │
                ┌─────────┴─────────┐
                │                   │
         ┌──────▼─────┐     ┌───────▼──────┐
         │    Blue    │     │    Green     │
         │Environment │     │ Environment  │
         └──────┬─────┘     └───────┬──────┘
                │                   │
                └─────────┬─────────┘
                          │
                    ┌─────▼─────┐
                    │ Database  │
                    └───────────┘
```

## Setup

### Prerequisites

- Docker and Docker Compose installed
- Git repository access
- Appropriate permissions to deploy to staging and production environments

### Initial Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/cryoprotect.git
   cd cryoprotect
   ```

2. Make the deployment scripts executable:
   ```bash
   ./scripts/make-scripts-executable.sh
   ```

3. Initialize the blue/green deployment environment with an initial version:
   ```bash
   ./scripts/init-blue-green.sh --version v1.0.0
   ```

This will:
- Set up the NGINX load balancer
- Deploy the initial version to the blue environment
- Configure traffic to route to the blue environment

## Deployment Process

### Deploying to Blue Environment

To deploy a new version to the blue environment:

```bash
./scripts/deploy-blue.sh --version v1.1.0
```

This will:
1. Pull the specified Docker image
2. Deploy it to the blue environment
3. Run health checks to ensure the deployment was successful
4. Apply any necessary database migrations

To automatically switch traffic to the blue environment after successful deployment:

```bash
./scripts/deploy-blue.sh --version v1.1.0 --auto-switch
```

### Deploying to Green Environment

To deploy a new version to the green environment:

```bash
./scripts/deploy-green.sh --version v1.1.0
```

To automatically switch traffic to the green environment after successful deployment:

```bash
./scripts/deploy-green.sh --version v1.1.0 --auto-switch
```

### Switching Traffic

To manually switch traffic between environments:

```bash
# Switch to blue
./scripts/switch-to-blue.sh

# Switch to green
./scripts/switch-to-green.sh
```

## Database Migrations

Database migrations are handled as part of the deployment process. The deployment scripts automatically run migrations after deploying a new version.

To manually run migrations:

```bash
./scripts/migrate-database.sh --environment blue
```

To run migrations in dry-run mode (without applying changes):

```bash
./scripts/migrate-database.sh --environment blue --dry-run
```

### Migration Safety

The migration process includes several safety measures:

1. **Database Backups**: A backup is created before running migrations
2. **Rollback Capability**: Failed migrations can be rolled back
3. **Dry Run Mode**: Migrations can be tested without applying changes

## Health Checking

Health checks are performed automatically during deployment to ensure the new version is functioning correctly.

To manually check the health of the environments:

```bash
./scripts/check-health.sh
```

To check a specific environment:

```bash
./scripts/check-health.sh --environment blue
```

For detailed health information:

```bash
./scripts/check-health.sh --verbose
```

### Health Check Criteria

The health check verifies:

1. Container status (running and healthy)
2. Application health endpoint response
3. NGINX configuration validity
4. Database connectivity

## Rollback Procedures

### Automated Rollback

If a deployment fails (e.g., health checks fail), the deployment script will automatically roll back to the previous version.

### Manual Rollback

To manually roll back to a previous version:

```bash
./scripts/rollback.sh --version v1.0.0 --environment blue
```

This will:
1. Stop the current environment
2. Deploy the specified version
3. Run health checks
4. Switch traffic to the rolled-back environment

### Emergency Rollback

In case of a critical failure:

1. Identify the last stable version
2. Run the rollback script with that version
3. Verify the rollback was successful using the health check script

## Troubleshooting

### Common Issues

#### Deployment Fails with Health Check Errors

1. Check the application logs:
   ```bash
   docker logs $(docker ps --filter "name=cryoprotect-blue" --quiet)
   ```

2. Verify the database connection:
   ```bash
   ./scripts/check-health.sh --environment blue --verbose
   ```

3. Check for migration errors:
   ```bash
   docker exec $(docker ps --filter "name=cryoprotect-blue" --quiet) /opt/conda/envs/cryoprotect/bin/python -m database.migrations status
   ```

#### NGINX Configuration Errors

1. Check NGINX configuration:
   ```bash
   docker exec $(docker ps --filter "name=nginx" --quiet) nginx -t
   ```

2. Check NGINX logs:
   ```bash
   docker logs $(docker ps --filter "name=nginx" --quiet)
   ```

#### Database Migration Failures

1. Check migration logs:
   ```bash
   docker exec $(docker ps --filter "name=cryoprotect-blue" --quiet) cat /app/logs/migration_log.log
   ```

2. Restore from backup if necessary:
   ```bash
   ./scripts/migrate-database.sh --environment blue
   # When prompted, choose to restore from backup
   ```

## Best Practices

1. **Always Test in Staging First**: Deploy to staging environment before production
2. **Use Semantic Versioning**: Follow semantic versioning for releases
3. **Monitor Deployments**: Watch logs and metrics during and after deployments
4. **Regular Backups**: Ensure database backups are created regularly
5. **Automate Everything**: Use CI/CD pipelines to automate the deployment process
6. **Gradual Rollouts**: Consider implementing canary deployments for critical updates
7. **Document Changes**: Keep a detailed changelog for each version

## Conclusion

The blue/green deployment system provides a robust, zero-downtime deployment process for CryoProtect v2. By following the procedures outlined in this document, you can safely deploy new versions, handle database migrations, and quickly roll back if issues arise.