# Deployment Guide

This document provides comprehensive instructions for deploying CryoProtect v2 to various environments.

## Prerequisites

- Docker and Docker Compose installed
- Access to GitHub Container Registry (ghcr.io)
- PostgreSQL database (local deployment) or Supabase credentials (cloud deployment)
- SSH access to deployment server (for manual deployment)

## Environment Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/organization/cryoprotect.git
   cd cryoprotect
   ```

2. Create environment configuration:
   ```bash
   # For development
   cp .env.template .env
   
   # For staging
   cp .env.staging .env
   
   # For production
   cp .env.production .env
   ```

3. Edit the `.env` file with appropriate values for your environment.

## Local Deployment

1. Initialize the local database:
   ```bash
   python database/init_local_db.py
   ```

2. Run the application using Docker Compose:
   ```bash
   docker-compose up -d
   ```

3. Verify the deployment:
   ```bash
   curl http://localhost:5000/api/health
   ```

## Cloud Deployment

### Automated Deployment via CI/CD

1. Push changes to the `master` branch to trigger the deployment pipeline.
2. The CI/CD pipeline will:
   - Run tests
   - Build and push a Docker image
   - Deploy to the staging environment
   - (Manual approval required for production deployment)

### Manual Deployment

1. Build the Docker image:
   ```bash
   docker build -t cryoprotect:latest .
   ```

2. Push the image to a registry:
   ```bash
   docker tag cryoprotect:latest ghcr.io/organization/cryoprotect:latest
   docker push ghcr.io/organization/cryoprotect:latest
   ```

3. SSH into the deployment server:
   ```bash
   ssh user@server
   ```

4. Deploy using Docker Compose:
   ```bash
   cd /path/to/deployment
   docker-compose pull
   docker-compose down
   docker-compose up -d
   ```

## Blue-Green Deployment

For zero-downtime production deployments, we use a blue-green deployment strategy:

1. Set up the environment:
   ```bash
   ./scripts/init-blue-green.sh
   ```

2. Deploy to the inactive environment:
   ```bash
   # If blue is active, deploy to green
   ./scripts/deploy-green.sh
   
   # If green is active, deploy to blue
   ./scripts/deploy-blue.sh
   ```

3. Switch traffic after verifying the deployment:
   ```bash
   # Switch to blue
   ./scripts/switch-to-blue.sh
   
   # Switch to green
   ./scripts/switch-to-green.sh
   ```

4. In case of issues, rollback:
   ```bash
   ./scripts/rollback.sh
   ```

## Database Migration

1. Apply migrations:
   ```bash
   ./scripts/migrate-database.sh
   ```

## Health Checks

1. Check application health:
   ```bash
   ./scripts/check-health.sh
   ```

## Troubleshooting

### Common Issues

1. **Database Connection Failure**
   - Check the database credentials in the `.env` file
   - Verify network connectivity to the database
   - Check if the database service is running

2. **Image Pull Failure**
   - Verify Docker registry credentials
   - Check network connectivity to the registry

3. **Application Not Starting**
   - Check application logs: `docker-compose logs app`
   - Verify environment variables
   - Check resource availability on the host

For additional help, refer to the [Troubleshooting Guide](./troubleshooting_guide.md).