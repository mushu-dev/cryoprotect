# CryoProtect v2 CI/CD Pipeline Documentation

This document provides an overview of the CI/CD pipeline for the CryoProtect v2 project.

## Pipeline Overview

The CI/CD pipeline is implemented using GitHub Actions and consists of the following stages:

1. **Testing**: Runs unit, integration, and security tests
2. **Docker Build**: Builds and tests the Docker image
3. **Version Determination**: Determines the version number for the build
4. **Deployment**: Deploys to development, staging, and production environments
5. **Blue/Green Deployment**: Implements zero-downtime deployments in production

## Pipeline Triggers

The pipeline can be triggered in the following ways:

- **Push to main branch**: Triggers deployment to staging and production
- **Push to staging branch**: Triggers deployment to staging only
- **Push to dev branch**: Triggers deployment to development only
- **Push of version tags** (v*.*.*)**: Triggers a release deployment
- **Manual trigger**: Using the workflow_dispatch event with environment selection

## Environment Configuration

The pipeline uses GitHub Environments for deployment. Before running the pipeline, ensure the following environments are configured in your GitHub repository settings:

1. **dev**: Development environment
2. **staging**: Staging environment
3. **production**: Production environment

For each environment, configure the appropriate secrets:

### Development Environment Secrets
- `DEV_SSH_HOST`: SSH host for development server
- `DEV_SSH_USER`: SSH username for development server
- `DEV_SSH_KEY`: SSH private key for development server
- `DEV_SUPABASE_URL`: Supabase URL for development
- `DEV_SUPABASE_KEY`: Supabase key for development
- `DEV_SECRET_KEY`: Secret key for development

### Staging Environment Secrets
- `STAGING_SSH_HOST`: SSH host for staging server
- `STAGING_SSH_USER`: SSH username for staging server
- `STAGING_SSH_KEY`: SSH private key for staging server
- `STAGING_SUPABASE_URL`: Supabase URL for staging
- `STAGING_SUPABASE_KEY`: Supabase key for staging
- `STAGING_SECRET_KEY`: Secret key for staging
- `REDIS_URL`: Redis URL for staging

### Production Environment Secrets
- `PRODUCTION_SSH_HOST`: SSH host for production server
- `PRODUCTION_SSH_USER`: SSH username for production server
- `PRODUCTION_SSH_KEY`: SSH private key for production server
- `PRODUCTION_SUPABASE_URL`: Supabase URL for production
- `PRODUCTION_SUPABASE_KEY`: Supabase key for production
- `PRODUCTION_SECRET_KEY`: Secret key for production
- `REDIS_URL`: Redis URL for production

### Notification Secrets
- `SLACK_WEBHOOK`: Webhook URL for Slack notifications
- `MAIL_SERVER`: SMTP server for email notifications
- `MAIL_PORT`: SMTP port for email notifications
- `MAIL_USERNAME`: SMTP username for email notifications
- `MAIL_PASSWORD`: SMTP password for email notifications
- `DEVOPS_EMAIL`: Email address for DevOps team
- `QA_EMAIL`: Email address for QA team
- `PRODUCT_EMAIL`: Email address for Product team

## Troubleshooting

### Common Pipeline Issues

1. **Environment Issues**:
   - Ensure all required secrets are configured in GitHub repository settings
   - Verify that GitHub environments (dev, staging, production) are properly configured

2. **SSH Connection Issues**:
   - Verify SSH keys are correctly formatted and have proper permissions
   - Ensure the SSH user has appropriate permissions on the target server

3. **Docker Build Issues**:
   - Check Dockerfile for errors
   - Verify that all required dependencies are available

4. **Deployment Issues**:
   - Verify that the target server has Docker and docker-compose installed
   - Check that the appropriate docker-compose files exist on the target server

## Adding New Deployment Targets

To add a new deployment target:

1. Create a new GitHub Environment in repository settings
2. Add the required secrets for the new environment
3. Create a new job in the workflow file following the pattern of existing deployment jobs
4. Configure the appropriate conditions for when this deployment should run
5. Update the workflow_dispatch inputs to include the new environment option

## Blue/Green Deployment

The production deployment uses a blue/green strategy for zero-downtime deployments:

1. The pipeline determines which environment (blue or green) is currently active
2. It deploys to the inactive environment
3. It runs health checks to verify the deployment
4. If successful, it switches traffic to the newly deployed environment

This requires the following scripts to be available on the production server:
- `scripts/deploy-blue.sh`
- `scripts/deploy-green.sh`
- `scripts/check-health.sh`
- `scripts/switch-to-blue.sh`
- `scripts/switch-to-green.sh`