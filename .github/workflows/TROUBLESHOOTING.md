# CryoProtect v2 CI/CD Pipeline Troubleshooting Guide

This document provides guidance for troubleshooting common issues with the CryoProtect v2 CI/CD pipeline.

## Environment Configuration Issues

### Invalid Environment Names

If you see errors like `Value 'dev' is not valid` or `Value 'staging' is not valid`, ensure that:

1. GitHub Environments are properly configured in your repository settings
2. Navigate to your repository → Settings → Environments
3. Create environments with the exact names:
   - `development`
   - `staging`
   - `production`

### Missing Secrets

If you see errors related to missing secrets:

1. Check that all required secrets are configured in each environment
2. Navigate to your repository → Settings → Environments → [Environment Name] → Environment secrets
3. Add the required secrets as listed in the README.md file

## SSH Connection Issues

If deployment fails with SSH connection errors:

1. Verify that the SSH keys are correctly formatted and have proper permissions
2. Ensure the SSH user has appropriate permissions on the target server
3. Check that the target server is reachable from GitHub Actions
4. Verify that the SSH host, user, and key are correctly configured in the environment secrets

## Docker Build Issues

If Docker build fails:

1. Check the Dockerfile for errors
2. Verify that all required dependencies are available
3. Check that the Docker build context is correct
4. Ensure Docker Buildx is properly configured

## Deployment Issues

### Docker Compose Errors

If deployment fails with Docker Compose errors:

1. Verify that Docker and docker-compose are installed on the target server
2. Check that the appropriate docker-compose files exist on the target server
3. Ensure the Docker Compose file is compatible with the installed Docker Compose version
4. Verify that the Docker Compose file references the correct image names and tags

### Docker Secrets Errors

If deployment fails with Docker secrets errors:

1. Verify that Docker Swarm is initialized on the target server
2. Check that the user has permissions to create Docker secrets
3. Ensure that the secrets are correctly formatted

## Blue/Green Deployment Issues

If blue/green deployment fails:

1. Verify that the required scripts exist on the production server:
   - `scripts/deploy-blue.sh`
   - `scripts/deploy-green.sh`
   - `scripts/check-health.sh`
   - `scripts/switch-to-blue.sh`
   - `scripts/switch-to-green.sh`
2. Check that the scripts have execute permissions
3. Ensure that the Nginx configuration is correctly set up for blue/green deployment
4. Verify that the health checks are correctly configured

## Notification Issues

### Slack Notification Errors

If Slack notifications fail:

1. Verify that the Slack webhook URL is correctly configured
2. Check that the Slack channel exists and the webhook has permissions to post to it
3. Ensure the Slack webhook is active and not revoked

### Email Notification Errors

If email notifications fail:

1. Verify that the SMTP server details are correctly configured
2. Check that the SMTP credentials are valid
3. Ensure the email addresses are correctly formatted
4. Check if the SMTP server requires additional authentication or has rate limits

## Version Determination Issues

If version determination fails:

1. Verify that Git is correctly configured in the GitHub Actions runner
2. Check that the repository has at least one tag for proper versioning
3. Ensure the Git history is properly fetched with `fetch-depth: 0`

## Testing Issues

If tests fail:

1. Check the test logs for specific error messages
2. Verify that all dependencies are correctly installed
3. Ensure the test environment is correctly configured
4. Check if tests are timing out due to performance issues

## GitHub Actions Runner Issues

If the GitHub Actions runner fails:

1. Check the GitHub Actions status page for any ongoing issues
2. Verify that the runner has sufficient resources (memory, disk space)
3. Check if the workflow is hitting GitHub Actions usage limits
4. Ensure the workflow is not running for too long (max 6 hours)

## Debugging Tips

1. Add debug steps to the workflow to print environment variables and configuration
2. Use the `actions/upload-artifact` action to save logs and other artifacts for inspection
3. Add timeouts to steps to prevent workflows from hanging indefinitely
4. Use the `continue-on-error` option for non-critical steps to prevent the entire workflow from failing