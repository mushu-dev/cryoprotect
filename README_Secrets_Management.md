# CryoProtect v2 - Secrets Management

This document outlines the secure secrets management approach for the CryoProtect v2 project, covering local development, staging, and production environments.

## Overview

Sensitive values (such as API keys, database credentials, and other secrets) should never be stored in plaintext in your codebase or deployment artifacts. This guide explains how we manage secrets securely across different environments.

## Environments and Secrets Management

### Local Development

For local development, we use a `.env` file that is **not** committed to the repository. This file contains environment-specific variables needed for development.

1. Copy the `.env.template` file to `.env`
2. Fill in your local development credentials
3. The `.env` file is listed in `.gitignore` to prevent accidental commits

```bash
# Example for local development
cp .env.template .env
# Edit .env with your credentials
```

### Staging and Production

For staging and production environments, we use Docker Secrets to securely manage sensitive values:

1. Docker Secrets are stored securely in the Docker Swarm infrastructure
2. Secrets are mounted as files in the container at runtime
3. Our application reads secrets from the mounted files
4. CI/CD pipelines create and manage these secrets automatically

## Docker Secrets

### How Docker Secrets Work

Docker Secrets provide a secure way to manage sensitive data:

1. Secrets are encrypted at rest in the Docker Swarm manager nodes
2. Secrets are only transmitted to nodes that run services that need them
3. Secrets are mounted as files in `/run/secrets/` inside the container
4. Our `docker-entrypoint.sh` script reads these files and sets environment variables

### Creating Docker Secrets Manually

If you need to create Docker Secrets manually:

```bash
# Create a secret
echo "my_secret_value" | docker secret create cryoprotect_supabase_key -

# Update an existing secret (must delete first)
docker secret rm cryoprotect_supabase_key
echo "new_secret_value" | docker secret create cryoprotect_supabase_key -

# List all secrets
docker secret ls
```

### Required Secrets for Production

The following secrets are required for production deployments:

| Secret Name | Description |
|-------------|-------------|
| cryoprotect_supabase_url | Supabase URL for the production environment |
| cryoprotect_supabase_key | Supabase service role key for production |
| cryoprotect_secret_key | Secret key for session encryption and JWT signing |
| cryoprotect_redis_url | Redis URL for caching (optional) |

## CI/CD Secrets Management

Our CI/CD pipelines (GitHub Actions) use GitHub Secrets to store sensitive values:

1. Secrets are stored securely in GitHub and are encrypted
2. Secrets are only exposed to authorized workflows
3. Our workflows use these secrets to create Docker Secrets on the deployment servers

### Required GitHub Secrets

The following GitHub Secrets are required for CI/CD:

| Secret Name | Description |
|-------------|-------------|
| STAGING_SUPABASE_URL | Supabase URL for staging |
| STAGING_SUPABASE_KEY | Supabase service role key for staging |
| STAGING_SECRET_KEY | Secret key for staging |
| STAGING_SSH_HOST | SSH host for staging deployment |
| STAGING_SSH_USER | SSH user for staging deployment |
| STAGING_SSH_KEY | SSH private key for staging deployment |
| PRODUCTION_SUPABASE_URL | Supabase URL for production |
| PRODUCTION_SUPABASE_KEY | Supabase service role key for production |
| PRODUCTION_SECRET_KEY | Secret key for production |
| PRODUCTION_SSH_HOST | SSH host for production deployment |
| PRODUCTION_SSH_USER | SSH user for production deployment |
| PRODUCTION_SSH_KEY | SSH private key for production deployment |
| REDIS_URL | Redis URL for production caching |

## Local Docker Development with Secrets

For local Docker development, you can use a development-specific docker-compose file that uses environment variables instead of secrets:

```bash
# Use the development docker-compose file
docker-compose -f docker-compose.dev.yml up
```

Or you can create local secrets for development:

```bash
# Create local secrets for development
echo "http://localhost:54321" | docker secret create cryoprotect_supabase_url -
echo "your-local-key" | docker secret create cryoprotect_supabase_key -
echo "your-local-secret-key" | docker secret create cryoprotect_secret_key -

# Run with secrets
USE_EXTERNAL_SECRETS=true docker-compose up
```

## Best Practices

1. **Never commit secrets** to the repository, even in configuration files
2. **Rotate secrets regularly**, especially for production environments
3. **Use different secrets** for each environment (development, staging, production)
4. **Limit access** to production secrets to only those who need them
5. **Audit secret usage** regularly to ensure they are being used properly
6. **Use the principle of least privilege** when creating service accounts and API keys

## Troubleshooting

### Checking if Secrets are Mounted

To check if secrets are properly mounted in a container:

```bash
docker exec -it container_name ls -la /run/secrets/
```

### Verifying Secret Values (for debugging only)

In development environments only, you can verify secret values:

```bash
docker exec -it container_name cat /run/secrets/cryoprotect_supabase_url
```

**Warning**: Never do this in production environments!

### Common Issues

1. **Secret not found**: Ensure the secret is created before starting the container
2. **Permission denied**: Check that the container user has permission to read the secret
3. **Secret not loaded**: Verify that the `docker-entrypoint.sh` script is properly loading the secret

## Migration Guide

If you're migrating from the old approach (environment variables) to Docker Secrets:

1. Create Docker Secrets for all sensitive values
2. Update your docker-compose.yml file to use secrets
3. Deploy the updated application
4. Verify that the application is working correctly with secrets
5. Remove any plaintext environment variables from your deployment scripts

## Additional Resources

- [Docker Secrets Documentation](https://docs.docker.com/engine/swarm/secrets/)
- [GitHub Actions Secrets Documentation](https://docs.github.com/en/actions/security-guides/encrypted-secrets)