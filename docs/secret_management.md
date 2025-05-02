# Secret Management in CryoProtect v2

This document describes the secure secret management system implemented in CryoProtect v2, which follows Docker and container security best practices.

## Overview

CryoProtect v2 uses Docker secrets as the primary mechanism for managing sensitive information. The system is designed to:

1. Keep secrets out of image layers and container environments
2. Support multiple secret providers (Docker Swarm, Kubernetes, AWS, HashiCorp Vault, Azure, GCP)
3. Handle environment-specific secrets (staging, production, blue/green)
4. Support secret rotation and validation
5. Provide robust error handling for missing secrets

## Secret Types

The system manages the following types of secrets:

### Required Secrets

These secrets are essential for the application to function. By default, the application will fail to start if any of these are missing:

- `SUPABASE_URL`: URL for the Supabase instance
- `SUPABASE_KEY`: API key for Supabase
- `SECRET_KEY`: Application secret key for Flask

### Optional Secrets

These secrets enhance functionality but are not strictly required:

- `SUPABASE_SERVICE_KEY`: Service key for Supabase
- `SUPABASE_SERVICE_ROLE_KEY`: Service role key for Supabase
- `SUPABASE_DB_PASSWORD`: Database password for Supabase
- `DATABASE_URL`: Direct database connection URL
- `REDIS_URL`: URL for Redis connection
- `NOTIFY_SMTP_PASSWORD`: Password for SMTP notifications
- `MAIL_PASSWORD`: Password for email service
- `API_KEY`: API key for external services
- `JWT_SECRET`: Secret for JWT token generation
- `ENCRYPTION_KEY`: Key for data encryption
- `BACKUP_ENCRYPTION_KEY`: Key for backup encryption
- `MONITORING_API_KEY`: API key for monitoring services
- `HEALTH_CHECK_TOKEN`: Token for health check endpoints

### Environment-Specific Secrets

Secrets can be prefixed with environment names to override the default values:

- `STAGING_*`: Used in staging environment
- `PRODUCTION_*`: Used in production environment
- `BLUE_*` / `GREEN_*`: Used in blue/green deployment environments

### SSH Keys

Special handling for SSH keys used in deployment:

- `STAGING_SSH_KEY`: SSH key for staging deployment
- `PRODUCTION_SSH_KEY`: SSH key for production deployment

## Using Docker Secrets

### In Docker Swarm

1. Create secrets in Docker Swarm:

```bash
echo "your-secret-value" | docker secret create cryoprotect_supabase_url -
```

2. Reference them in docker-compose.yml:

```yaml
secrets:
  supabase_url:
    external: true
    name: cryoprotect_supabase_url
```

3. Mount them in services:

```yaml
services:
  cryoprotect:
    secrets:
      - source: supabase_url
        target: SUPABASE_URL
        mode: 0400
```

### In Docker Compose (Development)

For development, you can use file-based secrets:

1. Create a `.secrets` directory (add to `.gitignore`)
2. Add secret files (e.g., `.secrets/SUPABASE_URL`)
3. Update docker-compose.yml:

```yaml
secrets:
  supabase_url:
    file: ./.secrets/SUPABASE_URL
```

## External Secret Providers

The system supports fetching secrets from external providers:

### AWS Secrets Manager

Set `EXTERNAL_SECRET_PROVIDER=aws` and provide:
- `/run/secrets/AWS_ACCESS_KEY_ID`
- `/run/secrets/AWS_SECRET_ACCESS_KEY`
- `/run/secrets/AWS_REGION` (optional)
- `AWS_SECRET_PREFIX` environment variable

### HashiCorp Vault

Set `EXTERNAL_SECRET_PROVIDER=vault` and provide:
- `/run/secrets/VAULT_TOKEN`
- `/run/secrets/VAULT_ADDR` (optional)
- `VAULT_SECRET_PATH` environment variable

### Azure Key Vault

Set `EXTERNAL_SECRET_PROVIDER=azure` and provide:
- `/run/secrets/AZURE_CLIENT_ID`
- `/run/secrets/AZURE_CLIENT_SECRET`
- `/run/secrets/AZURE_TENANT_ID`
- `AZURE_KEYVAULT_NAME` environment variable

### Google Cloud Secret Manager

Set `EXTERNAL_SECRET_PROVIDER=gcp` and provide:
- `/run/secrets/GOOGLE_APPLICATION_CREDENTIALS`
- `GCP_PROJECT_ID` environment variable
- `GCP_SECRET_PREFIX` environment variable

## Kubernetes Integration

When running in Kubernetes, the system can automatically fetch secrets from the Kubernetes API:

1. Set `K8S_SECRET_NAMESPACE` and `K8S_SECRET_NAME` environment variables
2. Ensure the pod has appropriate RBAC permissions to read the secret

## Secret Rotation

The system supports secret rotation through:

1. A timestamp file at `/run/secrets/SECRET_ROTATION_TIMESTAMP`
2. The `SECRET_ROTATION_THRESHOLD_DAYS` environment variable
3. The `STRICT_ROTATION_MODE` environment variable (set to `true` to enforce rotation)

## Configuration Options

The following environment variables control the secret management behavior:

- `STRICT_SECRET_MODE`: When `true` (default), the application will fail to start if required secrets are missing
- `STRICT_ENV_MODE`: When `true` (default), the application will fail to start if required environment variables are missing
- `STRICT_ROTATION_MODE`: When `true`, the application will fail to start if secrets are older than the threshold
- `SECRET_ROTATION_THRESHOLD_DAYS`: Number of days after which secrets should be rotated
- `EXTERNAL_SECRET_PROVIDER`: External secret provider to use (`aws`, `vault`, `azure`, `gcp`)
- `AWS_SECRET_PREFIX`: Prefix for AWS secrets
- `VAULT_SECRET_PATH`: Path for HashiCorp Vault secrets
- `AZURE_KEYVAULT_NAME`: Name of Azure Key Vault
- `GCP_PROJECT_ID`: Google Cloud project ID
- `GCP_SECRET_PREFIX`: Prefix for Google Cloud secrets
- `K8S_SECRET_NAMESPACE`: Kubernetes namespace for secrets
- `K8S_SECRET_NAME`: Kubernetes secret name

## Best Practices

1. **Never hardcode secrets** in Dockerfiles, source code, or environment variables
2. **Rotate secrets regularly** using the secret rotation mechanism
3. **Use least privilege** when setting up access to secret providers
4. **Audit secret access** through logging and monitoring
5. **Use different secrets** for different environments
6. **Encrypt secrets at rest** in external secret providers
7. **Limit secret scope** to only the services that need them

## Troubleshooting

If the application fails to start due to missing secrets:

1. Check the logs for specific error messages about which secrets are missing
2. Verify that the secrets are correctly defined in Docker Swarm or the external provider
3. Set `STRICT_SECRET_MODE=false` temporarily to debug (not recommended for production)
4. Ensure the container has the correct permissions to access the secrets

## Security Considerations

- The entrypoint script logs which secrets were loaded (but not their values)
- Secret files are read with minimal permissions (mode 0400)
- SSH keys are stored with proper permissions (mode 0600)
- The container runs as a non-root user
- The `/run/secrets` directory has restricted permissions (mode 0700)