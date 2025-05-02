# Guide for Adding New Deployment Targets to CryoProtect v2 CI/CD Pipeline

This document provides step-by-step instructions for adding new deployment targets to the CryoProtect v2 CI/CD pipeline.

## Overview

The CryoProtect v2 CI/CD pipeline currently supports three deployment environments:
- Development
- Staging
- Production

You may need to add additional deployment targets for various purposes, such as:
- Feature-specific environments
- QA testing environments
- Demo environments
- Customer-specific deployments

## Prerequisites

Before adding a new deployment target, ensure you have:

1. Access to the GitHub repository settings
2. SSH access to the target server
3. Docker and docker-compose installed on the target server
4. Necessary credentials for the new environment (Supabase, Redis, etc.)

## Step 1: Create a GitHub Environment

1. Navigate to your repository on GitHub
2. Go to Settings → Environments
3. Click "New environment"
4. Enter a name for your environment (e.g., "qa", "demo", "customer1")
5. Configure environment protection rules if needed (e.g., required reviewers)
6. Click "Configure environment"

## Step 2: Add Environment Secrets

Add the following secrets to your new environment:

1. `<ENV>_SSH_HOST`: SSH host for the target server
2. `<ENV>_SSH_USER`: SSH username for the target server
3. `<ENV>_SSH_KEY`: SSH private key for the target server
4. `<ENV>_SUPABASE_URL`: Supabase URL for the environment
5. `<ENV>_SUPABASE_KEY`: Supabase key for the environment
6. `<ENV>_SECRET_KEY`: Secret key for the environment
7. Any other environment-specific secrets needed

Replace `<ENV>` with an uppercase prefix for your environment (e.g., "QA", "DEMO", "CUSTOMER1").

## Step 3: Update Workflow Dispatch Inputs

Modify the workflow dispatch inputs in `.github/workflows/deploy.yml` to include your new environment:

```yaml
workflow_dispatch:
  inputs:
    environment:
      description: 'Environment to deploy to'
      required: true
      default: 'staging'
      type: choice
      options:
        - dev
        - staging
        - production
        - your-new-environment  # Add your new environment here
```

## Step 4: Add a New Deployment Job

Add a new job to the `.github/workflows/deploy.yml` file following this template:

```yaml
deploy-your-environment:
  name: Deploy to Your Environment
  runs-on: ubuntu-latest
  needs: [docker, version]
  if: github.event_name == 'push' && github.ref == 'refs/heads/your-branch' || (github.event_name == 'workflow_dispatch' && github.event.inputs.environment == 'your-environment')
  environment: your-environment-name
  
  steps:
    - name: Checkout code
      uses: actions/checkout@v3
    
    - name: Set up Docker Buildx
      uses: docker/setup-buildx-action@v2
    
    - name: Login to GitHub Container Registry
      uses: docker/login-action@v2
      with:
        registry: ghcr.io
        username: ${{ github.actor }}
        password: ${{ secrets.GITHUB_TOKEN }}
    
    - name: Build and push Docker image
      uses: docker/build-push-action@v4
      with:
        context: .
        push: true
        build-args: |
          FLASK_ENV=your-environment
          APP_VERSION=${{ needs.version.outputs.version }}
        tags: |
          ghcr.io/${{ github.repository }}/cryoprotect:your-environment
          ghcr.io/${{ github.repository }}/cryoprotect:${{ needs.version.outputs.version }}-your-environment
        cache-from: type=registry,ref=ghcr.io/${{ github.repository }}/cryoprotect:your-environment-cache
        cache-to: type=registry,ref=ghcr.io/${{ github.repository }}/cryoprotect:your-environment-cache,mode=max
    
    - name: Deploy to your environment
      env:
        YOUR_ENV_SSH_HOST: ${{ secrets.YOUR_ENV_SSH_HOST }}
        YOUR_ENV_SSH_USER: ${{ secrets.YOUR_ENV_SSH_USER }}
        YOUR_ENV_SSH_KEY: ${{ secrets.YOUR_ENV_SSH_KEY }}
        APP_VERSION: ${{ needs.version.outputs.version }}
      run: |
        echo "Deploying to your environment..."
        # Write SSH key to file
        mkdir -p ~/.ssh
        echo "$YOUR_ENV_SSH_KEY" > ~/.ssh/your_env_key
        chmod 600 ~/.ssh/your_env_key
        
        # Create Docker secrets on target server
        ssh -i ~/.ssh/your_env_key -o StrictHostKeyChecking=no $YOUR_ENV_SSH_USER@$YOUR_ENV_SSH_HOST "
          # Create Docker secrets for sensitive values
          echo '${{ secrets.YOUR_ENV_SUPABASE_URL }}' | docker secret create cryoprotect_supabase_url - || docker secret rm cryoprotect_supabase_url && echo '${{ secrets.YOUR_ENV_SUPABASE_URL }}' | docker secret create cryoprotect_supabase_url - &&
          echo '${{ secrets.YOUR_ENV_SUPABASE_KEY }}' | docker secret create cryoprotect_supabase_key - || docker secret rm cryoprotect_supabase_key && echo '${{ secrets.YOUR_ENV_SUPABASE_KEY }}' | docker secret create cryoprotect_supabase_key - &&
          echo '${{ secrets.YOUR_ENV_SECRET_KEY }}' | docker secret create cryoprotect_secret_key - || docker secret rm cryoprotect_secret_key && echo '${{ secrets.YOUR_ENV_SECRET_KEY }}' | docker secret create cryoprotect_secret_key - &&
          
          # Pull the latest image
          docker pull ghcr.io/${{ github.repository }}/cryoprotect:your-environment &&
          
          # Set non-sensitive environment variables
          export FLASK_ENV=your-environment &&
          export USE_EXTERNAL_SECRETS=true &&
          export APP_VERSION=$APP_VERSION &&
          
          # Deploy with docker-compose
          docker-compose -f docker-compose.your-env.yml up -d
        "
    
    - name: Notify deployment status
      if: always()
      uses: rtCamp/action-slack-notify@v2
      env:
        SLACK_WEBHOOK: ${{ secrets.SLACK_WEBHOOK }}
        SLACK_CHANNEL: deployments
        SLACK_COLOR: ${{ job.status }}
        SLACK_TITLE: Your Environment Deployment v${{ needs.version.outputs.version }}
        SLACK_MESSAGE: "Deployment to your environment ${{ job.status == 'success' && 'succeeded' || 'failed' }}! :rocket:"
        SLACK_FOOTER: 'CryoProtect v2 CI/CD Pipeline'
    
    - name: Send email notification
      if: always()
      uses: dawidd6/action-send-mail@v3
      with:
        server_address: ${{ secrets.MAIL_SERVER }}
        server_port: ${{ secrets.MAIL_PORT }}
        username: ${{ secrets.MAIL_USERNAME }}
        password: ${{ secrets.MAIL_PASSWORD }}
        subject: "CryoProtect v2 - Your Environment Deployment ${{ job.status }}"
        body: |
          Your Environment Deployment Status: ${{ job.status }}
          Version: ${{ needs.version.outputs.version }}
          Commit: ${{ github.sha }}
          Deployed by: ${{ github.actor }}
          
          See details: ${{ github.server_url }}/${{ github.repository }}/actions/runs/${{ github.run_id }}
        to: ${{ secrets.DEVOPS_EMAIL }},${{ secrets.YOUR_ENV_EMAIL }}
        from: CryoProtect CI/CD <${{ secrets.MAIL_USERNAME }}>
```

Replace all instances of:
- `your-environment` with your environment name (e.g., "qa", "demo", "customer1")
- `YOUR_ENV` with your environment prefix in uppercase (e.g., "QA", "DEMO", "CUSTOMER1")
- `your-branch` with the branch that should trigger deployment to this environment (if applicable)
- `docker-compose.your-env.yml` with the appropriate Docker Compose file for this environment

## Step 5: Create a Docker Compose File

Create a Docker Compose file for your new environment:

1. Create a file named `docker-compose.your-env.yml` in the root of your repository
2. Configure it appropriately for your environment
3. Ensure it uses the correct Docker image tag and environment variables

## Step 6: Update Deployment Dependencies

If your new environment should be part of the deployment pipeline (e.g., deploy to QA after staging but before production), update the `needs` parameter of the relevant jobs:

```yaml
deploy-production:
  name: Deploy to Production
  runs-on: ubuntu-latest
  needs: [deploy-your-environment, version]  # Add your environment here
  # ... rest of the job configuration
```

## Step 7: Test the New Deployment Target

1. Commit and push your changes
2. Go to the Actions tab in your GitHub repository
3. Select the "CryoProtect v2 Deployment Pipeline" workflow
4. Click "Run workflow"
5. Select your new environment from the dropdown
6. Click "Run workflow"
7. Monitor the workflow execution and verify that deployment succeeds

## Step 8: Update Documentation

Update the following documentation files to include your new environment:

1. `.github/workflows/README.md`: Add your environment to the Environment Configuration section
2. `.github/workflows/TROUBLESHOOTING.md`: Add any environment-specific troubleshooting tips

## Conclusion

You have successfully added a new deployment target to the CryoProtect v2 CI/CD pipeline. The pipeline can now deploy to your new environment either manually via workflow dispatch or automatically based on your configured triggers.