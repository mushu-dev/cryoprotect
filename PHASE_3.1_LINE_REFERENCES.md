# Phase 3.1: Deployment Infrastructure - Line Number References

This document provides specific line number references for key areas that need modification in larger files, allowing agents to focus precisely on relevant sections without having to read entire files.

## Critical File Line References

### Dockerfile

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Base Image | 1-5 | Update base image selection for security and performance |
| Dependency Installation | 10-25 | Replace with multi-stage build pattern for smaller image |
| User Permissions | 30-40 | Add non-root user for security |
| Entrypoint Configuration | 45-50 | Update entrypoint for proper initialization |
| Health Check | 55-60 | Add Docker health check directives |

### docker-compose.yml

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Service Definition | 5-20 | Update service configuration for production |
| Volume Configuration | 25-35 | Add proper volume mounts for persistence |
| Network Setup | 40-50 | Configure networks for security |
| Resource Limits | 55-65 | Add CPU and memory limits |
| Dependency Config | 70-85 | Configure service dependencies properly |

### .github/workflows/deploy.yml

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Trigger Configuration | 5-15 | Update event triggers for deployment |
| Testing Integration | 20-40 | Add comprehensive testing before deployment |
| Environment Selection | 45-60 | Add environment-specific deployment logic |
| Deployment Commands | 65-85 | Add actual deployment commands |
| Notification Config | 90-110 | Add notification steps for deployment status |

### config.py

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Configuration Class | 10-30 | Update base configuration class |
| Environment Variables | 35-50 | Add all required environment variables |
| Secret Management | 55-70 | Improve secret handling |
| Environment Selection | 75-90 | Update environment selection logic |

## Environment Configuration Files

### config_staging.py

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Class Definition | 5-15 | Create staging-specific configuration class |
| Database Settings | 20-30 | Configure staging database settings |
| Cache Settings | 35-45 | Configure staging cache settings |
| Logging Level | 50-60 | Set appropriate logging levels |

### config_production.py

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Class Definition | 5-15 | Create production-specific configuration class |
| Security Settings | 20-30 | Add production security settings |
| Performance Settings | 35-45 | Add production performance settings |
| Error Handling | 50-60 | Configure production error handling |

## Deployment Scripts

### scripts/deploy_blue_green.sh

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Environment Variables | 5-20 | Define all needed environment variables |
| Blue Environment | 25-45 | Setup blue environment deployment logic |
| Green Environment | 50-70 | Setup green environment deployment logic |
| Health Checking | 75-95 | Implement health check logic |
| Traffic Switching | 100-120 | Implement traffic switch logic |
| Rollback Procedure | 125-145 | Implement rollback logic |

## Implementation Priorities with Line References

1. Update `Dockerfile` multi-stage build (lines 10-25)
2. Improve `.github/workflows/deploy.yml` testing integration (lines 20-40)
3. Enhance `docker-compose.yml` service definition (lines 5-20)
4. Update `config.py` configuration class (lines 10-30)
5. Create `scripts/deploy_blue_green.sh` deployment logic (lines 25-70)

This line-specific guidance will help agents focus precisely on the areas that need work without wasting resources reading through entire files. When working on a specific issue, only request the relevant line range rather than the entire file.