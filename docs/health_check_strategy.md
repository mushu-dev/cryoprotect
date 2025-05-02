# CryoProtect v2 Health Check Strategy

This document outlines the comprehensive health check strategy implemented for CryoProtect v2 in Docker production environments.

## Overview

The health check strategy provides multiple levels of health monitoring:

1. **Container-level health checks** - Docker HEALTHCHECK directives
2. **Application-level health endpoints** - Flask API endpoints
3. **Service-level health checks** - Checks for supporting services
4. **Load balancer health monitoring** - NGINX health check integration
5. **Blue/Green deployment health verification** - Deployment-specific health checks

## Health Check Endpoints

### Main Application Health Endpoints

| Endpoint | Purpose | HTTP Status Codes |
|----------|---------|-------------------|
| `/health` | Overall system health | 200 (OK), 207 (Multi-Status), 500 (Error) |
| `/health/liveness` | Simple liveness probe | 200 (OK), 500 (Error) |
| `/health/readiness` | Readiness to serve traffic | 200 (OK), 503 (Not Ready) |
| `/health/startup` | Startup completion status | 200 (OK), 503 (Starting) |
| `/health/backup` | Backup system status | 200 (OK), 503 (Error) |

### NGINX Health Endpoints

| Endpoint | Purpose | HTTP Status Codes |
|----------|---------|-------------------|
| `:8080/health` | NGINX internal health | 200 (OK) |
| `:8080/health/blue` | Blue environment health | 200 (OK), 503 (Down) |
| `:8080/health/green` | Green environment health | 200 (OK), 503 (Down) |
| `:8080/health/blue/liveness` | Blue environment liveness | 200 (OK), 503 (Down) |
| `:8080/health/blue/readiness` | Blue environment readiness | 200 (OK), 503 (Down) |
| `:8080/health/green/liveness` | Green environment liveness | 200 (OK), 503 (Down) |
| `:8080/health/green/readiness` | Green environment readiness | 200 (OK), 503 (Down) |
| `:8080/health/dashboard` | Health check dashboard | 200 (OK) |

## Health Check Response Format

### Main Health Check (`/health`)

```json
{
  "status": "ok|degraded|error",
  "version": "2.0.0",
  "deployment": "blue|green|production",
  "timestamp": "2025-04-25T13:56:00.000Z",
  "services": {
    "database": "connected|unknown|error",
    "redis": "connected|unreachable|not_configured|error",
    "disk": "ok|low|error",
    "memory": "ok|high|not_available",
    "environment": "ok|missing: VAR1, VAR2"
  }
}
```

### Liveness Check (`/health/liveness`)

```json
{
  "status": "alive",
  "timestamp": "2025-04-25T13:56:00.000Z"
}
```

### Readiness Check (`/health/readiness`)

```json
{
  "status": "ready|not_ready",
  "reason": "Required services not available", // Only present if not_ready
  "timestamp": "2025-04-25T13:56:00.000Z"
}
```

## Docker Health Check Configuration

Health checks are configured in:

1. **Dockerfile** - HEALTHCHECK directive checks `/health`, `/health/liveness`, and `/health/readiness`
2. **docker-compose.yml** - Service-specific health checks for all services
3. **docker-compose.blue-green.yml** - Blue/Green deployment health checks

## Health Check Dashboard

A visual health check dashboard is available at `:8080/health/dashboard`. This dashboard provides:

- Current status of blue and green environments
- Active environment information
- Real-time health status updates

Access credentials:
- Username: `admin`
- Password: Contact system administrator

## Health Check Failure Actions

| Service | Failure Action | Recovery Action |
|---------|----------------|-----------------|
| Main application | Container restart | Automatic via Docker restart policy |
| Database | Connection retry, fallback | Automatic reconnection |
| Redis | Connection retry | Automatic reconnection |
| NGINX | Service restart | Automatic via Docker restart policy |
| Blue/Green | Traffic routing to healthy environment | Automatic via NGINX |

## Blue/Green Deployment Health Verification

The `cutover.sh` script verifies health before switching traffic:

```bash
./cutover.sh --target blue|green [--force] [--timeout 300] [--interval 5]
```

The script:
1. Checks target environment health via `:8080/health/blue` or `:8080/health/green`
2. Only performs cutover if health check passes (unless `--force` is used)
3. Verifies successful cutover

## Monitoring Integration

Health check data is integrated with:

1. **Prometheus** - Health metrics exposed via `/metrics`
2. **Grafana** - Health dashboards
3. **ELK Stack** - Health check logs

## Implementing Custom Health Checks

To add custom health checks:

1. Add check logic to the `/health` endpoint in `app.py`
2. Update the health response JSON structure
3. Configure appropriate failure thresholds

## Troubleshooting

If health checks are failing:

1. Check individual service health endpoints
2. Review logs for specific errors
3. Verify connectivity between services
4. Check resource utilization (memory, disk, CPU)
5. Verify environment variables and configuration