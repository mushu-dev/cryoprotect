# Docker Production Architecture Specification

## Overview
This document specifies the optimized Dockerfile, enhanced docker-compose files, and comprehensive container health check strategy for CryoProtect v2 production deployments. The design addresses image size, security, dependency minimization, build caching, resource allocation, secret management, health checks, CI/CD compatibility, and blue/green deployment.

---

## 1. Dockerfile Optimization

### Multi-Stage Build
- Use a builder stage to install dependencies and build the environment.
- Copy only necessary files to the final image.
- Clean up caches and temporary files in each stage.

### Image Size Minimization
- Use Miniconda as a base for Python/conda environments.
- Only copy production dependencies (`environment.prod.yml`, `requirements_updated.txt`).
- Remove build tools, caches, and unused packages after installation.
- Use `.dockerignore` to exclude unnecessary files from the build context.

### Security Hardening
- Run as a non-root user (`appuser`).
- Set strict permissions on `/app` and `/run/secrets`.
- Use Docker secrets for sensitive data.
- Avoid exposing unnecessary ports or files.
- Pin dependency versions in conda and pip.
- Regularly update base images and dependencies.

### Dependency Minimization
- Only install production dependencies.
- Use `conda clean -afy` and `pip cache purge` to remove caches.
- Remove build-only dependencies after build stage.

### Build Caching
- Order layers to maximize cache hits (copy requirements before app code).
- Separate environment and code copy steps.
- Use multi-stage to avoid leaking build tools into production image.

---

## 2. docker-compose.yml Enhancements

### Resource Allocation
- Use `deploy.resources.limits` and `reservations` for memory and CPU.
- Set conservative defaults, allow override via environment variables.

### Secret Management
- Use Docker secrets for all sensitive values.
- Externalize secrets for production (integrate with secret managers if possible).
- Load secrets at runtime via entrypoint script.

### Logging & Monitoring
- Centralized logging (ELK stack).
- Prometheus and Grafana for metrics.
- Healthchecks for supporting services (Elasticsearch, Logstash, Kibana).

### Volumes & Persistence
- Use named volumes for persistent data (logs, backups, database, etc.).
- Mount only required directories.

### Profiles & Environments
- Use Docker Compose profiles for prod, dev, staging, monitoring, backup, etc.
- Development service mounts source code for live reload.

---

## 3. docker-compose.blue-green.yml (Blue/Green Deployment)

### Load Balancing
- Use NGINX as a reverse proxy to route traffic to blue/green environments.

### Blue/Green Services
- Deploy two parallel environments (`cryoprotect-blue`, `cryoprotect-green`) with different image tags.
- Use environment variable `DEPLOYMENT_COLOR` for color awareness.

### Health Checks
- Each environment exposes a `/health` endpoint.
- Healthchecks in compose ensure only healthy containers receive traffic.

### Cutover Process
- Update NGINX config to switch traffic between blue and green.
- Rollback by switching back if healthchecks fail.

---

## 4. Entrypoint & Secret Management

- `docker-entrypoint.sh` loads secrets from `/run/secrets` into environment variables.
- Supports environment-specific and SSH secrets.
- No secrets are hardcoded or baked into images.

---

## 5. Health Check Strategy

- Dockerfile: `HEALTHCHECK` on port 5000 (root or `/health` endpoint).
- Compose: Service-level healthchecks for blue/green (`/health`).
- Supporting services (Elasticsearch, Logstash, Kibana) have their own healthchecks.
- Consider adding readiness/liveness endpoints for more granular checks.

---

## 6. CI/CD & Build/Deploy Integration

- Use build arguments and environment variables for dynamic configuration.
- Build and push tagged images for blue/green deployment.
- Inject secrets at deploy time (never in CI logs or images).
- Run tests in CI pipeline before deployment.
- Support for multi-arch builds if needed.

---

## 7. Security & Compliance

- Use non-root user, strict permissions, and Docker secrets.
- Pin all dependencies and scan images for vulnerabilities (e.g., Trivy, Snyk).
- Generate and store SBOM (Software Bill of Materials).
- Consider rootless Docker for additional isolation.
- Regularly update base images and dependencies.

---

## 8. Recommendations for Further Improvement

- Implement SBOM and automated image scanning in CI.
- Use multi-arch builds for ARM/x86 compatibility.
- Add readiness/liveness endpoints in the app for more robust healthchecks.
- Integrate with external secret managers (AWS Secrets Manager, HashiCorp Vault).
- Use rootless Docker where possible.
- Automate blue/green cutover and rollback in CI/CD pipeline.

---

## 9. Acceptance Criteria

- Dockerfile and compose files follow all above best practices.
- All secrets are managed via Docker secrets and never hardcoded.
- Healthchecks are robust and cover all critical services.
- Blue/green deployment is fully supported and documented.
- CI/CD pipeline can build, test, and deploy using these artifacts.
- Security scanning and SBOM generation are integrated or planned.

---

## References

- [Docker Best Practices](https://docs.docker.com/develop/develop-images/dockerfile_best-practices/)
- [Docker Compose Secrets](https://docs.docker.com/compose/use-secrets/)
- [Blue/Green Deployment with Docker](https://docs.docker.com/samples/blue-green-deployment/)
- [OWASP Docker Security Cheat Sheet](https://cheatsheetseries.owasp.org/cheatsheets/Docker_Security_Cheat_Sheet.html)