# --------
# Optimized Multi-stage Production Dockerfile for CryoProtect v2
# Hardened: non-root user, healthcheck, minimized image size, security best practices
# --------

# --- Builder Stage ---
FROM continuumio/miniconda3:4.12.0 AS builder

WORKDIR /app

# Install mamba, update OS, and clean up in a single layer
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y curl && \
    conda install -y -c conda-forge mamba=1.4.2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/* /var/tmp/*

# Copy only files needed for dependency installation
COPY environment.prod.yml ./environment.yml
COPY requirements_updated.txt ./requirements.txt

# Create conda environment with explicit version pinning and clean up caches
RUN mamba env create -f environment.yml && \
    conda clean -afy && \
    rm -rf /root/.cache /opt/conda/pkgs/* && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    find /opt/conda/ -follow -type f -name '*.pyc' -delete && \
    find /opt/conda/ -follow -type d -name '__pycache__' -exec rm -rf {} +

# Activate environment for subsequent RUN commands
SHELL ["conda", "run", "-n", "cryoprotect", "/bin/bash", "-c"]

# Install pip requirements with explicit version pinning and clean pip cache
RUN pip install --upgrade pip==23.1.2 && \
    pip install -r requirements.txt && \
    pip install --upgrade 'cryptography>=39.0.1' && \
    pip cache purge && \
    find /opt/conda/envs/cryoprotect -follow -type f -name '*.pyc' -delete && \
    find /opt/conda/envs/cryoprotect -follow -type d -name '__pycache__' -exec rm -rf {} +

# Copy only necessary application code
COPY app.py .
COPY config*.py .
COPY api/ ./api/
COPY static/ ./static/
COPY templates/ ./templates/
COPY migrations/models/ ./models/

# Pre-compile Python files to speed up startup
RUN python -m compileall .

# --- Production Image ---
FROM continuumio/miniconda3:4.12.0-slim AS production

# Set working directory
WORKDIR /app

# Install curl for healthcheck and security packages, then clean up
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
        curl \
        ca-certificates \
        tzdata \
        tini && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    # Set strict umask for security
    umask 027

# Copy only the cryoprotect environment from builder (optimized)
COPY --from=builder /opt/conda/envs/cryoprotect /opt/conda/envs/cryoprotect

# Link the environment as the default
RUN ln -s /opt/conda/envs/cryoprotect /opt/conda/envs/default && \
    ln -s /opt/conda/envs/cryoprotect /opt/conda/envs/base && \
    # Remove unnecessary files to reduce image size
    find /opt/conda/envs/cryoprotect -follow -type f -name '*.a' -delete && \
    find /opt/conda/envs/cryoprotect -follow -type f -name '*.js.map' -delete && \
    find /opt/conda/envs/cryoprotect -follow -type f -name 'test_*.py' -delete

# Create non-root user and set up directories with proper permissions
RUN useradd -m -s /bin/bash -u 1000 appuser && \
    mkdir -p /app /run/secrets /app/logs /app/cache && \
    chown -R appuser:appuser /app && \
    chown -R appuser:appuser /run/secrets && \
    chmod 750 /app && \
    chmod 700 /run/secrets

# Copy app code from builder
COPY --from=builder /app /app
RUN chown -R appuser:appuser /app

# Switch to non-root user
USER appuser

# Expose port
EXPOSE 5000

# Set environment variables
ENV FLASK_APP=app.py \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/opt/conda/envs/cryoprotect/bin:$PATH" \
    PYTHONOPTIMIZE=2 \
    PYTHONHASHSEED=random \
    PYTHONFAULTHANDLER=1
    
ARG FLASK_ENV=production
ENV FLASK_ENV=${FLASK_ENV}

# Add entrypoint script
COPY --chown=appuser:appuser docker-entrypoint.sh /app/
RUN chmod +x /app/docker-entrypoint.sh

# Add comprehensive health checks with faster startup detection
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
  CMD curl -f http://localhost:5000/health && \
      curl -f http://localhost:5000/health/liveness && \
      curl -f http://localhost:5000/health/readiness || exit 1

# Set security labels (if supported by the host)
LABEL org.opencontainers.image.vendor="CryoProtect" \
      org.opencontainers.image.title="CryoProtect v2" \
      org.opencontainers.image.description="CryoProtect v2 Production Image" \
      org.opencontainers.image.version="2.0.0" \
      org.opencontainers.image.created="2025-04-30" \
      org.opencontainers.image.source="https://github.com/cryoprotect/v2" \
      security.hardened="true" \
      org.opencontainers.image.licenses="Proprietary" \
      org.opencontainers.image.documentation="https://github.com/cryoprotect/v2/docs"

# Use tini as init to handle signals properly
ENTRYPOINT ["/usr/bin/tini", "--", "/app/docker-entrypoint.sh"]

# Run the application with optimized settings for faster startup
CMD ["/opt/conda/envs/cryoprotect/bin/python", "-O", "app.py"]