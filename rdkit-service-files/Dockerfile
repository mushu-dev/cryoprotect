FROM continuumio/miniconda3:4.12.0

WORKDIR /app

# Install mamba for faster conda installs and essential tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    ca-certificates && \
    conda install -y -c conda-forge mamba && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy environment file and requirements
COPY environment.yml .
COPY requirements.txt .

# Create the conda environment with RDKit
RUN mamba env create -f environment.yml && \
    conda clean -afy && \
    rm -rf /root/.cache
    
# Install Python requirements
RUN /opt/conda/envs/cryoprotect/bin/pip install -r requirements.txt

# Copy application code
COPY app/ ./app/
COPY docker-entrypoint-rdkit.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Expose port for the Flask application
EXPOSE 8080

# Set up environment activation in the entrypoint
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["gunicorn", "--bind", "0.0.0.0:8080", "--workers", "2", "app.app:app"]