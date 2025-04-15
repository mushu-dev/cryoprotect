FROM continuumio/miniconda3:latest

WORKDIR /app

# Copy environment files
COPY environment.yml .
COPY requirements_updated.txt requirements.txt

# Create conda environment
RUN conda env create -f environment.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "cryoprotect", "/bin/bash", "-c"]

# Install pip requirements
RUN pip install -r requirements.txt

# Copy application code
COPY . .

# Expose port
EXPOSE 5000

# Set environment variables
ENV FLASK_APP=app.py
ENV FLASK_ENV=development

# Run the application
CMD conda run -n cryoprotect python app.py