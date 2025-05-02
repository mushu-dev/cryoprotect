# Docker Image Security Scanning and SBOM Generation

This document describes the security scanning and Software Bill of Materials (SBOM) generation process for Docker images in the CryoProtect v2 project.

## Overview

The CryoProtect v2 project implements comprehensive security scanning and SBOM generation for Docker images as part of the CI/CD pipeline. This ensures that:

1. All Docker images are scanned for vulnerabilities before deployment
2. Critical vulnerabilities are remediated before release
3. A Software Bill of Materials (SBOM) is generated and stored for each release
4. Security artifacts are preserved for audit and compliance purposes

## Security Scanning

### Tools Used

- **Trivy**: A comprehensive vulnerability scanner for containers and other artifacts
- **Custom Scripts**: `scripts/scan_docker_image.sh` and `scripts/scan_docker_image.bat` for cross-platform scanning

### Scanning Process

1. **CI/CD Integration**: Security scanning is integrated into both GitHub Actions and GitLab CI pipelines
2. **Severity Levels**: Images are scanned for CRITICAL, HIGH, and MEDIUM severity vulnerabilities
3. **Build Failure**: The build will fail if CRITICAL vulnerabilities are found
4. **Artifact Storage**: Scan results are stored as artifacts in the CI/CD pipeline and on the production server

### Usage

To scan a Docker image manually:

```bash
# Linux/macOS
./scripts/scan_docker_image.sh --format json --output scan-results.json --severity CRITICAL,HIGH,MEDIUM --exit-code 1 --fail-on CRITICAL your-image:tag

# Windows
scripts\scan_docker_image.bat --format json --output scan-results.json --severity CRITICAL,HIGH,MEDIUM --exit-code 1 --fail-on CRITICAL your-image:tag
```

Options:
- `--format`: Output format (table, json, sarif, cyclonedx)
- `--output`: Output file name
- `--severity`: Severities to scan for (comma-separated)
- `--exit-code`: Exit code when vulnerabilities are found
- `--fail-on`: Fail on specific severity (CRITICAL, HIGH, MEDIUM, LOW, UNKNOWN)

## SBOM Generation

### What is an SBOM?

A Software Bill of Materials (SBOM) is a formal record containing the details and supply chain relationships of various components used in building software. It provides transparency about the components in a software product.

### Tools Used

- **Trivy**: Used to generate SBOMs in CycloneDX format
- **Custom Scripts**: `scripts/generate_sbom.sh` and `scripts/generate_sbom.bat` for cross-platform SBOM generation

### SBOM Generation Process

1. **CI/CD Integration**: SBOM generation is integrated into both GitHub Actions and GitLab CI pipelines
2. **Format**: SBOMs are generated in CycloneDX format, which is a standard for SBOMs
3. **Storage**: SBOMs are stored as artifacts in the CI/CD pipeline and on the production server
4. **Versioning**: SBOMs are versioned with the image tag and timestamp

### Usage

To generate an SBOM manually:

```bash
# Linux/macOS
./scripts/generate_sbom.sh --format cyclonedx --output sbom.json --store-dir ./sbom --version 1.0.0 your-image:tag

# Windows
scripts\generate_sbom.bat --format cyclonedx --output sbom.json --store-dir .\sbom --version 1.0.0 your-image:tag
```

Options:
- `--format`: Output format (cyclonedx, spdx, json)
- `--output`: Output file name
- `--store-dir`: Directory to store SBOM files
- `--version`: Version tag for the SBOM file

## CI/CD Integration

### GitHub Actions

The GitHub Actions workflow includes:

1. Security scanning in the `docker` job
2. SBOM generation in the `docker` and `blue-green-deployment` jobs
3. Storage of security artifacts on the production server

### GitLab CI

The GitLab CI pipeline includes:

1. Security scanning in the `security-scan` job
2. SBOM generation in the `generate-sbom` and `blue-green-deployment` jobs
3. Storage of security artifacts on the production server

## Storage and Retention

Security artifacts are stored in the following locations:

1. **CI/CD Artifacts**: Stored in the CI/CD pipeline for a specified retention period
2. **Production Server**: Stored in `/var/log/cryoprotect/sbom` and `/var/log/cryoprotect/security` directories
3. **Versioned Storage**: Each artifact is versioned with the image tag and timestamp

## Compliance and Auditing

The security scanning and SBOM generation process supports compliance with:

1. **Executive Order 14028**: Requires SBOMs for software used by federal agencies
2. **NIST Secure Software Development Framework (SSDF)**: Recommends security scanning and SBOMs
3. **ISO/IEC 27001**: Supports information security management requirements

## Troubleshooting

If you encounter issues with security scanning or SBOM generation:

1. Ensure Trivy is installed and up to date
2. Check that the Docker image is accessible
3. Verify that the scripts have execute permissions
4. Review the CI/CD logs for error messages

## References

- [Trivy Documentation](https://aquasecurity.github.io/trivy/latest/)
- [CycloneDX Specification](https://cyclonedx.org/specification/overview/)
- [SPDX Specification](https://spdx.dev/specifications/)
- [NIST Secure Software Development Framework](https://csrc.nist.gov/Projects/ssdf)