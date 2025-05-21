# CryoProtect Fedora OS Setup and Configuration Prompt

## Project Context

CryoProtect is a Flask-based molecular analysis application that integrates with Supabase for data storage and RDKit for molecular calculations. The project has been developed primarily on Windows and now needs to be properly set up on Fedora Linux. This prompt outlines the comprehensive steps needed to initialize and configure this project, ensuring all dependencies are properly installed and all Fedora-specific requirements are met.

## Required Tasks

### 1. Create Fedora-Specific Installation Script

Create a specialized Fedora installation script (`fedora_setup.sh`) that:
- Installs all required system dependencies via DNF
- Sets up PostgreSQL with proper SELinux configurations
- Configures the proper environment for RDKit
- Handles Fedora-specific requirements (firewall, SELinux, etc.)
- Automatically detects and adapts to different Fedora versions

### 2. Environment Configuration

Create a complete Fedora-ready `.env` file from the template:
- Set proper database connection parameters
- Configure Supabase integration
- Set security keys and authentication parameters
- Ensure logging paths are compatible with Fedora's filesystem structure
- Add any Fedora-specific environment variables

### 3. Update `.gitignore` for Fedora

Enhance the existing `.gitignore` file to include:
- Fedora-specific temporary files
- SELinux context files
- DNF/RPM related artifacts
- Fedora-specific IDE paths
- Any other Linux-specific files that should be excluded

### 4. Create Fedora RDKit Compatibility Layer

Develop solutions for RDKit compatibility on Fedora:
- Create a Fedora-specific RDKit installation script
- Update RDKit visualization modules for X11 compatibility on Fedora
- Provide fixes for any Fedora-specific font rendering issues
- Ensure proper library paths for RDKit dependencies

### 5. Database Setup and Migration

Create a comprehensive database setup script for Fedora that:
- Sets up PostgreSQL with proper SELinux contexts
- Configures PostgreSQL authentication for local development
- Applies all migrations properly
- Sets up default data with correct permissions
- Creates separate test and development databases

### 6. Docker Configuration for Fedora

Update Docker configuration to work seamlessly on Fedora:
- Ensure Docker installation script works with Fedora's security policies
- Update Docker volume paths for SELinux compatibility
- Configure Docker networking properly for Fedora's firewall
- Optimize Dockerfile to use Fedora-compatible base images when appropriate

### 7. SELinux Configuration

Create a comprehensive SELinux policy setup for the project:
- Write a script to set proper SELinux contexts for all project directories
- Create SELinux policy modules for the application
- Configure SELinux contexts for database connections
- Document SELinux troubleshooting steps specific to the project

### 8. Create Comprehensive Verification Tests

Develop a verification script specific to Fedora that checks:
- All system dependencies are correctly installed
- Database connection works properly
- RDKit is functioning with proper visualizations
- Docker containers run correctly
- SELinux permissions are correctly set
- Firewall configurations allow necessary connections
- X11 forwarding works for molecular visualization

### 9. Develop Fedora-Specific Documentation

Create comprehensive documentation for Fedora users:
- Installation guide specific to Fedora
- Troubleshooting guide for Fedora-specific issues
- Performance optimization recommendations for Fedora
- Security hardening specific to Fedora
- SELinux configuration guide for the project

### 10. Create Streamlined Workflow Scripts

Develop convenience scripts for Fedora development:
- Script to toggle SELinux enforcement for development
- Script to manage firewall rules for the application
- Script to quickly reset the development environment
- Script to check for Fedora updates that might affect the project
- Script to install only minimal dependencies for Fedora Minimal installations

## Technical Requirements

1. **Python Environment**:
   - Update to use Python 3.10 as specified in environment.yml
   - Ensure compatibility with Fedora's Python packages

2. **Database Setup**:
   - Configure for PostgreSQL 14+ (Fedora's default version)
   - Set up proper SELinux contexts for database directories

3. **System Dependencies**:
   - Use DNF package groups where appropriate
   - Install development headers for all required libraries
   - Include X11 libraries required by RDKit

4. **Security Considerations**:
   - Configure firewalld properly
   - Set up SELinux contexts correctly
   - Use systemd service files for production deployment
   - Implement Fedora-specific security best practices

5. **Performance Optimization**:
   - Optimize for Fedora's default filesystem
   - Configure proper resource limits
   - Set up systemd slice configurations
   - Optimize PostgreSQL for Fedora

## Implementation Guidelines

1. Prefer DNF module streams over manual package installation
2. Use Fedora's default locations for service files and configurations
3. Implement proper SELinux contexts throughout the project
4. Follow Fedora Packaging Guidelines where applicable
5. Use systemd for service management
6. Implement firewalld configurations rather than direct iptables
7. Prefer ostree/rpm-ostree compatible approaches when possible
8. Ensure compatibility with both Workstation and Server editions
9. Document differences between Fedora versions (36-39)
10. Include compatibility with Fedora Silverblue/Kinoite

## Success Criteria

1. All installation steps are fully automated
2. Setup can be completed with a single command
3. All verification tests pass
4. SELinux can remain in enforcing mode
5. Docker containers function correctly
6. Database connections work properly
7. RDKit visualizations display correctly
8. All project features function as expected
9. Documentation is comprehensive and Fedora-specific
10. Development workflow is smooth and intuitive

## Deliverables

1. `fedora_setup.sh` script
2. Fedora-specific `.env` template
3. Updated `.gitignore` for Fedora
4. SELinux policy module files
5. RDKit compatibility layer
6. Docker configuration updates
7. Verification test script
8. Fedora-specific documentation
9. Systemd service files
10. Workflow convenience scripts

This prompt provides comprehensive guidance to set up the CryoProtect project on Fedora OS, addressing all aspects of the migration and ensuring compatibility with Fedora's specific requirements and security features.