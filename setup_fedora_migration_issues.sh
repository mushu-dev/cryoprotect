#!/bin/bash
# Setup Fedora Migration Issues using RooRoo GitHub Issue Manager

set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

ROOROO_PATH=".github/rooroo/issue-manager.sh"

# Function to print section headers
print_section() {
    echo -e "\n${BLUE}======== $1 ========${NC}\n"
}

# Check if RooRoo exists
if [ ! -f "$ROOROO_PATH" ]; then
    echo -e "${RED}Error: RooRoo not found at $ROOROO_PATH${NC}"
    echo "Please run setup_rooroo.sh first."
    exit 1
fi

print_section "CryoProtect Fedora Migration - Issue Setup"

echo "This script will create GitHub issues for the CryoProtect Fedora migration project."
echo
echo -e "${YELLOW}Note: This will create actual GitHub issues in your repository.${NC}"
echo

read -p "Enter the repository owner/name (e.g., yourusername/CryoProtect): " REPO
if [ -z "$REPO" ]; then
    echo -e "${RED}Repository name is required. Exiting.${NC}"
    exit 1
fi

# Create Fedora-specific labels
print_section "Creating Fedora-specific labels"

$ROOROO_PATH create-labels \
  --repo "$REPO" \
  --labels "fedora-migration,selinux,database,docker,environment,testing,documentation,p0-blocker,p1-critical,p2-important,p3-enhancement"

echo -e "${GREEN}Labels created successfully!${NC}"

# Create Epic issues
print_section "Creating Epic Issues"

# Epic 1: Initial Setup and Analysis
EPIC1_NUMBER=$($ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Epic: Initial Setup and Analysis" \
  --body "# Epic: Initial Setup and Analysis

This epic covers the initial analysis, path normalization, and creation of essential scripts for cross-platform compatibility.

## Objectives
- Analyze the codebase for Windows-specific components
- Create Linux shell script equivalents for Windows batch files
- Normalize path separators across the codebase
- Create verification scripts for cross-platform compatibility

## Success Criteria
- All Windows-specific components identified
- Linux equivalents created for all essential batch files
- Path separators standardized throughout the codebase
- Verification scripts successfully test cross-platform functionality

## Related Epics
- Database Migration
- Fedora Core Setup" \
  --labels "fedora-migration,p0-blocker" \
  --output-number)

echo "Created Epic 1 with issue number: $EPIC1_NUMBER"

# Epic 2: Fedora Core Setup
EPIC2_NUMBER=$($ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Epic: Fedora Core Setup" \
  --body "# Epic: Fedora Core Setup

This epic covers the implementation of Fedora-specific configuration, including SELinux policies, firewall rules, and system dependencies.

## Objectives
- Create comprehensive Fedora setup script
- Configure SELinux contexts and policies
- Set up firewall rules for required services
- Install and configure system dependencies

## Success Criteria
- fedora_setup.sh script works end-to-end
- Application runs with SELinux in enforcing mode
- All necessary ports are properly configured in firewalld
- All system dependencies are properly installed and configured

## Related Epics
- Initial Setup and Analysis
- Database Migration
- Environment Configuration" \
  --labels "fedora-migration,selinux,p0-blocker" \
  --output-number)

echo "Created Epic 2 with issue number: $EPIC2_NUMBER"

# Epic 3: Database Migration
EPIC3_NUMBER=$($ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Epic: Database Migration" \
  --body "# Epic: Database Migration

This epic covers PostgreSQL setup on Fedora, including authentication, performance tuning, and data migration validation.

## Objectives
- Configure PostgreSQL installation on Fedora
- Set up proper authentication for local development
- Define SELinux contexts for database directories
- Optimize PostgreSQL configuration for Fedora

## Success Criteria
- PostgreSQL runs properly with Fedora security features
- Database authentication works for local and remote connections
- SELinux contexts allow proper database operation
- Performance is optimized for Fedora's filesystem and memory management

## Related Epics
- Fedora Core Setup
- Testing and Verification" \
  --labels "fedora-migration,database,selinux,p0-blocker" \
  --output-number)

echo "Created Epic 3 with issue number: $EPIC3_NUMBER"

# Epic 4: Environment Configuration
EPIC4_NUMBER=$($ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Epic: Environment Configuration" \
  --body "# Epic: Environment Configuration

This epic covers conda/Python environment setup, RDKit installation, and dependency management for Fedora.

## Objectives
- Create conda environment configuration for Fedora
- Resolve RDKit installation and X11 dependencies
- Update Python dependency management for Fedora
- Create environment variable templates for Fedora

## Success Criteria
- Conda environment is properly set up with all dependencies
- RDKit works properly with visualization functions
- All Python dependencies are correctly installed
- Environment variables properly configured for Fedora paths

## Related Epics
- Fedora Core Setup
- Testing and Verification" \
  --labels "fedora-migration,environment,p1-critical" \
  --output-number)

echo "Created Epic 4 with issue number: $EPIC4_NUMBER"

# Epic 5: Containerization
EPIC5_NUMBER=$($ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Epic: Containerization" \
  --body "# Epic: Containerization

This epic covers Docker configuration updates, volume mounting with SELinux contexts, and container security on Fedora.

## Objectives
- Update Dockerfile for Fedora compatibility
- Configure Docker volumes with proper SELinux contexts
- Update docker-compose configuration for Fedora
- Implement container security best practices for Fedora

## Success Criteria
- Docker containers build and run successfully on Fedora
- Volume mounts work properly with SELinux
- Docker Compose creates the entire environment correctly
- Container security follows Fedora best practices

## Related Epics
- Fedora Core Setup
- Testing and Verification" \
  --labels "fedora-migration,docker,selinux,p1-critical" \
  --output-number)

echo "Created Epic 5 with issue number: $EPIC5_NUMBER"

# Epic 6: Testing and Verification
EPIC6_NUMBER=$($ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Epic: Testing and Verification" \
  --body "# Epic: Testing and Verification

This epic covers the creation and execution of test scripts to validate the migration across all components.

## Objectives
- Create test scripts for SELinux configuration
- Create test scripts for database connectivity
- Create test scripts for RDKit functionality
- Create integration tests for the full application stack

## Success Criteria
- All tests pass on Fedora Linux
- SELinux remains in enforcing mode during tests
- Database operations work correctly
- RDKit visualization functions render properly
- Full application stack works end-to-end

## Related Epics
- All other epics" \
  --labels "fedora-migration,testing,p1-critical" \
  --output-number)

echo "Created Epic 6 with issue number: $EPIC6_NUMBER"

# Create some specific tasks under Epic 2: Fedora Core Setup
print_section "Creating Tasks under Epic 2: Fedora Core Setup"

# Task 1: Create fedora_setup.sh script
TASK1_NUMBER=$($ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Create fedora_setup.sh master script" \
  --body "# Task: Create fedora_setup.sh master script

## Description
Create a comprehensive setup script for Fedora Linux that handles all aspects of setting up the CryoProtect environment.

## Requirements
- Script should detect if running on Fedora
- Install all required system packages using DNF
- Set up PostgreSQL with proper configuration
- Configure SELinux contexts for application directories
- Set up firewall rules
- Create conda environment with all dependencies
- Verify the installation works correctly

## Success Criteria
- Script runs without errors on a fresh Fedora installation
- All components are properly installed and configured
- SELinux contexts are correctly set
- Firewall rules are properly configured
- Application starts successfully after running the script

## Related Issues
- Epic: Fedora Core Setup (#$EPIC2_NUMBER)" \
  --labels "fedora-migration,selinux,p0-blocker" \
  --output-number)

echo "Created Task 1 with issue number: $TASK1_NUMBER"

# Link Task 1 to Epic 2
$ROOROO_PATH add-sub-issue-to-issue \
  --repo "$REPO" \
  --parent-issue-number "$EPIC2_NUMBER" \
  --sub-issue-number "$TASK1_NUMBER"

echo "Linked Task 1 to Epic 2"

# Task 2: Configure SELinux for CryoProtect
TASK2_NUMBER=$($ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Configure SELinux contexts and policies for CryoProtect" \
  --body "# Task: Configure SELinux contexts and policies for CryoProtect

## Description
Create SELinux configurations that allow CryoProtect to run properly with SELinux in enforcing mode.

## Requirements
- Identify all directories that need special SELinux contexts
- Create a SELinux policy module for CryoProtect
- Set up proper contexts for application directories
- Configure contexts for database connections
- Document SELinux setup process

## Success Criteria
- Application runs with SELinux in enforcing mode
- No SELinux denials are generated during normal operation
- Database connections work properly with SELinux
- File operations work correctly with proper contexts

## Related Issues
- Epic: Fedora Core Setup (#$EPIC2_NUMBER)" \
  --labels "fedora-migration,selinux,p0-blocker" \
  --output-number)

echo "Created Task 2 with issue number: $TASK2_NUMBER"

# Link Task 2 to Epic 2
$ROOROO_PATH add-sub-issue-to-issue \
  --repo "$REPO" \
  --parent-issue-number "$EPIC2_NUMBER" \
  --sub-issue-number "$TASK2_NUMBER"

echo "Linked Task 2 to Epic 2"

# Create task under Epic 3: Database Migration
print_section "Creating Task under Epic 3: Database Migration"

# Task 3: Configure PostgreSQL for Fedora
TASK3_NUMBER=$($ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Configure PostgreSQL for Fedora with SELinux support" \
  --body "# Task: Configure PostgreSQL for Fedora with SELinux support

## Description
Set up PostgreSQL on Fedora with proper SELinux contexts and configuration.

## Requirements
- Install PostgreSQL packages for Fedora
- Configure PostgreSQL for local development
- Set up proper authentication methods
- Configure SELinux contexts for database directories
- Ensure database connections work with SELinux enforcing

## Success Criteria
- PostgreSQL starts successfully on Fedora
- Database connections work with SELinux in enforcing mode
- Authentication works for both local and remote connections
- SELinux contexts are properly set for all database files

## Related Issues
- Epic: Database Migration (#$EPIC3_NUMBER)" \
  --labels "fedora-migration,database,selinux,p0-blocker" \
  --output-number)

echo "Created Task 3 with issue number: $TASK3_NUMBER"

# Link Task 3 to Epic 3
$ROOROO_PATH add-sub-issue-to-issue \
  --repo "$REPO" \
  --parent-issue-number "$EPIC3_NUMBER" \
  --sub-issue-number "$TASK3_NUMBER"

echo "Linked Task 3 to Epic 3"

print_section "GitHub Issues Setup Complete"

echo -e "${GREEN}Successfully created the basic issue structure for the Fedora migration project!${NC}"
echo
echo "The following issues were created:"
echo "  Epic 1: Initial Setup and Analysis (#$EPIC1_NUMBER)"
echo "  Epic 2: Fedora Core Setup (#$EPIC2_NUMBER)"
echo "  Epic 3: Database Migration (#$EPIC3_NUMBER)"
echo "  Epic 4: Environment Configuration (#$EPIC4_NUMBER)"
echo "  Epic 5: Containerization (#$EPIC5_NUMBER)"
echo "  Epic 6: Testing and Verification (#$EPIC6_NUMBER)"
echo
echo "Tasks:"
echo "  Task 1: Create fedora_setup.sh master script (#$TASK1_NUMBER)"
echo "  Task 2: Configure SELinux contexts and policies for CryoProtect (#$TASK2_NUMBER)"
echo "  Task 3: Configure PostgreSQL for Fedora with SELinux support (#$TASK3_NUMBER)"
echo
echo "You can now begin working on these issues or create additional tasks as needed."
echo "Use the RooRoo issue manager to track progress and create additional tasks."