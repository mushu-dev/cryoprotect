#!/bin/bash

# Creates only the epic issues without closing any existing issues
# This will create 9 well-structured epics to organize the repository

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Creating Epic Issues for Repository Organization${NC}"
echo -e "${YELLOW}This will create 9 epic issues to organize the repository.${NC}"

# Track created epic numbers to save to a file
EPIC_NUMBERS=""

# Database epic
echo -e "${BLUE}Creating Database Implementation epic...${NC}"
DB_EPIC=$(gh issue create \
  --title "Epic: Database Implementation and Optimization" \
  --body "# Database Implementation and Optimization

Database schema, migrations, optimization, and data population tasks.

This epic consolidates all database-related issues, including schema design, migrations, optimization, and data population.

## Key Tasks
- Implement proper database schema
- Create migrations for schema changes
- Optimize database performance
- Set up data population scripts
- Implement RLS policies
- Handle connection pooling

## Consolidated From
Multiple database-related issues have been consolidated here for better tracking and organization." \
  --label "type:epic,status:planning,area:database,priority:high")
echo -e "${GREEN}Created: $DB_EPIC${NC}"
DB_EPIC_NUM=$(echo $DB_EPIC | grep -o '#[0-9]*' | tr -d '#')
EPIC_NUMBERS="$EPIC_NUMBERS $DB_EPIC_NUM"
sleep 2

# ChEMBL Integration epic
echo -e "${BLUE}Creating ChEMBL Integration epic...${NC}"
CHEMBL_EPIC=$(gh issue create \
  --title "Epic: ChEMBL Integration and Data Pipeline" \
  --body "# ChEMBL Integration and Data Pipeline

ChEMBL data import, processing, and integration tasks.

This epic consolidates all issues related to ChEMBL data integration, including data import, processing, and pipeline development.

## Key Tasks
- Implement ChEMBL data import scripts
- Process and validate imported data
- Create data pipelines for automation
- Set up verification and validation
- Integrate with existing data
- Set up ChEMBL data update mechanisms

## Consolidated From
Multiple ChEMBL-related issues have been consolidated here for better tracking and organization." \
  --label "type:epic,status:planning,area:database,priority:high")
echo -e "${GREEN}Created: $CHEMBL_EPIC${NC}"
CHEMBL_EPIC_NUM=$(echo $CHEMBL_EPIC | grep -o '#[0-9]*' | tr -d '#')
EPIC_NUMBERS="$EPIC_NUMBERS $CHEMBL_EPIC_NUM"
sleep 2

# PubChem Integration epic
echo -e "${BLUE}Creating PubChem Integration epic...${NC}"
PUBCHEM_EPIC=$(gh issue create \
  --title "Epic: PubChem Integration and Molecule Management" \
  --body "# PubChem Integration and Molecule Management

PubChem data import, molecule handling, and property management.

This epic consolidates all issues related to PubChem data integration, molecule management, and property handling.

## Key Tasks
- Implement PubChem data import mechanisms
- Process and validate molecule data
- Manage molecular properties
- Handle CID-based identifiers
- Set up molecule search and retrieval
- Implement molecule deduplication and consolidation

## Consolidated From
Multiple PubChem and molecule-related issues have been consolidated here for better tracking and organization." \
  --label "type:epic,status:planning,area:database,priority:high")
echo -e "${GREEN}Created: $PUBCHEM_EPIC${NC}"
PUBCHEM_EPIC_NUM=$(echo $PUBCHEM_EPIC | grep -o '#[0-9]*' | tr -d '#')
EPIC_NUMBERS="$EPIC_NUMBERS $PUBCHEM_EPIC_NUM"
sleep 2

# API Development epic
echo -e "${BLUE}Creating API Development epic...${NC}"
API_EPIC=$(gh issue create \
  --title "Epic: API Development and Standardization" \
  --body "# API Development and Standardization

API endpoints, standardization, documentation, and testing.

This epic consolidates all issues related to API development, including endpoint creation, standardization, documentation, and testing.

## Key Tasks
- Design RESTful API endpoints
- Standardize API responses and error handling
- Create API documentation
- Implement authentication and authorization
- Set up API testing
- Handle API versioning

## Consolidated From
Multiple API-related issues have been consolidated here for better tracking and organization." \
  --label "type:epic,status:planning,area:api,priority:high")
echo -e "${GREEN}Created: $API_EPIC${NC}"
API_EPIC_NUM=$(echo $API_EPIC | grep -o '#[0-9]*' | tr -d '#')
EPIC_NUMBERS="$EPIC_NUMBERS $API_EPIC_NUM"
sleep 2

# Frontend Implementation epic
echo -e "${BLUE}Creating Frontend Implementation epic...${NC}"
UI_EPIC=$(gh issue create \
  --title "Epic: Frontend Implementation and User Experience" \
  --body "# Frontend Implementation and User Experience

Frontend features, UI components, and user experience improvements.

This epic consolidates all issues related to frontend development, including UI components, user experience, and frontend features.

## Key Tasks
- Develop user interface components
- Implement frontend features
- Enhance user experience
- Set up frontend testing
- Handle frontend state management
- Implement responsive design

## Consolidated From
Multiple frontend-related issues have been consolidated here for better tracking and organization." \
  --label "type:epic,status:planning,area:ui,priority:high")
echo -e "${GREEN}Created: $UI_EPIC${NC}"
UI_EPIC_NUM=$(echo $UI_EPIC | grep -o '#[0-9]*' | tr -d '#')
EPIC_NUMBERS="$EPIC_NUMBERS $UI_EPIC_NUM"
sleep 2

# Authentication epic
echo -e "${BLUE}Creating Authentication epic...${NC}"
AUTH_EPIC=$(gh issue create \
  --title "Epic: Authentication and Security Implementation" \
  --body "# Authentication and Security Implementation

User authentication, authorization, RLS policies, and security.

This epic consolidates all issues related to authentication, authorization, and security, including RLS policies and user management.

## Key Tasks
- Implement user authentication
- Set up authorization mechanisms
- Develop RLS policies
- Handle JWT tokens
- Implement service roles
- Set up security testing

## Consolidated From
Multiple authentication and security-related issues have been consolidated here for better tracking and organization." \
  --label "type:epic,status:planning,area:auth,priority:high")
echo -e "${GREEN}Created: $AUTH_EPIC${NC}"
AUTH_EPIC_NUM=$(echo $AUTH_EPIC | grep -o '#[0-9]*' | tr -d '#')
EPIC_NUMBERS="$EPIC_NUMBERS $AUTH_EPIC_NUM"
sleep 2

# RDKit Integration epic
echo -e "${BLUE}Creating RDKit Integration epic...${NC}"
RDKIT_EPIC=$(gh issue create \
  --title "Epic: RDKit Integration and Chemical Functionality" \
  --body "# RDKit Integration and Chemical Functionality

RDKit features, chemical property calculations, and molecular analysis.

This epic consolidates all issues related to RDKit integration, chemical property calculations, and molecular analysis.

## Key Tasks
- Integrate RDKit into the application
- Implement chemical property calculations
- Set up molecular analysis functions
- Handle chemical data visualization
- Set up RDKit in containers
- Develop RDKit-based algorithms

## Consolidated From
Multiple RDKit-related issues have been consolidated here for better tracking and organization." \
  --label "type:epic,status:planning,area:chembl,priority:high")
echo -e "${GREEN}Created: $RDKIT_EPIC${NC}"
RDKIT_EPIC_NUM=$(echo $RDKIT_EPIC | grep -o '#[0-9]*' | tr -d '#')
EPIC_NUMBERS="$EPIC_NUMBERS $RDKIT_EPIC_NUM"
sleep 2

# Infrastructure epic
echo -e "${BLUE}Creating Infrastructure epic...${NC}"
INFRA_EPIC=$(gh issue create \
  --title "Epic: Infrastructure and Deployment" \
  --body "# Infrastructure and Deployment

Containerization, deployment, environment setup, and DevOps.

This epic consolidates all issues related to infrastructure, deployment, and environment setup, including containerization and DevOps.

## Key Tasks
- Set up containerization
- Implement deployment pipelines
- Configure environments
- Set up monitoring and logging
- Handle infrastructure as code
- Implement DevOps best practices

## Consolidated From
Multiple infrastructure-related issues have been consolidated here for better tracking and organization." \
  --label "type:epic,status:planning,area:devops,priority:high")
echo -e "${GREEN}Created: $INFRA_EPIC${NC}"
INFRA_EPIC_NUM=$(echo $INFRA_EPIC | grep -o '#[0-9]*' | tr -d '#')
EPIC_NUMBERS="$EPIC_NUMBERS $INFRA_EPIC_NUM"
sleep 2

# Testing epic
echo -e "${BLUE}Creating Testing epic...${NC}"
TESTING_EPIC=$(gh issue create \
  --title "Epic: Testing and Quality Assurance" \
  --body "# Testing and Quality Assurance

Test frameworks, validation procedures, and quality assurance processes.

This epic consolidates all issues related to testing, validation, and quality assurance, including test frameworks and procedures.

## Key Tasks
- Set up test frameworks
- Implement unit tests
- Develop integration tests
- Create validation procedures
- Handle test automation
- Set up CI/CD testing

## Consolidated From
Multiple testing-related issues have been consolidated here for better tracking and organization." \
  --label "type:epic,status:planning,area:testing,priority:high")
echo -e "${GREEN}Created: $TESTING_EPIC${NC}"
TESTING_EPIC_NUM=$(echo $TESTING_EPIC | grep -o '#[0-9]*' | tr -d '#')
EPIC_NUMBERS="$EPIC_NUMBERS $TESTING_EPIC_NUM"

# Save epic numbers to a file for use by other scripts
echo -e "${YELLOW}Saving epic issue numbers to /tmp/epic_issues.txt${NC}"
echo $EPIC_NUMBERS > /tmp/epic_issues.txt

echo -e "${GREEN}Successfully created 9 epic issues!${NC}"
echo -e "${BLUE}Epic issue numbers:${NC} $EPIC_NUMBERS"
echo ""
echo -e "${BLUE}Next step:${NC} Run the close-remaining-issues.sh script to close all issues except these epics and key issues."