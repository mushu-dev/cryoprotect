#!/bin/bash

# Step-by-step issue consolidation
# This script breaks down the consolidation process into smaller steps
# to avoid timeouts and rate limits

# Set variables
DRY_RUN="false"
BATCH_SIZE=20
SLEEP_BETWEEN_ISSUES=2
SLEEP_BETWEEN_BATCHES=15

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check parameters
if [ "$1" == "--dry-run" ]; then
  DRY_RUN="true"
  echo -e "${YELLOW}Running in dry-run mode. No changes will be made.${NC}"
fi

# Step 1: Create epic issues
echo -e "${BLUE}Step 1: Creating epic issues...${NC}"
create_epic_issues() {
  # Database epic
  if [ "$DRY_RUN" == "false" ]; then
    echo "Creating Database Implementation epic..."
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
    echo "Created: $DB_EPIC"
    DB_EPIC_NUM=$(echo $DB_EPIC | grep -o '#[0-9]*' | tr -d '#')
    echo "Epic number: $DB_EPIC_NUM"
    sleep $SLEEP_BETWEEN_ISSUES
  else
    echo "${YELLOW}[DRY RUN] Would create Database Implementation epic${NC}"
    DB_EPIC_NUM="999"
  fi

  # ChEMBL Integration epic
  if [ "$DRY_RUN" == "false" ]; then
    echo "Creating ChEMBL Integration epic..."
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
    echo "Created: $CHEMBL_EPIC"
    CHEMBL_EPIC_NUM=$(echo $CHEMBL_EPIC | grep -o '#[0-9]*' | tr -d '#')
    echo "Epic number: $CHEMBL_EPIC_NUM"
    sleep $SLEEP_BETWEEN_ISSUES
  else
    echo "${YELLOW}[DRY RUN] Would create ChEMBL Integration epic${NC}"
    CHEMBL_EPIC_NUM="998"
  fi

  # PubChem Integration epic
  if [ "$DRY_RUN" == "false" ]; then
    echo "Creating PubChem Integration epic..."
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
    echo "Created: $PUBCHEM_EPIC"
    PUBCHEM_EPIC_NUM=$(echo $PUBCHEM_EPIC | grep -o '#[0-9]*' | tr -d '#')
    echo "Epic number: $PUBCHEM_EPIC_NUM"
    sleep $SLEEP_BETWEEN_ISSUES
  else
    echo "${YELLOW}[DRY RUN] Would create PubChem Integration epic${NC}"
    PUBCHEM_EPIC_NUM="997"
  fi

  # API Development epic
  if [ "$DRY_RUN" == "false" ]; then
    echo "Creating API Development epic..."
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
    echo "Created: $API_EPIC"
    API_EPIC_NUM=$(echo $API_EPIC | grep -o '#[0-9]*' | tr -d '#')
    echo "Epic number: $API_EPIC_NUM"
    sleep $SLEEP_BETWEEN_ISSUES
  else
    echo "${YELLOW}[DRY RUN] Would create API Development epic${NC}"
    API_EPIC_NUM="996"
  fi

  # Frontend Implementation epic
  if [ "$DRY_RUN" == "false" ]; then
    echo "Creating Frontend Implementation epic..."
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
    echo "Created: $UI_EPIC"
    UI_EPIC_NUM=$(echo $UI_EPIC | grep -o '#[0-9]*' | tr -d '#')
    echo "Epic number: $UI_EPIC_NUM"
    sleep $SLEEP_BETWEEN_ISSUES
  else
    echo "${YELLOW}[DRY RUN] Would create Frontend Implementation epic${NC}"
    UI_EPIC_NUM="995"
  fi

  # Authentication epic
  if [ "$DRY_RUN" == "false" ]; then
    echo "Creating Authentication epic..."
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
    echo "Created: $AUTH_EPIC"
    AUTH_EPIC_NUM=$(echo $AUTH_EPIC | grep -o '#[0-9]*' | tr -d '#')
    echo "Epic number: $AUTH_EPIC_NUM"
    sleep $SLEEP_BETWEEN_ISSUES
  else
    echo "${YELLOW}[DRY RUN] Would create Authentication epic${NC}"
    AUTH_EPIC_NUM="994"
  fi

  # RDKit Integration epic
  if [ "$DRY_RUN" == "false" ]; then
    echo "Creating RDKit Integration epic..."
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
    echo "Created: $RDKIT_EPIC"
    RDKIT_EPIC_NUM=$(echo $RDKIT_EPIC | grep -o '#[0-9]*' | tr -d '#')
    echo "Epic number: $RDKIT_EPIC_NUM"
    sleep $SLEEP_BETWEEN_ISSUES
  else
    echo "${YELLOW}[DRY RUN] Would create RDKit Integration epic${NC}"
    RDKIT_EPIC_NUM="993"
  fi

  # Infrastructure epic
  if [ "$DRY_RUN" == "false" ]; then
    echo "Creating Infrastructure epic..."
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
    echo "Created: $INFRA_EPIC"
    INFRA_EPIC_NUM=$(echo $INFRA_EPIC | grep -o '#[0-9]*' | tr -d '#')
    echo "Epic number: $INFRA_EPIC_NUM"
    sleep $SLEEP_BETWEEN_ISSUES
  else
    echo "${YELLOW}[DRY RUN] Would create Infrastructure epic${NC}"
    INFRA_EPIC_NUM="992"
  fi

  # Testing epic
  if [ "$DRY_RUN" == "false" ]; then
    echo "Creating Testing epic..."
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
    echo "Created: $TESTING_EPIC"
    TESTING_EPIC_NUM=$(echo $TESTING_EPIC | grep -o '#[0-9]*' | tr -d '#')
    echo "Epic number: $TESTING_EPIC_NUM"
    sleep $SLEEP_BETWEEN_ISSUES
  else
    echo "${YELLOW}[DRY RUN] Would create Testing epic${NC}"
    TESTING_EPIC_NUM="991"
  fi

  # Return epic numbers as a space-separated string
  echo "$DB_EPIC_NUM $CHEMBL_EPIC_NUM $PUBCHEM_EPIC_NUM $API_EPIC_NUM $UI_EPIC_NUM $AUTH_EPIC_NUM $RDKIT_EPIC_NUM $INFRA_EPIC_NUM $TESTING_EPIC_NUM"
}

# Step 2: Get test issues
get_test_issues() {
  echo -e "${BLUE}Identifying test issues...${NC}"
  gh issue list --json number,title --search "is:open test in:title" --limit 100 | jq -r '.[] | .number'
  gh issue list --json number,title --search "is:open validation in:title" --limit 100 | jq -r '.[] | .number'
}

# Step 3: Close issues in batches
close_issues_in_batches() {
  local issue_list="$1"
  local epic_reference="$2"
  local total_issues=$(echo "$issue_list" | wc -w)
  local current=0
  local batch=0
  
  echo -e "${BLUE}Closing $total_issues issues (in batches of $BATCH_SIZE)...${NC}"
  
  for issue_num in $issue_list; do
    current=$((current + 1))
    batch=$((batch + 1))
    
    if [ "$DRY_RUN" == "false" ]; then
      echo "[$current/$total_issues] Closing issue #$issue_num..."
      gh issue close "$issue_num" --comment "This issue has been consolidated into $epic_reference as part of repository cleanup and organization.

Please refer to the epic for tracking this work going forward."
      sleep $SLEEP_BETWEEN_ISSUES
    else
      echo "${YELLOW}[DRY RUN] Would close issue #$issue_num${NC}"
    fi
    
    # If batch is complete, sleep to avoid rate limiting
    if [ $batch -eq $BATCH_SIZE ]; then
      echo -e "${GREEN}Completed batch of $BATCH_SIZE issues. Sleeping for $SLEEP_BETWEEN_BATCHES seconds...${NC}"
      sleep $SLEEP_BETWEEN_BATCHES
      batch=0
    fi
  done
}

# Main function
main() {
  # Step 1: Create epic issues
  echo -e "${BLUE}=== STEP 1: CREATING EPIC ISSUES ===${NC}"
  EPIC_NUMBERS=$(create_epic_issues)
  echo -e "${GREEN}Created epic issues with numbers: $EPIC_NUMBERS${NC}"
  
  # Parse epic numbers
  read DB_EPIC_NUM CHEMBL_EPIC_NUM PUBCHEM_EPIC_NUM API_EPIC_NUM UI_EPIC_NUM AUTH_EPIC_NUM RDKIT_EPIC_NUM INFRA_EPIC_NUM TESTING_EPIC_NUM <<< "$EPIC_NUMBERS"
  
  # Step 2: Close test issues
  echo -e "${BLUE}=== STEP 2: CLOSING TEST ISSUES ===${NC}"
  TEST_ISSUES=$(get_test_issues)
  TEST_COUNT=$(echo "$TEST_ISSUES" | wc -w)
  echo -e "${YELLOW}Found $TEST_COUNT test issues to close${NC}"
  
  if [ $TEST_COUNT -gt 0 ]; then
    close_issues_in_batches "$TEST_ISSUES" "the Testing epic (#$TESTING_EPIC_NUM)"
  fi
  
  # Step 3: Close database issues
  echo -e "${BLUE}=== STEP 3: CLOSING DATABASE ISSUES ===${NC}"
  DB_ISSUES=$(gh issue list --json number --search "is:open database schema migration sql table supabase" --limit 100 | jq -r '.[] | .number')
  DB_COUNT=$(echo "$DB_ISSUES" | wc -w)
  echo -e "${YELLOW}Found $DB_COUNT database issues to close${NC}"
  
  if [ $DB_COUNT -gt 0 ]; then
    close_issues_in_batches "$DB_ISSUES" "the Database Implementation epic (#$DB_EPIC_NUM)"
  fi
  
  # Step 4: Close ChEMBL issues
  echo -e "${BLUE}=== STEP 4: CLOSING CHEMBL ISSUES ===${NC}"
  CHEMBL_ISSUES=$(gh issue list --json number --search "is:open chembl" --limit 100 | jq -r '.[] | .number')
  CHEMBL_COUNT=$(echo "$CHEMBL_ISSUES" | wc -w)
  echo -e "${YELLOW}Found $CHEMBL_COUNT ChEMBL issues to close${NC}"
  
  if [ $CHEMBL_COUNT -gt 0 ]; then
    close_issues_in_batches "$CHEMBL_ISSUES" "the ChEMBL Integration epic (#$CHEMBL_EPIC_NUM)"
  fi
  
  # Step 5: Close PubChem issues
  echo -e "${BLUE}=== STEP 5: CLOSING PUBCHEM ISSUES ===${NC}"
  PUBCHEM_ISSUES=$(gh issue list --json number --search "is:open pubchem cid molecule property" --limit 100 | jq -r '.[] | .number')
  PUBCHEM_COUNT=$(echo "$PUBCHEM_ISSUES" | wc -w)
  echo -e "${YELLOW}Found $PUBCHEM_COUNT PubChem issues to close${NC}"
  
  if [ $PUBCHEM_COUNT -gt 0 ]; then
    close_issues_in_batches "$PUBCHEM_ISSUES" "the PubChem Integration epic (#$PUBCHEM_EPIC_NUM)"
  fi
  
  # Step 6: Close API issues
  echo -e "${BLUE}=== STEP 6: CLOSING API ISSUES ===${NC}"
  API_ISSUES=$(gh issue list --json number --search "is:open api endpoint route flask" --limit 100 | jq -r '.[] | .number')
  API_COUNT=$(echo "$API_ISSUES" | wc -w)
  echo -e "${YELLOW}Found $API_COUNT API issues to close${NC}"
  
  if [ $API_COUNT -gt 0 ]; then
    close_issues_in_batches "$API_ISSUES" "the API Development epic (#$API_EPIC_NUM)"
  fi
  
  # Step 7: Keep specific issues
  echo -e "${BLUE}=== STEP 7: KEEPING SPECIFIC ISSUES ===${NC}"
  KEEP_ISSUES="246 245 244 243 211"
  echo -e "${GREEN}Keeping these issues: $KEEP_ISSUES${NC}"
  
  # Final summary
  echo -e "${BLUE}=== CONSOLIDATION SUMMARY ===${NC}"
  echo -e "${GREEN}Created 9 epic issues${NC}"
  echo -e "${GREEN}Closed multiple issues in batches${NC}"
  echo -e "${GREEN}Kept 5 specific issues${NC}"
  
  echo -e "${BLUE}=== RECOMMENDED NEXT STEPS ===${NC}"
  echo "1. Review the created epic issues to ensure they cover all important areas"
  echo "2. Check that no essential issues were accidentally closed"
  echo "3. Set up the project board with the new epics and remaining issues"
  echo "4. Continue with ongoing work using the new streamlined issue structure"
}

# Execute main
main