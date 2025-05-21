#!/bin/bash

# GitHub Issue Management Setup Script
# This script sets up the complete GitHub issue management system for CryoProtect
# It creates all required labels, updates issue templates, and configures project settings

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if the GitHub CLI is installed and authenticated
check_gh() {
  if ! command -v gh &> /dev/null; then
    echo -e "${RED}Error: GitHub CLI (gh) is not installed.${NC}"
    echo "Please install it from https://cli.github.com/"
    exit 1
  fi
  
  # Check if authenticated
  if ! gh auth status &> /dev/null; then
    echo -e "${RED}Error: Not authenticated with GitHub CLI.${NC}"
    echo "Please run 'gh auth login' first."
    exit 1
  fi
}

# Create all the required labels
create_labels() {
  echo -e "${BLUE}Creating standard labels...${NC}"
  
  # Type labels
  gh label create "type:epic" --color "6E49CB" --description "Large initiative containing multiple issues"
  gh label create "type:feature" --color "b60205" --description "New feature implementation"
  gh label create "type:bugfix" --color "88E3AE" --description "Bug fix"
  gh label create "type:refactor" --color "F7C4A7" --description "Code refactoring"
  gh label create "type:documentation" --color "B01778" --description "Documentation updates"
  gh label create "type:validation" --color "6BA5AF" --description "Validation and testing tasks"
  
  # Status labels
  gh label create "status:planning" --color "FEF2C0" --description "In planning/design phase"
  gh label create "status:ready" --color "0E8A16" --description "Ready for implementation"
  gh label create "status:in-progress" --color "BFDADC" --description "Currently being worked on"
  gh label create "status:needs-review" --color "008F91" --description "Implementation complete, needs review"
  gh label create "status:completed" --color "806362" --description "Fully completed and validated"
  gh label create "status:blocked" --color "D93F0B" --description "Blocked by dependencies or issues"
  
  # Area labels
  gh label create "area:database" --color "0052CC" --description "Database-related issues"
  gh label create "area:api" --color "52CC00" --description "API-related issues"
  gh label create "area:ui" --color "CC5200" --description "User interface issues"
  gh label create "area:auth" --color "5200CC" --description "Authentication and security issues"
  gh label create "area:testing" --color "00CC52" --description "Testing and QA issues"
  gh label create "area:chembl" --color "CC0052" --description "ChEMBL integration issues"
  gh label create "area:docs" --color "00CCCC" --description "Documentation issues"
  gh label create "area:devops" --color "CCCC00" --description "Infrastructure and DevOps issues"
  
  # Priority labels
  gh label create "priority:high" --color "B60205" --description "High priority"
  gh label create "priority:medium" --color "F2CB05" --description "Medium priority" 
  gh label create "priority:low" --color "4CBB17" --description "Low priority"
  
  # Special labels
  gh label create "stale" --color "D3D3D3" --description "No activity for 30+ days"
  gh label create "duplicate" --color "EEEEEE" --description "Duplicate of another issue"
  gh label create "ai-generated" --color "FF00FF" --description "AI-generated test issue"
  
  echo -e "${GREEN}Labels created successfully!${NC}"
}

# Create the GitHub project board
create_project_board() {
  echo -e "${BLUE}Setting up GitHub project board...${NC}"
  
  # Create project board (if it doesn't exist)
  # Note: This is more complex with the GitHub CLI and may require GraphQL API
  echo -e "${YELLOW}Project board setup requires manual configuration in the GitHub UI${NC}"
  echo "Please create a project board with the following columns:"
  echo "  - Todo"
  echo "  - Ready"
  echo "  - In Progress"
  echo "  - Review"
  echo "  - Done"
  echo "  - Blocked"
}

# Create Epic issues
create_epics() {
  echo -e "${BLUE}Creating epic issues...${NC}"
  
  # Define the 9 epic categories
  declare -A epics
  epics["database"]="Database Implementation and Optimization"
  epics["chembl"]="ChEMBL Integration and Data Pipeline"
  epics["pubchem"]="PubChem Integration and Molecule Management"
  epics["api"]="API Development and Standardization"
  epics["frontend"]="Frontend Implementation and User Experience"
  epics["auth"]="Authentication and Security Implementation"
  epics["rdkit"]="RDKit Integration and Chemical Functionality"
  epics["infrastructure"]="Infrastructure and Deployment"
  epics["testing"]="Testing and Quality Assurance"
  
  # Create each epic
  for key in "${!epics[@]}"; do
    echo "Creating epic: ${epics[$key]}"
    
    # Determine appropriate area label
    area_label=""
    case $key in
      database|chembl|pubchem)
        area_label="area:database"
        ;;
      api)
        area_label="area:api"
        ;;
      frontend)
        area_label="area:ui"
        ;;
      auth)
        area_label="area:auth"
        ;;
      rdkit)
        area_label="area:chembl"
        ;;
      infrastructure)
        area_label="area:devops"
        ;;
      testing)
        area_label="area:testing"
        ;;
    esac
    
    # Create the epic issue
    gh issue create \
      --title "Epic: ${epics[$key]}" \
      --label "type:epic" \
      --label "status:planning" \
      --label "priority:high" \
      --label "$area_label" \
      --body "# ${epics[$key]}

## Description
This epic tracks all work related to ${epics[$key]} for the CryoProtect project.

## Goals
- Add specific goals for this epic

## User Stories/Requirements
- Add key requirements here

## Technical Scope
- Define the technical scope here

## Tasks
- [ ] Task 1
- [ ] Task 2

## Definition of Done
- [ ] All subtasks are completed
- [ ] Documentation is updated
- [ ] Tests are passing
- [ ] Code review is completed
- [ ] Validation is completed"
    
    # Avoid rate limiting
    sleep 2
  done
  
  echo -e "${GREEN}Epic issues created successfully!${NC}"
}

# Update workflow files
update_workflows() {
  echo -e "${BLUE}Updating GitHub Actions workflow files...${NC}"
  
  # Update project automation workflow file if needed
  # This is a placeholder - actual implementation would depend on current workflows
  echo -e "${YELLOW}Workflow updates need to be customized to your current CI/CD setup${NC}"
  echo "Please review and adjust the workflow files in .github/workflows/ directory"
}

# Main function
main() {
  echo -e "${BLUE}================================${NC}"
  echo -e "${BLUE}GitHub Issue Management Setup${NC}"
  echo -e "${BLUE}================================${NC}"
  
  # Check GitHub CLI installation
  check_gh
  
  # Show options menu
  echo -e "\n${YELLOW}Select an option:${NC}"
  echo "1) Create all labels"
  echo "2) Set up project board (guidance)"
  echo "3) Create epic issues"
  echo "4) Update workflow files (guidance)"
  echo "5) Run complete setup (all of the above)"
  echo "q) Quit"
  
  read -p "Enter your choice: " choice
  
  case $choice in
    1)
      create_labels
      ;;
    2)
      create_project_board
      ;;
    3)
      create_epics
      ;;
    4)
      update_workflows
      ;;
    5)
      create_labels
      create_project_board
      create_epics
      update_workflows
      ;;
    q|Q)
      echo "Exiting..."
      exit 0
      ;;
    *)
      echo -e "${RED}Invalid option${NC}"
      exit 1
      ;;
  esac
  
  echo -e "\n${GREEN}Setup completed!${NC}"
}

# Run the main function
main