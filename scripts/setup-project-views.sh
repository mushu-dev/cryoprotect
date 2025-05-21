#!/bin/bash
# Script to set up project views for the GitHub project board

# Set colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Set up variables
LOG_FILE="project_views_setup_$(date +%Y%m%d_%H%M%S).log"

echo -e "${BLUE}CryoProtect Project Views Setup Tool${NC}"
echo "========================================"
echo "This script will set up standard views for your project board"
echo "Log will be saved to $LOG_FILE"

# Log start time
echo "Setup started at $(date)" > "$LOG_FILE"
echo "-----------------------------------------------" >> "$LOG_FILE"

# Ask user for the project number
read -p "Enter the GitHub project number (you can find this in the project URL): " PROJECT_NUMBER

if [ -z "$PROJECT_NUMBER" ]; then
  echo -e "${RED}No project number provided. Exiting.${NC}"
  echo "No project number provided. Exiting." >> "$LOG_FILE"
  exit 1
fi

# Verify project exists
echo -e "${BLUE}Checking project #$PROJECT_NUMBER...${NC}"
PROJECT_INFO=$(gh project view "$PROJECT_NUMBER" --owner mushu-dev --json title,url 2>/dev/null)

if [ -z "$PROJECT_INFO" ]; then
  echo -e "${RED}Project #$PROJECT_NUMBER not found or access denied. Exiting.${NC}"
  echo "Project #$PROJECT_NUMBER not found or access denied. Exiting." >> "$LOG_FILE"
  exit 1
fi

PROJECT_TITLE=$(echo "$PROJECT_INFO" | jq -r '.title')
PROJECT_URL=$(echo "$PROJECT_INFO" | jq -r '.url')

echo -e "${GREEN}Found project:${NC} $PROJECT_TITLE"
echo -e "${GREEN}URL:${NC} $PROJECT_URL"
echo "Working with project: $PROJECT_TITLE ($PROJECT_URL)" >> "$LOG_FILE"

# Create Kanban board view
echo -e "${BLUE}Creating Kanban board view...${NC}"
echo "Creating Kanban board view with status field" >> "$LOG_FILE"

VIEW_RESPONSE=$(gh project field-list "$PROJECT_NUMBER" --owner mushu-dev --format json)
STATUS_FIELD_ID=$(echo "$VIEW_RESPONSE" | jq -r '.[] | select(.name=="Status") | .id')

if [ -z "$STATUS_FIELD_ID" ] || [ "$STATUS_FIELD_ID" == "null" ]; then
  echo -e "${YELLOW}Warning: Status field not found. Creating the field...${NC}"
  echo "Warning: Status field not found. Creating the field..." >> "$LOG_FILE"
  
  # Create status field
  FIELD_RESPONSE=$(gh project field create "$PROJECT_NUMBER" --data-type single-select --name "Status" --owner mushu-dev --format json)
  STATUS_FIELD_ID=$(echo "$FIELD_RESPONSE" | jq -r '.id')
  
  # Add options to status field
  gh project field-option add "$PROJECT_NUMBER" Status --owner mushu-dev --option "Todo" --color "#E99695"
  gh project field-option add "$PROJECT_NUMBER" Status --owner mushu-dev --option "In Progress" --color "#F4D03F"
  gh project field-option add "$PROJECT_NUMBER" Status --owner mushu-dev --option "In Review" --color "#6AB04C"
  gh project field-option add "$PROJECT_NUMBER" Status --owner mushu-dev --option "Done" --color "#2ECC71"
  gh project field-option add "$PROJECT_NUMBER" Status --owner mushu-dev --option "Blocked" --color "#E74C3C"
fi

echo -e "${YELLOW}Note: GitHub API does not support creating views through the CLI.${NC}"
echo -e "${YELLOW}Please follow these manual steps to set up your views:${NC}"
echo ""
echo -e "${GREEN}1. Create Kanban Board View:${NC}"
echo "   - Go to: $PROJECT_URL"
echo "   - Click on the '+' button next to the 'Table' tab"
echo "   - Select 'Board'"
echo "   - Choose 'Status' as the field to group by"
echo "   - Click 'Create'"
echo ""
echo -e "${GREEN}2. Create Priority View:${NC}"
echo "   - Click on the '+' button next to the tabs"
echo "   - Select 'Board'"
echo "   - Choose 'Priority' as the field to group by"
echo "   - Click 'Create'"
echo ""
echo -e "${GREEN}3. Create Area View:${NC}"
echo "   - Click on the '+' button next to the tabs"
echo "   - Select 'Board'"
echo "   - Choose 'Area' as the field to group by"
echo "   - Click 'Create'"
echo ""
echo -e "${GREEN}4. Set up Automation:${NC}"
echo "   - Click on the project menu (three dots in the upper right)"
echo "   - Select 'Workflows'"
echo "   - Set up automation rules to update Status based on Issue status changes"
echo ""

# Provide instructions for bulk editing items
echo -e "${BLUE}Instructions for bulk editing project items:${NC}"
echo "1. In the Table view, select multiple items using checkboxes"
echo "2. Click on the field you want to edit (e.g., Status)"
echo "3. Set the value to apply to all selected items"
echo "4. This is useful for assigning priority or area to multiple issues at once"

# Log completion
echo "-----------------------------------------------" >> "$LOG_FILE"
echo "Setup instructions provided at $(date)" >> "$LOG_FILE"

echo -e "\n${GREEN}Project view setup instructions provided!${NC}"
echo -e "${YELLOW}Follow the manual steps to complete the setup.${NC}"
echo -e "${YELLOW}Log saved to:${NC} $LOG_FILE"

echo -e "\n${YELLOW}Project URL:${NC} $PROJECT_URL"