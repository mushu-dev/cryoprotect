#!/bin/bash
# Script to add all issues to the GitHub project board

# Set colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Set up variables
REPO="mushu-dev/cryoprotect"
LOG_FILE="project_board_population_$(date +%Y%m%d_%H%M%S).log"

echo -e "${BLUE}CryoProtect Project Board Population Tool${NC}"
echo "========================================"
echo "This script will add all issues to the project board"
echo "Log will be saved to $LOG_FILE"

# Log start time
echo "Population started at $(date)" > "$LOG_FILE"
echo "Repository: $REPO" >> "$LOG_FILE"
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

# Get all issues from repository
echo -e "${BLUE}Fetching issues from $REPO...${NC}"
ISSUES=$(gh issue list -s all -L 500 --json number,title,url -R "$REPO")

# Check if we got any issues
if [ -z "$ISSUES" ]; then
  echo -e "${RED}No issues found in $REPO or error fetching issues${NC}"
  echo "No issues found in $REPO or error fetching issues" >> "$LOG_FILE"
  exit 1
fi

# Count issues
TOTAL_ISSUES=$(echo "$ISSUES" | jq '. | length')
echo -e "${BLUE}Found $TOTAL_ISSUES issues to add to project board${NC}"
echo "Found $TOTAL_ISSUES issues to add to project board" >> "$LOG_FILE"

# Confirm with user
read -p "Do you want to add all $TOTAL_ISSUES issues to the project board? (y/N): " CONFIRM

if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
  echo "Operation cancelled by user."
  echo "Operation cancelled by user at $(date)" >> "$LOG_FILE"
  exit 0
fi

# Process each issue
ADDED_COUNT=0
FAILED_COUNT=0

echo "$ISSUES" | jq -c '.[]' | while read -r issue; do
  ISSUE_NUMBER=$(echo "$issue" | jq -r '.number')
  ISSUE_TITLE=$(echo "$issue" | jq -r '.title')
  ISSUE_URL=$(echo "$issue" | jq -r '.url')
  
  echo -e "${BLUE}Adding issue #$ISSUE_NUMBER to project...${NC}"
  
  # Add to project
  if gh project item-add "$PROJECT_NUMBER" --owner mushu-dev --url "$ISSUE_URL"; then
    echo "[SUCCESS] Added issue #$ISSUE_NUMBER to project: $ISSUE_TITLE" >> "$LOG_FILE"
    ADDED_COUNT=$((ADDED_COUNT + 1))
  else
    echo "[FAILED] Could not add issue #$ISSUE_NUMBER to project: $ISSUE_TITLE" >> "$LOG_FILE"
    FAILED_COUNT=$((FAILED_COUNT + 1))
  fi
  
  echo -e "${GREEN}Progress:${NC} $ADDED_COUNT / $TOTAL_ISSUES added"
  
  # Sleep to avoid rate limiting
  sleep 1
done

# Log completion
echo "-----------------------------------------------" >> "$LOG_FILE"
echo "Population completed at $(date)" >> "$LOG_FILE"
echo "Total issues: $TOTAL_ISSUES" >> "$LOG_FILE"
echo "Successfully added: $ADDED_COUNT" >> "$LOG_FILE"
echo "Failed additions: $FAILED_COUNT" >> "$LOG_FILE"

echo -e "\n${GREEN}Project board population complete!${NC}"
echo -e "${YELLOW}Summary:${NC}"
echo "- Successfully added: $ADDED_COUNT issues"
echo "- Failed additions: $FAILED_COUNT issues"
echo "- Log saved to $LOG_FILE"

echo -e "\n${BLUE}Next steps:${NC}"
echo "1. Visit the project URL to configure views (Kanban, etc.)"
echo "2. Set up automation rules for status changes"
echo "3. Map status field to issue labels for better synchronization"
echo "4. Create custom filtered views based on priorities or areas"

echo -e "\n${YELLOW}Project URL:${NC} $PROJECT_URL"