#!/bin/bash

# Script to create well-structured GitHub task issues within the proper epic
# Usage: ./create-task-issue.sh [--title "Title"] [--epic database|chembl|pubchem|api|frontend|auth|rdkit|infrastructure|testing] [--priority high|medium|low] [--description "Description"]

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
TITLE=""
EPIC="uncategorized"
PRIORITY="medium"
DESCRIPTION=""

# Parse command line options
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --title=*) TITLE="${1#*=}" ;;
    --title) TITLE="$2"; shift ;;
    --epic=*) EPIC="${1#*=}" ;;
    --epic) EPIC="$2"; shift ;;
    --priority=*) PRIORITY="${1#*=}" ;;
    --priority) PRIORITY="$2"; shift ;;
    --description=*) DESCRIPTION="${1#*=}" ;;
    --description) DESCRIPTION="$2"; shift ;;
    *) echo -e "${RED}Unknown parameter: $1${NC}"; exit 1 ;;
  esac
  shift
done

# Validate input
if [ -z "$TITLE" ]; then
  echo -e "${RED}Error: Title is required.${NC}"
  echo "Usage: ./create-task-issue.sh --title \"Title\" [--epic database|chembl|pubchem|api|frontend|auth|rdkit|infrastructure|testing] [--priority high|medium|low] [--description \"Description\"]"
  exit 1
fi

# Map epic to area label
case $EPIC in
  database) AREA="area:database" ;;
  chembl) AREA="area:database" ;;
  pubchem) AREA="area:database" ;;
  api) AREA="area:api" ;;
  frontend) AREA="area:ui" ;;
  auth) AREA="area:auth" ;;
  rdkit) AREA="area:chembl" ;;
  infrastructure) AREA="area:devops" ;;
  testing) AREA="area:testing" ;;
  *) AREA="" ;;
esac

# Map epic to related epic number
# Note: These should be updated with your actual epic issue numbers
EPIC_NUMBERS=(
  ["database"]="get-epic-number database"
  ["chembl"]="get-epic-number chembl"
  ["pubchem"]="get-epic-number pubchem"
  ["api"]="get-epic-number api"
  ["frontend"]="get-epic-number frontend"
  ["auth"]="get-epic-number auth"
  ["rdkit"]="get-epic-number rdkit"
  ["infrastructure"]="get-epic-number infrastructure"
  ["testing"]="get-epic-number testing"
)

# Function to get epic issue number (needs to be adapted for your repository)
get_epic_number() {
  local category="$1"
  # This is a placeholder - replace with actual code to get your epic issue numbers
  # For example, you could store them in a file or hardcode them here
  echo "UNKNOWN"
}

# Build label list
LABELS="type:feature,status:planning,priority:$PRIORITY"
if [ ! -z "$AREA" ]; then
  LABELS="$LABELS,$AREA"
fi

# Build issue body
BODY="# Task Description
$DESCRIPTION

## Epic
This task is part of the $EPIC epic.

## Acceptance Criteria
- [ ] Criterion 1
- [ ] Criterion 2
- [ ] Criterion 3

## Implementation Notes
Add implementation notes here.

## Dependencies
List any dependencies here.
"

echo -e "${BLUE}Creating new task issue...${NC}"
echo -e "${YELLOW}Title:${NC} $TITLE"
echo -e "${YELLOW}Epic:${NC} $EPIC"
echo -e "${YELLOW}Priority:${NC} $PRIORITY"
echo -e "${YELLOW}Labels:${NC} $LABELS"

# Create the issue
ISSUE=$(gh issue create --title "$TITLE" --body "$BODY" --label "$LABELS")
ISSUE_NUMBER=$(echo $ISSUE | grep -o '#[0-9]*' | tr -d '#')

echo -e "${GREEN}Created issue: $ISSUE${NC}"

# Add comment linking to epic if we know the epic number
if [ "$EPIC" != "uncategorized" ]; then
  EPIC_NUMBER=$(eval ${EPIC_NUMBERS[$EPIC]})
  if [ "$EPIC_NUMBER" != "UNKNOWN" ]; then
    echo -e "${BLUE}Linking to epic #$EPIC_NUMBER...${NC}"
    gh issue comment "$ISSUE_NUMBER" --body "This issue is part of epic #$EPIC_NUMBER"
  fi
fi

echo -e "${GREEN}Task issue created successfully!${NC}"