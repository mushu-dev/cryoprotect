#!/bin/bash

# Issue Cleanup Script
# This script performs a complete GitHub issue cleanup
# - Closes test issues
# - Identifies and labels stale issues
# - Finds potential duplicates
# - Detects AI-generated spam
# - Updates project board

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

echo -e "${BLUE}====================================${NC}"
echo -e "${BLUE}GitHub Issue Comprehensive Cleanup${NC}"
echo -e "${BLUE}====================================${NC}"
echo ""

# Check if the user wants to run in dry run mode
DRY_RUN="--dry-run"
read -p "Run in dry-run mode (no changes will be made)? (Y/n): " DRY_RUN_RESPONSE
if [[ $DRY_RUN_RESPONSE =~ ^[Nn]$ ]]; then
  DRY_RUN="--no-dry-run"
  echo -e "${YELLOW}WARNING: Changes will be applied!${NC}"
else
  echo -e "${GREEN}Running in dry-run mode (no changes will be made).${NC}"
fi

echo ""
echo -e "${BLUE}Step 1: Running issue consolidation script...${NC}"
"$SCRIPT_DIR/consolidate-issues.sh" $DRY_RUN --close-test-issues --find-duplicates --stale-days=30

echo ""
echo -e "${BLUE}Step 2: Detecting AI-generated issues...${NC}"
"$SCRIPT_DIR/detect-ai-spam.py" --threshold=0.75 $([[ $DRY_RUN == "--no-dry-run" ]] && echo "--close-confirmed")

echo ""
echo -e "${BLUE}Step 3: Synchronizing issues with project board...${NC}"
if [[ $DRY_RUN == "--no-dry-run" ]]; then
  echo -e "Syncing project board..."
  
  # Get project ID
  PROJECT_ID=$(gh project list --owner "$(gh repo view --json owner -q .owner.login)" --format json | jq -r '.[0].id')
  
  if [ -z "$PROJECT_ID" ]; then
    echo -e "${RED}Error: Could not find project ID.${NC}"
  else
    echo -e "Found project ID: $PROJECT_ID"
    
    # Add recent issues to the project
    echo -e "Adding recent issues to project..."
    RECENT_ISSUES=$(gh issue list --json number --limit 20 --search "created:>2025-05-01" | jq -r '.[].number')
    
    for ISSUE_NUM in $RECENT_ISSUES; do
      echo -e "  Checking issue #${ISSUE_NUM}..."
      
      # Check if issue is already in project (this is a simplified check)
      ALREADY_IN_PROJECT=$(gh project item-list "$PROJECT_ID" --owner "$(gh repo view --json owner -q .owner.login)" --format json | jq -r '.items[] | select(.content.number == '$ISSUE_NUM') | .id')
      
      if [ -z "$ALREADY_IN_PROJECT" ]; then
        echo -e "  Adding issue #${ISSUE_NUM} to project..."
        gh project item-add "$PROJECT_ID" --owner "$(gh repo view --json owner -q .owner.login)" --url "$(gh repo view --json url -q .url)/issues/$ISSUE_NUM"
      else
        echo -e "  Issue #${ISSUE_NUM} is already in project."
      fi
    done
  fi
else
  echo -e "${YELLOW}Skipping project board sync in dry-run mode.${NC}"
fi

echo ""
echo -e "${GREEN}Cleanup process completed.${NC}"
echo -e "${BLUE}For complete details on issue management best practices, see GITHUB_ISSUE_MANAGEMENT_GUIDE.md${NC}"