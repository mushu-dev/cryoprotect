#!/bin/bash

# Script to close all remaining issues except for a specific set of kept issues
# Usage: ./close-remaining-issues.sh [--dry-run]

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
DRY_RUN=false
SLEEP_BETWEEN_ISSUES=2
BATCH_SIZE=10
SLEEP_BETWEEN_BATCHES=15

# Process arguments
if [ "$1" == "--dry-run" ]; then
  DRY_RUN=true
  echo -e "${YELLOW}Running in dry-run mode. No changes will be made.${NC}"
fi

# Issues to keep - MODIFY THIS LIST AS NEEDED
# Our epic issues (high-level categories) and a few essential specific issues
KEPT_ISSUES_STRING="246 245 244 243 211"

# Check if we should exclude epic issues
if [ -f "/tmp/epic_issues.txt" ]; then
  EPIC_ISSUES=$(cat /tmp/epic_issues.txt)
  KEPT_ISSUES_STRING="$KEPT_ISSUES_STRING $EPIC_ISSUES"
  echo -e "${GREEN}Found epic issues to keep: $EPIC_ISSUES${NC}"
else
  echo -e "${YELLOW}No epic issues file found. Only keeping predefined issues.${NC}"
fi

# Convert to array for easier checking
read -ra KEPT_ISSUES <<< "$KEPT_ISSUES_STRING"
echo -e "${BLUE}Keeping these issues:${NC} ${KEPT_ISSUES[@]}"

# Get all open issues
echo -e "${BLUE}Getting all open issues...${NC}"
ALL_ISSUES=$(gh issue list --json number --limit 200 | jq -r '.[].number')
ISSUE_COUNT=$(echo "$ALL_ISSUES" | wc -w)
echo -e "${YELLOW}Found $ISSUE_COUNT open issues${NC}"

# Identify issues to close
ISSUES_TO_CLOSE=()
for issue in $ALL_ISSUES; do
  KEEP=false
  for kept in "${KEPT_ISSUES[@]}"; do
    if [ "$issue" == "$kept" ]; then
      KEEP=true
      break
    fi
  done
  
  if [ "$KEEP" == "false" ]; then
    ISSUES_TO_CLOSE+=("$issue")
  fi
done

CLOSE_COUNT=${#ISSUES_TO_CLOSE[@]}
echo -e "${YELLOW}Found $CLOSE_COUNT issues to close${NC}"

# Close issues in batches
if [ $CLOSE_COUNT -gt 0 ]; then
  echo -e "${BLUE}Closing issues in batches of $BATCH_SIZE...${NC}"
  
  CURRENT=0
  BATCH=0
  
  for issue in "${ISSUES_TO_CLOSE[@]}"; do
    CURRENT=$((CURRENT + 1))
    BATCH=$((BATCH + 1))
    
    if [ "$DRY_RUN" == "true" ]; then
      echo -e "${YELLOW}[DRY RUN] [$CURRENT/$CLOSE_COUNT] Would close issue #$issue${NC}"
    else
      echo -e "${GREEN}[$CURRENT/$CLOSE_COUNT] Closing issue #$issue...${NC}"
      gh issue close "$issue" --comment "This issue has been closed as part of repository cleanup and consolidation.

The repository is being reorganized to use a streamlined issue structure with epics for better tracking and organization."
      sleep $SLEEP_BETWEEN_ISSUES
    fi
    
    # If batch is complete, sleep to avoid rate limiting
    if [ $BATCH -eq $BATCH_SIZE ]; then
      echo -e "${BLUE}Completed batch of $BATCH_SIZE issues. Sleeping for $SLEEP_BETWEEN_BATCHES seconds...${NC}"
      sleep $SLEEP_BETWEEN_BATCHES
      BATCH=0
    fi
  done
  
  echo -e "${GREEN}Completed closing $CLOSE_COUNT issues${NC}"
else
  echo -e "${GREEN}No issues to close${NC}"
fi

echo -e "${BLUE}=== SUMMARY ===${NC}"
echo -e "${GREEN}Kept ${#KEPT_ISSUES[@]} issues${NC}"
echo -e "${GREEN}Closed $CLOSE_COUNT issues${NC}"
echo -e "${YELLOW}Total remaining: ${#KEPT_ISSUES[@]}${NC}"

echo -e "${BLUE}=== RECOMMENDED NEXT STEPS ===${NC}"
echo "1. Review the remaining issues to ensure they cover all important areas"
echo "2. Set up the project board with the remaining issues"
echo "3. Continue with ongoing work using the new streamlined issue structure"