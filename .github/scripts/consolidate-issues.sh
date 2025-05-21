#!/bin/bash

# Consolidate Issues Script
# This script helps identify and manage duplicate/stale GitHub issues
# Usage: ./consolidate-issues.sh [--dry-run] [--close-test-issues] [--find-duplicates]

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default options
DRY_RUN=true
CLOSE_TEST_ISSUES=false
FIND_DUPLICATES=false
STALE_DAYS=30

# Parse command line options
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --dry-run)
      DRY_RUN=true
      ;;
    --no-dry-run)
      DRY_RUN=false
      ;;
    --close-test-issues)
      CLOSE_TEST_ISSUES=true
      ;;
    --find-duplicates)
      FIND_DUPLICATES=true
      ;;
    --stale-days=*)
      STALE_DAYS="${1#*=}"
      ;;
    *)
      echo -e "${RED}Unknown parameter: $1${NC}"
      exit 1
      ;;
  esac
  shift
done

echo -e "${BLUE}================================${NC}"
echo -e "${BLUE}GitHub Issue Consolidation Tool${NC}"
echo -e "${BLUE}================================${NC}"
echo ""
echo -e "Dry run mode: ${YELLOW}$DRY_RUN${NC}"
echo -e "Close test issues: ${YELLOW}$CLOSE_TEST_ISSUES${NC}"
echo -e "Find duplicates: ${YELLOW}$FIND_DUPLICATES${NC}"
echo -e "Stale threshold: ${YELLOW}$STALE_DAYS days${NC}"
echo ""

# Function to check if gh CLI is installed
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

# Function to close test issues
close_test_issues() {
  echo -e "${BLUE}Identifying test issues to close...${NC}"
  
  # Find issues with "Test Issue" in the title
  TEST_ISSUES=$(gh issue list --search "Test Issue in:title is:open" --json number,title,createdAt --limit 100)
  
  if [ -z "$TEST_ISSUES" ]; then
    echo -e "${GREEN}No test issues found.${NC}"
    return
  fi
  
  # Count issues
  NUM_ISSUES=$(echo "$TEST_ISSUES" | jq '. | length')
  echo -e "${YELLOW}Found $NUM_ISSUES test issues.${NC}"
  
  # Show list of issues
  echo "$TEST_ISSUES" | jq -r '.[] | "  #\(.number) - \(.title) (Created: \(.createdAt))"'
  
  # Ask for confirmation if not in dry run mode
  if [ "$DRY_RUN" = false ]; then
    echo ""
    read -p "Close these issues? (y/n): " CONFIRM
    
    if [[ $CONFIRM =~ ^[Yy]$ ]]; then
      echo -e "${BLUE}Closing test issues...${NC}"
      echo "$TEST_ISSUES" | jq -r '.[] | .number' | while read -r ISSUE_NUM; do
        echo -e "  Closing issue #${ISSUE_NUM}..."
        gh issue close "$ISSUE_NUM" --comment "Closing test issue as part of repository cleanup."
      done
      echo -e "${GREEN}Closed $NUM_ISSUES test issues.${NC}"
    else
      echo -e "${YELLOW}Operation cancelled.${NC}"
    fi
  else
    echo -e "${YELLOW}DRY RUN: Issues would be closed with the above confirmation.${NC}"
  fi
}

# Function to find potential duplicate issues
find_duplicate_issues() {
  echo -e "${BLUE}Looking for potential duplicate issues...${NC}"
  
  # Get all open issues
  ISSUES=$(gh issue list --state open --json number,title --limit 200)
  
  if [ -z "$ISSUES" ]; then
    echo -e "${GREEN}No issues found.${NC}"
    return
  fi
  
  # Count issues
  NUM_ISSUES=$(echo "$ISSUES" | jq '. | length')
  echo -e "${YELLOW}Analyzing $NUM_ISSUES issues for duplicates...${NC}"
  
  # Extract titles and numbers for comparison
  echo "$ISSUES" | jq -c '.[] | {number: .number, title: .title}' > /tmp/issue_titles.json
  
  # Find issues with similar titles (a simple approach)
  echo -e "${YELLOW}Potential duplicates:${NC}"
  cat /tmp/issue_titles.json | while read -r ISSUE; do
    ISSUE_NUM=$(echo "$ISSUE" | jq -r '.number')
    ISSUE_TITLE=$(echo "$ISSUE" | jq -r '.title')
    
    # Clean title for comparison (lowercase, remove special chars)
    CLEAN_TITLE=$(echo "$ISSUE_TITLE" | tr '[:upper:]' '[:lower:]' | sed 's/[^a-z0-9 ]//g')
    
    # Look for similar titles
    SIMILAR_ISSUES=$(cat /tmp/issue_titles.json | jq -r '. | select(.number != '$ISSUE_NUM') | {number: .number, title: .title}' | 
                     while read -r OTHER_ISSUE; do
                       OTHER_NUM=$(echo "$OTHER_ISSUE" | jq -r '.number')
                       OTHER_TITLE=$(echo "$OTHER_ISSUE" | jq -r '.title')
                       
                       # Clean other title
                       CLEAN_OTHER=$(echo "$OTHER_TITLE" | tr '[:upper:]' '[:lower:]' | sed 's/[^a-z0-9 ]//g')
                       
                       # Simple check: If 3 consecutive words match
                       if echo "$CLEAN_TITLE $CLEAN_OTHER" | grep -q -E "([a-z0-9]+ [a-z0-9]+ [a-z0-9]+).*\1"; then
                         echo "$OTHER_NUM"
                       fi
                     done)
    
    if [ ! -z "$SIMILAR_ISSUES" ]; then
      echo -e "  #$ISSUE_NUM - $ISSUE_TITLE"
      echo -e "  Similar to:"
      for SIM_NUM in $SIMILAR_ISSUES; do
        SIM_TITLE=$(cat /tmp/issue_titles.json | jq -r '. | select(.number == '$SIM_NUM') | .title')
        echo -e "    #$SIM_NUM - $SIM_TITLE"
      done
      echo ""
    fi
  done
  
  # Clean up
  rm -f /tmp/issue_titles.json
}

# Function to find stale issues
find_stale_issues() {
  echo -e "${BLUE}Looking for stale issues (not updated in $STALE_DAYS days)...${NC}"
  
  # Calculate date threshold
  THRESHOLD_DATE=$(date -d "-$STALE_DAYS days" +%Y-%m-%d)
  
  # Find stale issues
  STALE_ISSUES=$(gh issue list --search "updated:<$THRESHOLD_DATE is:open" --json number,title,updatedAt,labels --limit 100)
  
  if [ -z "$STALE_ISSUES" ] || [ "$(echo "$STALE_ISSUES" | jq '. | length')" -eq 0 ]; then
    echo -e "${GREEN}No stale issues found.${NC}"
    return
  fi
  
  # Count issues
  NUM_ISSUES=$(echo "$STALE_ISSUES" | jq '. | length')
  echo -e "${YELLOW}Found $NUM_ISSUES stale issues.${NC}"
  
  # Show list of issues
  echo "$STALE_ISSUES" | jq -r '.[] | "#\(.number) - \(.title) (Last updated: \(.updatedAt)) - Labels: \(.labels | map(.name) | join(", "))"'
  
  # Ask for labeling if not in dry run mode
  if [ "$DRY_RUN" = false ]; then
    echo ""
    read -p "Label these issues as 'stale'? (y/n): " CONFIRM
    
    if [[ $CONFIRM =~ ^[Yy]$ ]]; then
      echo -e "${BLUE}Labeling stale issues...${NC}"
      echo "$STALE_ISSUES" | jq -r '.[] | .number' | while read -r ISSUE_NUM; do
        # Check if issue already has stale label
        HAS_STALE=$(echo "$STALE_ISSUES" | jq -r '.[] | select(.number == '$ISSUE_NUM') | .labels | map(.name) | contains(["stale"]) | tostring')
        
        if [ "$HAS_STALE" = "false" ]; then
          echo -e "  Labeling issue #${ISSUE_NUM} as stale..."
          gh issue edit "$ISSUE_NUM" --add-label "stale" \
            --body-file <(echo -e "This issue has been automatically marked as stale because it has not had recent activity.\n\n$(gh issue view "$ISSUE_NUM" --json body -q .body)")
        else
          echo -e "  Issue #${ISSUE_NUM} is already labeled as stale."
        fi
      done
      echo -e "${GREEN}Labeled stale issues.${NC}"
    else
      echo -e "${YELLOW}Operation cancelled.${NC}"
    fi
  else
    echo -e "${YELLOW}DRY RUN: Issues would be labeled as stale with the above confirmation.${NC}"
  fi
}

# Main function
main() {
  check_gh
  
  if [ "$CLOSE_TEST_ISSUES" = true ]; then
    close_test_issues
    echo ""
  fi
  
  if [ "$FIND_DUPLICATES" = true ]; then
    find_duplicate_issues
    echo ""
  fi
  
  # Always find stale issues
  find_stale_issues
  
  echo -e "${GREEN}Issue consolidation process completed.${NC}"
}

# Execute main function
main