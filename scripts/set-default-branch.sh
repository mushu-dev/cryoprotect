#!/bin/bash
# Script to set the default branch for the GitHub repository

# Set colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

REPO="mushu-dev/cryoprotect"
BRANCH="master"

echo -e "${BLUE}Setting default branch for $REPO to $BRANCH${NC}"
echo "========================================"

# Check if the branch exists
if ! gh api repos/$REPO/branches/$BRANCH &>/dev/null; then
  echo -e "${RED}Error: Branch '$BRANCH' does not exist in repository '$REPO'${NC}"
  exit 1
fi

# Set the default branch
echo -e "${BLUE}Changing default branch to $BRANCH...${NC}"
if gh api -X PATCH repos/$REPO -f default_branch="$BRANCH"; then
  echo -e "${GREEN}Default branch successfully changed to $BRANCH!${NC}"
else
  echo -e "${RED}Failed to change default branch.${NC}"
  exit 1
fi

# Verify the change
echo -e "${BLUE}Verifying default branch...${NC}"
DEFAULT_BRANCH=$(gh repo view $REPO --json defaultBranchRef --jq '.defaultBranchRef.name')

if [ "$DEFAULT_BRANCH" == "$BRANCH" ]; then
  echo -e "${GREEN}Verified: Default branch is now $DEFAULT_BRANCH${NC}"
else
  echo -e "${RED}Verification failed: Default branch is $DEFAULT_BRANCH, not $BRANCH${NC}"
  exit 1
fi

echo -e "\n${YELLOW}Important next steps:${NC}"
echo "1. Update any protection rules on the new default branch"
echo "2. Update CI/CD workflows to target the new default branch"
echo "3. Make sure team members know about the change"
echo "4. Update deployment platforms to point to the new default branch"
echo "5. Update any documentation that references the branch"