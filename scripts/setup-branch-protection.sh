#!/bin/bash
# Script to set up branch protection rules for the master branch

# Set colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

REPO="mushu-dev/cryoprotect"
BRANCH="master"

echo -e "${BLUE}Setting up branch protection rules for $BRANCH in $REPO${NC}"
echo "========================================"

# Set up branch protection using GitHub API
echo -e "${BLUE}Creating branch protection rules...${NC}"

# Basic protection rules to prevent direct pushes and require review before merging
gh api -X PUT "repos/$REPO/branches/$BRANCH/protection" \
  -f required_status_checks=null \
  -f enforce_admins=false \
  -f required_pull_request_reviews='{"required_approving_review_count":1,"dismiss_stale_reviews":true}' \
  -f restrictions=null

if [ $? -eq 0 ]; then
  echo -e "${GREEN}Branch protection rules successfully set up for $BRANCH!${NC}"
else
  echo -e "${RED}Failed to set up branch protection rules.${NC}"
  exit 1
fi

# Verify the branch protection rules
echo -e "${BLUE}Verifying branch protection rules...${NC}"
PROTECTION_INFO=$(gh api "repos/$REPO/branches/$BRANCH/protection")

if [ -n "$PROTECTION_INFO" ]; then
  echo -e "${GREEN}Branch protection rules verified!${NC}"
  echo -e "${BLUE}Protection Rules:${NC}"
  echo "$PROTECTION_INFO" | jq -r '.required_pull_request_reviews.required_approving_review_count + " review(s) required before merging"'
  if echo "$PROTECTION_INFO" | jq -e '.required_pull_request_reviews.dismiss_stale_reviews' > /dev/null; then
    echo "- Stale review dismissal: Enabled"
  fi
else
  echo -e "${RED}Failed to verify branch protection rules.${NC}"
  exit 1
fi

echo -e "\n${GREEN}Branch protection setup complete!${NC}"
echo -e "${YELLOW}Protection rules summary:${NC}"
echo "- Direct pushes to master are prevented"
echo "- Pull requests require at least 1 review before merging"
echo "- Stale reviews are dismissed when new commits are pushed"