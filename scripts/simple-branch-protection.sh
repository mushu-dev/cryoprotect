#!/bin/bash
# Simple script to set up branch protection for master branch

# Set colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

REPO="mushu-dev/cryoprotect"
BRANCH="master"

echo -e "${BLUE}Setting up simple branch protection for $BRANCH in $REPO${NC}"
echo "========================================"

# Set up basic branch protection using GitHub API
echo -e "${BLUE}Creating branch protection...${NC}"

# Create temporary JSON file for the request body
TEMP_JSON=$(mktemp)

cat > "$TEMP_JSON" << EOF
{
  "required_status_checks": null,
  "enforce_admins": false,
  "required_pull_request_reviews": {
    "dismissal_restrictions": {},
    "dismiss_stale_reviews": true,
    "require_code_owner_reviews": false,
    "required_approving_review_count": 1
  },
  "restrictions": null
}
EOF

# Send the request
curl -s -X PUT \
  -H "Accept: application/vnd.github+json" \
  -H "Authorization: token $(gh auth token)" \
  -H "X-GitHub-Api-Version: 2022-11-28" \
  "https://api.github.com/repos/$REPO/branches/$BRANCH/protection" \
  -d @"$TEMP_JSON" > /dev/null

if [ $? -eq 0 ]; then
  echo -e "${GREEN}Branch protection enabled for $BRANCH!${NC}"
else
  echo -e "${RED}Failed to set up branch protection.${NC}"
  rm "$TEMP_JSON"
  exit 1
fi

# Clean up
rm "$TEMP_JSON"

echo -e "\n${GREEN}Branch protection setup complete!${NC}"
echo -e "${YELLOW}Summary:${NC}"
echo "- Pull requests are now required for the master branch"
echo "- At least 1 review is required before merging"
echo "- Stale reviews are dismissed when new commits are pushed"