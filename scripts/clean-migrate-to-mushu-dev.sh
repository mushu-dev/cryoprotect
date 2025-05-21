#!/bin/bash
# Clean migration script from blueprint-house to mushu-dev (open issues only)

# Set colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Repository Clean Migration Tool${NC}"
echo "========================================"

# Step 1: Delete all existing issues in the new repository
echo -e "${BLUE}Step 1: Deleting existing issues in new repository...${NC}"
./scripts/delete-all-issues.sh

# Step 2: Create labels in the new repository
echo -e "\n${BLUE}Step 2: Creating labels in new repository...${NC}"
.github/scripts/create-labels.sh mushu-dev/cryoprotect

# Step 3: Create milestones in the new repository
echo -e "\n${BLUE}Step 3: Creating milestones in new repository...${NC}"
echo "Creating milestones from blueprint-house/cryoprotect to mushu-dev/cryoprotect..."

# Get milestones from the old repository
MILESTONES=$(gh api repos/blueprint-house/cryoprotect/milestones --paginate)

# Create each milestone in the new repository
echo "$MILESTONES" | jq -c '.[]' | while read -r milestone; do
  TITLE=$(echo "$milestone" | jq -r '.title')
  DESCRIPTION=$(echo "$milestone" | jq -r '.description')
  STATE=$(echo "$milestone" | jq -r '.state')
  DUE_DATE=$(echo "$milestone" | jq -r '.due_on')
  
  # Create the milestone in the new repository
  echo "Creating milestone: $TITLE"
  if [ "$DUE_DATE" != "null" ]; then
    gh api repos/mushu-dev/cryoprotect/milestones -f title="$TITLE" -f description="$DESCRIPTION" -f state="$STATE" -f due_on="$DUE_DATE" > /dev/null
  else
    gh api repos/mushu-dev/cryoprotect/milestones -f title="$TITLE" -f description="$DESCRIPTION" -f state="$STATE" > /dev/null
  fi
done

# Step 4: Recreate ONLY OPEN issues from the old repository to the new one
echo -e "\n${BLUE}Step 4: Recreating ONLY OPEN issues from old repository to new one...${NC}"
./scripts/recreate-open-issues.sh

# Step 5: Set up GitHub project board
echo -e "\n${BLUE}Step 5: Setting up GitHub project board...${NC}"
./scripts/setup-github-project.sh

# Step 6: Update local repository configuration
echo -e "\n${BLUE}Step 6: Verifying local repository configuration...${NC}"
REMOTE_URL=$(git config --get remote.origin.url)
echo "Current remote URL: $REMOTE_URL"

if [[ "$REMOTE_URL" != *"mushu-dev"* ]]; then
  echo "Updating remote URL to point to mushu-dev..."
  git remote set-url origin https://github.com/mushu-dev/cryoprotect.git
  echo "New remote URL: $(git config --get remote.origin.url)"
fi

# Verify branch configuration
BRANCH_REMOTE=$(git config --get branch.master.remote)
BRANCH_MERGE=$(git config --get branch.master.merge)

echo "Branch remote: $BRANCH_REMOTE"
echo "Branch merge: $BRANCH_MERGE"

if [ "$BRANCH_REMOTE" != "origin" ] || [ "$BRANCH_MERGE" != "refs/heads/master" ]; then
  echo "Updating branch configuration..."
  git config --local branch.master.remote origin
  git config --local branch.master.merge refs/heads/master
  echo "Updated branch configuration."
fi

# Step 7: Summary
echo -e "\n${GREEN}Repository clean migration completed!${NC}"
echo -e "${YELLOW}Summary:${NC}"
echo "- All existing issues deleted in mushu-dev/cryoprotect"
echo "- Labels created in mushu-dev/cryoprotect"
echo "- Milestones migrated from blueprint-house/cryoprotect to mushu-dev/cryoprotect"
echo "- OPEN issues recreated from blueprint-house/cryoprotect to mushu-dev/cryoprotect"
echo "- GitHub project board set up for mushu-dev/cryoprotect"
echo "- Local repository configuration verified and updated if needed"

echo -e "\n${BLUE}Next steps:${NC}"
echo "1. Verify the recreated issues in the new repository"
echo "2. Verify the project board setup and configure views"
echo "3. Set up branch protection rules in the new repository"
echo "4. Update deployment integrations (Vercel, Heroku, etc.)"
echo "5. Communicate the repository change to team members"

echo -e "\nRemember to update any references to the old repository URL in:"
echo "- Documentation files"
echo "- Configuration files"
echo "- CI/CD workflows" 
echo "- Deployment platforms"