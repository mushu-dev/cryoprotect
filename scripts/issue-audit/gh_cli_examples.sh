#!/bin/bash
# GitHub CLI Examples for Issue Management
# This script demonstrates common GitHub CLI commands for issue management

# Colors for terminal output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}GitHub CLI Examples for Issue Management${NC}"
echo "========================================="

# Check if gh CLI is installed and authenticated
if ! command -v gh &> /dev/null; then
    echo -e "${RED}Error: GitHub CLI (gh) is not installed or not in PATH${NC}"
    echo "Please install it from: https://cli.github.com/"
    exit 1
fi

# Check authentication status
auth_status=$(gh auth status 2>&1)
if echo "$auth_status" | grep -q "not logged"; then
    echo -e "${RED}Error: GitHub CLI is not authenticated${NC}"
    echo "Please run 'gh auth login' to authenticate"
    exit 1
fi

# Display current repository
current_repo=$(gh repo view --json name,owner | jq -r '"\(.owner.login)/\(.name)"')
echo -e "${BLUE}Current repository:${NC} $current_repo"
echo ""

# List examples by category
echo -e "${BLUE}GitHub CLI Examples by Category:${NC}"
echo ""

echo -e "${GREEN}1. Basic Issue Management${NC}"
echo "# List all open issues"
echo "gh issue list"
echo ""
echo "# List closed issues"
echo "gh issue list --state closed"
echo ""
echo "# List all issues with a specific label"
echo "gh issue list --label bug"
echo ""
echo "# View a specific issue"
echo "gh issue view 123"
echo ""
echo "# View a specific issue in JSON format"
echo "gh issue view 123 --json number,title,body,labels,assignees,milestone"
echo ""

echo -e "${GREEN}2. Creating and Editing Issues${NC}"
echo "# Create a new issue"
echo "gh issue create --title \"New Issue Title\" --body \"Issue description goes here\""
echo ""
echo "# Create an issue from a file"
echo "gh issue create --title \"New Issue Title\" --body-file issue_template.md"
echo ""
echo "# Edit an issue"
echo "gh issue edit 123 --title \"Updated Title\" --body \"Updated description\""
echo ""
echo "# Add a label to an issue"
echo "gh issue edit 123 --add-label bug"
echo ""
echo "# Remove a label from an issue"
echo "gh issue edit 123 --remove-label bug"
echo ""
echo "# Set a milestone"
echo "gh issue edit 123 --milestone \"Phase 1\""
echo ""
echo "# Close an issue"
echo "gh issue close 123"
echo ""
echo "# Reopen an issue"
echo "gh issue reopen 123"
echo ""

echo -e "${GREEN}3. Commenting on Issues${NC}"
echo "# Add a comment to an issue"
echo "gh issue comment 123 --body \"This is a comment\""
echo ""
echo "# List comments on an issue"
echo "gh issue view 123 --comments"
echo ""

echo -e "${GREEN}4. Advanced Issue Searching${NC}"
echo "# Search for issues with complex queries"
echo "gh issue list --search \"label:bug state:open\""
echo ""
echo "# Search for issues created by a specific user"
echo "gh issue list --search \"author:username\""
echo ""
echo "# Search for issues mentioning specific text"
echo "gh issue list --search \"database in:title,body\""
echo ""
echo "# Search for issues updated after a specific date"
echo "gh issue list --search \"updated:>2023-01-01\""
echo ""

echo -e "${GREEN}5. Bulk Operations with JQ${NC}"
echo "# Get all issues as JSON"
echo "gh issue list --limit 100 --json number,title,state,labels > issues.json"
echo ""
echo "# Filter issues with jq"
echo "cat issues.json | jq '.[] | select(.labels[].name == \"bug\")'"
echo ""
echo "# Process issues in a loop"
echo "gh issue list --json number,title --limit 10 | jq -c '.[]' | while read -r issue; do"
echo "  number=\$(echo \"\$issue\" | jq -r '.number')"
echo "  title=\$(echo \"\$issue\" | jq -r '.title')"
echo "  echo \"Processing issue #\$number: \$title\""
echo "done"
echo ""

echo -e "${GREEN}6. Label Management${NC}"
echo "# List all labels"
echo "gh api repos/OWNER/REPO/labels --paginate"
echo ""
echo "# Create a new label"
echo "gh api repos/OWNER/REPO/labels -X POST -f name=\"bug\" -f color=\"ff0000\" -f description=\"Bug reports\""
echo ""
echo "# Update a label"
echo "gh api repos/OWNER/REPO/labels/bug -X PATCH -f color=\"cc0000\" -f description=\"Updated description\""
echo ""
echo "# Delete a label"
echo "gh api repos/OWNER/REPO/labels/bug -X DELETE"
echo ""

echo -e "${GREEN}7. Milestone Management${NC}"
echo "# List all milestones"
echo "gh api repos/OWNER/REPO/milestones"
echo ""
echo "# Create a new milestone"
echo "gh api repos/OWNER/REPO/milestones -X POST -f title=\"Phase 1\" -f state=\"open\" -f description=\"First phase\""
echo ""
echo "# Update a milestone"
echo "gh api repos/OWNER/REPO/milestones/1 -X PATCH -f description=\"Updated description\""
echo ""
echo "# Close a milestone"
echo "gh api repos/OWNER/REPO/milestones/1 -X PATCH -f state=\"closed\""
echo ""

echo -e "${GREEN}8. Automated Issue Templates${NC}"
echo "# Create a template directory"
echo "mkdir -p .github/ISSUE_TEMPLATE"
echo ""
echo "# Create a template file (feature_request.md)"
echo "cat > .github/ISSUE_TEMPLATE/feature_request.md << EOL"
echo "---"
echo "name: Feature Request"
echo "about: Suggest a new feature"
echo "title: \"[FEATURE] \""
echo "labels: enhancement"
echo "---"
echo ""
echo "## Description"
echo "<!-- Describe the feature -->"
echo ""
echo "## Requirements"
echo "<!-- List specific requirements -->"
echo "- [ ] Requirement 1"
echo "- [ ] Requirement 2"
echo "EOL"
echo ""

echo -e "${GREEN}9. Useful Scripts${NC}"
echo "# Close all issues with a specific label"
echo "gh issue list --label \"wontfix\" --json number | jq -r '.[].number' | xargs -I{} gh issue close {}"
echo ""
echo "# Add a label to all open issues without labels"
echo "gh issue list --json number,labels --search \"no:label\" | jq -r '.[] | .number' | xargs -I{} gh issue edit {} --add-label \"needs-triage\""
echo ""
echo "# Move all issues from one milestone to another"
echo "gh issue list --milestone \"Old Milestone\" --json number | jq -r '.[].number' | xargs -I{} gh issue edit {} --milestone \"New Milestone\""
echo ""

echo -e "${BLUE}End of Examples${NC}"
echo "For more details, see: https://cli.github.com/manual/"