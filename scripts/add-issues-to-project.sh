#!/bin/bash
# Simple script to add all issues to GitHub project using basic gh commands

# Set up variables
REPO="mushu-dev/cryoprotect"
PROJECT_NUMBER="1"
PROJECT_ID="PVT_kwHOBMCBws4A4xSP"

echo "Adding all issues from $REPO to project #$PROJECT_NUMBER"
echo "============================================="

# Get all issues
echo "Fetching issues..."
ISSUES=$(gh issue list -R "$REPO" -s all -L 100 --json number,title,url)

if [ -z "$ISSUES" ]; then
  echo "No issues found or error fetching issues"
  exit 1
fi

# Count issues
TOTAL_ISSUES=$(echo "$ISSUES" | jq '. | length')
echo "Found $TOTAL_ISSUES issues to add"

# Add each issue to project
COUNT=0
echo "$ISSUES" | jq -c '.[]' | while read -r issue; do
  ISSUE_NUMBER=$(echo "$issue" | jq -r '.number')
  ISSUE_TITLE=$(echo "$issue" | jq -r '.title')
  ISSUE_URL=$(echo "$issue" | jq -r '.url')
  
  echo "Adding issue #$ISSUE_NUMBER: $ISSUE_TITLE"
  
  # Add to project using the gh api command
  gh api graphql -f query='
    mutation($project:ID!, $issue:ID!) {
      addProjectV2ItemById(input: {projectId: $project, contentId: $issue}) {
        item {
          id
        }
      }
    }
  ' -f project="$PROJECT_ID" -F issue="$(gh api repos/$REPO/issues/$ISSUE_NUMBER --jq '.node_id')" > /dev/null
  
  COUNT=$((COUNT + 1))
  echo "Progress: $COUNT / $TOTAL_ISSUES"
  
  # Sleep to avoid rate limiting
  sleep 1
done

echo "============================================="
echo "Completed! Added $COUNT issues to project #$PROJECT_NUMBER"