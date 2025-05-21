#!/bin/bash
# Script to help update GitHub issue descriptions with more comprehensive content

ISSUE_NUMBER=$1

if [ -z "$ISSUE_NUMBER" ]; then
  echo "Usage: $0 <issue_number>"
  exit 1
fi

# Get current issue details
gh issue view $ISSUE_NUMBER --json number,title,body > issue_details.json

# Extract the title and current body
TITLE=$(jq -r '.title' issue_details.json)
BODY=$(jq -r '.body' issue_details.json)

# Determine the issue type from the title
TYPE="feature"
if [[ "$TITLE" == *"Bug"* || "$TITLE" == *"Fix"* ]]; then
  TYPE="bugfix"
elif [[ "$TITLE" == *"Validation"* || "$TITLE" == *"Test"* ]]; then
  TYPE="validation"
elif [[ "$TITLE" == *"Documentation"* || "$TITLE" == *"Doc"* ]]; then
  TYPE="documentation"
fi

# Create a temporary file for editing
TMP_FILE=$(mktemp)

# Put existing content (if not just "# Task Description") as a starting point
if [[ "$BODY" == "# Task Description"* && ${#BODY} -lt 20 ]]; then
  # Create a new template
  cat > $TMP_FILE << EOF
# ${TYPE^}: $TITLE

## Description
[Add a detailed description of the task or issue here]

## Acceptance Criteria
- [ ] [Criterion 1]
- [ ] [Criterion 2]
- [ ] [Criterion 3]

## Dependencies
[List any issues this depends on, e.g., "Depends on #XX"]

## Technical Details
[Add any technical information or implementation notes here]
EOF
else
  # Use existing content as a starting point
  echo "$BODY" > $TMP_FILE
fi

# Open the temporary file in the default editor
${EDITOR:-vim} $TMP_FILE

# Check if the file was modified
if cmp -s "$TMP_FILE" <(echo "$BODY"); then
  echo "No changes made. Issue description not updated."
  rm $TMP_FILE
  exit 0
fi

# Update the issue with the new description
echo "Updating issue #$ISSUE_NUMBER with the new description..."
gh issue edit $ISSUE_NUMBER --body-file $TMP_FILE

echo "Issue #$ISSUE_NUMBER updated successfully."

# Clean up
rm $TMP_FILE