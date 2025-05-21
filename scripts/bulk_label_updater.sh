#!/bin/bash
# Script to help update GitHub issue labels in bulk

echo "GitHub Issue Bulk Label Updater"
echo "==============================="

# Check if issue numbers were provided
if [ "$#" -eq 0 ]; then
  # If no parameters, list issues with minimal descriptions for selection
  echo "Finding issues with minimal descriptions..."
  gh issue list --limit 100 --json number,title,body,labels > issues.json
  
  # Extract issues with minimal descriptions
  MINIMAL_ISSUES=$(jq -r '.[] | select(.body | test("^# Task Description") or length < 100) | "\(.number):\t\(.title)"' issues.json)
  
  echo -e "\nIssues with minimal descriptions:"
  echo "$MINIMAL_ISSUES"
  
  echo -e "\nEnter issue numbers to update (space-separated):"
  read -a ISSUE_NUMBERS
else
  # Use provided issue numbers
  ISSUE_NUMBERS=("$@")
fi

if [ ${#ISSUE_NUMBERS[@]} -eq 0 ]; then
  echo "No issues selected. Exiting."
  exit 0
fi

# Get available labels
echo -e "\nFetching available labels..."
LABELS=$(gh api repos/blueprint-house/CryoProtect/labels --paginate | jq -r '.[].name')

# Display available labels by category
echo -e "\nAvailable labels:"
echo -e "\nStatus labels:"
echo "$LABELS" | grep "status:" | sort
echo -e "\nType labels:"
echo "$LABELS" | grep "type:" | sort
echo -e "\nOther labels:"
echo "$LABELS" | grep -v "status:" | grep -v "type:" | sort

# Prompt for status label
echo -e "\nEnter status label to apply (e.g., status:planning):"
read STATUS_LABEL

# Prompt for type label
echo -e "\nEnter type label to apply (e.g., type:feature):"
read TYPE_LABEL

# Prompt for additional labels
echo -e "\nEnter any additional labels (space-separated):"
read -a ADDITIONAL_LABELS

# Process each issue
for ISSUE in "${ISSUE_NUMBERS[@]}"; do
  echo -e "\nUpdating issue #$ISSUE..."
  
  # Get current labels
  CURRENT_LABELS=$(gh issue view $ISSUE --json labels | jq -r '.labels[].name')
  
  # Remove existing status and type labels if new ones are provided
  if [ ! -z "$STATUS_LABEL" ]; then
    for LABEL in $CURRENT_LABELS; do
      if [[ $LABEL == status:* ]]; then
        gh issue edit $ISSUE --remove-label "$LABEL"
      fi
    done
    gh issue edit $ISSUE --add-label "$STATUS_LABEL"
  fi
  
  if [ ! -z "$TYPE_LABEL" ]; then
    for LABEL in $CURRENT_LABELS; do
      if [[ $LABEL == type:* ]]; then
        gh issue edit $ISSUE --remove-label "$LABEL"
      fi
    done
    gh issue edit $ISSUE --add-label "$TYPE_LABEL"
  fi
  
  # Add additional labels
  for LABEL in "${ADDITIONAL_LABELS[@]}"; do
    if [ ! -z "$LABEL" ]; then
      gh issue edit $ISSUE --add-label "$LABEL"
    fi
  done
  
  echo "Issue #$ISSUE updated successfully."
done

echo -e "\nAll issues updated successfully."