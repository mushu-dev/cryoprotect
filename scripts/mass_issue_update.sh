#!/bin/bash
# Mass issue update script to systematically fix all issues

echo "CryoProtect GitHub Issue Mass Update Tool"
echo "========================================="

# Set up variables
BATCH_SIZE=5
PROCESSED_FILE="processed_issues.txt"
SKIP_FILE="skipped_issues.txt"
MIN_DESCRIPTION_LENGTH=100
START_ISSUE=${1:-1}

# Create or clear processed and skipped files
> $PROCESSED_FILE
> $SKIP_FILE

# Initialize area labels if they don't exist
echo "Checking for required labels..."
gh api repos/blueprint-house/CryoProtect/labels --paginate > existing_labels.json

# Function to create label if it doesn't exist
create_label_if_missing() {
  local name=$1
  local color=$2
  local description=$3

  # Check if label exists
  if ! jq -e ".[] | select(.name == \"$name\")" existing_labels.json > /dev/null; then
    echo "Creating label: $name"
    gh api repos/blueprint-house/CryoProtect/labels --method POST \
      -f name="$name" -f color="$color" -f description="$description"
  else
    echo "Label already exists: $name"
  fi
}

# Create area labels
create_label_if_missing "area:database" "0052cc" "Database-related issues"
create_label_if_missing "area:api" "5319e7" "API-related issues"
create_label_if_missing "area:documentation" "0e8a16" "Documentation-related issues"
create_label_if_missing "area:testing" "fbca04" "Testing infrastructure issues"
create_label_if_missing "area:ui" "d93f0b" "User interface issues"
create_label_if_missing "area:devops" "b60205" "CI/CD and deployment issues"

# Create priority labels
create_label_if_missing "priority:critical" "b60205" "Critical priority issues"
create_label_if_missing "priority:high" "d93f0b" "High priority issues"
create_label_if_missing "priority:medium" "fbca04" "Medium priority issues"
create_label_if_missing "priority:low" "0e8a16" "Low priority issues"

# Get all issues and save to file
echo "Fetching all issues..."
gh issue list --limit 300 --json number,title,body,labels > all_issues.json

# Count total issues
TOTAL_ISSUES=$(jq '. | length' all_issues.json)
echo "Found $TOTAL_ISSUES issues to process"

# Get all issues with minimal descriptions
echo "Finding issues with minimal descriptions..."
MINIMAL_ISSUES=$(jq -r '.[] | select(.body | length < '$MIN_DESCRIPTION_LENGTH' or test("^# Task Description")) | .number' all_issues.json | sort -n)
MINIMAL_COUNT=$(echo "$MINIMAL_ISSUES" | wc -l)
echo "Found $MINIMAL_COUNT issues with minimal descriptions"

# Function to determine issue area based on title
determine_area() {
  local title=$(echo "$1" | tr '[:upper:]' '[:lower:]')
  
  if [[ "$title" == *"database"* || "$title" == *"supabase"* || "$title" == *"db"* || "$title" == *"sql"* || "$title" == *"schema"* || "$title" == *"table"* || "$title" == *"chembl"* || "$title" == *"pubchem"* ]]; then
    echo "area:database"
  elif [[ "$title" == *"api"* || "$title" == *"endpoint"* || "$title" == *"rest"* || "$title" == *"route"* ]]; then
    echo "area:api"
  elif [[ "$title" == *"documentation"* || "$title" == *"doc"* || "$title" == *"guide"* || "$title" == *"readme"* ]]; then
    echo "area:documentation"
  elif [[ "$title" == *"test"* || "$title" == *"validation"* || "$title" == *"verify"* ]]; then
    echo "area:testing"
  elif [[ "$title" == *"ui"* || "$title" == *"interface"* || "$title" == *"component"* || "$title" == *"layout"* ]]; then
    echo "area:ui"
  elif [[ "$title" == *"ci"* || "$title" == *"cd"* || "$title" == *"deploy"* || "$title" == *"build"* || "$title" == *"container"* ]]; then
    echo "area:devops"
  else
    echo ""
  fi
}

# Function to determine issue type based on title
determine_type() {
  local title=$(echo "$1" | tr '[:upper:]' '[:lower:]')
  
  if [[ "$title" == *"bug"* || "$title" == *"fix"* || "$title" == *"issue"* || "$title" == *"error"* ]]; then
    echo "type:bugfix"
  elif [[ "$title" == *"feature"* || "$title" == *"implement"* || "$title" == *"add"* || "$title" == *"new"* ]]; then
    echo "type:feature"
  elif [[ "$title" == *"documentation"* || "$title" == *"doc"* ]]; then
    echo "type:documentation"
  elif [[ "$title" == *"refactor"* || "$title" == *"improve"* || "$title" == *"enhance"* ]]; then
    echo "type:refactor"
  elif [[ "$title" == *"test"* || "$title" == *"validation"* || "$title" == *"verify"* ]]; then
    echo "type:validation"
  else
    echo "type:feature"  # Default to feature if can't determine
  fi
}

# Function to determine status based on existing labels
determine_status() {
  local labels="$1"
  
  if [[ "$labels" == *"status:Done"* || "$labels" == *"status:Implemented"* || "$labels" == *"status:Validated"* ]]; then
    echo "status:completed"
  elif [[ "$labels" == *"status:Running"* ]]; then
    echo "status:in-progress"
  elif [[ "$labels" == *"status:Failed"* || "$labels" == *"status:Error"* ]]; then
    echo "status:blocked"
  elif [[ "$labels" == *"status:Pending"* ]]; then
    echo "status:ready"
  elif [[ "$labels" == *"status:Blocked"* ]]; then
    echo "status:blocked"
  else
    echo "status:planning"  # Default to planning if no status
  fi
}

# Function to update issue description
update_description() {
  local issue_number=$1
  local title=$2
  local body=$3
  local issue_type=$4
  
  echo "Updating description for issue #$issue_number"
  
  # Skip issues that already have good descriptions
  if [[ ${#body} -gt $MIN_DESCRIPTION_LENGTH && "$body" != "# Task Description"* ]]; then
    echo "Issue #$issue_number already has a good description, skipping"
    echo "$issue_number" >> $SKIP_FILE
    return 0
  fi
  
  # Create template based on issue type
  local template=""
  if [[ "$issue_type" == "type:bugfix" ]]; then
    template="# Bug Report: $title

## Description
This issue addresses a bug related to [describe the affected functionality].

## Steps to Reproduce
1. [First Step]
2. [Second Step]
3. [Third Step]

## Expected Behavior
[What should happen]

## Current Behavior
[What actually happens]

## Fix Approach
[Suggested approach to fix the issue]

## Acceptance Criteria
- [ ] The bug has been fixed
- [ ] Unit tests have been added/updated
- [ ] The fix has been verified in [relevant environment]"
  elif [[ "$issue_type" == "type:validation" ]]; then
    template="# Validation: $title

## Description
This issue is for validating the changes made in [reference issue number].

## Validation Steps
1. [First Step]
2. [Second Step]
3. [Third Step]

## Expected Results
[What should be observed when validation is successful]

## Acceptance Criteria
- [ ] All validation steps executed successfully
- [ ] Results documented and shared
- [ ] Any issues found have been reported"
  else
    template="# Feature: $title

## Description
[Detailed description of the feature or task]

## Implementation Details
[Technical details about how this should be implemented]

## Dependencies
[List any issues this depends on]

## Acceptance Criteria
- [ ] [Criterion 1]
- [ ] [Criterion 2]
- [ ] [Criterion 3]

## Notes
[Any additional information, references, or context]"
  fi
  
  # Write template to temporary file
  echo "$template" > temp_description.txt
  
  # Update issue with new description
  gh issue edit $issue_number --body-file temp_description.txt
  
  # Add to processed file
  echo "$issue_number" >> $PROCESSED_FILE
  
  # Small delay to avoid rate limiting
  sleep 1
}

# Function to update issue labels
update_labels() {
  local issue_number=$1
  local title=$2
  local current_labels=$3
  
  echo "Updating labels for issue #$issue_number"
  
  # Determine area, type, and status
  local area=$(determine_area "$title")
  local type=$(determine_type "$title")
  local status=$(determine_status "$current_labels")
  
  # Add area label if determined
  if [[ ! -z "$area" ]]; then
    gh issue edit $issue_number --add-label "$area"
  fi
  
  # Update type label
  # First remove any existing type: labels
  for label in $(echo "$current_labels" | grep -o "type:[^ ,]*"); do
    gh issue edit $issue_number --remove-label "$label"
  done
  # Add new type label
  gh issue edit $issue_number --add-label "$type"
  
  # Update status label
  # First remove any existing status: labels
  for label in $(echo "$current_labels" | grep -o "status:[^ ,]*"); do
    gh issue edit $issue_number --remove-label "$label"
  done
  # Add new status label
  gh issue edit $issue_number --add-label "$status"
  
  # Small delay to avoid rate limiting
  sleep 1
}

# Process issues in batches
echo -e "\nBeginning issue processing from issue #$START_ISSUE..."
COUNTER=0
BATCH_COUNTER=0

# Extract just the issue numbers from the minimal issues
IFS=$'\n' read -d '' -a ISSUE_NUMBERS <<< "$MINIMAL_ISSUES"

for ISSUE in "${ISSUE_NUMBERS[@]}"; do
  # Skip issues that are below the starting point
  if [[ $ISSUE -lt $START_ISSUE ]]; then
    continue
  fi
  
  # Get issue details
  echo -e "\nGetting details for issue #$ISSUE"
  ISSUE_DATA=$(jq -r ".[] | select(.number == $ISSUE)" all_issues.json)
  
  if [[ -z "$ISSUE_DATA" ]]; then
    echo "Could not find data for issue #$ISSUE, skipping"
    continue
  fi
  
  TITLE=$(echo "$ISSUE_DATA" | jq -r '.title')
  BODY=$(echo "$ISSUE_DATA" | jq -r '.body')
  LABELS=$(echo "$ISSUE_DATA" | jq -r '.labels[] | .name' | tr '\n' ',')
  TYPE=$(determine_type "$TITLE")
  
  echo "Issue #$ISSUE: $TITLE"
  echo "Current labels: $LABELS"
  
  # Update issue description and labels
  update_description $ISSUE "$TITLE" "$BODY" "$TYPE"
  update_labels $ISSUE "$TITLE" "$LABELS"
  
  # Increment counters
  COUNTER=$((COUNTER+1))
  BATCH_COUNTER=$((BATCH_COUNTER+1))
  
  # Check if we need to take a break between batches
  if [[ $BATCH_COUNTER -ge $BATCH_SIZE ]]; then
    echo -e "\nCompleted batch of $BATCH_SIZE issues"
    echo "Processed $COUNTER issues in total"
    echo "Taking a short break to avoid rate limiting..."
    sleep 5
    BATCH_COUNTER=0
    
    # Ask if the user wants to continue
    read -p "Continue to next batch? (y/n): " CONTINUE
    if [[ "$CONTINUE" != "y" ]]; then
      echo "Stopping after processing $COUNTER issues"
      break
    fi
  fi
done

echo -e "\nIssue update process complete."
echo "Total issues processed: $COUNTER"
echo "Processed issues saved to $PROCESSED_FILE"
echo "Skipped issues saved to $SKIP_FILE"

# Clean up temporary files
rm -f temp_description.txt
echo "Temporary files cleaned up"