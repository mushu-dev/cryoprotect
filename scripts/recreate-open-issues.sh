#!/bin/bash
# This script recreates ONLY OPEN issues from blueprint-house/cryoprotect to mushu-dev/cryoprotect

# Set up variables
SOURCE_REPO="blueprint-house/cryoprotect"
DEST_REPO="mushu-dev/cryoprotect"
LOG_FILE="open_issue_recreation_$(date +%Y%m%d_%H%M%S).log"
ISSUE_COUNT=0
SUCCESS_COUNT=0
FAIL_COUNT=0

echo "Starting OPEN issue recreation from $SOURCE_REPO to $DEST_REPO"
echo "Recreation log will be saved to $LOG_FILE"
echo "-----------------------------------------------"

# Log start time
echo "Recreation started at $(date)" > "$LOG_FILE"
echo "Source: $SOURCE_REPO" >> "$LOG_FILE"
echo "Destination: $DEST_REPO" >> "$LOG_FILE"
echo "-----------------------------------------------" >> "$LOG_FILE"

# Get ONLY OPEN issues from source repository
echo "Fetching OPEN issues from $SOURCE_REPO..."
ISSUES=$(gh issue list -s open -L 500 --json number,title,body,state,labels,assignees,milestone -R "$SOURCE_REPO")

# Check if we got any issues
if [ -z "$ISSUES" ]; then
  echo "No open issues found in $SOURCE_REPO or error fetching issues"
  echo "No open issues found or error fetching issues" >> "$LOG_FILE"
  exit 1
fi

# Count issues
TOTAL_ISSUES=$(echo "$ISSUES" | jq '. | length')
echo "Found $TOTAL_ISSUES open issues to recreate"
echo "Total open issues to recreate: $TOTAL_ISSUES" >> "$LOG_FILE"

# Create issue mapping file
MAPPING_FILE="open_issue_mapping_$(date +%Y%m%d_%H%M%S).txt"
echo "Issue mapping will be saved to $MAPPING_FILE"
echo "Mapping of Open Issues from $SOURCE_REPO to $DEST_REPO" > "$MAPPING_FILE"
echo "Created at $(date)" >> "$MAPPING_FILE"
echo "=======================================" >> "$MAPPING_FILE"
echo "" >> "$MAPPING_FILE"

# Process each issue
echo "$ISSUES" | jq -c '.[]' | while read -r issue; do
  ISSUE_NUMBER=$(echo "$issue" | jq -r '.number')
  ISSUE_TITLE=$(echo "$issue" | jq -r '.title')
  ISSUE_BODY=$(echo "$issue" | jq -r '.body')
  ISSUE_STATE=$(echo "$issue" | jq -r '.state')
  
  # Create a mapping file entry
  echo "Original Issue #$ISSUE_NUMBER in $SOURCE_REPO" >> "$MAPPING_FILE"
  
  # Handle labels
  LABELS=""
  LABEL_JSON=$(echo "$issue" | jq -r '.labels')
  if [ "$LABEL_JSON" != "null" ] && [ "$LABEL_JSON" != "[]" ]; then
    # Extract label names and create a comma-separated list
    LABELS=$(echo "$LABEL_JSON" | jq -r '.[].name' | tr '\n' ',' | sed 's/,$//')
  fi
  
  # Add reference to original issue in body
  UPDATED_BODY="${ISSUE_BODY}

---
*Recreated from original issue blueprint-house/cryoprotect#${ISSUE_NUMBER}*"
  
  echo "Recreating issue #$ISSUE_NUMBER: $ISSUE_TITLE ($ISSUE_STATE)..."
  
  # Create a temporary file for the body
  TEMP_BODY_FILE=$(mktemp)
  echo "$UPDATED_BODY" > "$TEMP_BODY_FILE"
  
  # Process labels individually to avoid errors with non-existent labels
  NEW_ISSUE_URL=$(gh issue create --title "$ISSUE_TITLE" --body-file "$TEMP_BODY_FILE" -R "$DEST_REPO")
  
  # Extract the new issue number from the URL
  NEW_ISSUE_NUMBER=$(echo "$NEW_ISSUE_URL" | grep -oE '[0-9]+$')
  
  # If issue was created successfully and we have labels
  if [ ! -z "$NEW_ISSUE_URL" ] && [ ! -z "$LABELS" ]; then
    # Add each label individually
    for LABEL in $(echo "$LABELS" | tr ',' '\n'); do
      echo "  Adding label '$LABEL' to issue #$NEW_ISSUE_NUMBER..."
      gh issue edit "$NEW_ISSUE_NUMBER" --add-label "$LABEL" -R "$DEST_REPO" || echo "  Failed to add label '$LABEL'"
    done
  fi
  
  # Add to mapping file
  echo "  → New Issue #$NEW_ISSUE_NUMBER in $DEST_REPO" >> "$MAPPING_FILE"
  echo "  → $NEW_ISSUE_URL" >> "$MAPPING_FILE"
  echo "" >> "$MAPPING_FILE"
  
  # Clean up temp file
  rm "$TEMP_BODY_FILE"
  
  if [ ! -z "$NEW_ISSUE_URL" ]; then
    SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    echo "[SUCCESS] Recreated issue #$ISSUE_NUMBER as #$NEW_ISSUE_NUMBER: $ISSUE_TITLE" >> "$LOG_FILE"
    echo "  $NEW_ISSUE_URL" >> "$LOG_FILE"
  else
    echo "[FAILED] Could not recreate issue #$ISSUE_NUMBER: $ISSUE_TITLE" >> "$LOG_FILE"
    FAIL_COUNT=$((FAIL_COUNT + 1))
  fi
  
  ISSUE_COUNT=$((ISSUE_COUNT + 1))
  echo "Progress: $ISSUE_COUNT / $TOTAL_ISSUES"
  
  # Sleep to avoid rate limiting
  sleep 2
done

# Log completion
echo "-----------------------------------------------" >> "$LOG_FILE"
echo "Recreation completed at $(date)" >> "$LOG_FILE"
echo "Total open issues: $TOTAL_ISSUES" >> "$LOG_FILE"
echo "Successfully recreated: $SUCCESS_COUNT" >> "$LOG_FILE"
echo "Failed recreations: $FAIL_COUNT" >> "$LOG_FILE"

echo "-----------------------------------------------"
echo "Open issue recreation complete!"
echo "Successfully recreated: $SUCCESS_COUNT issues"
echo "Failed recreations: $FAIL_COUNT issues"
echo "Log saved to $LOG_FILE"
echo "Issue mapping saved to $MAPPING_FILE"