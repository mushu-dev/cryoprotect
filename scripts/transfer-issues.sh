#!/bin/bash
# This script transfers all issues from blueprint-house/cryoprotect to mushu-dev/cryoprotect

# Set up variables
SOURCE_REPO="blueprint-house/cryoprotect"
DEST_REPO="mushu-dev/cryoprotect"
LOG_FILE="issue_transfer_$(date +%Y%m%d_%H%M%S).log"
ISSUE_COUNT=0
SUCCESS_COUNT=0
FAIL_COUNT=0

echo "Starting issue transfer from $SOURCE_REPO to $DEST_REPO"
echo "Transfer log will be saved to $LOG_FILE"
echo "-----------------------------------------------"

# Log start time
echo "Transfer started at $(date)" > "$LOG_FILE"
echo "Source: $SOURCE_REPO" >> "$LOG_FILE"
echo "Destination: $DEST_REPO" >> "$LOG_FILE"
echo "-----------------------------------------------" >> "$LOG_FILE"

# Get all issues from source repository
echo "Fetching issues from $SOURCE_REPO..."
ISSUES=$(gh issue list -s all -L 500 --json number,title,state,labels -R "$SOURCE_REPO")

# Check if we got any issues
if [ -z "$ISSUES" ]; then
  echo "No issues found in $SOURCE_REPO or error fetching issues"
  echo "No issues found or error fetching issues" >> "$LOG_FILE"
  exit 1
fi

# Count issues
TOTAL_ISSUES=$(echo "$ISSUES" | jq '. | length')
echo "Found $TOTAL_ISSUES issues to transfer"
echo "Total issues to transfer: $TOTAL_ISSUES" >> "$LOG_FILE"

# Process each issue
echo "$ISSUES" | jq -c '.[]' | while read -r issue; do
  ISSUE_NUMBER=$(echo "$issue" | jq -r '.number')
  ISSUE_TITLE=$(echo "$issue" | jq -r '.title')
  ISSUE_STATE=$(echo "$issue" | jq -r '.state')
  
  echo "Transferring issue #$ISSUE_NUMBER: $ISSUE_TITLE ($ISSUE_STATE)..."
  
  # Try to transfer the issue
  if gh issue transfer "$ISSUE_NUMBER" "$DEST_REPO" -R "$SOURCE_REPO"; then
    echo "[SUCCESS] Transferred issue #$ISSUE_NUMBER: $ISSUE_TITLE" >> "$LOG_FILE"
    SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
  else
    echo "[FAILED] Could not transfer issue #$ISSUE_NUMBER: $ISSUE_TITLE" >> "$LOG_FILE"
    FAIL_COUNT=$((FAIL_COUNT + 1))
  fi
  
  ISSUE_COUNT=$((ISSUE_COUNT + 1))
  echo "Progress: $ISSUE_COUNT / $TOTAL_ISSUES"
  
  # Sleep to avoid rate limiting
  sleep 3
done

# Log completion
echo "-----------------------------------------------" >> "$LOG_FILE"
echo "Transfer completed at $(date)" >> "$LOG_FILE"
echo "Total issues: $TOTAL_ISSUES" >> "$LOG_FILE"
echo "Successfully transferred: $SUCCESS_COUNT" >> "$LOG_FILE"
echo "Failed transfers: $FAIL_COUNT" >> "$LOG_FILE"

echo "-----------------------------------------------"
echo "Issue transfer complete!"
echo "Successfully transferred: $SUCCESS_COUNT issues"
echo "Failed transfers: $FAIL_COUNT issues"
echo "Log saved to $LOG_FILE"