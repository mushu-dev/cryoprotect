#!/bin/bash
# Script to delete all issues in the new repository for a clean start

# Set up variables
REPO="mushu-dev/cryoprotect"
LOG_FILE="issue_deletion_$(date +%Y%m%d_%H%M%S).log"

echo "Starting issue deletion for $REPO"
echo "Deletion log will be saved to $LOG_FILE"
echo "-----------------------------------------------"

# Log start time
echo "Deletion started at $(date)" > "$LOG_FILE"
echo "Repository: $REPO" >> "$LOG_FILE"
echo "-----------------------------------------------" >> "$LOG_FILE"

# Get confirmation from user
read -p "Are you sure you want to delete ALL issues in $REPO? This cannot be undone! (y/N): " CONFIRM

if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
  echo "Operation cancelled by user."
  echo "Operation cancelled by user at $(date)" >> "$LOG_FILE"
  exit 0
fi

# Get all issues from the repository
echo "Fetching issues from $REPO..."
ISSUES=$(gh issue list -s all -L 500 --json number,title -R "$REPO")

# Check if we got any issues
if [ -z "$ISSUES" ]; then
  echo "No issues found in $REPO or error fetching issues"
  echo "No issues found or error fetching issues" >> "$LOG_FILE"
  exit 1
fi

# Count issues
TOTAL_ISSUES=$(echo "$ISSUES" | jq '. | length')
echo "Found $TOTAL_ISSUES issues to delete"
echo "Total issues to delete: $TOTAL_ISSUES" >> "$LOG_FILE"

# Process each issue
DELETED_COUNT=0
FAILED_COUNT=0

echo "$ISSUES" | jq -c '.[]' | while read -r issue; do
  ISSUE_NUMBER=$(echo "$issue" | jq -r '.number')
  ISSUE_TITLE=$(echo "$issue" | jq -r '.title')
  
  echo "Deleting issue #$ISSUE_NUMBER: $ISSUE_TITLE..."
  
  # Attempt to delete the issue (GitHub CLI doesn't support direct deletion, so we use API)
  if gh api -X DELETE "repos/$REPO/issues/$ISSUE_NUMBER" 2>/dev/null; then
    echo "[SUCCESS] Deleted issue #$ISSUE_NUMBER: $ISSUE_TITLE" >> "$LOG_FILE"
    DELETED_COUNT=$((DELETED_COUNT + 1))
  else
    # If direct deletion fails, try to close and lock it
    if gh issue close "$ISSUE_NUMBER" -R "$REPO" && \
       gh issue lock "$ISSUE_NUMBER" -R "$REPO" --reason resolved; then
      echo "[PARTIAL] Closed and locked issue #$ISSUE_NUMBER: $ISSUE_TITLE (Cannot delete)" >> "$LOG_FILE"
      DELETED_COUNT=$((DELETED_COUNT + 1))
    else
      echo "[FAILED] Could not delete or close issue #$ISSUE_NUMBER: $ISSUE_TITLE" >> "$LOG_FILE"
      FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
  fi
  
  echo "Progress: $DELETED_COUNT / $TOTAL_ISSUES"
  
  # Sleep to avoid rate limiting
  sleep 1
done

# Log completion
echo "-----------------------------------------------" >> "$LOG_FILE"
echo "Deletion completed at $(date)" >> "$LOG_FILE"
echo "Total issues: $TOTAL_ISSUES" >> "$LOG_FILE"
echo "Successfully deleted/closed: $DELETED_COUNT" >> "$LOG_FILE"
echo "Failed operations: $FAILED_COUNT" >> "$LOG_FILE"

echo "-----------------------------------------------"
echo "Issue deletion complete!"
echo "Successfully deleted/closed: $DELETED_COUNT issues"
echo "Failed operations: $FAILED_COUNT issues"
echo "Log saved to $LOG_FILE"