#!/bin/bash
#
# GitHub Stale Issue Cleanup Script
#
# This script identifies and manages stale issues in the GitHub repository.
# It adds the 'status:stale' label to issues that haven't been updated in 30+ days,
# and optionally closes issues that have been stale for 60+ days.
#
# Usage:
#   ./cleanup_stale_issues.sh [--dry-run] [--close-after=DAYS]
#
# Options:
#   --dry-run         Show what would be done without making changes
#   --close-after=N   Close issues that have been stale for N days (default: 60)
#   --help            Show this help message

set -e

# Default options
DRY_RUN=false
STALE_THRESHOLD=30  # days
CLOSE_THRESHOLD=60  # days
SHOULD_CLOSE=false

# Parse command line arguments
for arg in "$@"; do
  case $arg in
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    --close-after=*)
      CLOSE_THRESHOLD="${arg#*=}"
      SHOULD_CLOSE=true
      shift
      ;;
    --help)
      echo "Usage: ./cleanup_stale_issues.sh [--dry-run] [--close-after=DAYS]"
      echo ""
      echo "Options:"
      echo "  --dry-run         Show what would be done without making changes"
      echo "  --close-after=N   Close issues that have been stale for N days (default: 60)"
      echo "  --help            Show this help message"
      exit 0
      ;;
    *)
      # Unknown option
      echo "Unknown option: $arg"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

# Check if GitHub CLI is installed
if ! command -v gh &> /dev/null; then
    echo "Error: GitHub CLI (gh) is not installed."
    echo "Please install it from https://cli.github.com/"
    exit 1
fi

# Check if GitHub CLI is authenticated
if ! gh auth status &> /dev/null; then
    echo "Error: GitHub CLI is not authenticated."
    echo "Please run 'gh auth login' to authenticate."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p output

# Get current date
CURRENT_DATE=$(date +%s)

# Get all open issues in JSON format
echo "Fetching open issues..."
ISSUES=$(gh issue list --state open --json number,title,updatedAt,labels,url -L 1000)

# Count issues
ISSUE_COUNT=$(echo "$ISSUES" | jq '. | length')
echo "Found $ISSUE_COUNT open issues"

# Initialize counters
STALE_COUNT=0
ALREADY_STALE_COUNT=0
TO_CLOSE_COUNT=0

# Create arrays to store actions
declare -a STALE_ISSUES
declare -a CLOSE_ISSUES

# Process each issue
echo "Analyzing issues for staleness..."
while read -r number title updated_at labels url; do
    # Skip if not valid JSON
    if [[ -z "$number" || "$number" == "null" ]]; then
        continue
    fi
    
    # Check if already has stale label
    if echo "$labels" | grep -q "status:stale"; then
        ALREADY_STALE_COUNT=$((ALREADY_STALE_COUNT + 1))
        
        # Get the date the stale label was added
        # This requires an API call to get issue events
        if [[ "$SHOULD_CLOSE" == true ]]; then
            echo "Checking stale duration for issue #$number..."
            
            # Get the label events for this issue
            EVENTS=$(gh api repos/:owner/:repo/issues/$number/events)
            
            # Find when the stale label was added
            STALE_DATE=$(echo "$EVENTS" | jq -r '.[] | select(.event=="labeled" and .label.name=="status:stale") | .created_at' | sort | head -n1)
            
            if [[ -n "$STALE_DATE" && "$STALE_DATE" != "null" ]]; then
                # Convert to timestamp
                STALE_TIMESTAMP=$(date -d "$STALE_DATE" +%s)
                DAYS_STALE=$(( (CURRENT_DATE - STALE_TIMESTAMP) / 86400 ))
                
                if [[ $DAYS_STALE -ge $CLOSE_THRESHOLD ]]; then
                    echo "  Issue #$number has been stale for $DAYS_STALE days, marking for closure"
                    CLOSE_ISSUES+=("$number|$title|$url")
                    TO_CLOSE_COUNT=$((TO_CLOSE_COUNT + 1))
                fi
            fi
        fi
        
        continue
    fi
    
    # Convert updated_at to timestamp
    UPDATED_TIMESTAMP=$(date -d "$updated_at" +%s)
    
    # Calculate days since last update
    DAYS_SINCE_UPDATE=$(( (CURRENT_DATE - UPDATED_TIMESTAMP) / 86400 ))
    
    if [[ $DAYS_SINCE_UPDATE -ge $STALE_THRESHOLD ]]; then
        echo "  Issue #$number is stale ($DAYS_SINCE_UPDATE days without update)"
        STALE_ISSUES+=("$number|$title|$url")
        STALE_COUNT=$((STALE_COUNT + 1))
    fi
done < <(echo "$ISSUES" | jq -r '.[] | "\(.number) \(.title) \(.updatedAt) \(.labels) \(.url)"')

# Report findings
echo ""
echo "Summary:"
echo "  Total open issues: $ISSUE_COUNT"
echo "  Already marked as stale: $ALREADY_STALE_COUNT"
echo "  Newly identified stale issues: $STALE_COUNT"
if [[ "$SHOULD_CLOSE" == true ]]; then
    echo "  Issues to close (stale for $CLOSE_THRESHOLD+ days): $TO_CLOSE_COUNT"
fi

# Generate report file
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
REPORT_FILE="output/stale_issues_report_${TIMESTAMP}.md"
ACTIONS_FILE="output/stale_issues_actions_${TIMESTAMP}.sh"

echo "# Stale Issues Report" > "$REPORT_FILE"
echo "" >> "$REPORT_FILE"
echo "Generated on: $(date)" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"
echo "## Summary" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"
echo "- Total open issues: $ISSUE_COUNT" >> "$REPORT_FILE"
echo "- Already marked as stale: $ALREADY_STALE_COUNT" >> "$REPORT_FILE"
echo "- Newly identified stale issues: $STALE_COUNT" >> "$REPORT_FILE"
if [[ "$SHOULD_CLOSE" == true ]]; then
    echo "- Issues to close (stale for $CLOSE_THRESHOLD+ days): $TO_CLOSE_COUNT" >> "$REPORT_FILE"
fi
echo "" >> "$REPORT_FILE"

# Create action script header
echo "#!/bin/bash" > "$ACTIONS_FILE"
echo "#" >> "$ACTIONS_FILE"
echo "# Stale Issues Actions" >> "$ACTIONS_FILE"
echo "# Generated on: $(date)" >> "$ACTIONS_FILE"
echo "#" >> "$ACTIONS_FILE"
echo "" >> "$ACTIONS_FILE"

# Add newly stale issues to report
if [[ $STALE_COUNT -gt 0 ]]; then
    echo "## Newly Identified Stale Issues" >> "$REPORT_FILE"
    echo "" >> "$REPORT_FILE"
    
    echo "# Mark issues as stale" >> "$ACTIONS_FILE"
    
    for issue in "${STALE_ISSUES[@]}"; do
        IFS='|' read -r number title url <<< "$issue"
        echo "- #$number: [$title]($url) - No updates in $DAYS_SINCE_UPDATE days" >> "$REPORT_FILE"
        
        # Add command to action script
        echo "echo \"Marking issue #$number as stale...\"" >> "$ACTIONS_FILE"
        echo "gh issue edit $number --add-label \"status:stale\"" >> "$ACTIONS_FILE"
        echo "gh issue comment $number --body \"This issue has been automatically marked as stale because it has not had recent activity. It will be closed if no further activity occurs within $CLOSE_THRESHOLD days. Thank you for your contributions.\"" >> "$ACTIONS_FILE"
        echo "" >> "$ACTIONS_FILE"
    done
    echo "" >> "$REPORT_FILE"
fi

# Add issues to close to report
if [[ $TO_CLOSE_COUNT -gt 0 ]]; then
    echo "## Issues to Close" >> "$REPORT_FILE"
    echo "" >> "$REPORT_FILE"
    
    echo "# Close stale issues" >> "$ACTIONS_FILE"
    
    for issue in "${CLOSE_ISSUES[@]}"; do
        IFS='|' read -r number title url <<< "$issue"
        echo "- #$number: [$title]($url) - Stale for $DAYS_STALE days" >> "$REPORT_FILE"
        
        # Add command to action script
        echo "echo \"Closing issue #$number...\"" >> "$ACTIONS_FILE"
        echo "gh issue comment $number --body \"This issue has been automatically closed because it has been marked as stale for $DAYS_STALE days with no activity. If this issue is still relevant, please reopen it or create a new issue with updated information.\"" >> "$ACTIONS_FILE"
        echo "gh issue close $number --reason not_planned" >> "$ACTIONS_FILE"
        echo "" >> "$ACTIONS_FILE"
    done
    echo "" >> "$REPORT_FILE"
fi

echo "chmod +x \"$ACTIONS_FILE\"" >> "$ACTIONS_FILE"

# Make the actions file executable
chmod +x "$ACTIONS_FILE"

echo ""
echo "Report generated: $REPORT_FILE"
echo "Actions script generated: $ACTIONS_FILE"

# Execute actions if not in dry-run mode
if [[ "$DRY_RUN" == false ]]; then
    echo ""
    echo "Do you want to execute the generated actions? [y/N]"
    read -r EXECUTE
    
    if [[ "$EXECUTE" =~ ^[Yy]$ ]]; then
        echo "Executing actions..."
        bash "$ACTIONS_FILE"
        echo "Actions completed."
    else
        echo "Actions not executed. You can run them later with:"
        echo "  $ACTIONS_FILE"
    fi
else
    echo ""
    echo "Dry run complete. No changes were made."
    echo "To execute the actions, run:"
    echo "  $ACTIONS_FILE"
fi