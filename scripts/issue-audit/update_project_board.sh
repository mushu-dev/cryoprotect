#!/bin/bash
#
# GitHub Project Board Update Script
#
# This script updates the project board by adding missing issues,
# ensuring proper column placement based on status labels, and
# generating a comprehensive report.
#
# Usage:
#   ./update_project_board.sh [--dry-run] [--project=PROJECT_NUMBER]
#
# Options:
#   --dry-run            Show what would be done without making changes
#   --project=NUMBER     Specify project number (default: looks for first project)
#   --help               Show this help message

set -e

# Default options
DRY_RUN=false
PROJECT_NUMBER=""

# Parse command line arguments
for arg in "$@"; do
  case $arg in
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    --project=*)
      PROJECT_NUMBER="${arg#*=}"
      shift
      ;;
    --help)
      echo "Usage: ./update_project_board.sh [--dry-run] [--project=PROJECT_NUMBER]"
      echo ""
      echo "Options:"
      echo "  --dry-run            Show what would be done without making changes"
      echo "  --project=NUMBER     Specify project number (default: looks for first project)"
      echo "  --help               Show this help message"
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

# Get current date for report filenames
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
REPORT_FILE="output/project_board_update_${TIMESTAMP}.md"
ACTIONS_FILE="output/project_board_actions_${TIMESTAMP}.sh"

# Initialize report file
echo "# GitHub Project Board Update Report" > "$REPORT_FILE"
echo "" >> "$REPORT_FILE"
echo "Generated on: $(date)" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"

# Initialize actions script
echo "#!/bin/bash" > "$ACTIONS_FILE"
echo "#" >> "$ACTIONS_FILE"
echo "# GitHub Project Board Update Actions" >> "$ACTIONS_FILE"
echo "# Generated on: $(date)" >> "$ACTIONS_FILE"
echo "#" >> "$ACTIONS_FILE"
echo "" >> "$ACTIONS_FILE"

# Determine the project ID if not specified
if [[ -z "$PROJECT_NUMBER" ]]; then
    echo "Project number not specified, looking for projects..."
    
    # Get all projects
    PROJECTS=$(gh project list --limit 10 --json number,title,url)
    
    # Check if any projects exist
    if [[ $(echo "$PROJECTS" | jq '. | length') -eq 0 ]]; then
        echo "Error: No projects found."
        exit 1
    fi
    
    # Use the first project
    PROJECT_NUMBER=$(echo "$PROJECTS" | jq -r '.[0].number')
    PROJECT_TITLE=$(echo "$PROJECTS" | jq -r '.[0].title')
    
    echo "Using project #$PROJECT_NUMBER: $PROJECT_TITLE"
else
    # Verify the project exists
    if ! gh project view "$PROJECT_NUMBER" --json title &> /dev/null; then
        echo "Error: Project #$PROJECT_NUMBER not found."
        exit 1
    fi
    
    PROJECT_TITLE=$(gh project view "$PROJECT_NUMBER" --json title | jq -r '.title')
    echo "Using project #$PROJECT_NUMBER: $PROJECT_TITLE"
fi

# Add project info to report
echo "## Project Information" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"
echo "- Project Number: $PROJECT_NUMBER" >> "$REPORT_FILE"
echo "- Project Title: $PROJECT_TITLE" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"

# Get project fields to determine column field name
echo "Getting project field information..."
PROJECT_FIELDS=$(gh project field-list "$PROJECT_NUMBER" --json name,id,dataType)

# Find the status field
STATUS_FIELD_ID=$(echo "$PROJECT_FIELDS" | jq -r '.[] | select(.name == "Status") | .id')
if [[ -z "$STATUS_FIELD_ID" || "$STATUS_FIELD_ID" == "null" ]]; then
    echo "Error: Status field not found in project. Please create a Status field in your project."
    exit 1
fi

echo "Found Status field: $STATUS_FIELD_ID"

# Get status field options
STATUS_OPTIONS=$(gh project field-options "$PROJECT_NUMBER" "$STATUS_FIELD_ID" --json name,id)

# Extract option mappings (name to ID)
echo "Getting status field options..."
TODO_OPTION_ID=$(echo "$STATUS_OPTIONS" | jq -r '.[] | select(.name == "To Do") | .id')
IN_PROGRESS_OPTION_ID=$(echo "$STATUS_OPTIONS" | jq -r '.[] | select(.name == "In Progress") | .id')
REVIEW_OPTION_ID=$(echo "$STATUS_OPTIONS" | jq -r '.[] | select(.name == "Review") | .id')
DONE_OPTION_ID=$(echo "$STATUS_OPTIONS" | jq -r '.[] | select(.name == "Done") | .id')

# Verify we found all required options
if [[ -z "$TODO_OPTION_ID" || "$TODO_OPTION_ID" == "null" ||
      -z "$IN_PROGRESS_OPTION_ID" || "$IN_PROGRESS_OPTION_ID" == "null" ||
      -z "$REVIEW_OPTION_ID" || "$REVIEW_OPTION_ID" == "null" ||
      -z "$DONE_OPTION_ID" || "$DONE_OPTION_ID" == "null" ]]; then
    echo "Error: Not all required status options found (To Do, In Progress, Review, Done)."
    echo "Please ensure your project has these status options."
    exit 1
fi

echo "Found status options:"
echo "  To Do: $TODO_OPTION_ID"
echo "  In Progress: $IN_PROGRESS_OPTION_ID"
echo "  Review: $REVIEW_OPTION_ID"
echo "  Done: $DONE_OPTION_ID"

# Get all open issues
echo "Fetching open issues..."
OPEN_ISSUES=$(gh issue list --state open --json number,title,labels,url --limit 1000)
OPEN_ISSUE_COUNT=$(echo "$OPEN_ISSUES" | jq '. | length')
echo "Found $OPEN_ISSUE_COUNT open issues"

# Get all project items
echo "Fetching project items..."
PROJECT_ITEMS=$(gh project item-list "$PROJECT_NUMBER" --json id,content,fieldValues --limit 1000)
PROJECT_ITEM_COUNT=$(echo "$PROJECT_ITEMS" | jq '. | length')
echo "Found $PROJECT_ITEM_COUNT items in project"

# Extract issue numbers already in the project
declare -A PROJECT_ISSUE_MAP
while read -r item_id issue_number status_id; do
    if [[ -n "$issue_number" && "$issue_number" != "null" ]]; then
        PROJECT_ISSUE_MAP["$issue_number"]="$item_id|$status_id"
    fi
done < <(echo "$PROJECT_ITEMS" | jq -r '.[] | 
    .id + " " + 
    (.content.number | tostring) + " " + 
    (.fieldValues[] | select(.field.id == "'"$STATUS_FIELD_ID"'") | .optionId)')

echo "Found ${#PROJECT_ISSUE_MAP[@]} issues already in project"

# Add issues that are missing from the project
echo "Identifying issues missing from project..."

ISSUES_TO_ADD=()
ISSUES_TO_UPDATE=()

while read -r issue_number title labels; do
    # Skip if not valid JSON
    if [[ -z "$issue_number" || "$issue_number" == "null" ]]; then
        continue
    fi
    
    # Skip if already in project
    if [[ -n "${PROJECT_ISSUE_MAP[$issue_number]}" ]]; then
        # Check if status needs updating
        IFS='|' read -r item_id current_status_id <<< "${PROJECT_ISSUE_MAP[$issue_number]}"
        
        # Determine correct status based on issue labels
        CORRECT_STATUS_ID="$TODO_OPTION_ID" # Default to "To Do"
        
        # Parse labels
        if echo "$labels" | grep -q '"name":"status:in-progress"'; then
            CORRECT_STATUS_ID="$IN_PROGRESS_OPTION_ID"
        elif echo "$labels" | grep -q '"name":"status:review"'; then
            CORRECT_STATUS_ID="$REVIEW_OPTION_ID"
        elif echo "$labels" | grep -q '"name":"status:completed"'; then
            CORRECT_STATUS_ID="$DONE_OPTION_ID"
        fi
        
        # Check if status needs updating
        if [[ "$current_status_id" != "$CORRECT_STATUS_ID" ]]; then
            ISSUES_TO_UPDATE+=("$issue_number|$item_id|$CORRECT_STATUS_ID")
        fi
        
        continue
    fi
    
    # Determine status for new issue
    STATUS_ID="$TODO_OPTION_ID" # Default to "To Do"
    
    # Parse labels to set appropriate status
    if echo "$labels" | grep -q '"name":"status:in-progress"'; then
        STATUS_ID="$IN_PROGRESS_OPTION_ID"
    elif echo "$labels" | grep -q '"name":"status:review"'; then
        STATUS_ID="$REVIEW_OPTION_ID"
    elif echo "$labels" | grep -q '"name":"status:completed"'; then
        STATUS_ID="$DONE_OPTION_ID"
    fi
    
    ISSUES_TO_ADD+=("$issue_number|$title|$STATUS_ID")
    
done < <(echo "$OPEN_ISSUES" | jq -r '.[] | "\(.number) \(.title) \(.labels)"')

# Report findings
echo "" >> "$REPORT_FILE"
echo "## Summary" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"
echo "- Total open issues: $OPEN_ISSUE_COUNT" >> "$REPORT_FILE"
echo "- Issues already in project: ${#PROJECT_ISSUE_MAP[@]}" >> "$REPORT_FILE"
echo "- Issues to add to project: ${#ISSUES_TO_ADD[@]}" >> "$REPORT_FILE"
echo "- Issues to update status: ${#ISSUES_TO_UPDATE[@]}" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"

# Report issues to add
if [[ ${#ISSUES_TO_ADD[@]} -gt 0 ]]; then
    echo "" >> "$REPORT_FILE"
    echo "## Issues to Add to Project" >> "$REPORT_FILE"
    echo "" >> "$REPORT_FILE"
    
    echo "# Add issues to project" >> "$ACTIONS_FILE"
    echo "" >> "$ACTIONS_FILE"
    
    for issue_info in "${ISSUES_TO_ADD[@]}"; do
        IFS='|' read -r number title status_id <<< "$issue_info"
        
        # Determine status name for display
        STATUS_NAME="To Do"
        if [[ "$status_id" == "$IN_PROGRESS_OPTION_ID" ]]; then
            STATUS_NAME="In Progress"
        elif [[ "$status_id" == "$REVIEW_OPTION_ID" ]]; then
            STATUS_NAME="Review"
        elif [[ "$status_id" == "$DONE_OPTION_ID" ]]; then
            STATUS_NAME="Done"
        fi
        
        echo "- Issue #$number: $title (Status: $STATUS_NAME)" >> "$REPORT_FILE"
        
        # Add command to action script
        echo "echo \"Adding issue #$number to project...\"" >> "$ACTIONS_FILE"
        echo "gh project item-add $PROJECT_NUMBER --id \"#$number\"" >> "$ACTIONS_FILE"
        echo "sleep 1" >> "$ACTIONS_FILE"
        
        # Get the newly created item ID
        echo "ITEM_ID=\$(gh project item-list $PROJECT_NUMBER --json id,content --limit 100 | jq -r '.[] | select(.content.number == $number) | .id')" >> "$ACTIONS_FILE"
        echo "if [[ -n \"\$ITEM_ID\" ]]; then" >> "$ACTIONS_FILE"
        echo "  echo \"Setting status for item \$ITEM_ID...\"" >> "$ACTIONS_FILE"
        echo "  gh project item-edit \$ITEM_ID --field-id $STATUS_FIELD_ID --option-id $status_id" >> "$ACTIONS_FILE"
        echo "else" >> "$ACTIONS_FILE"
        echo "  echo \"Error: Could not find item ID for issue #$number\"" >> "$ACTIONS_FILE"
        echo "fi" >> "$ACTIONS_FILE"
        echo "" >> "$ACTIONS_FILE"
    done
    
    echo "" >> "$REPORT_FILE"
fi

# Report issues to update
if [[ ${#ISSUES_TO_UPDATE[@]} -gt 0 ]]; then
    echo "" >> "$REPORT_FILE"
    echo "## Issues to Update Status" >> "$REPORT_FILE"
    echo "" >> "$REPORT_FILE"
    
    echo "# Update issue status" >> "$ACTIONS_FILE"
    echo "" >> "$ACTIONS_FILE"
    
    for issue_info in "${ISSUES_TO_UPDATE[@]}"; do
        IFS='|' read -r number item_id new_status_id <<< "$issue_info"
        
        # Get issue title from open issues
        TITLE=$(echo "$OPEN_ISSUES" | jq -r ".[] | select(.number == $number) | .title")
        
        # Determine status name for display
        NEW_STATUS_NAME="To Do"
        if [[ "$new_status_id" == "$IN_PROGRESS_OPTION_ID" ]]; then
            NEW_STATUS_NAME="In Progress"
        elif [[ "$new_status_id" == "$REVIEW_OPTION_ID" ]]; then
            NEW_STATUS_NAME="Review"
        elif [[ "$new_status_id" == "$DONE_OPTION_ID" ]]; then
            NEW_STATUS_NAME="Done"
        fi
        
        echo "- Issue #$number: $TITLE (New Status: $NEW_STATUS_NAME)" >> "$REPORT_FILE"
        
        # Add command to action script
        echo "echo \"Updating status for issue #$number...\"" >> "$ACTIONS_FILE"
        echo "gh project item-edit $item_id --field-id $STATUS_FIELD_ID --option-id $new_status_id" >> "$ACTIONS_FILE"
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