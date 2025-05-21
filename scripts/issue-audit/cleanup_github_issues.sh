#!/bin/bash
# GitHub Issue Cleanup Script
# This script processes GitHub issues based on audit results

# Colors for terminal output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

INPUT_DIR="/home/mushu/Projects/CryoProtect/scripts/issue-audit/output"
LOGS_DIR="/home/mushu/Projects/CryoProtect/scripts/issue-audit/logs"
mkdir -p "$LOGS_DIR"

echo -e "${BLUE}CryoProtect GitHub Issue Cleanup${NC}"
echo "=================================="

# Check if gh CLI is installed
if ! command -v gh &> /dev/null; then
    echo -e "${RED}Error: GitHub CLI (gh) is not installed or not in PATH${NC}"
    echo "Please install it from: https://cli.github.com/"
    exit 1
fi

# Check if jq is installed
if ! command -v jq &> /dev/null; then
    echo -e "${RED}Error: jq is not installed or not in PATH${NC}"
    echo "Please install it with: sudo apt-get install jq"
    exit 1
fi

# Check if input files exist
if [ ! -f "$INPUT_DIR/all_issues.json" ]; then
    echo -e "${RED}Error: Input files not found. Please run audit_github_issues.sh first.${NC}"
    exit 1
fi

# Load issue data
all_issues=$(cat "$INPUT_DIR/all_issues.json")
keep_issues=$(cat "$INPUT_DIR/keep_issues.json")
improve_issues=$(cat "$INPUT_DIR/improve_issues.json")
review_issues=$(cat "$INPUT_DIR/review_issues.json")

# Available actions for each issue category
echo -e "${BLUE}Available Actions:${NC}"
echo "1. Archive issues (add 'archived' label and close)"
echo "2. Standardize issues (add template structure to issues)"
echo "3. Fix issues (apply milestone and label corrections)"
echo "4. Close empty issues (close issues with no meaningful content)"
echo "5. Assign issues to milestones"
echo "6. Batch label issues"
echo "7. Generate issue report"
echo "8. Apply standard template to issues"

# Ask for action to perform
read -p "Enter action number (1-8): " action

# Function to process a list of issues
process_issues() {
    local issues=$1
    local process_type=$2
    local issue_count=$(echo "$issues" | jq '. | length')
    
    echo -e "${BLUE}Processing $issue_count issues for: $process_type${NC}"
    
    # Ask for confirmation
    read -p "Are you sure you want to process these issues? (y/n): " confirm
    if [ "$confirm" != "y" ]; then
        echo "Operation cancelled."
        return
    fi
    
    local log_file="$LOGS_DIR/${process_type}_$(date +%Y%m%d_%H%M%S).log"
    echo "Starting $process_type at $(date)" > "$log_file"
    
    # Process each issue
    for i in $(seq 0 $((issue_count-1))); do
        issue_number=$(echo "$issues" | jq -r ".[$i]")
        issue_data=$(echo "$all_issues" | jq ".[] | select(.number == $issue_number)")
        title=$(echo "$issue_data" | jq -r '.title')
        
        echo -e "${YELLOW}Processing issue #$issue_number: $title${NC}"
        echo "Processing issue #$issue_number: $title" >> "$log_file"
        
        case $process_type in
            "archive")
                # Add archived label and close
                gh issue edit "$issue_number" --add-label "archived" --state closed
                gh issue comment "$issue_number" --body "This issue has been archived during repository cleanup."
                ;;
            "standardize")
                # Add standard structure if missing
                body=$(echo "$issue_data" | jq -r '.body')
                if [ -z "$body" ] || [ "$body" == "null" ]; then
                    # Create standard template
                    new_body="## Description
                    
*[Describe the issue]*

## Requirements

- [ ] *[Requirement 1]*
- [ ] *[Requirement 2]*

## Technical Context

*[Provide technical context or background information]*

## Acceptance Criteria

- [ ] *[Criterion 1]*
- [ ] *[Criterion 2]*

---
Note: This issue was standardized during repository cleanup."
                    
                    gh issue edit "$issue_number" --body "$new_body"
                fi
                ;;
            "fix")
                # Add component label if missing
                labels=$(echo "$issue_data" | jq '.labels[].name')
                has_component=$(echo "$labels" | grep -c "component:")
                
                if [ "$has_component" -eq 0 ]; then
                    # Determine component based on title
                    title_lower=$(echo "$title" | tr '[:upper:]' '[:lower:]')
                    
                    if [[ "$title_lower" == *"database"* ]] || [[ "$title_lower" == *"db"* ]]; then
                        gh issue edit "$issue_number" --add-label "component:database"
                    elif [[ "$title_lower" == *"api"* ]]; then
                        gh issue edit "$issue_number" --add-label "component:api"
                    elif [[ "$title_lower" == *"ui"* ]] || [[ "$title_lower" == *"interface"* ]]; then
                        gh issue edit "$issue_number" --add-label "component:ui"
                    elif [[ "$title_lower" == *"auth"* ]]; then
                        gh issue edit "$issue_number" --add-label "component:auth"
                    elif [[ "$title_lower" == *"test"* ]]; then
                        gh issue edit "$issue_number" --add-label "component:testing"
                    elif [[ "$title_lower" == *"doc"* ]]; then
                        gh issue edit "$issue_number" --add-label "component:docs"
                    else
                        gh issue edit "$issue_number" --add-label "component:other"
                    fi
                fi
                
                # Add priority label if missing
                has_priority=$(echo "$labels" | grep -c "priority:")
                if [ "$has_priority" -eq 0 ]; then
                    gh issue edit "$issue_number" --add-label "priority:medium"
                fi
                
                # Add type label if missing
                has_type=$(echo "$labels" | grep -c "type:")
                if [ "$has_type" -eq 0 ]; then
                    gh issue edit "$issue_number" --add-label "type:task"
                fi
                ;;
            "close")
                # Close empty issues
                body=$(echo "$issue_data" | jq -r '.body')
                body_length=$(echo "$body" | wc -c)
                
                if [ -z "$body" ] || [ "$body" == "null" ] || [ "$body_length" -lt 30 ]; then
                    gh issue edit "$issue_number" --add-label "empty" --state closed
                    gh issue comment "$issue_number" --body "This issue has been closed during repository cleanup as it contained insufficient information."
                fi
                ;;
            "milestones")
                # Prompt for milestone
                read -p "Enter milestone title for issue #$issue_number: " milestone
                if [ -n "$milestone" ]; then
                    gh issue edit "$issue_number" --milestone "$milestone"
                fi
                ;;
            "label")
                # Prompt for label
                read -p "Enter comma-separated labels to add to issue #$issue_number: " labels
                if [ -n "$labels" ]; then
                    IFS=',' read -ra LABEL_ARRAY <<< "$labels"
                    for label in "${LABEL_ARRAY[@]}"; do
                        gh issue edit "$issue_number" --add-label "$(echo "$label" | xargs)"
                    done
                fi
                ;;
            "template")
                # Apply full template to issue
                read -p "Apply full template to issue #$issue_number? This will replace existing content (y/n): " confirm_template
                if [ "$confirm_template" == "y" ]; then
                    new_body="## Description

*[Describe the issue or feature request in detail]*

## Requirements

- [ ] *[Specific requirement 1]*
- [ ] *[Specific requirement 2]*
- [ ] *[Specific requirement 3]*

## Technical Context

*[Provide information about the existing architecture, systems, or code that relates to this issue]*

Key files:
- *[file1.py (purpose)]*
- *[file2.py (purpose)]*

## Constraints

- *[Constraint 1]*
- *[Constraint 2]*

## Acceptance Criteria

- [ ] *[Specific criterion 1]*
- [ ] *[Specific criterion 2]*
- [ ] *[Specific criterion 3]*

## Additional Notes

*[Any extra information, references, or context]*

---
*This issue was standardized using the CryoProtect issue template.*"
                    
                    gh issue edit "$issue_number" --body "$new_body"
                fi
                ;;
        esac
        
        echo "Completed processing issue #$issue_number" >> "$log_file"
        sleep 1  # Add delay to avoid API rate limiting
    done
    
    echo -e "${GREEN}Completed $process_type process. Log saved to $log_file${NC}"
}

# Execute selected action
case $action in
    1)
        # Archive issues (review category)
        process_issues "$review_issues" "archive"
        ;;
    2)
        # Standardize issues (improve category)
        process_issues "$improve_issues" "standardize"
        ;;
    3)
        # Fix issues (all categories)
        echo -e "${BLUE}Select issues to fix:${NC}"
        echo "1. High quality issues"
        echo "2. Medium quality issues"
        echo "3. Low quality issues"
        read -p "Enter selection (1-3): " fix_selection
        
        case $fix_selection in
            1) process_issues "$keep_issues" "fix" ;;
            2) process_issues "$improve_issues" "fix" ;;
            3) process_issues "$review_issues" "fix" ;;
            *) echo "Invalid selection" ;;
        esac
        ;;
    4)
        # Close empty issues
        blank_issues=$(cat "$INPUT_DIR/blank_issues.json" | jq '[.[].number]')
        process_issues "$blank_issues" "close"
        ;;
    5)
        # Assign issues to milestones
        echo -e "${BLUE}Select issues to assign milestones:${NC}"
        echo "1. High quality issues"
        echo "2. Medium quality issues"
        echo "3. Low quality issues"
        read -p "Enter selection (1-3): " milestone_selection
        
        case $milestone_selection in
            1) process_issues "$keep_issues" "milestones" ;;
            2) process_issues "$improve_issues" "milestones" ;;
            3) process_issues "$review_issues" "milestones" ;;
            *) echo "Invalid selection" ;;
        esac
        ;;
    6)
        # Batch label issues
        echo -e "${BLUE}Select issues to batch label:${NC}"
        echo "1. High quality issues"
        echo "2. Medium quality issues"
        echo "3. Low quality issues"
        read -p "Enter selection (1-3): " label_selection
        
        case $label_selection in
            1) process_issues "$keep_issues" "label" ;;
            2) process_issues "$improve_issues" "label" ;;
            3) process_issues "$review_issues" "label" ;;
            *) echo "Invalid selection" ;;
        esac
        ;;
    7)
        # Generate issue report
        report_file="$LOGS_DIR/issue_report_$(date +%Y%m%d_%H%M%S).md"
        echo "# CryoProtect GitHub Issue Report" > "$report_file"
        echo "" >> "$report_file"
        echo "Generated on: $(date)" >> "$report_file"
        echo "" >> "$report_file"
        
        echo "## Issue Categories" >> "$report_file"
        echo "" >> "$report_file"
        echo "| Category | Count |" >> "$report_file"
        echo "|----------|-------|" >> "$report_file"
        echo "| High Quality | $(echo "$keep_issues" | jq '. | length') |" >> "$report_file"
        echo "| Medium Quality | $(echo "$improve_issues" | jq '. | length') |" >> "$report_file"
        echo "| Low Quality | $(echo "$review_issues" | jq '. | length') |" >> "$report_file"
        echo "" >> "$report_file"
        
        # Add high quality issues
        echo "## High Quality Issues" >> "$report_file"
        echo "" >> "$report_file"
        echo "| Number | Title |" >> "$report_file"
        echo "|--------|-------|" >> "$report_file"
        
        high_quality_count=$(echo "$keep_issues" | jq '. | length')
        for i in $(seq 0 $((high_quality_count-1))); do
            issue_number=$(echo "$keep_issues" | jq -r ".[$i]")
            title=$(echo "$all_issues" | jq -r ".[] | select(.number == $issue_number) | .title")
            echo "| #$issue_number | $title |" >> "$report_file"
        done
        
        echo "" >> "$report_file"
        
        # Add medium quality issues
        echo "## Medium Quality Issues" >> "$report_file"
        echo "" >> "$report_file"
        echo "| Number | Title |" >> "$report_file"
        echo "|--------|-------|" >> "$report_file"
        
        medium_quality_count=$(echo "$improve_issues" | jq '. | length')
        for i in $(seq 0 $((medium_quality_count-1))); do
            issue_number=$(echo "$improve_issues" | jq -r ".[$i]")
            title=$(echo "$all_issues" | jq -r ".[] | select(.number == $issue_number) | .title")
            echo "| #$issue_number | $title |" >> "$report_file"
        done
        
        echo "" >> "$report_file"
        
        # Add low quality issues
        echo "## Low Quality Issues" >> "$report_file"
        echo "" >> "$report_file"
        echo "| Number | Title |" >> "$report_file"
        echo "|--------|-------|" >> "$report_file"
        
        low_quality_count=$(echo "$review_issues" | jq '. | length')
        for i in $(seq 0 $((low_quality_count-1))); do
            issue_number=$(echo "$review_issues" | jq -r ".[$i]")
            title=$(echo "$all_issues" | jq -r ".[] | select(.number == $issue_number) | .title")
            echo "| #$issue_number | $title |" >> "$report_file"
        done
        
        echo -e "${GREEN}Report generated at: $report_file${NC}"
        ;;
    8)
        # Apply standard template
        echo -e "${BLUE}Select issues to apply template:${NC}"
        echo "1. High quality issues"
        echo "2. Medium quality issues"
        echo "3. Low quality issues"
        read -p "Enter selection (1-3): " template_selection
        
        case $template_selection in
            1) process_issues "$keep_issues" "template" ;;
            2) process_issues "$improve_issues" "template" ;;
            3) process_issues "$review_issues" "template" ;;
            *) echo "Invalid selection" ;;
        esac
        ;;
    *)
        echo -e "${RED}Invalid action selected${NC}"
        ;;
esac

echo -e "${GREEN}Cleanup operations completed!${NC}"