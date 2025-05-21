#!/bin/bash
# GitHub Repository Standardization Script
# Sets up standard labels, milestones, and issue templates for the CryoProtect project

# Colors for terminal output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

OUTPUT_DIR="/home/mushu/Projects/CryoProtect/scripts/issue-audit/output"
mkdir -p "$OUTPUT_DIR"

echo -e "${BLUE}CryoProtect GitHub Repository Standardization${NC}"
echo "============================================="

# Check if gh CLI is installed
if ! command -v gh &> /dev/null; then
    echo -e "${RED}Error: GitHub CLI (gh) is not installed or not in PATH${NC}"
    echo "Please install it from: https://cli.github.com/"
    exit 1
fi

# Get repository information
repo_info=$(gh repo view --json name,owner,isPrivate,description)
repo_name=$(echo "$repo_info" | jq -r '.name')
repo_owner=$(echo "$repo_info" | jq -r '.owner.login')

echo -e "${BLUE}Repository:${NC} $repo_owner/$repo_name"
echo ""

# Display menu of actions
echo -e "${BLUE}Available Actions:${NC}"
echo "1. Create standard labels"
echo "2. Create project milestones"
echo "3. Create issue template files"
echo "4. Clean up existing labels"
echo "5. Check existing issue templates"
echo "6. Create standard workflow files"
echo "7. Export current project structure to JSON"
echo "8. Delete unnecessary labels"
echo "9. Apply metadata to all issues"
echo "0. Exit"

read -p "Select action (0-9): " action

# Standard label definitions
create_standard_labels() {
    echo -e "${BLUE}Creating standard labels...${NC}"
    
    # Component labels (blue)
    gh label create "component:api" --color "0366d6" --description "API-related issues"
    gh label create "component:database" --color "0366d6" --description "Database-related issues"
    gh label create "component:ui" --color "0366d6" --description "User interface issues"
    gh label create "component:auth" --color "0366d6" --description "Authentication and authorization"
    gh label create "component:docs" --color "0366d6" --description "Documentation-related issues"
    gh label create "component:infrastructure" --color "0366d6" --description "Infrastructure and deployment"
    gh label create "component:testing" --color "0366d6" --description "Testing-related issues"
    
    # Priority labels (red/orange/yellow)
    gh label create "priority:critical" --color "b60205" --description "Critical priority issues that must be fixed ASAP"
    gh label create "priority:high" --color "d93f0b" --description "High priority issues"
    gh label create "priority:medium" --color "fbca04" --description "Medium priority issues"
    gh label create "priority:low" --color "fef2c0" --description "Low priority issues"
    
    # Type labels (purple)
    gh label create "type:bug" --color "6f42c1" --description "Something isn't working as expected"
    gh label create "type:feature" --color "6f42c1" --description "New feature request"
    gh label create "type:enhancement" --color "6f42c1" --description "Enhancement to existing functionality"
    gh label create "type:task" --color "6f42c1" --description "General task or work item"
    gh label create "type:refactor" --color "6f42c1" --description "Code refactoring"
    gh label create "type:optimization" --color "6f42c1" --description "Performance optimization"
    gh label create "type:documentation" --color "6f42c1" --description "Documentation updates"
    gh label create "type:question" --color "6f42c1" --description "Question or request for clarification"
    
    # Status labels (green)
    gh label create "status:ready" --color "0e8a16" --description "Ready for implementation"
    gh label create "status:in-progress" --color "1d76db" --description "Work in progress"
    gh label create "status:blocked" --color "b60205" --description "Blocked by dependency or issue"
    gh label create "status:needs-review" --color "5319e7" --description "Needs review"
    gh label create "status:needs-info" --color "fbca04" --description "Needs more information"
    
    # Other labels
    gh label create "good-first-issue" --color "7057ff" --description "Good for newcomers"
    gh label create "help-wanted" --color "008672" --description "Extra attention is needed"
    gh label create "wontfix" --color "ffffff" --description "This will not be worked on"
    gh label create "duplicate" --color "cccccc" --description "This issue already exists"
    gh label create "invalid" --color "e4e669" --description "This doesn't seem right"
    gh label create "archived" --color "cccccc" --description "Archived issues from cleanup"
    
    echo -e "${GREEN}Standard labels created successfully!${NC}"
}

# Create project milestones based on project phases
create_project_milestones() {
    echo -e "${BLUE}Creating project milestones...${NC}"
    
    # Project Timeline
    current_date=$(date +"%Y-%m-%d")
    phase1_date=$(date -d "+1 month" +"%Y-%m-%d")
    phase2_date=$(date -d "+2 months" +"%Y-%m-%d")
    phase3_date=$(date -d "+3 months" +"%Y-%m-%d")
    phase4_date=$(date -d "+4 months" +"%Y-%m-%d")
    
    # Phase 1: Technical Foundation
    gh api repos/$repo_owner/$repo_name/milestones -X POST -f title="Phase 1: Technical Foundation" \
      -f state="open" \
      -f description="Database architecture, authentication system, and core infrastructure" \
      -f due_on="$phase1_date"
    
    # Phase 2: Feature Completion
    gh api repos/$repo_owner/$repo_name/milestones -X POST -f title="Phase 2: Feature Completion" \
      -f state="open" \
      -f description="API layer completion, core functionality, and user interface" \
      -f due_on="$phase2_date"
    
    # Phase 3: Production Readiness  
    gh api repos/$repo_owner/$repo_name/milestones -X POST -f title="Phase 3: Production Readiness" \
      -f state="open" \
      -f description="Deployment infrastructure, monitoring, security enhancements" \
      -f due_on="$phase3_date"
    
    # Phase 4: Documentation & Knowledge Transfer
    gh api repos/$repo_owner/$repo_name/milestones -X POST -f title="Phase 4: Documentation & Knowledge Transfer" \
      -f state="open" \
      -f description="Documentation, knowledge transfer, and final delivery" \
      -f due_on="$phase4_date"
    
    echo -e "${GREEN}Project milestones created successfully!${NC}"
}

# Create issue templates in .github/ISSUE_TEMPLATE directory
create_issue_templates() {
    echo -e "${BLUE}Creating issue templates...${NC}"
    
    # Ensure directory exists
    mkdir -p .github/ISSUE_TEMPLATE
    
    # Create config.yml for the issue template chooser
    cat > .github/ISSUE_TEMPLATE/config.yml << 'EOL'
blank_issues_enabled: false
contact_links:
  - name: CryoProtect Documentation
    url: https://github.com/username/CryoProtect/tree/master/docs
    about: Check the documentation before opening an issue
EOL
    
    # Feature request template
    cat > .github/ISSUE_TEMPLATE/feature_request.md << 'EOL'
---
name: Feature request
about: Suggest a new feature or enhancement for CryoProtect
title: "[FEATURE] "
labels: type:feature
assignees: ''
---

## Description
<!-- Provide a clear and concise description of the feature you're requesting -->

## Requirements
<!-- List specific requirements or functionality this feature should have -->
- [ ] Requirement 1
- [ ] Requirement 2
- [ ] Requirement 3

## Technical Context
<!-- Provide information about the existing architecture, systems, or code that relates to this feature -->

Key files:
- file1.py (purpose)
- file2.py (purpose)

## Constraints
<!-- List any constraints or limitations to consider -->
- Constraint 1
- Constraint 2

## Acceptance Criteria
<!-- Define clear criteria for when this feature would be considered complete -->
- [ ] Criterion 1
- [ ] Criterion 2
- [ ] Criterion 3

## Additional Notes
<!-- Add any other context, mockups, or examples about the feature request here -->
EOL
    
    # Bug report template
    cat > .github/ISSUE_TEMPLATE/bug_report.md << 'EOL'
---
name: Bug report
about: Report a bug or unexpected behavior in CryoProtect
title: "[BUG] "
labels: type:bug
assignees: ''
---

## Description
<!-- Provide a clear and concise description of the bug -->

## Steps to Reproduce
<!-- Detailed steps to reproduce the behavior -->
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

## Expected Behavior
<!-- A clear and concise description of what you expected to happen -->

## Actual Behavior
<!-- What actually happened (include error messages, screenshots if applicable) -->

## Technical Details
<!-- Information about your environment and relevant configuration -->
- Version: [e.g., commit hash or version number]
- Environment: [e.g., local development, production]
- Browser/OS: [if applicable]

## Possible Solution
<!-- If you have suggestions on how to fix the bug -->

## Additional Context
<!-- Add any other context about the problem here -->
EOL

    # Task template
    cat > .github/ISSUE_TEMPLATE/task.md << 'EOL'
---
name: Task
about: Create a general development task for CryoProtect
title: "[TASK] "
labels: type:task
assignees: ''
---

## Description
<!-- Describe the task in detail -->

## Requirements
<!-- List specific requirements or steps needed to complete this task -->
- [ ] Requirement 1
- [ ] Requirement 2
- [ ] Requirement 3

## Technical Context
<!-- Provide information about the existing architecture, systems, or code that relates to this task -->

Key files:
- file1.py (purpose)
- file2.py (purpose)

## Constraints
<!-- List any constraints or limitations to consider -->
- Constraint 1
- Constraint 2

## Acceptance Criteria
<!-- Define clear criteria for when this task would be considered complete -->
- [ ] Criterion 1
- [ ] Criterion 2
- [ ] Criterion 3

## Additional Notes
<!-- Add any other relevant information or context about the task -->
EOL

    # Validation task template
    cat > .github/ISSUE_TEMPLATE/validation_task.md << 'EOL'
---
name: Validation Task
about: Create a verification/validation task for CryoProtect
title: "[VALIDATION] "
labels: type:task, component:testing
assignees: ''
---

## Description
<!-- Describe what needs to be validated/verified -->

## Validation Steps
<!-- List specific steps to validate the functionality -->
- [ ] Step 1
- [ ] Step 2
- [ ] Step 3

## Expected Results
<!-- Describe the expected outcomes from each validation step -->
- Expected result 1
- Expected result 2
- Expected result 3

## Test Data
<!-- Specify any test data or setup needed -->

## Validation Environment
<!-- Specify which environment(s) should be used for validation -->
- [ ] Local development
- [ ] Staging
- [ ] Production

## Additional Notes
<!-- Add any other relevant information or context -->
EOL

    echo -e "${GREEN}Issue templates created successfully!${NC}"
}

# Clean up existing labels
clean_up_labels() {
    echo -e "${BLUE}Retrieving existing labels...${NC}"
    
    # Get existing labels
    existing_labels=$(gh api repos/$repo_owner/$repo_name/labels --paginate | jq -c '.[] | {name, color, description}')
    echo "$existing_labels" > "$OUTPUT_DIR/existing_labels.json"
    
    label_count=$(echo "$existing_labels" | wc -l)
    echo -e "${BLUE}Found ${label_count} labels in the repository${NC}"
    
    # Display labels
    echo -e "${YELLOW}Existing labels:${NC}"
    echo "$existing_labels" | jq -r '.name'
    
    echo -e "${BLUE}Labels saved to $OUTPUT_DIR/existing_labels.json${NC}"
    echo -e "${YELLOW}Review this file and run the 'Delete unnecessary labels' action to clean up${NC}"
}

# Check existing issue templates
check_issue_templates() {
    echo -e "${BLUE}Checking existing issue templates...${NC}"
    
    # Check if .github/ISSUE_TEMPLATE directory exists
    if [ -d ".github/ISSUE_TEMPLATE" ]; then
        echo -e "${GREEN}Issue template directory exists${NC}"
        
        # List existing templates
        echo -e "${BLUE}Existing templates:${NC}"
        ls -la .github/ISSUE_TEMPLATE/
    else
        echo -e "${YELLOW}No issue template directory found${NC}"
        echo "Run 'Create issue template files' action to set up templates"
    fi
}

# Create standard workflow files
create_workflow_files() {
    echo -e "${BLUE}Creating standard workflow files...${NC}"
    
    # Ensure workflow directory exists
    mkdir -p .github/workflows
    
    # Create issue-automation.yml
    cat > .github/workflows/issue-automation.yml << 'EOL'
name: Issue Automation

on:
  issues:
    types: [opened, labeled, unlabeled, edited]

jobs:
  triage-issue:
    runs-on: ubuntu-latest
    steps:
      - name: Initial issue labeling
        if: github.event.action == 'opened' && !contains(github.event.issue.labels.*.name, 'priority:') && !contains(github.event.issue.labels.*.name, 'component:')
        uses: actions/github-script@v6
        with:
          github-token: ${{secrets.GITHUB_TOKEN}}
          script: |
            github.rest.issues.addLabels({
              issue_number: context.issue.number,
              owner: context.repo.owner,
              repo: context.repo.repo,
              labels: ['status:needs-triage', 'priority:medium']
            })

      - name: Welcome new contributor
        if: github.event.action == 'opened'
        uses: actions/github-script@v6
        with:
          github-token: ${{secrets.GITHUB_TOKEN}}
          script: |
            const creator = context.payload.issue.user.login;
            const repo = context.repo.repo;
            
            const response = await github.rest.issues.listForRepo({
              owner: context.repo.owner,
              repo: context.repo.repo,
              creator: creator,
              state: 'all'
            });
            
            // If this is their first issue, welcome them
            if (response.data.length === 1) {
              github.rest.issues.createComment({
                issue_number: context.issue.number,
                owner: context.repo.owner,
                repo: context.repo.repo,
                body: `Thanks for opening your first issue in the ${repo} project! We'll review this issue soon and provide feedback.`
              });
            }
EOL

    # Create stale-issues.yml
    cat > .github/workflows/stale-issues.yml << 'EOL'
name: Handle Stale Issues

on:
  schedule:
    - cron: '0 0 * * *' # Run daily at midnight UTC

jobs:
  stale:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/stale@v8
        with:
          repo-token: ${{ secrets.GITHUB_TOKEN }}
          stale-issue-message: 'This issue has been automatically marked as stale because it has not had recent activity. It will be closed in 7 days if no further activity occurs.'
          close-issue-message: 'This issue has been automatically closed due to inactivity. Please reopen if this issue is still relevant.'
          days-before-stale: 30
          days-before-close: 7
          exempt-issue-labels: 'priority:critical,priority:high,status:in-progress'
          stale-issue-label: 'status:stale'
EOL

    echo -e "${GREEN}Standard workflow files created successfully!${NC}"
}

# Export current project structure to JSON
export_project_structure() {
    echo -e "${BLUE}Exporting current project structure...${NC}"
    
    # Get milestone data
    milestones=$(gh api repos/$repo_owner/$repo_name/milestones --paginate | jq -c '.')
    echo "$milestones" > "$OUTPUT_DIR/milestones.json"
    
    # Get label data
    labels=$(gh api repos/$repo_owner/$repo_name/labels --paginate | jq -c '.')
    echo "$labels" > "$OUTPUT_DIR/labels.json"
    
    # Get issues (open and closed)
    open_issues=$(gh api repos/$repo_owner/$repo_name/issues --paginate -X GET -f state=open -f per_page=100 | jq -c '.')
    closed_issues=$(gh api repos/$repo_owner/$repo_name/issues --paginate -X GET -f state=closed -f per_page=100 | jq -c '.')
    
    echo "$open_issues" > "$OUTPUT_DIR/open_issues.json"
    echo "$closed_issues" > "$OUTPUT_DIR/closed_issues.json"
    
    # Create a summary report
    repo_info=$(gh api repos/$repo_owner/$repo_name)
    
    cat > "$OUTPUT_DIR/repository_summary.md" << EOL
# Repository Summary: $repo_owner/$repo_name

## Overview
- **Description:** $(echo "$repo_info" | jq -r '.description // "No description"')
- **Created:** $(echo "$repo_info" | jq -r '.created_at')
- **Last Updated:** $(echo "$repo_info" | jq -r '.updated_at')
- **Default Branch:** $(echo "$repo_info" | jq -r '.default_branch')
- **Open Issues Count:** $(echo "$repo_info" | jq -r '.open_issues_count')
- **Language:** $(echo "$repo_info" | jq -r '.language // "Not specified"')

## Statistics
- **Stars:** $(echo "$repo_info" | jq -r '.stargazers_count')
- **Forks:** $(echo "$repo_info" | jq -r '.forks_count')
- **Watchers:** $(echo "$repo_info" | jq -r '.subscribers_count // "Unknown"')

## Milestones
$(echo "$milestones" | jq -r '.[] | "- **" + .title + ":** " + (.description // "No description") + " (Due: " + (.due_on // "No due date") + ")"')

## Labels
$(echo "$labels" | jq -r '.[] | "- **" + .name + ":** " + (.description // "No description")')

## Issues
- **Open Issues:** $(echo "$open_issues" | jq '. | length')
- **Closed Issues:** $(echo "$closed_issues" | jq '. | length')

## Repository Structure

Repository data has been exported to the following files:
- milestones.json
- labels.json
- open_issues.json
- closed_issues.json

Generated on: $(date)
EOL
    
    echo -e "${GREEN}Project structure exported successfully!${NC}"
    echo -e "${BLUE}Summary report: $OUTPUT_DIR/repository_summary.md${NC}"
}

# Delete unnecessary labels
delete_labels() {
    echo -e "${BLUE}Delete unnecessary labels...${NC}"
    
    # Check if we have the list of existing labels
    if [ ! -f "$OUTPUT_DIR/existing_labels.json" ]; then
        echo -e "${RED}Labels file not found. Run 'Clean up existing labels' action first.${NC}"
        return
    fi
    
    # Load existing labels
    existing_labels=$(cat "$OUTPUT_DIR/existing_labels.json")
    
    # Display labels for selection
    echo -e "${BLUE}Available labels to delete:${NC}"
    echo "$existing_labels" | jq -r '.name'
    
    # Ask for labels to delete
    read -p "Enter comma-separated labels to delete (or 'all' to see deletion menu): " labels_input
    
    if [ "$labels_input" == "all" ]; then
        echo -e "${BLUE}Select label category to delete:${NC}"
        echo "1. Delete default GitHub labels (bug, documentation, duplicate, etc.)"
        echo "2. Delete all custom labels (non-standard)"
        echo "3. Delete all labels"
        echo "4. Cancel"
        
        read -p "Enter option (1-4): " delete_option
        
        case $delete_option in
            1)
                # Default GitHub labels
                default_labels=("bug" "documentation" "duplicate" "enhancement" "good first issue" "help wanted" "invalid" "question" "wontfix")
                for label in "${default_labels[@]}"; do
                    echo -e "${YELLOW}Deleting label: $label${NC}"
                    gh api -X DELETE repos/$repo_owner/$repo_name/labels/"$label" || echo "Label '$label' not found"
                done
                ;;
            2)
                # Delete all custom labels (non-standard)
                default_labels=("bug" "documentation" "duplicate" "enhancement" "good first issue" "help wanted" "invalid" "question" "wontfix")
                echo "$existing_labels" | jq -r '.name' | while read -r label; do
                    if [[ ! " ${default_labels[@]} " =~ " ${label} " ]]; then
                        echo -e "${YELLOW}Deleting label: $label${NC}"
                        gh api -X DELETE repos/$repo_owner/$repo_name/labels/"$label" || echo "Label '$label' not found"
                    fi
                done
                ;;
            3)
                # Delete all labels
                echo "$existing_labels" | jq -r '.name' | while read -r label; do
                    echo -e "${YELLOW}Deleting label: $label${NC}"
                    gh api -X DELETE repos/$repo_owner/$repo_name/labels/"$label" || echo "Label '$label' not found"
                done
                ;;
            4)
                echo "Operation cancelled."
                ;;
            *)
                echo -e "${RED}Invalid option${NC}"
                ;;
        esac
    else
        # Delete specific labels
        IFS=',' read -ra LABELS <<< "$labels_input"
        for label in "${LABELS[@]}"; do
            label=$(echo "$label" | xargs)  # Trim whitespace
            echo -e "${YELLOW}Deleting label: $label${NC}"
            gh api -X DELETE repos/$repo_owner/$repo_name/labels/"$label" || echo "Label '$label' not found"
        done
    fi
    
    echo -e "${GREEN}Label deletion completed!${NC}"
    echo -e "${YELLOW}Run 'Clean up existing labels' again to see the updated label list${NC}"
}

# Apply metadata to all issues
apply_metadata() {
    echo -e "${BLUE}Apply metadata to all issues...${NC}"
    
    # Get all issues
    all_issues=$(gh issue list --limit 500 --state all --json number,title,state,labels)
    
    # Metadata types
    echo -e "${BLUE}Select type of metadata to apply:${NC}"
    echo "1. Apply component labels based on title"
    echo "2. Apply priority labels to all unlabeled issues"
    echo "3. Add 'needs-triage' label to issues without labels"
    echo "4. Apply milestone to issues without milestone"
    echo "5. Cancel"
    
    read -p "Enter option (1-5): " metadata_option
    
    case $metadata_option in
        1)
            # Apply component labels based on title
            echo "$all_issues" | jq -c '.[]' | while read -r issue; do
                issue_number=$(echo "$issue" | jq -r '.number')
                issue_title=$(echo "$issue" | jq -r '.title')
                issue_labels=$(echo "$issue" | jq -r '.labels[].name')
                
                # Check if issue already has a component label
                if echo "$issue_labels" | grep -q "component:"; then
                    echo -e "${YELLOW}Issue #$issue_number already has a component label. Skipping.${NC}"
                    continue
                fi
                
                # Determine component based on title
                title_lower=$(echo "$issue_title" | tr '[:upper:]' '[:lower:]')
                
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
                elif [[ "$title_lower" == *"deploy"* ]] || [[ "$title_lower" == *"infrastructure"* ]]; then
                    gh issue edit "$issue_number" --add-label "component:infrastructure"
                else
                    gh issue edit "$issue_number" --add-label "component:other"
                fi
                
                echo -e "${GREEN}Applied component label to issue #$issue_number${NC}"
            done
            ;;
        2)
            # Apply priority labels to all unlabeled issues
            echo "$all_issues" | jq -c '.[]' | while read -r issue; do
                issue_number=$(echo "$issue" | jq -r '.number')
                issue_labels=$(echo "$issue" | jq -r '.labels[].name')
                
                # Check if issue already has a priority label
                if echo "$issue_labels" | grep -q "priority:"; then
                    echo -e "${YELLOW}Issue #$issue_number already has a priority label. Skipping.${NC}"
                    continue
                fi
                
                # Apply medium priority by default
                gh issue edit "$issue_number" --add-label "priority:medium"
                echo -e "${GREEN}Applied 'priority:medium' label to issue #$issue_number${NC}"
            done
            ;;
        3)
            # Add 'needs-triage' label to issues without labels
            echo "$all_issues" | jq -c '.[]' | while read -r issue; do
                issue_number=$(echo "$issue" | jq -r '.number')
                labels_count=$(echo "$issue" | jq '.labels | length')
                
                if [ "$labels_count" -eq 0 ]; then
                    gh issue edit "$issue_number" --add-label "status:needs-triage"
                    echo -e "${GREEN}Applied 'status:needs-triage' label to issue #$issue_number${NC}"
                fi
            done
            ;;
        4)
            # Apply milestone to issues without milestone
            echo -e "${BLUE}Available milestones:${NC}"
            milestones=$(gh api repos/$repo_owner/$repo_name/milestones)
            echo "$milestones" | jq -r '.[] | .number, .title, ""'
            
            read -p "Enter milestone number to apply: " milestone_number
            
            echo "$all_issues" | jq -c '.[]' | while read -r issue; do
                issue_number=$(echo "$issue" | jq -r '.number')
                issue_milestone=$(echo "$issue" | jq -r '.milestone')
                
                if [ "$issue_milestone" == "null" ]; then
                    gh issue edit "$issue_number" --milestone "$milestone_number"
                    echo -e "${GREEN}Applied milestone #$milestone_number to issue #$issue_number${NC}"
                fi
            done
            ;;
        5)
            echo "Operation cancelled."
            ;;
        *)
            echo -e "${RED}Invalid option${NC}"
            ;;
    esac
    
    echo -e "${GREEN}Metadata application completed!${NC}"
}

# Execute selected action
case $action in
    1) create_standard_labels ;;
    2) create_project_milestones ;;
    3) create_issue_templates ;;
    4) clean_up_labels ;;
    5) check_issue_templates ;;
    6) create_workflow_files ;;
    7) export_project_structure ;;
    8) delete_labels ;;
    9) apply_metadata ;;
    0) echo "Exiting..." && exit 0 ;;
    *) echo -e "${RED}Invalid option selected${NC}" ;;
esac