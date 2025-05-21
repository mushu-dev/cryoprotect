#!/bin/bash
# github-prd.sh - Helper script for managing Product Requirements Documents (PRDs)
# and integrating with TaskMaster AI
# This command is automatically run when invoked with /github-prd

# Colors for terminal output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Ensure scripts directory exists
mkdir -p scripts

# Display help information
show_help() {
  echo -e "${BLUE}GitHub PRD Manager - TaskMaster Integration Helper${NC}"
  echo ""
  echo "Usage: $0 [command] [options]"
  echo ""
  echo "Commands:"
  echo "  create [title] [issue-number]  - Create a new PRD from a GitHub issue and save to scripts/"
  echo "  extract [issue-number]         - Extract PRD sections from an existing GitHub issue"
  echo "  generate-tasks [prd-file]      - Generate TaskMaster tasks from a PRD file"
  echo "  view [prd-file]                - Display a PRD file with syntax highlighting"
  echo "  list                           - List all PRD files in scripts/ directory"
  echo "  template [title]               - Create a blank PRD template"
  echo "  link [prd-file] [issue-number] - Link a PRD file to a GitHub issue for tracking"
  echo "  update-mapping [prd] [issue#]  - Update GitHub to TaskMaster task mapping"
  echo "  mapping [issue#]               - Show GitHub issue to TaskMaster task mapping"
  echo "  help                           - Show this help message"
  echo ""
  echo "Examples:"
  echo "  $0 create \"Authentication System Redesign\" 42"
  echo "  $0 extract 42"
  echo "  $0 generate-tasks scripts/prd_auth_system.txt"
  echo "  $0 view scripts/prd_auth_system.txt"
  echo ""
}

# Create a new PRD from a GitHub issue
create_prd() {
  if [ -z "$1" ] || [ -z "$2" ]; then
    echo -e "${RED}Error: Missing title or issue number${NC}"
    echo "Usage: $0 create \"Title\" issue-number"
    exit 1
  fi

  title="$1"
  issue_number="$2"
  
  # Create filename from title
  filename=$(echo "$title" | tr '[:upper:]' '[:lower:]' | tr ' ' '_')
  prd_file="scripts/prd_${filename}.txt"
  
  echo -e "${BLUE}Creating PRD from issue #${issue_number}...${NC}"
  
  # Get issue details using GitHub CLI
  issue_json=$(gh issue view "$issue_number" --json title,body,labels)
  
  if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to retrieve issue #${issue_number}${NC}"
    exit 1
  fi
  
  # Extract issue body
  issue_body=$(echo "$issue_json" | jq -r '.body')
  issue_title=$(echo "$issue_json" | jq -r '.title')
  
  # Create PRD with standardized format
  cat > "$prd_file" << EOL
# Product Requirements Document: $title
# Related Issue: #$issue_number
# Created: $(date +"%Y-%m-%d")

## Overview
$(echo "$issue_body" | grep -Pzo "(?s)## Overview.*?(?=##|\Z)" | tail -n +2)

## Requirements
$(echo "$issue_body" | grep -Pzo "(?s)## Requirements.*?(?=##|\Z)" | tail -n +2)

## Technical Context
$(echo "$issue_body" | grep -Pzo "(?s)## Technical Context.*?(?=##|\Z)" | tail -n +2)

## Constraints
$(echo "$issue_body" | grep -Pzo "(?s)## Constraints.*?(?=##|\Z)" | tail -n +2)

## Acceptance Criteria
$(echo "$issue_body" | grep -Pzo "(?s)## Acceptance Criteria.*?(?=##|\Z)" | tail -n +2)

## Additional Notes
- Generated with github-prd.sh from issue #$issue_number
- Reference this PRD in code comments and commits
- Update this PRD as implementation details evolve
EOL

  echo -e "${GREEN}PRD created: ${prd_file}${NC}"
  echo -e "${YELLOW}Next steps:${NC}"
  echo -e "  - Review and edit the PRD as needed"
  echo -e "  - Generate tasks with: github-prd generate-tasks ${prd_file}"
}

# Extract PRD sections from a GitHub issue
extract_prd() {
  if [ -z "$1" ]; then
    echo -e "${RED}Error: Missing issue number${NC}"
    echo "Usage: $0 extract issue-number"
    exit 1
  fi

  issue_number="$1"
  
  echo -e "${BLUE}Extracting PRD sections from issue #${issue_number}...${NC}"
  
  # Get issue details using GitHub CLI
  issue_json=$(gh issue view "$issue_number" --json title,body,labels)
  
  if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to retrieve issue #${issue_number}${NC}"
    exit 1
  fi
  
  # Extract issue body and title
  issue_body=$(echo "$issue_json" | jq -r '.body')
  issue_title=$(echo "$issue_json" | jq -r '.title')
  
  # Check for PRD sections
  if ! echo "$issue_body" | grep -q "## Overview"; then
    echo -e "${YELLOW}Warning: Issue doesn't appear to have PRD sections${NC}"
    echo "Would you like to create a template? (y/n)"
    read answer
    if [[ "$answer" == "y" ]]; then
      title="$issue_title"
      filename=$(echo "$title" | tr '[:upper:]' '[:lower:]' | tr ' ' '_')
      template_prd "scripts/prd_${filename}.txt" "$title" "$issue_number"
      echo -e "${GREEN}Template created: scripts/prd_${filename}.txt${NC}"
      echo -e "${YELLOW}Please edit the template and update the issue${NC}"
    fi
    exit 0
  fi
  
  # Display PRD sections
  echo -e "${BLUE}Overview:${NC}"
  echo "$issue_body" | grep -Pzo "(?s)## Overview.*?(?=##|\Z)" | tail -n +2
  
  echo -e "${BLUE}Requirements:${NC}"
  echo "$issue_body" | grep -Pzo "(?s)## Requirements.*?(?=##|\Z)" | tail -n +2
  
  echo -e "${BLUE}Technical Context:${NC}"
  echo "$issue_body" | grep -Pzo "(?s)## Technical Context.*?(?=##|\Z)" | tail -n +2
  
  echo -e "${BLUE}Constraints:${NC}"
  echo "$issue_body" | grep -Pzo "(?s)## Constraints.*?(?=##|\Z)" | tail -n +2
  
  echo -e "${BLUE}Acceptance Criteria:${NC}"
  echo "$issue_body" | grep -Pzo "(?s)## Acceptance Criteria.*?(?=##|\Z)" | tail -n +2
}

# Generate TaskMaster tasks from a PRD file
generate_tasks() {
  if [ -z "$1" ]; then
    echo -e "${RED}Error: Missing PRD file${NC}"
    echo "Usage: $0 generate-tasks prd-file"
    exit 1
  fi

  prd_file="$1"
  
  if [ ! -f "$prd_file" ]; then
    echo -e "${RED}Error: PRD file not found: ${prd_file}${NC}"
    exit 1
  fi
  
  echo -e "${BLUE}Generating TaskMaster tasks from ${prd_file}...${NC}"
  
  # Check if task-master is installed
  if ! command -v npx &> /dev/null; then
    echo -e "${RED}Error: npx not found. Please install Node.js${NC}"
    exit 1
  fi
  
  # Generate tasks using TaskMaster
  npx -y task-master-ai parse-prd "$prd_file"
  
  if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to generate tasks${NC}"
    echo "Please make sure TaskMaster is properly installed and configured"
    exit 1
  fi
  
  echo -e "${GREEN}Tasks generated successfully${NC}"
  echo -e "${YELLOW}Next steps:${NC}"
  echo -e "  - Review tasks with: npx task-master-ai list"
  echo -e "  - Analyze complexity with: npx task-master-ai analyze-complexity"
  echo -e "  - Generate task files with: npx task-master-ai generate"
}

# View a PRD file with syntax highlighting
view_prd() {
  if [ -z "$1" ]; then
    echo -e "${RED}Error: Missing PRD file${NC}"
    echo "Usage: $0 view prd-file"
    exit 1
  fi

  prd_file="$1"
  
  if [ ! -f "$prd_file" ]; then
    echo -e "${RED}Error: PRD file not found: ${prd_file}${NC}"
    exit 1
  fi
  
  # Check if bat (syntax highlighter) is available, otherwise use cat
  if command -v bat &> /dev/null; then
    bat --style=plain "$prd_file"
  else
    cat "$prd_file"
  fi
}

# List all PRD files in scripts/ directory
list_prds() {
  echo -e "${BLUE}Available PRD files:${NC}"
  
  if ! ls scripts/prd_*.txt 2>/dev/null; then
    echo -e "${YELLOW}No PRD files found in scripts/ directory${NC}"
    echo "Create a PRD with: $0 create \"Title\" issue-number"
    echo "Or create a template with: $0 template \"Title\""
  fi
}

# Create a blank PRD template
template_prd() {
  local file="$1"
  local title="$2"
  local issue="$3"
  
  cat > "$file" << EOL
# Product Requirements Document: $title
# Related Issue: #$issue
# Created: $(date +"%Y-%m-%d")

## Overview
[Provide a high-level description of the feature or component]

## Requirements
- [Requirement 1]
- [Requirement 2]
- [Requirement 3]

## Technical Context
- **Current Architecture**: [Describe current system architecture related to this PRD]
- **Integration Points**: [List systems/components this will interact with]
- **Existing Code**: [Reference relevant existing code]

## Constraints
- **Performance**: [Performance requirements]
- **Security**: [Security requirements]
- **Compatibility**: [Compatibility requirements]
- **Resource Limitations**: [Any resource constraints]

## Acceptance Criteria
- [ ] [Criterion 1]
- [ ] [Criterion 2]
- [ ] [Criterion 3]

## Additional Notes
- [Any additional context, references, or notes]
- [Links to similar features or implementations]
- [Known challenges]
EOL
}

# Create a blank PRD template with a given title
create_template() {
  if [ -z "$1" ]; then
    echo -e "${RED}Error: Missing title${NC}"
    echo "Usage: $0 template \"Title\""
    exit 1
  fi

  title="$1"
  issue="N/A"
  if [ ! -z "$2" ]; then
    issue="$2"
  fi
  
  # Create filename from title
  filename=$(echo "$title" | tr '[:upper:]' '[:lower:]' | tr ' ' '_')
  prd_file="scripts/prd_${filename}.txt"
  
  # Create PRD template
  template_prd "$prd_file" "$title" "$issue"
  
  echo -e "${GREEN}Template created: ${prd_file}${NC}"
  echo -e "${YELLOW}Next steps:${NC}"
  echo -e "  - Edit the template with your requirements"
  echo -e "  - Generate tasks with: github-prd generate-tasks ${prd_file}"
}

# Link a PRD file to a GitHub issue for tracking
link_prd() {
  if [ -z "$1" ] || [ -z "$2" ]; then
    echo -e "${RED}Error: Missing PRD file or issue number${NC}"
    echo "Usage: $0 link prd-file issue-number"
    exit 1
  fi

  prd_file="$1"
  issue_number="$2"
  
  if [ ! -f "$prd_file" ]; then
    echo -e "${RED}Error: PRD file not found: ${prd_file}${NC}"
    exit 1
  fi
  
  echo -e "${BLUE}Linking ${prd_file} to issue #${issue_number}...${NC}"
  
  # Get PRD content
  prd_content=$(cat "$prd_file")
  
  # Create comment with PRD content
  gh issue comment "$issue_number" --body "## Linked PRD\n\`\`\`\n${prd_content}\n\`\`\`"
  
  if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to link PRD to issue #${issue_number}${NC}"
    exit 1
  fi
  
  echo -e "${GREEN}PRD successfully linked to issue #${issue_number}${NC}"
}

# Create task mapping from GitHub issue to TaskMaster
update_task_mapping() {
  if [ -z "$1" ] || [ -z "$2" ]; then
    echo -e "${RED}Error: Missing PRD file or issue number${NC}"
    echo "Usage: $0 update-mapping prd-file issue-number"
    exit 1
  fi

  prd_file="$1"
  issue_number="$2"
  
  if [ ! -f "$prd_file" ]; then
    echo -e "${RED}Error: PRD file not found: ${prd_file}${NC}"
    exit 1
  fi
  
  # Ensure tasks/ directory exists
  mkdir -p tasks
  
  # Create or update mapping file
  mapping_file="tasks/github_mapping.json"
  
  if [ ! -f "$mapping_file" ]; then
    # Create initial mapping file
    echo '{
  "mappings": [],
  "last_updated": "'"$(date -Iseconds)"'"
}' > "$mapping_file"
  fi
  
  # Get task IDs if TaskMaster tasks.json exists
  if [ -f "tasks/tasks.json" ]; then
    task_ids=$(jq -r '.tasks[].id' tasks/tasks.json 2>/dev/null | jq -Rs 'split("\n") | map(select(length > 0))')
    
    # Get PRD title
    prd_title=$(grep "^# " "$prd_file" | sed 's/^# Product Requirements Document: //' | sed 's/^# //')
    
    # Update mapping file
    jq --arg issue "$issue_number" \
       --arg title "$prd_title" \
       --arg prd "$prd_file" \
       --argjson tasks "$task_ids" \
       --arg date "$(date -Iseconds)" \
       '.mappings = (.mappings | map(select(.github_issue != ($issue | tonumber)))) + [{
         "github_issue": ($issue | tonumber),
         "title": $title,
         "taskmaster_tasks": $tasks,
         "prd_path": $prd,
         "last_synced": $date
       }] | .last_updated = $date' "$mapping_file" > tmp.json && mv tmp.json "$mapping_file"
    
    echo -e "${GREEN}Task mapping updated for GitHub issue #${issue_number}${NC}"
    echo -e "${BLUE}Linked TaskMaster tasks:${NC}"
    jq -r '.mappings[] | select(.github_issue == '"$issue_number"') | .taskmaster_tasks[] | tostring' "$mapping_file" | sed 's/^/- /'
  else
    echo -e "${YELLOW}Warning: tasks.json not found. Run task-master parse-prd first${NC}"
  fi
}

# Show mapping between GitHub issues and TaskMaster tasks
show_mapping() {
  mapping_file="tasks/github_mapping.json"
  
  if [ ! -f "$mapping_file" ]; then
    echo -e "${YELLOW}No GitHub-TaskMaster mapping found${NC}"
    echo "Create mappings with: $0 update-mapping prd-file issue-number"
    return
  fi
  
  if [ -z "$1" ]; then
    # Show all mappings
    echo -e "${BLUE}GitHub Issues to TaskMaster Task Mappings:${NC}"
    jq -r '.mappings[] | "Issue #\(.github_issue): \(.title) → \(.taskmaster_tasks | length) tasks"' "$mapping_file"
  else
    # Show specific mapping
    issue_number="$1"
    echo -e "${BLUE}Mapping for GitHub Issue #${issue_number}:${NC}"
    if jq -e '.mappings[] | select(.github_issue == '"$issue_number"')' "$mapping_file" > /dev/null; then
      jq -r '.mappings[] | select(.github_issue == '"$issue_number"') | "Title: \(.title)\nPRD: \(.prd_path)\nTasks:\n" + (.taskmaster_tasks | map("- " + (. | tostring)) | join("\n"))' "$mapping_file"
    else
      echo -e "${YELLOW}No mapping found for GitHub issue #${issue_number}${NC}"
    fi
  fi
}

# Main command processing
case "$1" in
  create)
    create_prd "$2" "$3"
    ;;
  extract)
    extract_prd "$2"
    ;;
  generate-tasks)
    generate_tasks "$2"
    # Auto-update mapping if this is a PRD linked to an issue
    issue_number=$(grep -o "Related Issue: #[0-9]\+" "$2" | grep -o "[0-9]\+")
    if [ -n "$issue_number" ]; then
      update_task_mapping "$2" "$issue_number"
    fi
    ;;
  view)
    view_prd "$2"
    ;;
  list)
    list_prds
    ;;
  template)
    create_template "$2" "$3"
    ;;
  link)
    link_prd "$2" "$3"
    ;;
  update-mapping)
    update_task_mapping "$2" "$3"
    ;;
  mapping)
    show_mapping "$2"
    ;;
  help|*)
    show_help
    ;;
esac