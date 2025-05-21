#!/bin/bash
# Setup GitHub Project Board for the new repository

# Set up variables
REPO="mushu-dev/cryoprotect"
BOARD_NAME="CryoProtect Development"
LOG_FILE="github_project_setup_$(date +%Y%m%d_%H%M%S).log"

echo "Setting up GitHub Project for $REPO"
echo "Log will be saved to $LOG_FILE"
echo "-----------------------------------------------"

# Log start time
echo "Setup started at $(date)" > "$LOG_FILE"
echo "Repository: $REPO" >> "$LOG_FILE"
echo "-----------------------------------------------" >> "$LOG_FILE"

# Create the project board
echo "Creating project board '$BOARD_NAME'..."
PROJECT_URL=$(gh project create "$BOARD_NAME" --owner mushu-dev --public)
PROJECT_NUMBER=$(echo "$PROJECT_URL" | grep -oE '[0-9]+$')

if [ -z "$PROJECT_NUMBER" ]; then
  echo "Failed to create project board"
  echo "Failed to create project board" >> "$LOG_FILE"
  exit 1
fi

echo "Created project board #$PROJECT_NUMBER: $PROJECT_URL"
echo "Created project board #$PROJECT_NUMBER: $PROJECT_URL" >> "$LOG_FILE"

# Create default fields
echo "Creating default fields..."

# Status field with standard values
echo "Adding custom field 'Status'..."
gh project field create "$PROJECT_NUMBER" --data-type single-select --name "Status" --owner mushu-dev >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Status --owner mushu-dev --option "Todo" --color "#E99695" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Status --owner mushu-dev --option "In Progress" --color "#F4D03F" >> "$LOG_FILE" 
gh project field-option add "$PROJECT_NUMBER" Status --owner mushu-dev --option "In Review" --color "#6AB04C" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Status --owner mushu-dev --option "Done" --color "#2ECC71" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Status --owner mushu-dev --option "Blocked" --color "#E74C3C" >> "$LOG_FILE"

# Priority field
echo "Adding custom field 'Priority'..."
gh project field create "$PROJECT_NUMBER" --data-type single-select --name "Priority" --owner mushu-dev >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Priority --owner mushu-dev --option "Critical" --color "#FF0000" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Priority --owner mushu-dev --option "High" --color "#FFA500" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Priority --owner mushu-dev --option "Medium" --color "#FFFF00" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Priority --owner mushu-dev --option "Low" --color "#90EE90" >> "$LOG_FILE"

# Area field
echo "Adding custom field 'Area'..."
gh project field create "$PROJECT_NUMBER" --data-type single-select --name "Area" --owner mushu-dev >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Area --owner mushu-dev --option "Database" --color "#3498DB" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Area --owner mushu-dev --option "API" --color "#9B59B6" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Area --owner mushu-dev --option "Frontend" --color "#2ECC71" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Area --owner mushu-dev --option "DevOps" --color "#F1C40F" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Area --owner mushu-dev --option "Security" --color "#E74C3C" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Area --owner mushu-dev --option "Testing" --color "#1ABC9C" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Area --owner mushu-dev --option "Documentation" --color "#95A5A6" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Area --owner mushu-dev --option "RDKit" --color "#D35400" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Area --owner mushu-dev --option "ChEMBL" --color "#8E44AD" >> "$LOG_FILE"
gh project field-option add "$PROJECT_NUMBER" Area --owner mushu-dev --option "PubChem" --color "#16A085" >> "$LOG_FILE"

# Story points field
echo "Adding custom field 'Story Points'..."
gh project field create "$PROJECT_NUMBER" --data-type number --name "Story Points" --owner mushu-dev >> "$LOG_FILE"

# Add existing issues to the project
echo "Adding issues to the project board..."
gh issue list -s all -L 500 --json number -R "$REPO" | jq -r '.[] | .number' | while read -r issue; do
  echo "Adding issue #$issue to project..."
  gh project item-add "$PROJECT_NUMBER" --owner mushu-dev --url "https://github.com/$REPO/issues/$issue" >> "$LOG_FILE"
done

# Log completion
echo "-----------------------------------------------" >> "$LOG_FILE"
echo "Setup completed at $(date)" >> "$LOG_FILE"

echo "-----------------------------------------------"
echo "GitHub Project setup complete!"
echo "Project board: $PROJECT_URL"
echo "Log saved to $LOG_FILE"
echo ""
echo "Next steps:"
echo "1. Visit the project URL to create custom views (Kanban, etc.)"
echo "2. Set up automation rules for status changes"
echo "3. Connect the project to repository workflows"