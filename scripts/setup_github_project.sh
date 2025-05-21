#\!/bin/bash
# Script to set up a GitHub project board with proper columns and automation

# Get the repository owner and name
REPO_FULL_NAME=$(gh repo view --json nameWithOwner -q .nameWithOwner)
REPO_OWNER=$(echo $REPO_FULL_NAME  < /dev/null |  cut -d'/' -f1)
REPO_NAME=$(echo $REPO_FULL_NAME | cut -d'/' -f2)

echo "Setting up GitHub project board for $REPO_FULL_NAME"

# Create a new project
echo "Creating new project board..."
PROJECT_ID=$(gh api graphql -f query='
  mutation {
    createProjectV2(input: {ownerId: "'"$REPO_OWNER"'", title: "CryoProtect Development", repositoryId: "'"$REPO_NAME"'"}) {
      projectV2 {
        id
        number
      }
    }
  }
' --jq '.data.createProjectV2.projectV2.id')

echo "Created project with ID: $PROJECT_ID"

# Create columns (views)
echo "Creating project columns..."

# Status columns to create
COLUMNS=("Backlog" "Ready" "In Progress" "Needs Review" "Blocked" "Completed")
COLUMN_FIELDS=()

for COLUMN in "${COLUMNS[@]}"; do
  echo "Creating column: $COLUMN"
  
  # Create the column (view)
  VIEW_ID=$(gh api graphql -f query='
    mutation {
      createProjectV2View(input: {projectId: "'"$PROJECT_ID"'", name: "'"$COLUMN"'"}) {
        projectV2View {
          id
        }
      }
    }
  ' --jq '.data.createProjectV2View.projectV2View.id')
  
  echo "Created column $COLUMN with ID: $VIEW_ID"
  COLUMN_FIELDS+=("$COLUMN:$VIEW_ID")
  
  # Add filter to the view based on status
  STATUS_LABEL=""
  case "$COLUMN" in
    "Backlog")
      STATUS_LABEL="status:planning"
      ;;
    "Ready")
      STATUS_LABEL="status:ready"
      ;;
    "In Progress")
      STATUS_LABEL="status:in-progress"
      ;;
    "Needs Review")
      STATUS_LABEL="status:needs-review"
      ;;
    "Blocked")
      STATUS_LABEL="status:blocked"
      ;;
    "Completed")
      STATUS_LABEL="status:completed"
      ;;
  esac
  
  if [ -n "$STATUS_LABEL" ]; then
    echo "Setting filter for $COLUMN to label: $STATUS_LABEL"
    
    # Set up the filter for this view
    gh api graphql -f query='
      mutation {
        updateProjectV2View(input: {
          projectId: "'"$PROJECT_ID"'",
          viewId: "'"$VIEW_ID"'",
          filteredBy: [{
            fieldName: "labels",
            operator: EQ,
            value: "'"$STATUS_LABEL"'"
          }]
        }) {
          projectV2View {
            id
          }
        }
      }
    '
  fi
done

echo "Project board setup complete."
echo "Project board is available at: https://github.com/$REPO_FULL_NAME/projects"
echo ""
echo "Next steps:"
echo "1. Go to the project board URL above"
echo "2. Add all existing issues to the board using the 'Add items' button"
echo "3. Set up additional views for different areas (database, API, etc.)"
echo "4. Configure automation rules in the project settings"
