# GitHub Project Board Setup Guide

Since the automatic setup script requires the 'project' scope which is not yet available in the current token, this guide provides instructions for manually setting up a GitHub project board.

## Step 1: Create a New Project

1. Navigate to the repository: https://github.com/blueprint-house/CryoProtect
2. Click on the "Projects" tab
3. Click "New project"
4. Select "Table" as the template
5. Enter "CryoProtect Development" as the project name
6. Click "Create"

## Step 2: Set Up Status Columns

1. In the project board, click "+ Add fields" in the top right
2. Select "Single select" to create a custom field
3. Name the field "Status"
4. Add the following status values (with their recommended colors):
   - Planning (Yellow)
   - Ready (Green)
   - In Progress (Blue)
   - Needs Review (Purple)
   - Blocked (Red)
   - Completed (Gray)
5. Click "Save"

## Step 3: Create Status Views

1. Click on "Views" in the sidebar
2. Click the "+" button to create a new view
3. Select "Board" as the view type
4. Select "Status" as the field to group by
5. Name the view "Status Board"

Repeat this process to create additional filtered views:

1. **Database Tasks View**:
   - Create a new "Board" view
   - Set up a filter for "Label contains area:database"
   - Group by "Status"
   - Name it "Database Tasks"

2. **API Tasks View**:
   - Create a new "Board" view
   - Set up a filter for "Label contains area:api"
   - Group by "Status"
   - Name it "API Tasks"

3. **Testing Tasks View**:
   - Create a new "Board" view
   - Set up a filter for "Label contains area:testing"
   - Group by "Status"
   - Name it "Testing Tasks"

## Step 4: Add Issues to the Project

1. Click "Add items" in the toolbar
2. Select all relevant issues from the repository
3. Click "Add selected items"

## Step 5: Set Up Automation Rules

1. Click on the "..." menu (three dots) in the top right
2. Select "Workflows"
3. Configure the following automation rules:

**Status Label Sync**:
- When: Issues with the label `status:planning` are added to the project
- Action: Set Status field value to "Planning"

Repeat for all status labels:
- `status:ready` → "Ready"
- `status:in-progress` → "In Progress"
- `status:needs-review` → "Needs Review"
- `status:blocked` → "Blocked"
- `status:completed` → "Completed"

**Milestone Tracking**:
- Add a "Milestone" field to track the project milestones
- Create an automation to set milestone values based on issue milestones

## Step 6: Configure Project Settings

1. Click on the "..." menu (three dots) in the top right
2. Select "Settings"
3. Configure the following settings:
   - Set the project description
   - Adjust permissions if needed
   - Enable/disable features as appropriate

## Best Practices

1. **Keep Status In Sync**: Ensure that issue status labels and project board status remain synchronized
2. **Regular Updates**: Move cards to appropriate columns as work progresses
3. **Use Filters**: Create filtered views for different teams or components
4. **Milestone Tracking**: Group issues by milestone for sprint planning
5. **Priority Indicators**: Consider adding priority indicators to high-priority items

## Alternative: Script-Based Setup with New Token

If you prefer to use the script-based setup:

1. Create a new GitHub personal access token with 'repo' and 'project' scopes
2. Log in with the new token: `gh auth login --with-token < token.txt`
3. Run the setup script: `./scripts/setup_github_project.sh`
</content>