# Adding Issues to GitHub Project Manually

If the automated script doesn't work, you can use these step-by-step commands to add issues to your project board.

## Prerequisites

Make sure GitHub CLI is authenticated:

```bash
gh auth status
```

## Step 1: Get the Project ID

Your project number is 1, but you need the internal project ID to use the GraphQL API.
This command has already been run and your project ID is: `PVT_kwHOBMCBws4A4xSP`

If you need to check project IDs in the future:

```bash
gh project list --owner mushu-dev
```

## Step 2: Get All Issue Numbers

This command will list all issue numbers:

```bash
gh issue list -R mushu-dev/cryoprotect -s all -L 100 --json number | jq -r '.[].number'
```

## Step 3: Add Each Issue to the Project

For each issue number, run the following command (replace ISSUE_NUMBER with the actual number):

```bash
PROJECT_ID="PVT_kwHOBMCBws4A4xSP"
ISSUE_NUMBER="10"  # Replace with actual issue number
REPO="mushu-dev/cryoprotect"

gh api graphql -f query='
  mutation($project:ID!, $issue:ID!) {
    addProjectV2ItemById(input: {projectId: $project, contentId: $issue}) {
      item {
        id
      }
    }
  }
' -f project="$PROJECT_ID" -F issue="$(gh api repos/$REPO/issues/$ISSUE_NUMBER --jq '.node_id')"
```

## Step 4: Add All Issues with a One-liner

This command will add all issues to the project in one go:

```bash
PROJECT_ID="PVT_kwHOBMCBws4A4xSP"
REPO="mushu-dev/cryoprotect"

gh issue list -R $REPO -s all -L 100 --json number | jq -r '.[].number' | xargs -I{} sh -c 'echo "Adding issue #{}" && gh api graphql -f query='\''mutation($project:ID!, $issue:ID!) {addProjectV2ItemById(input: {projectId: $project, contentId: $issue}) {item {id}}}'\'' -f project="$0" -F issue="$(gh api repos/$1/issues/{} --jq '\''.node_id'\'')" && sleep 1' "$PROJECT_ID" "$REPO"
```

## Step 5: Verify Issues Were Added

You can view your project board to verify all issues are now included:

```bash
open https://github.com/orgs/mushu-dev/projects/1
```

Or with the gh CLI:

```bash
gh project view 1 --owner mushu-dev
```