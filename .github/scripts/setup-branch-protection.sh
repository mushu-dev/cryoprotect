#!/bin/bash
# Script to set up branch protection rules for CryoProtect repository
# Requires GitHub CLI and appropriate permissions

# Check if gh is installed
if ! command -v gh &> /dev/null; then
    echo "GitHub CLI (gh) is not installed. Please install it first."
    exit 1
fi

# Check if user is authenticated
if ! gh auth status &> /dev/null; then
    echo "You are not authenticated with GitHub CLI. Please run 'gh auth login' first."
    exit 1
fi

# Repository details
REPO="blueprint-house/CryoProtect"
BRANCH="master"

echo "Setting up branch protection rules for $REPO"

# Create branch protection rule for master branch
echo "Creating branch protection rule for $BRANCH branch..."

gh api \
  --method PUT \
  -H "Accept: application/vnd.github+json" \
  -H "X-GitHub-Api-Version: 2022-11-28" \
  "/repos/$REPO/branches/$BRANCH/protection" \
  -f required_status_checks='{
    "strict": true,
    "contexts": [
      "Run Tests (python)",
      "Run Tests (node)",
      "Security Scan",
      "deploy",
      "build"
    ]
  }' \
  -f enforce_admins=true \
  -f required_pull_request_reviews='{
    "dismissal_restrictions": {},
    "dismiss_stale_reviews": true,
    "require_code_owner_reviews": false,
    "required_approving_review_count": 1,
    "bypass_pull_request_allowances": {}
  }' \
  -f restrictions=null \
  -f required_linear_history=true \
  -f allow_force_pushes=false \
  -f allow_deletions=false \
  -f block_creations=false \
  -f required_conversation_resolution=true \
  -f lock_branch=false \
  -f allow_fork_syncing=true

# Check if the command was successful
if [ $? -eq 0 ]; then
    echo "Branch protection rules for $BRANCH successfully set up!"
else
    echo "Failed to set up branch protection rules. Check your permissions."
    exit 1
fi

# Set up tag protection
echo "Setting up tag protection rule for v*.*.* pattern..."

gh api \
  --method POST \
  -H "Accept: application/vnd.github+json" \
  -H "X-GitHub-Api-Version: 2022-11-28" \
  "/repos/$REPO/tags/protection" \
  -f pattern="v*.*.*" \
  -f allows_deletions=false \
  -f allows_force_pushes=false

# Check if the command was successful
if [ $? -eq 0 ]; then
    echo "Tag protection rule successfully set up!"
else
    echo "Failed to set up tag protection rule. Check your permissions."
    exit 1
fi

echo "All protection rules have been set up successfully!"
echo "Don't forget to add the required secrets for GitHub Actions workflows:"
echo "- VERCEL_TOKEN, VERCEL_ORG_ID, VERCEL_PROJECT_ID"
echo "- HEROKU_API_KEY, HEROKU_APP_NAME, HEROKU_EMAIL, HEROKU_PIPELINE_ID"