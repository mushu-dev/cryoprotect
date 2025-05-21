#!/bin/bash
# Script to create a standardized workflow GitHub issue template

# Set terminal colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

function show_help {
  echo -e "${BLUE}GitHub Workflow Issue Template Generator${NC}"
  echo -e "Usage: /github-workflow-template [title] [component] [priority]"
  echo
  echo -e "Creates a new GitHub issue with a standardized template for complex workflows."
  echo -e "Parameters:"
  echo -e "  ${GREEN}title${NC}      Issue title (in quotes)"
  echo -e "  ${GREEN}component${NC}  Component (database, api, ui, chembl, authentication, testing)"
  echo -e "  ${GREEN}priority${NC}   Priority (high, medium, low)"
  echo
  echo -e "Example:"
  echo -e "  /github-workflow-template \"Configure SELinux for Podman\" database high"
}

# Ensure GitHub CLI is available
if ! command -v gh &> /dev/null; then
  if [ -f "/home/mushu/Projects/CryoProtect/gh-launcher.sh" ]; then
    GH_CMD="/home/mushu/Projects/CryoProtect/gh-launcher.sh"
  else
    echo -e "${RED}Error: GitHub CLI not found. Please install it first.${NC}"
    exit 1
  fi
else
  GH_CMD="gh"
fi

# Check required arguments
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
  show_help
  exit 1
fi

TITLE="$1"
COMPONENT="$2"
PRIORITY="$3"

# Validate component and priority
VALID_COMPONENTS=("database" "api" "ui" "chembl" "authentication" "testing")
VALID_PRIORITIES=("high" "medium" "low")

COMPONENT_VALID=false
for comp in "${VALID_COMPONENTS[@]}"; do
  if [ "$comp" = "$COMPONENT" ]; then
    COMPONENT_VALID=true
    break
  fi
done

PRIORITY_VALID=false
for prio in "${VALID_PRIORITIES[@]}"; do
  if [ "$prio" = "$PRIORITY" ]; then
    PRIORITY_VALID=true
    break
  fi
done

if [ "$COMPONENT_VALID" = false ]; then
  echo -e "${RED}Error: Invalid component \"$COMPONENT\".${NC}"
  echo -e "Valid components: ${YELLOW}${VALID_COMPONENTS[@]}${NC}"
  exit 1
fi

if [ "$PRIORITY_VALID" = false ]; then
  echo -e "${RED}Error: Invalid priority \"$PRIORITY\".${NC}"
  echo -e "Valid priorities: ${YELLOW}${VALID_PRIORITIES[@]}${NC}"
  exit 1
fi

# Create issue template content
ISSUE_TEMPLATE="# $TITLE

## Overview
*Brief description of the task and its purpose...*

## Prerequisites
- [ ] Required tool 1
- [ ] Required tool 2
- [ ] Required access or permissions

## Environment Configuration

### Local Environment
\`\`\`bash
# Configuration commands for local environment
\`\`\`

### Supabase Database
\`\`\`sql
-- Required database configuration or queries
\`\`\`

### Vercel Deployment
*Configuration steps for Vercel if applicable...*

## Step-by-Step Process

### Phase 1: Preparation
- [ ] Step 1
- [ ] Step 2
- [ ] Step 3

### Phase 2: Implementation
- [ ] Step 1
- [ ] Step 2
- [ ] Step 3

### Phase 3: Verification
- [ ] Step 1
- [ ] Step 2
- [ ] Step 3

## Verification Steps
- [ ] Verify requirement 1
- [ ] Verify requirement 2
- [ ] Verify requirement 3

## Rollback Procedure
\`\`\`bash
# Commands to undo changes if needed
\`\`\`

## Working Notes
*Use this section as a scratchpad for progress, observations, and intermediate results...*

## References
- [Related documentation link 1]()
- [Related documentation link 2]()
- Related issues: #
"

# Create the issue
echo -e "${BLUE}Creating workflow issue: ${YELLOW}$TITLE${NC}"
echo -e "${BLUE}Component: ${YELLOW}$COMPONENT${NC}, Priority: ${YELLOW}$PRIORITY${NC}"

$GH_CMD issue create \
  --title "$TITLE" \
  --body "$ISSUE_TEMPLATE" \
  --label "component:$COMPONENT,priority:$PRIORITY,type:workflow"

echo -e "${GREEN}Workflow issue created successfully!${NC}"