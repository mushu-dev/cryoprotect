#!/bin/bash
# Script to consolidate and organize GitHub issues

# Set terminal colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

function show_help {
  echo -e "${BLUE}GitHub Issue Consolidation Tool${NC}"
  echo -e "Usage: /github-consolidate [command] [options]"
  echo
  echo -e "Commands:"
  echo -e "  ${GREEN}setup-milestones${NC}     Create standard project milestones"
  echo -e "  ${GREEN}setup-labels${NC}         Create standard project labels"
  echo -e "  ${GREEN}link-issues${NC} [id1] [id2] Link two related issues"
  echo -e "  ${GREEN}categorize${NC} [id] [component] [priority] Categorize an issue"
  echo -e "  ${GREEN}assign-milestone${NC} [id] [milestone] Assign issue to milestone"
  echo -e "  ${GREEN}status${NC}              Show current repo organization status"
  echo
  echo -e "Examples:"
  echo -e "  /github-consolidate setup-milestones"
  echo -e "  /github-consolidate link-issues 123 124"
  echo -e "  /github-consolidate categorize 125 database high"
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

# Process command
CMD=${1:-status}

case $CMD in
  "setup-milestones")
    echo -e "${BLUE}Creating standard project milestones...${NC}"
    $GH_CMD api repos/:owner/:repo/milestones -X POST -F title="Phase 1: Technical Foundation" -F state="open" -F description="Database Architecture & Authentication System"
    $GH_CMD api repos/:owner/:repo/milestones -X POST -F title="Phase 2: Feature Completion" -F state="open" -F description="API Layer, Core Functionality, User Interface"
    $GH_CMD api repos/:owner/:repo/milestones -X POST -F title="Phase 3: Production Readiness" -F state="open" -F description="Deployment, Monitoring, Security"
    $GH_CMD api repos/:owner/:repo/milestones -X POST -F title="Phase 4: Documentation" -F state="open" -F description="Documentation and Knowledge Transfer"
    $GH_CMD api repos/:owner/:repo/milestones -X POST -F title="ChEMBL Integration" -F state="open" -F description="ChEMBL Database Integration"
    echo -e "${GREEN}Milestones created successfully!${NC}"
    ;;
    
  "setup-labels")
    echo -e "${BLUE}Creating standard project labels...${NC}"
    $GH_CMD label create "priority:high" --color "#FF0000" --description "High priority task"
    $GH_CMD label create "priority:medium" --color "#FFA500" --description "Medium priority task"
    $GH_CMD label create "priority:low" --color "#FFFF00" --description "Low priority task"
    $GH_CMD label create "component:database" --color "#0075CA" --description "Related to database functionality"
    $GH_CMD label create "component:api" --color "#A2EEEF" --description "Related to API functionality"
    $GH_CMD label create "component:ui" --color "#7057FF" --description "Related to UI functionality"
    $GH_CMD label create "component:chembl" --color "#D876E3" --description "Related to ChEMBL integration"
    $GH_CMD label create "component:authentication" --color "#008672" --description "Related to authentication system"
    $GH_CMD label create "component:testing" --color "#FBCA04" --description "Related to testing infrastructure"
    $GH_CMD label create "type:feature" --color "#C5DEF5" --description "New feature implementation"
    $GH_CMD label create "type:bugfix" --color "#D93F0B" --description "Bug fix"
    $GH_CMD label create "type:refactor" --color "#BFDADC" --description "Code refactoring"
    $GH_CMD label create "type:documentation" --color "#0E8A16" --description "Documentation update"
    $GH_CMD label create "status:implemented" --color "#C5DEF5" --description "Feature implemented"
    $GH_CMD label create "status:validated" --color "#0E8A16" --description "Implementation validated"
    echo -e "${GREEN}Labels created successfully!${NC}"
    ;;
    
  "link-issues")
    if [ -z "$2" ] || [ -z "$3" ]; then
      echo -e "${RED}Error: Two issue IDs required.${NC}"
      echo "Usage: /github-consolidate link-issues [issue-id1] [issue-id2]"
      exit 1
    fi
    echo -e "${BLUE}Linking issues #$2 and #$3...${NC}"
    $GH_CMD issue comment $2 -b "This issue is linked to #$3. Both issues will be tracked together."
    $GH_CMD issue comment $3 -b "This issue is linked to #$2. Both issues will be tracked together."
    echo -e "${GREEN}Issues linked successfully!${NC}"
    ;;
    
  "categorize")
    if [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
      echo -e "${RED}Error: Missing parameters.${NC}"
      echo "Usage: /github-consolidate categorize [issue-id] [component] [priority]"
      exit 1
    fi
    echo -e "${BLUE}Categorizing issue #$2...${NC}"
    $GH_CMD issue edit $2 --add-label "component:$3,priority:$4"
    echo -e "${GREEN}Issue #$2 categorized successfully!${NC}"
    ;;
    
  "assign-milestone")
    if [ -z "$2" ] || [ -z "$3" ]; then
      echo -e "${RED}Error: Missing parameters.${NC}"
      echo "Usage: /github-consolidate assign-milestone [issue-id] [milestone-id]"
      exit 1
    fi
    echo -e "${BLUE}Assigning issue #$2 to milestone #$3...${NC}"
    $GH_CMD issue edit $2 --milestone $3
    echo -e "${GREEN}Issue #$2 assigned to milestone #$3 successfully!${NC}"
    ;;
    
  "status")
    echo -e "${BLUE}Current repository organization status...${NC}"
    echo -e "${YELLOW}Milestones:${NC}"
    $GH_CMD api repos/:owner/:repo/milestones --jq '.[] | [.number, .title, .open_issues] | @tsv' | \
      awk -F'\t' '{printf "%-6s %-40s %s open issues\n", $1, $2, $3}'
      
    echo -e "\n${YELLOW}Labels:${NC}"
    $GH_CMD label list --limit 30
    
    echo -e "\n${YELLOW}Open issues by milestone:${NC}"
    $GH_CMD issue list --limit 20
    ;;
    
  "help"|*)
    show_help
    ;;
esac