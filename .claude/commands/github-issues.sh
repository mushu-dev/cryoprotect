#!/bin/bash
# Script to manage GitHub issues for CryoProtect

# Set terminal colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

function show_help {
  echo -e "${BLUE}GitHub Issue Management Commands${NC}"
  echo -e "Usage: /github-issues [command] [options]"
  echo
  echo -e "Commands:"
  echo -e "  ${GREEN}list${NC}           List open issues (default)"
  echo -e "  ${GREEN}list-all${NC}       List all issues (open and closed)"
  echo -e "  ${GREEN}list-mine${NC}      List issues assigned to you"
  echo -e "  ${GREEN}list-milestone${NC} List issues in a specific milestone"
  echo -e "  ${GREEN}view${NC} [id]      View details of a specific issue"
  echo -e "  ${GREEN}create${NC}         Create a new issue"
  echo -e "  ${GREEN}edit${NC} [id]      Edit an existing issue"
  echo -e "  ${GREEN}close${NC} [id]     Close an issue"
  echo -e "  ${GREEN}reopen${NC} [id]    Reopen a closed issue"
  echo -e "  ${GREEN}comment${NC} [id]   Add a comment to an issue"
  echo -e "  ${GREEN}milestones${NC}     List all milestones"
  echo -e "  ${GREEN}labels${NC}         List all labels"
  echo
  echo -e "Examples:"
  echo -e "  /github-issues list"
  echo -e "  /github-issues view 123"
  echo -e "  /github-issues list-milestone \"Phase 1: Technical Foundation\""
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
CMD=${1:-list}

case $CMD in
  "list")
    echo -e "${BLUE}Listing open issues...${NC}"
    $GH_CMD issue list
    ;;
    
  "list-all")
    echo -e "${BLUE}Listing all issues...${NC}"
    $GH_CMD issue list --state all
    ;;
    
  "list-mine")
    echo -e "${BLUE}Listing issues assigned to you...${NC}"
    $GH_CMD issue list --assignee @me
    ;;
    
  "list-milestone")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No milestone specified.${NC}"
      echo "Usage: /github-issues list-milestone \"Milestone Name\""
      exit 1
    fi
    echo -e "${BLUE}Listing issues in milestone '$2'...${NC}"
    $GH_CMD issue list --milestone "$2"
    ;;
    
  "view")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No issue ID specified.${NC}"
      echo "Usage: /github-issues view [issue-id]"
      exit 1
    fi
    echo -e "${BLUE}Viewing issue #$2...${NC}"
    $GH_CMD issue view $2
    ;;
    
  "create")
    echo -e "${BLUE}Creating a new issue...${NC}"
    echo -e "${YELLOW}You'll be prompted to enter the issue details.${NC}"
    $GH_CMD issue create
    ;;
    
  "edit")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No issue ID specified.${NC}"
      echo "Usage: /github-issues edit [issue-id]"
      exit 1
    fi
    echo -e "${BLUE}Editing issue #$2...${NC}"
    $GH_CMD issue edit $2
    ;;
    
  "close")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No issue ID specified.${NC}"
      echo "Usage: /github-issues close [issue-id]"
      exit 1
    fi
    echo -e "${BLUE}Closing issue #$2...${NC}"
    $GH_CMD issue close $2
    ;;
    
  "reopen")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No issue ID specified.${NC}"
      echo "Usage: /github-issues reopen [issue-id]"
      exit 1
    fi
    echo -e "${BLUE}Reopening issue #$2...${NC}"
    $GH_CMD issue reopen $2
    ;;
    
  "comment")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No issue ID specified.${NC}"
      echo "Usage: /github-issues comment [issue-id]"
      exit 1
    fi
    echo -e "${BLUE}Adding comment to issue #$2...${NC}"
    $GH_CMD issue comment $2
    ;;
    
  "milestones")
    echo -e "${BLUE}Listing milestones...${NC}"
    $GH_CMD api repos/:owner/:repo/milestones --jq '.[] | [.number, .title, .open_issues] | @tsv' | \
      awk -F'\t' '{printf "%-6s %-40s %s open issues\n", $1, $2, $3}'
    ;;
    
  "labels")
    echo -e "${BLUE}Listing labels...${NC}"
    $GH_CMD label list
    ;;
    
  "help"|*)
    show_help
    ;;
esac