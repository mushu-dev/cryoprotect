#!/bin/bash
# Script to manage GitHub pull requests for CryoProtect

# Set terminal colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

function show_help {
  echo -e "${BLUE}GitHub Pull Request Management Commands${NC}"
  echo -e "Usage: /github-prs [command] [options]"
  echo
  echo -e "Commands:"
  echo -e "  ${GREEN}list${NC}           List open pull requests (default)"
  echo -e "  ${GREEN}list-all${NC}       List all pull requests (open and closed)"
  echo -e "  ${GREEN}list-mine${NC}      List pull requests created by you"
  echo -e "  ${GREEN}view${NC} [id]      View details of a specific pull request"
  echo -e "  ${GREEN}create${NC}         Create a new pull request"
  echo -e "  ${GREEN}checkout${NC} [id]  Checkout the branch of a pull request"
  echo -e "  ${GREEN}checks${NC} [id]    View CI status checks for a pull request"
  echo -e "  ${GREEN}diff${NC} [id]      View changes in a pull request"
  echo -e "  ${GREEN}review${NC} [id]    Add a review to a pull request"
  echo -e "  ${GREEN}merge${NC} [id]     Merge a pull request"
  echo -e "  ${GREEN}close${NC} [id]     Close a pull request without merging"
  echo -e "  ${GREEN}reopen${NC} [id]    Reopen a closed pull request"
  echo
  echo -e "Examples:"
  echo -e "  /github-prs list"
  echo -e "  /github-prs view 123"
  echo -e "  /github-prs diff 123"
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
    echo -e "${BLUE}Listing open pull requests...${NC}"
    $GH_CMD pr list
    ;;
    
  "list-all")
    echo -e "${BLUE}Listing all pull requests...${NC}"
    $GH_CMD pr list --state all
    ;;
    
  "list-mine")
    echo -e "${BLUE}Listing pull requests created by you...${NC}"
    $GH_CMD pr list --author @me
    ;;
    
  "view")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No PR ID specified.${NC}"
      echo "Usage: /github-prs view [pr-id]"
      exit 1
    fi
    echo -e "${BLUE}Viewing pull request #$2...${NC}"
    $GH_CMD pr view $2
    ;;
    
  "create")
    echo -e "${BLUE}Creating a new pull request...${NC}"
    echo -e "${YELLOW}You'll be prompted to enter the PR details.${NC}"
    $GH_CMD pr create
    ;;
    
  "checkout")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No PR ID specified.${NC}"
      echo "Usage: /github-prs checkout [pr-id]"
      exit 1
    fi
    echo -e "${BLUE}Checking out pull request #$2...${NC}"
    $GH_CMD pr checkout $2
    ;;
    
  "checks")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No PR ID specified.${NC}"
      echo "Usage: /github-prs checks [pr-id]"
      exit 1
    fi
    echo -e "${BLUE}Viewing CI checks for pull request #$2...${NC}"
    $GH_CMD pr checks $2
    ;;
    
  "diff")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No PR ID specified.${NC}"
      echo "Usage: /github-prs diff [pr-id]"
      exit 1
    fi
    echo -e "${BLUE}Viewing diff for pull request #$2...${NC}"
    $GH_CMD pr diff $2
    ;;
    
  "review")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No PR ID specified.${NC}"
      echo "Usage: /github-prs review [pr-id]"
      exit 1
    fi
    echo -e "${BLUE}Adding review to pull request #$2...${NC}"
    $GH_CMD pr review $2
    ;;
    
  "merge")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No PR ID specified.${NC}"
      echo "Usage: /github-prs merge [pr-id]"
      exit 1
    fi
    echo -e "${BLUE}Merging pull request #$2...${NC}"
    $GH_CMD pr merge $2
    ;;
    
  "close")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No PR ID specified.${NC}"
      echo "Usage: /github-prs close [pr-id]"
      exit 1
    fi
    echo -e "${BLUE}Closing pull request #$2 without merging...${NC}"
    $GH_CMD pr close $2
    ;;
    
  "reopen")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No PR ID specified.${NC}"
      echo "Usage: /github-prs reopen [pr-id]"
      exit 1
    fi
    echo -e "${BLUE}Reopening pull request #$2...${NC}"
    $GH_CMD pr reopen $2
    ;;
    
  "help"|*)
    show_help
    ;;
esac