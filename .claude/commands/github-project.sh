#!/bin/bash
# Script to manage GitHub project status and provide project-wide insights

# Set terminal colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

function show_help {
  echo -e "${BLUE}GitHub Project Status Commands${NC}"
  echo -e "Usage: /github-project [command] [options]"
  echo
  echo -e "Commands:"
  echo -e "  ${GREEN}status${NC}         Show overall project status (default)"
  echo -e "  ${GREEN}milestone-progress${NC} [milestone-id] Show progress of a milestone"
  echo -e "  ${GREEN}component-status${NC} [component] Show status of a component"
  echo -e "  ${GREEN}sprint-summary${NC}  Summarize current sprint status"
  echo -e "  ${GREEN}burndown${NC}       Generate a burndown report"
  echo -e "  ${GREEN}priorities${NC}     List high priority outstanding issues"
  echo -e "  ${GREEN}health${NC}         Check repository health"
  echo
  echo -e "Examples:"
  echo -e "  /github-project status"
  echo -e "  /github-project milestone-progress 1"
  echo -e "  /github-project component-status database"
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
  "status")
    echo -e "${BLUE}Project Status Overview${NC}"
    
    # Show milestone summary
    echo -e "\n${YELLOW}Milestone Progress:${NC}"
    $GH_CMD api repos/:owner/:repo/milestones --jq '.[] | [.number, .title, .open_issues, .closed_issues] | @tsv' | \
      awk -F'\t' '{
        if ($3 + $4 > 0) {
          percent = $4 / ($3 + $4) * 100;
        } else {
          percent = 0;
        }
        printf "%-6s %-40s %d/%d (%.1f%% complete)\n", $1, $2, $4, ($3 + $4), percent;
      }'
    
    # Show component summary
    echo -e "\n${YELLOW}Component Status:${NC}"
    COMPONENTS=("database" "api" "ui" "chembl" "authentication" "testing")
    for component in "${COMPONENTS[@]}"; do
      open_count=$($GH_CMD issue list --label "component:$component" --state open --json number | jq '. | length')
      closed_count=$($GH_CMD issue list --label "component:$component" --state closed --json number | jq '. | length')
      total=$((open_count + closed_count))
      if [ $total -gt 0 ]; then
        percent=$((closed_count * 100 / total))
      else
        percent=0
      fi
      printf "%-15s %d/%d (%d%% complete)\n" "$component" "$closed_count" "$total" "$percent"
    done
    
    # Show priority summary
    echo -e "\n${YELLOW}Priority Summary:${NC}"
    high_count=$($GH_CMD issue list --label "priority:high" --state open --json number | jq '. | length')
    medium_count=$($GH_CMD issue list --label "priority:medium" --state open --json number | jq '. | length')
    low_count=$($GH_CMD issue list --label "priority:low" --state open --json number | jq '. | length')
    echo "High Priority: $high_count open issues"
    echo "Medium Priority: $medium_count open issues"
    echo "Low Priority: $low_count open issues"
    
    # Show recent activity
    echo -e "\n${YELLOW}Recent Activity (last 7 days):${NC}"
    $GH_CMD issue list --search "updated:>$(date -d '7 days ago' '+%Y-%m-%d')" --limit 10
    ;;
    
  "milestone-progress")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No milestone ID specified.${NC}"
      echo "Usage: /github-project milestone-progress [milestone-id]"
      exit 1
    fi
    
    echo -e "${BLUE}Progress for milestone #$2${NC}"
    
    # Get milestone details
    milestone=$($GH_CMD api repos/:owner/:repo/milestones/$2 --jq '.title, .description, .open_issues, .closed_issues')
    read -r title description open_issues closed_issues <<< "$milestone"
    
    total_issues=$((open_issues + closed_issues))
    if [ $total_issues -gt 0 ]; then
      percent=$((closed_issues * 100 / total_issues))
    else
      percent=0
    fi
    
    echo -e "${YELLOW}$title${NC}"
    echo "$description"
    echo -e "Progress: $closed_issues/$total_issues ($percent% complete)"
    
    # Show open issues in milestone
    echo -e "\n${YELLOW}Open issues in this milestone:${NC}"
    $GH_CMD issue list --milestone $2 --state open
    ;;
    
  "component-status")
    if [ -z "$2" ]; then
      echo -e "${RED}Error: No component specified.${NC}"
      echo "Usage: /github-project component-status [component]"
      exit 1
    fi
    
    echo -e "${BLUE}Status for component: $2${NC}"
    
    # Get open issues for the component
    echo -e "\n${YELLOW}Open issues:${NC}"
    $GH_CMD issue list --label "component:$2" --state open
    
    # Get recently closed issues for the component
    echo -e "\n${YELLOW}Recently closed issues (last 14 days):${NC}"
    $GH_CMD issue list --label "component:$2" --state closed --search "closed:>$(date -d '14 days ago' '+%Y-%m-%d')"
    ;;
    
  "sprint-summary")
    echo -e "${BLUE}Current Sprint Summary${NC}"
    
    # Get sprint issues (assuming current milestone is the sprint)
    current_milestone=$($GH_CMD api repos/:owner/:repo/milestones --jq '.[0].number')
    
    echo -e "\n${YELLOW}Sprint progress:${NC}"
    $GH_CMD api repos/:owner/:repo/milestones/$current_milestone --jq '[.title, .open_issues, .closed_issues] | @tsv' | \
      awk -F'\t' '{
        if ($2 + $3 > 0) {
          percent = $3 / ($2 + $3) * 100;
        } else {
          percent = 0;
        }
        printf "%s: %d/%d (%.1f%% complete)\n", $1, $3, ($2 + $3), percent;
      }'
    
    # Show remaining issues
    echo -e "\n${YELLOW}Remaining issues in sprint:${NC}"
    $GH_CMD issue list --milestone $current_milestone --state open
    ;;
    
  "burndown")
    echo -e "${BLUE}Burndown Report${NC}"
    
    # Get open and closed issues counts for the last 14 days
    echo -e "\n${YELLOW}Issue closure in the last 14 days:${NC}"
    for i in {14..0}; do
      date=$(date -d "$i days ago" '+%Y-%m-%d')
      closed=$($GH_CMD issue list --search "closed:$date" --state closed --json number | jq '. | length')
      created=$($GH_CMD issue list --search "created:$date" --json number | jq '. | length')
      echo "$date: +$created created, -$closed closed"
    done
    ;;
    
  "priorities")
    echo -e "${BLUE}High Priority Outstanding Issues${NC}"
    
    # Get high priority open issues
    echo -e "\n${YELLOW}High priority issues:${NC}"
    $GH_CMD issue list --label "priority:high" --state open
    ;;
    
  "health")
    echo -e "${BLUE}Repository Health Check${NC}"
    
    # Check PRs needing review
    echo -e "\n${YELLOW}PRs needing review:${NC}"
    $GH_CMD pr list --search "review:required"
    
    # Check stale issues
    echo -e "\n${YELLOW}Stale issues (no activity for 30+ days):${NC}"
    $GH_CMD issue list --search "updated:<$(date -d '30 days ago' '+%Y-%m-%d')" --state open --limit 10
    
    # Check rate of issue closure
    closed_last_month=$($GH_CMD issue list --search "closed:>$(date -d '30 days ago' '+%Y-%m-%d')" --state closed --json number | jq '. | length')
    created_last_month=$($GH_CMD issue list --search "created:>$(date -d '30 days ago' '+%Y-%m-%d')" --json number | jq '. | length')
    echo -e "\n${YELLOW}Issue flow (last 30 days):${NC}"
    echo "Created: $created_last_month issues"
    echo "Closed: $closed_last_month issues"
    net=$((created_last_month - closed_last_month))
    if [ $net -gt 0 ]; then
      echo "Net: +$net issues (increasing backlog)"
    else
      echo "Net: $net issues (reducing backlog)"
    fi
    ;;
    
  "help"|*)
    show_help
    ;;
esac