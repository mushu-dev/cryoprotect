#!/bin/bash
# sync-tasks.sh - Helper script for synchronizing TaskMaster tasks between environments
# This command is automatically run when invoked with /sync-tasks

# Colors for terminal output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Ensure tasks directory exists
mkdir -p tasks

# Display help information
show_help() {
  echo -e "${BLUE}Task Synchronization Helper${NC}"
  echo ""
  echo "Usage: $0 [command]"
  echo ""
  echo "Commands:"
  echo "  status               - Show task synchronization status"
  echo "  pull                 - Pull latest task changes from git"
  echo "  push                 - Commit and push task changes to git"
  echo "  list                 - List all tasks (shortcut to task-master list)"
  echo "  move [id] [cursor|claude] - Move a task to Cursor or Claude environment"
  echo "  commit               - Commit task files to git with automatic message"
  echo "  help                 - Show this help message"
  echo ""
  echo "Examples:"
  echo "  $0 status"
  echo "  $0 push"
  echo "  $0 commit"
  echo ""
}

# Show task synchronization status
show_status() {
  echo -e "${BLUE}Task Synchronization Status${NC}"
  echo ""
  
  # Check if tasks directory exists
  if [ ! -d "tasks" ]; then
    echo -e "${YELLOW}Warning: No tasks directory found${NC}"
    echo "You need to initialize TaskMaster first with: task-master init"
    return 1
  fi
  
  # Check if tasks.json exists
  if [ ! -f "tasks/tasks.json" ]; then
    echo -e "${YELLOW}Warning: No tasks.json file found${NC}"
    echo "You need to create tasks first with: task-master parse-prd"
    return 1
  fi
  
  # Count tasks and task files
  task_count=$(jq '.tasks | length' tasks/tasks.json 2>/dev/null || echo "0")
  file_count=$(find tasks -name "task_*.txt" | wc -l)
  
  echo -e "Tasks in tasks.json: ${GREEN}${task_count}${NC}"
  echo -e "Task files found: ${GREEN}${file_count}${NC}"
  
  # Check git status for task files
  echo -e "\n${BLUE}Git Status for Task Files:${NC}"
  git status --porcelain tasks/
  
  # Show last modification time
  echo -e "\n${BLUE}Last Modified:${NC}"
  stat -c "tasks.json: %y" tasks/tasks.json 2>/dev/null || echo "tasks.json: Not found"
  
  # Show next task recommendation
  echo -e "\n${BLUE}Next Task Recommendation:${NC}"
  task-master next | grep -A 1 "Next Task to Work On:" | head -2
  
  echo -e "\n${YELLOW}Use 'sync-tasks commit' to save task changes before switching environments${NC}"
}

# Pull task changes from git
pull_tasks() {
  echo -e "${BLUE}Pulling latest task changes from git...${NC}"
  
  # Stash any local changes first
  git stash push -m "Auto-stash before pulling task changes" -- tasks/
  
  # Pull from git
  git pull
  
  # Pop stash if there was anything stashed
  if git stash list | grep -q "Auto-stash before pulling task changes"; then
    git stash pop
  fi
  
  echo -e "${GREEN}Task files synchronized from git${NC}"
  echo -e "${YELLOW}Run 'task-master list' to see current tasks${NC}"
}

# Commit and push task changes to git
push_tasks() {
  echo -e "${BLUE}Committing and pushing task changes to git...${NC}"
  
  # Check if there are changes to commit
  if [ -z "$(git status --porcelain tasks/)" ]; then
    echo -e "${YELLOW}No task changes to commit${NC}"
    return 0
  fi
  
  # Generate commit message
  task_count=$(jq '.tasks | length' tasks/tasks.json 2>/dev/null || echo "0")
  commit_msg="Update TaskMaster tasks (${task_count} tasks total)"
  
  # Add and commit changes
  git add tasks/
  git commit -m "$commit_msg"
  
  # Push to remote if there's a remote branch
  if git remote -v | grep -q origin; then
    git push
    echo -e "${GREEN}Task changes committed and pushed to git${NC}"
  else
    echo -e "${GREEN}Task changes committed locally${NC}"
    echo -e "${YELLOW}No remote repository found for pushing${NC}"
  fi
}

# List all tasks (wrapper for task-master list)
list_tasks() {
  task-master list
}

# Move a task between environments
move_task() {
  if [ -z "$1" ] || [ -z "$2" ]; then
    echo -e "${RED}Error: Missing task ID or destination${NC}"
    echo "Usage: $0 move [id] [cursor|claude]"
    return 1
  fi
  
  task_id="$1"
  destination="$2"
  
  if [ "$destination" != "cursor" ] && [ "$destination" != "claude" ]; then
    echo -e "${RED}Error: Destination must be 'cursor' or 'claude'${NC}"
    return 1
  fi
  
  # Check if task exists
  if ! jq -e ".tasks[] | select(.id == $task_id)" tasks/tasks.json &>/dev/null; then
    echo -e "${RED}Error: Task #${task_id} not found${NC}"
    return 1
  fi
  
  # For now, this is just a tagging mechanism since both environments use the same files
  # We could expand this in the future to add environment-specific metadata
  echo -e "${GREEN}Task #${task_id} moved to ${destination} environment${NC}"
  echo -e "${YELLOW}Both environments share the same task files, so this is just a logical move${NC}"
  echo -e "To work on this task in ${destination}, run: task-master show ${task_id}"
}

# Commit task changes to git
commit_tasks() {
  echo -e "${BLUE}Committing task changes to git...${NC}"
  
  # Check if there are changes to commit
  if [ -z "$(git status --porcelain tasks/)" ]; then
    echo -e "${YELLOW}No task changes to commit${NC}"
    return 0
  fi
  
  # Get task stats for commit message
  task_count=$(jq '.tasks | length' tasks/tasks.json 2>/dev/null || echo "0")
  pending_count=$(jq '.tasks | map(select(.status == "pending")) | length' tasks/tasks.json 2>/dev/null || echo "0")
  in_progress_count=$(jq '.tasks | map(select(.status == "in-progress")) | length' tasks/tasks.json 2>/dev/null || echo "0")
  done_count=$(jq '.tasks | map(select(.status == "done")) | length' tasks/tasks.json 2>/dev/null || echo "0")
  
  # Generate commit message
  commit_msg="Update TaskMaster tasks: ${task_count} total (${pending_count} pending, ${in_progress_count} in progress, ${done_count} done)"
  
  # Add and commit changes
  git add tasks/
  git commit -m "$commit_msg"
  
  echo -e "${GREEN}Task changes committed to git${NC}"
  echo -e "${YELLOW}Run 'sync-tasks push' to push the changes to the remote repository${NC}"
}

# Main command processing
case "$1" in
  status)
    show_status
    ;;
  pull)
    pull_tasks
    ;;
  push)
    push_tasks
    ;;
  list)
    list_tasks
    ;;
  move)
    move_task "$2" "$3"
    ;;
  commit)
    commit_tasks
    ;;
  help|*)
    show_help
    ;;
esac