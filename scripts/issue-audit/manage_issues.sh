#!/bin/bash
# GitHub Issue Management Hub Script
# Entry point for all GitHub issue management tasks

# Colors for terminal output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m' # No Color

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if gh CLI is installed
if ! command -v gh &> /dev/null; then
    echo -e "${RED}Error: GitHub CLI (gh) is not installed or not in PATH${NC}"
    echo "Please install it from: https://cli.github.com/"
    exit 1
fi

# Check if jq is installed
if ! command -v jq &> /dev/null; then
    echo -e "${RED}Error: jq is not installed or not in PATH${NC}"
    echo "Please install it with: sudo apt-get install jq"
    exit 1
fi

# Clear screen and show header
clear
echo -e "${BOLD}${BLUE}=======================================================${NC}"
echo -e "${BOLD}${BLUE}          CryoProtect GitHub Issue Manager${NC}"
echo -e "${BOLD}${BLUE}=======================================================${NC}"
echo ""

# Get repository information
repo_info=$(gh repo view --json name,owner,isPrivate,description,url)
repo_name=$(echo "$repo_info" | jq -r '.name')
repo_owner=$(echo "$repo_info" | jq -r '.owner.login')
repo_url=$(echo "$repo_info" | jq -r '.url')
repo_description=$(echo "$repo_info" | jq -r '.description // "No description"')
repo_type=$(echo "$repo_info" | jq -r 'if .isPrivate then "Private" else "Public" end')

echo -e "${BLUE}Repository:${NC} $repo_owner/$repo_name ($repo_type)"
echo -e "${BLUE}Description:${NC} $repo_description"
echo -e "${BLUE}URL:${NC} $repo_url"
echo ""

# Get issue statistics
open_issues_count=$(gh issue list --json number --state open | jq '. | length')
closed_issues_count=$(gh issue list --json number --state closed --limit 500 | jq '. | length')
total_issues_count=$((open_issues_count + closed_issues_count))

echo -e "${BLUE}Issue Statistics:${NC}"
echo -e "- Open Issues: ${GREEN}$open_issues_count${NC}"
echo -e "- Closed Issues: ${GREEN}$closed_issues_count${NC}"
echo -e "- Total Issues: ${GREEN}$total_issues_count${NC}"
echo ""

# Display main menu
show_main_menu() {
    echo -e "${BOLD}${BLUE}Main Menu:${NC}"
    echo "1. Audit Issues (analyze quality and identify problems)"
    echo "2. Standardize Repository (labels, milestones, templates)"
    echo "3. Clean Up Issues (based on audit results)"
    echo "4. Consolidate Issues (merge, move, bulk update)"
    echo "5. View GitHub CLI Examples"
    echo "6. TaskMaster Integration"
    echo "7. View Reports"
    echo "0. Exit"
    echo ""
    read -p "Select an option (0-7): " option
    
    case $option in
        1) 
            echo -e "${YELLOW}Launching Issue Audit...${NC}"
            "$SCRIPT_DIR/audit_github_issues.sh"
            ;;
        2) 
            echo -e "${YELLOW}Launching Repository Standardization...${NC}"
            "$SCRIPT_DIR/standardize_github_repository.sh"
            ;;
        3) 
            echo -e "${YELLOW}Launching Issue Cleanup...${NC}"
            "$SCRIPT_DIR/cleanup_github_issues.sh"
            ;;
        4) 
            echo -e "${YELLOW}Launching Issue Consolidation...${NC}"
            "$SCRIPT_DIR/consolidate_github_issues.sh"
            ;;
        5)
            echo -e "${YELLOW}Displaying GitHub CLI Examples...${NC}"
            "$SCRIPT_DIR/gh_cli_examples.sh" | less -R
            show_main_menu
            ;;
        6)
            show_taskmaster_menu
            ;;
        7)
            show_reports_menu
            ;;
        0)
            echo -e "${GREEN}Exiting Issue Manager. Goodbye!${NC}"
            exit 0
            ;;
        *)
            echo -e "${RED}Invalid option. Please try again.${NC}"
            show_main_menu
            ;;
    esac
}

# TaskMaster integration menu
show_taskmaster_menu() {
    clear
    echo -e "${BOLD}${BLUE}TaskMaster Integration:${NC}"
    echo "1. Update GitHub to TaskMaster mapping"
    echo "2. Verify PRD to GitHub issue consistency"
    echo "3. Synchronize Tasks (/sync-tasks)"
    echo "4. Generate Tasks from PRD"
    echo "5. List TaskMaster Tasks"
    echo "6. Back to Main Menu"
    echo ""
    read -p "Select an option (1-6): " option
    
    case $option in
        1)
            read -p "Enter PRD file path: " prd_file
            read -p "Enter GitHub issue number: " issue_number
            if [ -f "$prd_file" ] && [ -n "$issue_number" ]; then
                echo -e "${YELLOW}Updating GitHub issue #$issue_number to TaskMaster mapping...${NC}"
                /home/mushu/Projects/CryoProtect/.claude/commands/github-prd.sh update-mapping "$prd_file" "$issue_number"
            else
                echo -e "${RED}Invalid PRD file or issue number${NC}"
            fi
            read -p "Press Enter to continue..."
            show_taskmaster_menu
            ;;
        2)
            echo -e "${YELLOW}Verifying PRD to GitHub issue consistency...${NC}"
            /home/mushu/Projects/CryoProtect/scripts/verify_github_taskmaster_consistency.sh
            read -p "Press Enter to continue..."
            show_taskmaster_menu
            ;;
        3)
            read -p "Enter sync-tasks command (status|pull|push|commit): " sync_cmd
            if [ -n "$sync_cmd" ]; then
                echo -e "${YELLOW}Running 'sync-tasks $sync_cmd'...${NC}"
                /home/mushu/Projects/CryoProtect/.claude/commands/sync-tasks.sh "$sync_cmd"
            else
                echo -e "${RED}Invalid command${NC}"
            fi
            read -p "Press Enter to continue..."
            show_taskmaster_menu
            ;;
        4)
            read -p "Enter PRD file path: " prd_file
            if [ -f "$prd_file" ]; then
                echo -e "${YELLOW}Generating tasks from PRD '$prd_file'...${NC}"
                /home/mushu/Projects/CryoProtect/.claude/commands/github-prd.sh generate-tasks "$prd_file"
            else
                echo -e "${RED}Invalid PRD file${NC}"
            fi
            read -p "Press Enter to continue..."
            show_taskmaster_menu
            ;;
        5)
            echo -e "${YELLOW}Listing TaskMaster tasks...${NC}"
            if command -v task-master &> /dev/null; then
                task-master list
            else
                echo -e "${RED}TaskMaster not found. Make sure it's installed and in your PATH${NC}"
            fi
            read -p "Press Enter to continue..."
            show_taskmaster_menu
            ;;
        6)
            show_main_menu
            ;;
        *)
            echo -e "${RED}Invalid option. Please try again.${NC}"
            show_taskmaster_menu
            ;;
    esac
}

# Reports menu
show_reports_menu() {
    clear
    echo -e "${BOLD}${BLUE}Issue Management Reports:${NC}"
    echo "1. List Available Reports"
    echo "2. View Latest Audit Report"
    echo "3. View Latest Consolidation Report"
    echo "4. View GitHub-TaskMaster Mapping"
    echo "5. Back to Main Menu"
    echo ""
    read -p "Select an option (1-5): " option
    
    OUTPUT_DIR="$SCRIPT_DIR/output"
    
    case $option in
        1)
            echo -e "${YELLOW}Available Reports:${NC}"
            find "$OUTPUT_DIR" -name "*.md" -o -name "*.json" | sort -r
            read -p "Enter report file to view (or press Enter to go back): " report_file
            if [ -n "$report_file" ] && [ -f "$report_file" ]; then
                if [[ "$report_file" == *.md ]]; then
                    # Display markdown files
                    cat "$report_file" | less
                elif [[ "$report_file" == *.json ]]; then
                    # Display JSON files with formatting
                    jq '.' "$report_file" | less
                fi
            fi
            read -p "Press Enter to continue..."
            show_reports_menu
            ;;
        2)
            latest_audit=$(find "$OUTPUT_DIR" -name "issue_audit_summary.md" | sort -r | head -1)
            if [ -n "$latest_audit" ] && [ -f "$latest_audit" ]; then
                cat "$latest_audit" | less
            else
                echo -e "${RED}No audit reports found${NC}"
            fi
            read -p "Press Enter to continue..."
            show_reports_menu
            ;;
        3)
            latest_consolidation=$(find "$OUTPUT_DIR" -name "consolidation_report_*.md" | sort -r | head -1)
            if [ -n "$latest_consolidation" ] && [ -f "$latest_consolidation" ]; then
                cat "$latest_consolidation" | less
            else
                echo -e "${RED}No consolidation reports found${NC}"
            fi
            read -p "Press Enter to continue..."
            show_reports_menu
            ;;
        4)
            github_mapping="/home/mushu/Projects/CryoProtect/tasks/github_mapping.json"
            if [ -f "$github_mapping" ]; then
                echo -e "${YELLOW}GitHub to TaskMaster Mapping:${NC}"
                jq '.' "$github_mapping" | less
            else
                echo -e "${RED}GitHub mapping file not found${NC}"
            fi
            read -p "Press Enter to continue..."
            show_reports_menu
            ;;
        5)
            show_main_menu
            ;;
        *)
            echo -e "${RED}Invalid option. Please try again.${NC}"
            show_reports_menu
            ;;
    esac
}

# Run the main menu
show_main_menu