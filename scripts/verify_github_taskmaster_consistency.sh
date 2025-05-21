#!/bin/bash
# GitHub Issue to TaskMaster Consistency Verification
# This script checks the consistency between GitHub issues, PRDs, and TaskMaster tasks

# Colors for terminal output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Key directories
PRD_DIR="/home/mushu/Projects/CryoProtect/scripts"
TASK_DIR="/home/mushu/Projects/CryoProtect/tasks"
MAPPING_FILE="${TASK_DIR}/github_mapping.json"
OUTPUT_DIR="/home/mushu/Projects/CryoProtect/scripts/issue-audit/output"
mkdir -p "$OUTPUT_DIR"

echo -e "${BLUE}GitHub-TaskMaster Consistency Verification${NC}"
echo "==========================================="

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

# Check for mapping file
if [ ! -f "$MAPPING_FILE" ]; then
    echo -e "${YELLOW}Warning: GitHub to TaskMaster mapping file not found${NC}"
    echo "Creating an empty mapping file"
    echo '{
  "mappings": [],
  "last_updated": "'"$(date -Iseconds)"'"
}' > "$MAPPING_FILE"
fi

# Get repository information
repo_info=$(gh repo view --json name,owner,isPrivate,description)
repo_name=$(echo "$repo_info" | jq -r '.name')
repo_owner=$(echo "$repo_info" | jq -r '.owner.login')

echo -e "${BLUE}Repository:${NC} $repo_owner/$repo_name"
echo ""

# Get GitHub issues
echo -e "${BLUE}Retrieving GitHub issues...${NC}"
github_issues=$(gh issue list --limit 500 --state all --json number,title,url,state)
issue_count=$(echo "$github_issues" | jq '. | length')
echo -e "${GREEN}Found $issue_count GitHub issues${NC}"

# Get PRD files
echo -e "${BLUE}Scanning for PRD files...${NC}"
prd_files=$(find "$PRD_DIR" -name "prd_*.txt")
prd_count=$(echo "$prd_files" | wc -l)
echo -e "${GREEN}Found $prd_count PRD files${NC}"

# Get TaskMaster tasks
echo -e "${BLUE}Retrieving TaskMaster tasks...${NC}"
if [ -f "${TASK_DIR}/tasks.json" ]; then
    taskmaster_tasks=$(cat "${TASK_DIR}/tasks.json")
    task_count=$(echo "$taskmaster_tasks" | jq '.tasks | length')
    echo -e "${GREEN}Found $task_count TaskMaster tasks${NC}"
else
    echo -e "${YELLOW}Warning: TaskMaster tasks.json not found${NC}"
    taskmaster_tasks='{}'
    task_count=0
fi

# Get GitHub-TaskMaster mapping
echo -e "${BLUE}Checking GitHub-TaskMaster mapping...${NC}"
github_mappings=$(cat "$MAPPING_FILE" | jq '.mappings')
mapping_count=$(echo "$github_mappings" | jq '. | length')
echo -e "${GREEN}Found $mapping_count GitHub-TaskMaster mappings${NC}"

echo ""
echo -e "${BLUE}Checking consistency...${NC}"

# Initialize counters
orphaned_prds=0
missing_github_issues=0
orphaned_tasks=0
invalid_mappings=0
inconsistent_github_issues=0

# Initialize output files
verification_report="$OUTPUT_DIR/github_taskmaster_verification_$(date +%Y%m%d_%H%M%S).md"
issues_file="$OUTPUT_DIR/inconsistent_issues_$(date +%Y%m%d_%H%M%S).json"
prds_file="$OUTPUT_DIR/orphaned_prds_$(date +%Y%m%d_%H%M%S).json"
tasks_file="$OUTPUT_DIR/orphaned_tasks_$(date +%Y%m%d_%H%M%S).json"
mappings_file="$OUTPUT_DIR/invalid_mappings_$(date +%Y%m%d_%H%M%S).json"

# Create report header
cat > "$verification_report" << EOL
# GitHub-TaskMaster Consistency Verification Report

Repository: **$repo_owner/$repo_name**
Generated on: **$(date)**

## Overview

- **GitHub Issues:** $issue_count
- **PRD Files:** $prd_count
- **TaskMaster Tasks:** $task_count
- **GitHub-TaskMaster Mappings:** $mapping_count

## Consistency Check Results

EOL

# Create initial JSON arrays
echo "[]" > "$issues_file"
echo "[]" > "$prds_file"
echo "[]" > "$tasks_file"
echo "[]" > "$mappings_file"

# Check PRDs for GitHub issue references
echo -e "\n${BLUE}Checking PRDs for GitHub issue references...${NC}"
for prd_file in $prd_files; do
    prd_basename=$(basename "$prd_file")
    
    # Check for GitHub issue reference
    issue_number=$(grep -o "GitHub issue #[0-9]\+" "$prd_file" | grep -o "[0-9]\+")
    related_issue=$(grep -o "Related Issue: #[0-9]\+" "$prd_file" | grep -o "[0-9]\+")
    
    # Use issue_number if found, otherwise use related_issue
    check_issue=${issue_number:-$related_issue}
    
    if [ -z "$check_issue" ]; then
        echo -e "${YELLOW}Warning: PRD $prd_basename has no GitHub issue reference${NC}"
        orphaned_prds=$((orphaned_prds + 1))
        
        # Get PRD title
        prd_title=$(grep -m 1 "^# " "$prd_file" | sed 's/^# //')
        
        # Add to JSON
        temp_file=$(mktemp)
        jq --arg file "$prd_basename" \
           --arg path "$prd_file" \
           --arg title "$prd_title" \
           '. + [{"file": $file, "path": $path, "title": $title, "issue": null}]' \
           "$prds_file" > "$temp_file"
        mv "$temp_file" "$prds_file"
        continue
    fi
    
    # Check if referenced issue exists
    if ! echo "$github_issues" | jq -e ".[] | select(.number == $check_issue)" > /dev/null; then
        echo -e "${YELLOW}Warning: PRD $prd_basename references non-existent issue #$check_issue${NC}"
        missing_github_issues=$((missing_github_issues + 1))
        
        # Get PRD title
        prd_title=$(grep -m 1 "^# " "$prd_file" | sed 's/^# //')
        
        # Add to JSON
        temp_file=$(mktemp)
        jq --arg file "$prd_basename" \
           --arg path "$prd_file" \
           --arg issue "$check_issue" \
           --arg title "$prd_title" \
           '. + [{"file": $file, "path": $path, "title": $title, "issue": $issue}]' \
           "$prds_file" > "$temp_file"
        mv "$temp_file" "$prds_file"
    fi
done

# Check TaskMaster tasks that are not in any mapping
if [ "$task_count" -gt 0 ]; then
    echo -e "\n${BLUE}Checking for orphaned TaskMaster tasks...${NC}"
    
    # Get all task IDs from TaskMaster
    task_ids=$(echo "$taskmaster_tasks" | jq -r '.tasks[].id')
    
    # Get all mapped task IDs
    mapped_task_ids=$(echo "$github_mappings" | jq -r '.[].taskmaster_tasks[]? | tostring')
    
    # Find orphaned tasks
    for task_id in $task_ids; do
        if ! echo "$mapped_task_ids" | grep -q "^$task_id$"; then
            echo -e "${YELLOW}Warning: TaskMaster task $task_id is not mapped to any GitHub issue${NC}"
            orphaned_tasks=$((orphaned_tasks + 1))
            
            # Get task details
            task_title=$(echo "$taskmaster_tasks" | jq -r ".tasks[] | select(.id == $task_id) | .title")
            
            # Add to JSON
            temp_file=$(mktemp)
            jq --arg id "$task_id" \
               --arg title "$task_title" \
               '. + [{"id": $id, "title": $title}]' \
               "$tasks_file" > "$temp_file"
            mv "$temp_file" "$tasks_file"
        fi
    done
fi

# Check GitHub-TaskMaster mappings
echo -e "\n${BLUE}Checking GitHub-TaskMaster mappings...${NC}"
echo "$github_mappings" | jq -c '.[]' | while read -r mapping; do
    issue_number=$(echo "$mapping" | jq -r '.github_issue')
    prd_path=$(echo "$mapping" | jq -r '.prd_path')
    task_ids=$(echo "$mapping" | jq -r '.taskmaster_tasks[]? | tostring' 2>/dev/null)
    
    mapping_valid=true
    mapping_issues=()
    
    # Check if GitHub issue exists
    if ! echo "$github_issues" | jq -e ".[] | select(.number == $issue_number)" > /dev/null; then
        echo -e "${YELLOW}Warning: Mapping references non-existent GitHub issue #$issue_number${NC}"
        mapping_valid=false
        mapping_issues+=("GitHub issue #$issue_number not found")
    fi
    
    # Check if PRD file exists
    if [ ! -f "$prd_path" ]; then
        echo -e "${YELLOW}Warning: Mapping references non-existent PRD file $prd_path${NC}"
        mapping_valid=false
        mapping_issues+=("PRD file $prd_path not found")
    fi
    
    # Check if TaskMaster tasks exist
    if [ -n "$task_ids" ]; then
        for task_id in $task_ids; do
            if [ "$task_count" -gt 0 ]; then
                if ! echo "$taskmaster_tasks" | jq -e ".tasks[] | select(.id == $task_id)" > /dev/null; then
                    echo -e "${YELLOW}Warning: Mapping references non-existent TaskMaster task $task_id${NC}"
                    mapping_valid=false
                    mapping_issues+=("TaskMaster task $task_id not found")
                fi
            fi
        done
    fi
    
    if [ "$mapping_valid" = false ]; then
        invalid_mappings=$((invalid_mappings + 1))
        
        # Add to JSON
        temp_file=$(mktemp)
        jq --argjson mapping "$mapping" \
           --argjson issues "$(echo "${mapping_issues[@]}" | jq -R -s 'split("\n") | map(select(length > 0))')" \
           '. + [{"mapping": $mapping, "issues": $issues}]' \
           "$mappings_file" > "$temp_file"
        mv "$temp_file" "$mappings_file"
    fi
done

# Check GitHub issues with PRD references for consistency
echo -e "\n${BLUE}Checking GitHub issues with PRD references for consistency...${NC}"
echo "$github_issues" | jq -c '.[]' | while read -r issue; do
    issue_number=$(echo "$issue" | jq -r '.number')
    issue_title=$(echo "$issue" | jq -r '.title')
    issue_url=$(echo "$issue" | jq -r '.url')
    
    # Get full issue body to check for PRD references
    issue_body=$(gh issue view "$issue_number" --json body | jq -r '.body')
    
    # Check for PRD references in body
    if echo "$issue_body" | grep -q "PRD:" || echo "$issue_body" | grep -q "Product Requirements Document"; then
        # Issue mentions a PRD, check if there's a mapping
        if ! echo "$github_mappings" | jq -e ".[] | select(.github_issue == $issue_number)" > /dev/null; then
            echo -e "${YELLOW}Warning: GitHub issue #$issue_number mentions a PRD but has no mapping${NC}"
            inconsistent_github_issues=$((inconsistent_github_issues + 1))
            
            # Add to JSON
            temp_file=$(mktemp)
            jq --arg number "$issue_number" \
               --arg title "$issue_title" \
               --arg url "$issue_url" \
               --arg reason "mentions PRD but has no mapping" \
               '. + [{"number": $number, "title": $title, "url": $url, "reason": $reason}]' \
               "$issues_file" > "$temp_file"
            mv "$temp_file" "$issues_file"
        fi
    fi
done

# Add results to report
cat >> "$verification_report" << EOL
### Summary

- **Orphaned PRDs:** $orphaned_prds (PRDs without GitHub issue references)
- **Missing GitHub Issues:** $missing_github_issues (PRDs reference non-existent issues)
- **Orphaned TaskMaster Tasks:** $orphaned_tasks (Tasks not mapped to any GitHub issue)
- **Invalid Mappings:** $invalid_mappings (Mappings with missing components)
- **Inconsistent GitHub Issues:** $inconsistent_github_issues (Issues mention PRDs but lack mappings)

### Orphaned PRDs

PRDs without valid GitHub issue references:

$(cat "$prds_file" | jq -r '.[] | "- **" + .file + "**: " + .title')

### Orphaned TaskMaster Tasks

TaskMaster tasks not mapped to any GitHub issue:

$(cat "$tasks_file" | jq -r '.[] | "- **Task " + .id + "**: " + .title')

### Invalid Mappings

Mappings with issues:

$(cat "$mappings_file" | jq -r '.[] | "- **GitHub Issue #" + (.mapping.github_issue | tostring) + "**: " + (.issues | join(", "))')

### Inconsistent GitHub Issues

GitHub issues that need mapping updates:

$(cat "$issues_file" | jq -r '.[] | "- **Issue #" + .number + "**: " + .title + " (" + .reason + ")"')

## Recommendations

1. **Fix Orphaned PRDs**: Add GitHub issue references to PRDs that don't have them
2. **Fix Missing GitHub Issues**: Create missing GitHub issues or update PRD references
3. **Map Orphaned Tasks**: Add TaskMaster tasks to appropriate GitHub issue mappings
4. **Fix Invalid Mappings**: Update mappings with correct GitHub issues, PRDs, and TaskMaster tasks
5. **Create Missing Mappings**: Add mappings for GitHub issues that mention PRDs

## Next Steps

Run these commands to fix the identified issues:

1. Add GitHub issue reference to PRD:
   - Edit the PRD file and add: "## GitHub Issue Reference: This PRD is based on GitHub issue #XXX"

2. Create a new mapping:
   - \`/github-prd update-mapping [prd-file] [issue-number]\`

3. Update TaskMaster tasks:
   - \`task-master set-reference --id=[task-id] --reference="GH-[issue-number]"\`

4. Verify fixes:
   - Run this verification script again
EOL

# Display summary
echo -e "\n${BLUE}Verification Summary:${NC}"
echo "- Orphaned PRDs: $orphaned_prds"
echo "- Missing GitHub Issues: $missing_github_issues"
echo "- Orphaned TaskMaster Tasks: $orphaned_tasks"
echo "- Invalid Mappings: $invalid_mappings"
echo "- Inconsistent GitHub Issues: $inconsistent_github_issues"

echo -e "\n${GREEN}Verification report generated: $verification_report${NC}"
echo -e "${YELLOW}Review the report for details and recommended actions${NC}"

# Cleanup empty JSON files
for file in "$issues_file" "$prds_file" "$tasks_file" "$mappings_file"; do
    if [ "$(jq '. | length' "$file")" -eq 0 ]; then
        rm "$file"
    fi
done