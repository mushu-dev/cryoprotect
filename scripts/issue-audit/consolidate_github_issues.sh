#!/bin/bash
# GitHub Issue Consolidation Script
# This script helps with consolidating duplicate or related issues

# Colors for terminal output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

OUTPUT_DIR="/home/mushu/Projects/CryoProtect/scripts/issue-audit/output"
LOGS_DIR="/home/mushu/Projects/CryoProtect/scripts/issue-audit/logs"
mkdir -p "$OUTPUT_DIR" "$LOGS_DIR"

echo -e "${BLUE}CryoProtect GitHub Issue Consolidation${NC}"
echo "========================================"

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

# Get repository information
repo_info=$(gh repo view --json name,owner,isPrivate,description)
repo_name=$(echo "$repo_info" | jq -r '.name')
repo_owner=$(echo "$repo_info" | jq -r '.owner.login')

echo -e "${BLUE}Repository:${NC} $repo_owner/$repo_name"
echo ""

# Display menu of actions
echo -e "${BLUE}Available Actions:${NC}"
echo "1. Find duplicate issues by title similarity"
echo "2. Find related issues by keyword"
echo "3. Consolidate issues (merge content and close duplicates)"
echo "4. Move issues between milestones"
echo "5. Mass close issues with a comment"
echo "6. Generate consolidation report"
echo "7. Find issues by content similarity"
echo "8. Find orphaned PRD references"
echo "9. Bulk update issue titles"
echo "0. Exit"

read -p "Select action (0-9): " action

# Function to find duplicate issues by title similarity
find_duplicate_issues() {
    echo -e "${BLUE}Finding issues with similar titles...${NC}"
    
    # Get all issues
    all_issues=$(gh issue list --limit 500 --state all --json number,title,url,state)
    
    # Create a file for similarity results
    similarity_file="$OUTPUT_DIR/similar_issues_$(date +%Y%m%d_%H%M%S).json"
    
    # Extract titles for comparison
    echo "$all_issues" | jq -c '.[]' | while read -r issue1; do
        issue1_number=$(echo "$issue1" | jq -r '.number')
        issue1_title=$(echo "$issue1" | jq -r '.title')
        issue1_url=$(echo "$issue1" | jq -r '.url')
        issue1_state=$(echo "$issue1" | jq -r '.state')
        
        # Skip empty titles
        if [ -z "$issue1_title" ] || [ "$issue1_title" == "null" ]; then
            continue
        fi
        
        # Compare with all other issues
        echo "$all_issues" | jq -c '.[]' | while read -r issue2; do
            issue2_number=$(echo "$issue2" | jq -r '.number')
            issue2_title=$(echo "$issue2" | jq -r '.title')
            issue2_url=$(echo "$issue2" | jq -r '.url')
            issue2_state=$(echo "$issue2" | jq -r '.state')
            
            # Skip comparing the same issue or empty titles
            if [ "$issue1_number" -eq "$issue2_number" ] || [ -z "$issue2_title" ] || [ "$issue2_title" == "null" ]; then
                continue
            fi
            
            # Calculate similarity using Levenshtein distance
            # This is a simplified approach - more sophisticated analysis would require additional tools
            similarity=$(python3 -c "
import difflib
import sys
title1 = '''$issue1_title'''.lower()
title2 = '''$issue2_title'''.lower()
similarity = difflib.SequenceMatcher(None, title1, title2).ratio()
print(f'{similarity:.2f}')
")
            
            # If similarity is above threshold, record it
            if (( $(echo "$similarity > 0.6" | bc -l) )); then
                echo "{\"issue1\": {\"number\": $issue1_number, \"title\": \"$issue1_title\", \"url\": \"$issue1_url\", \"state\": \"$issue1_state\"}, \"issue2\": {\"number\": $issue2_number, \"title\": \"$issue2_title\", \"url\": \"$issue2_url\", \"state\": \"$issue2_state\"}, \"similarity\": $similarity}" >> "$similarity_file"
            fi
        done
    done
    
    # Format the similarity file as proper JSON
    if [ -f "$similarity_file" ]; then
        # Count lines in file
        line_count=$(wc -l < "$similarity_file")
        
        if [ "$line_count" -gt 0 ]; then
            # Create a proper JSON array
            echo "[" > "$similarity_file.tmp"
            cat "$similarity_file" | sed '$!s/$/,/' >> "$similarity_file.tmp"
            echo "]" >> "$similarity_file.tmp"
            mv "$similarity_file.tmp" "$similarity_file"
            
            # Display results
            echo -e "${GREEN}Found $(jq '. | length' "$similarity_file") potential duplicate issue pairs${NC}"
            echo -e "${BLUE}Top similar pairs:${NC}"
            jq -r 'sort_by(.similarity) | reverse | .[0:5] | .[] | "Issue #\(.issue1.number) and #\(.issue2.number) - Similarity: \(.similarity) - \(.issue1.title) / \(.issue2.title)"' "$similarity_file"
            echo -e "${YELLOW}Full results saved to $similarity_file${NC}"
        else
            echo -e "${YELLOW}No similar issues found${NC}"
            rm "$similarity_file"
        fi
    else
        echo -e "${YELLOW}No similar issues found${NC}"
    fi
}

# Function to find related issues by keyword
find_related_issues() {
    echo -e "${BLUE}Find related issues by keyword...${NC}"
    
    # Ask for keyword
    read -p "Enter keyword to search for: " keyword
    
    if [ -z "$keyword" ]; then
        echo -e "${RED}No keyword provided. Exiting.${NC}"
        return
    fi
    
    # Get all issues that contain the keyword in title or body
    echo -e "${BLUE}Searching for issues containing '$keyword'...${NC}"
    
    # Search in titles
    title_results=$(gh issue list --limit 500 --state all --json number,title,url,state --search "in:title $keyword")
    
    # Search in body
    body_results=$(gh issue list --limit 500 --state all --json number,title,url,state --search "in:body $keyword")
    
    # Combine results and remove duplicates
    all_results=$(echo "$title_results $body_results" | jq -s 'add | unique_by(.number)')
    
    # Save results
    results_file="$OUTPUT_DIR/keyword_issues_${keyword}_$(date +%Y%m%d_%H%M%S).json"
    echo "$all_results" > "$results_file"
    
    # Display results
    result_count=$(echo "$all_results" | jq '. | length')
    
    if [ "$result_count" -gt 0 ]; then
        echo -e "${GREEN}Found $result_count issues containing '$keyword'${NC}"
        echo -e "${BLUE}Issues:${NC}"
        echo "$all_results" | jq -r '.[] | "Issue #\(.number) - \(.title) - \(.state)"'
        echo -e "${YELLOW}Full results saved to $results_file${NC}"
    else
        echo -e "${YELLOW}No issues found containing '$keyword'${NC}"
    fi
}

# Function to consolidate issues
consolidate_issues() {
    echo -e "${BLUE}Consolidate issues (merge content and close duplicates)...${NC}"
    
    # Ask for primary issue
    read -p "Enter primary issue number (issue to keep): " primary_issue
    
    if [ -z "$primary_issue" ]; then
        echo -e "${RED}No primary issue provided. Exiting.${NC}"
        return
    fi
    
    # Ask for issues to consolidate
    read -p "Enter comma-separated issue numbers to merge into primary issue: " duplicate_issues
    
    if [ -z "$duplicate_issues" ]; then
        echo -e "${RED}No duplicate issues provided. Exiting.${NC}"
        return
    fi
    
    # Verify with user
    echo -e "${YELLOW}You are about to consolidate the following issues into #$primary_issue:${NC}"
    
    IFS=',' read -ra ISSUES <<< "$duplicate_issues"
    for issue in "${ISSUES[@]}"; do
        issue=$(echo "$issue" | xargs)  # Trim whitespace
        echo -e "Issue #$issue"
    done
    
    read -p "Proceed? (y/n): " confirm
    
    if [ "$confirm" != "y" ]; then
        echo -e "${YELLOW}Operation cancelled.${NC}"
        return
    fi
    
    # Get primary issue details
    primary_issue_data=$(gh issue view "$primary_issue" --json number,title,body,url)
    primary_title=$(echo "$primary_issue_data" | jq -r '.title')
    primary_body=$(echo "$primary_issue_data" | jq -r '.body')
    primary_url=$(echo "$primary_issue_data" | jq -r '.url')
    
    # Create a consolidated comment for the primary issue
    consolidated_comment="## Consolidated Issues\n\nThe following issues have been consolidated into this one:\n\n"
    
    # Process each duplicate issue
    for issue in "${ISSUES[@]}"; do
        issue=$(echo "$issue" | xargs)  # Trim whitespace
        echo -e "${BLUE}Processing issue #$issue...${NC}"
        
        # Get issue details
        issue_data=$(gh issue view "$issue" --json number,title,body,url)
        issue_title=$(echo "$issue_data" | jq -r '.title')
        issue_body=$(echo "$issue_data" | jq -r '.body')
        issue_url=$(echo "$issue_data" | jq -r '.url')
        
        # Create consolidation comment for the primary issue
        consolidated_comment="${consolidated_comment}- #$issue: $issue_title\n"
        consolidated_comment="${consolidated_comment}  - ${issue_url}\n\n"
        consolidated_comment="${consolidated_comment}  <details>\n  <summary>Original Description</summary>\n\n  $(echo "$issue_body" | sed 's/^/  /')\n  </details>\n\n"
        
        # Add comment to duplicate issue
        gh issue comment "$issue" --body "This issue has been consolidated into #$primary_issue ($primary_title).\n\nPlease refer to that issue for further updates."
        
        # Close duplicate issue with appropriate label
        gh issue edit "$issue" --add-label "duplicate" --state closed
        
        echo -e "${GREEN}Issue #$issue consolidated and closed${NC}"
    done
    
    # Add consolidated comment to primary issue
    gh issue comment "$primary_issue" --body "$consolidated_comment"
    
    echo -e "${GREEN}Consolidation complete! All issues have been merged into #$primary_issue${NC}"
}

# Function to move issues between milestones
move_issues() {
    echo -e "${BLUE}Move issues between milestones...${NC}"
    
    # Get available milestones
    echo -e "${BLUE}Available milestones:${NC}"
    milestones=$(gh api repos/$repo_owner/$repo_name/milestones)
    echo "$milestones" | jq -r '.[] | "#\(.number): \(.title)"'
    
    # Ask for source milestone
    read -p "Enter source milestone number (or 'none' for issues without milestone): " source_milestone
    
    # Ask for target milestone
    read -p "Enter target milestone number: " target_milestone
    
    if [ -z "$target_milestone" ]; then
        echo -e "${RED}No target milestone provided. Exiting.${NC}"
        return
    fi
    
    # Verify target milestone exists
    if ! echo "$milestones" | jq -e ".[] | select(.number == $target_milestone)" > /dev/null; then
        echo -e "${RED}Target milestone #$target_milestone not found${NC}"
        return
    fi
    
    # Get issues in source milestone
    if [ "$source_milestone" == "none" ]; then
        # Get issues without milestone
        source_issues=$(gh issue list --limit 500 --state all --json number,title,milestone --milestone none)
    else
        # Get issues in specified milestone
        source_issues=$(gh issue list --limit 500 --state all --json number,title,milestone --milestone "$source_milestone")
    fi
    
    # Count issues
    issue_count=$(echo "$source_issues" | jq '. | length')
    
    if [ "$issue_count" -eq 0 ]; then
        echo -e "${YELLOW}No issues found in the source milestone${NC}"
        return
    fi
    
    # Display issues
    echo -e "${BLUE}Found $issue_count issues to move:${NC}"
    echo "$source_issues" | jq -r '.[] | "#\(.number): \(.title)"'
    
    # Confirm with user
    read -p "Move these issues to milestone #$target_milestone? (y/n): " confirm
    
    if [ "$confirm" != "y" ]; then
        echo -e "${YELLOW}Operation cancelled.${NC}"
        return
    fi
    
    # Move issues
    target_milestone_title=$(echo "$milestones" | jq -r ".[] | select(.number == $target_milestone) | .title")
    log_file="$LOGS_DIR/move_issues_$(date +%Y%m%d_%H%M%S).log"
    
    echo "Moving issues to milestone #$target_milestone ($target_milestone_title)" > "$log_file"
    echo "Started at $(date)" >> "$log_file"
    echo "" >> "$log_file"
    
    echo "$source_issues" | jq -c '.[]' | while read -r issue; do
        issue_number=$(echo "$issue" | jq -r '.number')
        issue_title=$(echo "$issue" | jq -r '.title')
        
        echo -e "${BLUE}Moving issue #$issue_number: $issue_title${NC}"
        gh issue edit "$issue_number" --milestone "$target_milestone"
        
        echo "Moved issue #$issue_number: $issue_title" >> "$log_file"
    done
    
    echo "" >> "$log_file"
    echo "Completed at $(date)" >> "$log_file"
    
    echo -e "${GREEN}All issues moved to milestone #$target_milestone ($target_milestone_title)${NC}"
    echo -e "${YELLOW}Log saved to $log_file${NC}"
}

# Function to mass close issues with a comment
mass_close_issues() {
    echo -e "${BLUE}Mass close issues with a comment...${NC}"
    
    # Ask for query
    echo -e "${YELLOW}Enter a search query to select issues (e.g., label:bug state:open)${NC}"
    read -p "Query: " query
    
    if [ -z "$query" ]; then
        echo -e "${RED}No query provided. Exiting.${NC}"
        return
    fi
    
    # Get matching issues
    echo -e "${BLUE}Searching for issues matching: $query${NC}"
    matching_issues=$(gh issue list --limit 500 --json number,title,state --search "$query")
    
    # Count issues
    issue_count=$(echo "$matching_issues" | jq '. | length')
    
    if [ "$issue_count" -eq 0 ]; then
        echo -e "${YELLOW}No issues found matching the query${NC}"
        return
    fi
    
    # Display issues
    echo -e "${BLUE}Found $issue_count issues:${NC}"
    echo "$matching_issues" | jq -r '.[] | "#\(.number): \(.title) (\(.state))"'
    
    # Ask for confirmation
    read -p "Are you sure you want to close these $issue_count issues? (y/n): " confirm
    
    if [ "$confirm" != "y" ]; then
        echo -e "${YELLOW}Operation cancelled.${NC}"
        return
    fi
    
    # Ask for closing comment
    echo -e "${YELLOW}Enter a comment to add when closing these issues:${NC}"
    read -p "Comment: " comment
    
    if [ -z "$comment" ]; then
        comment="This issue has been closed during repository cleanup."
    fi
    
    # Ask for label
    read -p "Add a label when closing? (leave empty for none): " label
    
    # Close issues
    log_file="$LOGS_DIR/mass_close_$(date +%Y%m%d_%H%M%S).log"
    
    echo "Closing issues matching query: $query" > "$log_file"
    echo "Started at $(date)" >> "$log_file"
    echo "Closing comment: $comment" >> "$log_file"
    if [ -n "$label" ]; then
        echo "Adding label: $label" >> "$log_file"
    fi
    echo "" >> "$log_file"
    
    echo "$matching_issues" | jq -c '.[]' | while read -r issue; do
        issue_number=$(echo "$issue" | jq -r '.number')
        issue_title=$(echo "$issue" | jq -r '.title')
        issue_state=$(echo "$issue" | jq -r '.state')
        
        # Skip already closed issues
        if [ "$issue_state" == "CLOSED" ]; then
            echo -e "${YELLOW}Issue #$issue_number already closed, skipping${NC}"
            echo "Skipped issue #$issue_number (already closed): $issue_title" >> "$log_file"
            continue
        fi
        
        echo -e "${BLUE}Closing issue #$issue_number: $issue_title${NC}"
        
        # Add comment
        gh issue comment "$issue_number" --body "$comment"
        
        # Close issue
        if [ -n "$label" ]; then
            gh issue edit "$issue_number" --add-label "$label" --state closed
        else
            gh issue edit "$issue_number" --state closed
        fi
        
        echo "Closed issue #$issue_number: $issue_title" >> "$log_file"
    done
    
    echo "" >> "$log_file"
    echo "Completed at $(date)" >> "$log_file"
    
    echo -e "${GREEN}All matching issues have been closed${NC}"
    echo -e "${YELLOW}Log saved to $log_file${NC}"
}

# Function to generate consolidation report
generate_report() {
    echo -e "${BLUE}Generating consolidation report...${NC}"
    
    # Get repository statistics
    repo_info=$(gh api repos/$repo_owner/$repo_name)
    
    # Get issues by state
    open_issues=$(gh issue list --limit 500 --state open --json number,title,milestone,labels)
    closed_issues=$(gh issue list --limit 500 --state closed --json number,title,milestone,labels)
    
    # Get milestones
    milestones=$(gh api repos/$repo_owner/$repo_name/milestones)
    
    # Get labels
    labels=$(gh api repos/$repo_owner/$repo_name/labels --paginate)
    
    # Create report file
    report_file="$OUTPUT_DIR/consolidation_report_$(date +%Y%m%d_%H%M%S).md"
    
    # Write report
    cat > "$report_file" << EOL
# GitHub Issue Consolidation Report

Repository: **$repo_owner/$repo_name**
Generated on: **$(date)**

## Overview

- **Open Issues:** $(echo "$open_issues" | jq '. | length')
- **Closed Issues:** $(echo "$closed_issues" | jq '. | length')
- **Total Issues:** $(( $(echo "$open_issues" | jq '. | length') + $(echo "$closed_issues" | jq '. | length') ))
- **Milestones:** $(echo "$milestones" | jq '. | length')
- **Labels:** $(echo "$labels" | jq '. | length')

## Issues by Milestone

$(echo "$milestones" | jq -r '.[] | "### " + .title + "\n\n- **Due date:** " + (.due_on // "No due date") + "\n- **Description:** " + (.description // "No description") + "\n- **Open issues:** " + (.open_issues | tostring) + "\n- **Closed issues:** " + (.closed_issues | tostring) + "\n"')

## Issues by Component

$(echo "$open_issues" | jq '[.[] | .labels[].name] | map(select(startswith("component:"))) | group_by(.) | map({component: .[0], count: length}) | sort_by(.count) | reverse[]' | jq -r '. | "- **" + (.component | gsub("component:"; "")) + ":** " + (.count | tostring) + " issues"')

## Issues by Priority

$(echo "$open_issues" | jq '[.[] | .labels[].name] | map(select(startswith("priority:"))) | group_by(.) | map({priority: .[0], count: length}) | sort_by(.count) | reverse[]' | jq -r '. | "- **" + (.priority | gsub("priority:"; "")) + ":** " + (.count | tostring) + " issues"')

## Issues by Type

$(echo "$open_issues" | jq '[.[] | .labels[].name] | map(select(startswith("type:"))) | group_by(.) | map({type: .[0], count: length}) | sort_by(.count) | reverse[]' | jq -r '. | "- **" + (.type | gsub("type:"; "")) + ":** " + (.count | tostring) + " issues"')

## Issues Without Labels

$(echo "$open_issues" | jq '[.[] | select(.labels | length == 0)]' | jq -r '.[] | "- #" + (.number | tostring) + ": " + .title')

## Issues Without Milestone

$(echo "$open_issues" | jq '[.[] | select(.milestone == null)]' | jq -r '.[] | "- #" + (.number | tostring) + ": " + .title')

## Recommendations

Based on the analysis of the repository issues, here are some recommendations:

1. **Labeling**: Ensure all issues have appropriate component, priority, and type labels
2. **Milestones**: Assign all issues to relevant milestones
3. **Duplicate Detection**: Review potential duplicate issues identified during analysis
4. **Cleanup**: Close or archive issues that are no longer relevant
5. **Templates**: Use standardized issue templates for future issues

## Recent Consolidation Actions

*This section should be manually updated with recent consolidation actions taken:*

- Consolidated issues #X, #Y into issue #Z
- Moved N issues from milestone A to milestone B
- Closed X stale issues
- Standardized Y issue titles and descriptions

EOL
    
    echo -e "${GREEN}Consolidation report generated: $report_file${NC}"
}

# Function to find issues by content similarity
find_content_similarity() {
    echo -e "${BLUE}Find issues by content similarity...${NC}"
    
    # This is a simplified implementation - for more advanced content similarity,
    # consider using external tools or a more sophisticated algorithm
    
    # Ask for minimum body length to consider
    read -p "Enter minimum body length to consider (recommended: 100): " min_length
    
    if [ -z "$min_length" ]; then
        min_length=100
    fi
    
    # Get all issues with substantial bodies
    echo -e "${BLUE}Retrieving issues with substantial content...${NC}"
    all_issues=$(gh issue list --limit 500 --state all --json number,title,body,url,state)
    
    # Filter issues with sufficient body length
    echo -e "${BLUE}Filtering issues with content length >= $min_length...${NC}"
    filtered_issues=$(echo "$all_issues" | jq --arg min "$min_length" '[.[] | select(.body != null and .body != "" and (.body | length) >= ($min | tonumber))]')
    filtered_count=$(echo "$filtered_issues" | jq '. | length')
    
    echo -e "${BLUE}Found $filtered_count issues with substantial content${NC}"
    
    if [ "$filtered_count" -eq 0 ]; then
        echo -e "${YELLOW}No issues with substantial content found${NC}"
        return
    fi
    
    # Create similarity output file
    similarity_file="$OUTPUT_DIR/content_similarity_$(date +%Y%m%d_%H%M%S).json"
    
    # Compare issue contents (this is a simplified comparison)
    # For large repositories, this could be slow and memory-intensive
    echo -e "${BLUE}Comparing issue contents (this may take some time)...${NC}"
    
    echo "$filtered_issues" | jq -c '.[]' | while read -r issue1; do
        issue1_number=$(echo "$issue1" | jq -r '.number')
        issue1_title=$(echo "$issue1" | jq -r '.title')
        issue1_body=$(echo "$issue1" | jq -r '.body')
        issue1_url=$(echo "$issue1" | jq -r '.url')
        issue1_state=$(echo "$issue1" | jq -r '.state')
        
        echo "$filtered_issues" | jq -c '.[]' | while read -r issue2; do
            issue2_number=$(echo "$issue2" | jq -r '.number')
            
            # Skip comparing the same issue or already compared pairs
            if [ "$issue1_number" -ge "$issue2_number" ]; then
                continue
            fi
            
            issue2_title=$(echo "$issue2" | jq -r '.title')
            issue2_body=$(echo "$issue2" | jq -r '.body')
            issue2_url=$(echo "$issue2" | jq -r '.url')
            issue2_state=$(echo "$issue2" | jq -r '.state')
            
            # Calculate similarity using difflib
            similarity=$(python3 -c "
import difflib
import sys
body1 = '''$issue1_body'''.lower()
body2 = '''$issue2_body'''.lower()
similarity = difflib.SequenceMatcher(None, body1, body2).ratio()
print(f'{similarity:.2f}')
")
            
            # If similarity is above threshold, record it
            if (( $(echo "$similarity > 0.4" | bc -l) )); then
                echo "{\"issue1\": {\"number\": $issue1_number, \"title\": \"$issue1_title\", \"url\": \"$issue1_url\", \"state\": \"$issue1_state\"}, \"issue2\": {\"number\": $issue2_number, \"title\": \"$issue2_title\", \"url\": \"$issue2_url\", \"state\": \"$issue2_state\"}, \"similarity\": $similarity}" >> "$similarity_file"
            fi
        done
    done
    
    # Format the similarity file as proper JSON
    if [ -f "$similarity_file" ]; then
        # Count lines in file
        line_count=$(wc -l < "$similarity_file")
        
        if [ "$line_count" -gt 0 ]; then
            # Create a proper JSON array
            echo "[" > "$similarity_file.tmp"
            cat "$similarity_file" | sed '$!s/$/,/' >> "$similarity_file.tmp"
            echo "]" >> "$similarity_file.tmp"
            mv "$similarity_file.tmp" "$similarity_file"
            
            # Display results
            echo -e "${GREEN}Found $(jq '. | length' "$similarity_file") potential similar issue pairs${NC}"
            echo -e "${BLUE}Top similar pairs:${NC}"
            jq -r 'sort_by(.similarity) | reverse | .[0:5] | .[] | "Issue #\(.issue1.number) and #\(.issue2.number) - Similarity: \(.similarity) - \(.issue1.title) / \(.issue2.title)"' "$similarity_file"
            echo -e "${YELLOW}Full results saved to $similarity_file${NC}"
        else
            echo -e "${YELLOW}No similar issues found${NC}"
            rm "$similarity_file"
        fi
    else
        echo -e "${YELLOW}No similar issues found${NC}"
    fi
}

# Function to find orphaned PRD references
find_orphaned_prds() {
    echo -e "${BLUE}Finding orphaned PRD references...${NC}"
    
    # Get all PRD files
    echo -e "${BLUE}Searching for PRD files...${NC}"
    prd_files=$(find /home/mushu/Projects/CryoProtect/scripts -name "prd_*.txt")
    
    # Get all GitHub issues
    echo -e "${BLUE}Retrieving GitHub issues...${NC}"
    all_issues=$(gh issue list --limit 500 --state all --json number,title,url)
    
    # Check each PRD for GitHub issue references
    orphaned_prds=0
    orphaned_list="$OUTPUT_DIR/orphaned_prds_$(date +%Y%m%d_%H%M%S).json"
    echo "[]" > "$orphaned_list"
    
    for prd_file in $prd_files; do
        prd_basename=$(basename "$prd_file")
        
        # Check for GitHub issue reference
        issue_number=$(grep -o "GitHub issue #[0-9]\+" "$prd_file" | grep -o "[0-9]\+")
        related_issue=$(grep -o "Related Issue: #[0-9]\+" "$prd_file" | grep -o "[0-9]\+")
        
        if [ -z "$issue_number" ] && [ -z "$related_issue" ]; then
            echo -e "${YELLOW}PRD file $prd_basename has no GitHub issue reference${NC}"
            orphaned_prds=$((orphaned_prds + 1))
            
            # Add to orphaned list
            temp_file=$(mktemp)
            jq --arg file "$prd_basename" --arg path "$prd_file" \
               '. + [{"file": $file, "path": $path, "title": "Unknown", "type": "no_reference"}]' \
               "$orphaned_list" > "$temp_file"
            mv "$temp_file" "$orphaned_list"
            continue
        fi
        
        # Use issue_number if found, otherwise use related_issue
        check_issue=${issue_number:-$related_issue}
        
        # Check if the referenced issue exists
        if ! echo "$all_issues" | jq -e ".[] | select(.number == $check_issue)" > /dev/null; then
            echo -e "${YELLOW}PRD file $prd_basename references non-existent issue #$check_issue${NC}"
            orphaned_prds=$((orphaned_prds + 1))
            
            # Get PRD title
            prd_title=$(grep -m 1 "^# " "$prd_file" | sed 's/^# //')
            
            # Add to orphaned list
            temp_file=$(mktemp)
            jq --arg file "$prd_basename" --arg path "$prd_file" --arg issue "$check_issue" --arg title "$prd_title" \
               '. + [{"file": $file, "path": $path, "issue": $issue, "title": $title, "type": "missing_issue"}]' \
               "$orphaned_list" > "$temp_file"
            mv "$temp_file" "$orphaned_list"
        fi
    done
    
    if [ "$orphaned_prds" -eq 0 ]; then
        echo -e "${GREEN}No orphaned PRD references found${NC}"
        rm "$orphaned_list"
    else
        echo -e "${YELLOW}Found $orphaned_prds orphaned PRD references${NC}"
        echo -e "${BLUE}Orphaned PRDs:${NC}"
        jq -r '.[] | "- \(.file) (\(.type) - \(.title))"' "$orphaned_list"
        echo -e "${YELLOW}Full results saved to $orphaned_list${NC}"
    fi
}

# Function to bulk update issue titles
bulk_update_titles() {
    echo -e "${BLUE}Bulk update issue titles...${NC}"
    
    # Ask for query
    echo -e "${YELLOW}Enter a search query to select issues (e.g., label:bug state:open)${NC}"
    read -p "Query: " query
    
    if [ -z "$query" ]; then
        echo -e "${RED}No query provided. Exiting.${NC}"
        return
    fi
    
    # Get matching issues
    echo -e "${BLUE}Searching for issues matching: $query${NC}"
    matching_issues=$(gh issue list --limit 500 --json number,title,state --search "$query")
    
    # Count issues
    issue_count=$(echo "$matching_issues" | jq '. | length')
    
    if [ "$issue_count" -eq 0 ]; then
        echo -e "${YELLOW}No issues found matching the query${NC}"
        return
    fi
    
    # Display issues
    echo -e "${BLUE}Found $issue_count issues:${NC}"
    echo "$matching_issues" | jq -r '.[] | "#\(.number): \(.title) (\(.state))"'
    
    # Ask for title pattern and replacement
    echo -e "${YELLOW}Enter search pattern for title text:${NC}"
    read -p "Pattern: " pattern
    
    if [ -z "$pattern" ]; then
        echo -e "${RED}No pattern provided. Exiting.${NC}"
        return
    fi
    
    echo -e "${YELLOW}Enter replacement text:${NC}"
    read -p "Replacement: " replacement
    
    # Ask for confirmation
    read -p "Are you sure you want to update titles for these $issue_count issues? (y/n): " confirm
    
    if [ "$confirm" != "y" ]; then
        echo -e "${YELLOW}Operation cancelled.${NC}"
        return
    fi
    
    # Update titles
    log_file="$LOGS_DIR/title_updates_$(date +%Y%m%d_%H%M%S).log"
    
    echo "Updating titles for issues matching query: $query" > "$log_file"
    echo "Pattern: $pattern" >> "$log_file"
    echo "Replacement: $replacement" >> "$log_file"
    echo "Started at $(date)" >> "$log_file"
    echo "" >> "$log_file"
    
    echo "$matching_issues" | jq -c '.[]' | while read -r issue; do
        issue_number=$(echo "$issue" | jq -r '.number')
        issue_title=$(echo "$issue" | jq -r '.title')
        
        # Generate new title
        new_title=$(echo "$issue_title" | sed "s/$pattern/$replacement/g")
        
        # Skip if title didn't change
        if [ "$new_title" == "$issue_title" ]; then
            echo -e "${YELLOW}Issue #$issue_number title unchanged, skipping${NC}"
            echo "Skipped issue #$issue_number (unchanged): $issue_title" >> "$log_file"
            continue
        fi
        
        echo -e "${BLUE}Updating issue #$issue_number${NC}"
        echo -e "${YELLOW}Old: $issue_title${NC}"
        echo -e "${GREEN}New: $new_title${NC}"
        
        # Update title
        gh issue edit "$issue_number" --title "$new_title"
        
        echo "Updated issue #$issue_number" >> "$log_file"
        echo "  Old: $issue_title" >> "$log_file"
        echo "  New: $new_title" >> "$log_file"
        echo "" >> "$log_file"
    done
    
    echo "Completed at $(date)" >> "$log_file"
    
    echo -e "${GREEN}Title updates completed!${NC}"
    echo -e "${YELLOW}Log saved to $log_file${NC}"
}

# Execute selected action
case $action in
    1) find_duplicate_issues ;;
    2) find_related_issues ;;
    3) consolidate_issues ;;
    4) move_issues ;;
    5) mass_close_issues ;;
    6) generate_report ;;
    7) find_content_similarity ;;
    8) find_orphaned_prds ;;
    9) bulk_update_titles ;;
    0) echo "Exiting..." && exit 0 ;;
    *) echo -e "${RED}Invalid option selected${NC}" ;;
esac