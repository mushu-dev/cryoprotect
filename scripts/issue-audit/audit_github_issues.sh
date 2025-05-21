#!/bin/bash
# GitHub Issue Audit Script
# This script analyzes GitHub issues and generates a report to help clean up the repository

# Colors for terminal output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

OUTPUT_DIR="/home/mushu/Projects/CryoProtect/scripts/issue-audit/output"
mkdir -p "$OUTPUT_DIR"

echo -e "${BLUE}CryoProtect GitHub Issue Audit${NC}"
echo "==============================="
echo "Running comprehensive issue analysis..."

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

# Get all issues (both open and closed)
echo -e "\n${BLUE}Fetching all issues...${NC}"
all_issues=$(gh issue list --limit 500 --state all --json number,title,state,author,createdAt,updatedAt,labels,assignees,milestone,body)

# Save raw issue data
echo "$all_issues" > "$OUTPUT_DIR/all_issues.json"

# Count issues
total_issues=$(echo "$all_issues" | jq '. | length')
open_issues=$(echo "$all_issues" | jq '[.[] | select(.state == "OPEN")] | length')
closed_issues=$(echo "$all_issues" | jq '[.[] | select(.state == "CLOSED")] | length')

echo -e "${BLUE}Total issues:${NC} $total_issues"
echo -e "${BLUE}Open issues:${NC} $open_issues"
echo -e "${BLUE}Closed issues:${NC} $closed_issues"

# Identify potentially blank/meaningless issues
blank_issues=$(echo "$all_issues" | jq '[.[] | select(.body == null or .body == "" or (.body | length) < 50)]')
blank_count=$(echo "$blank_issues" | jq '. | length')
echo -e "${YELLOW}Potentially blank or short issues:${NC} $blank_count"
echo "$blank_issues" > "$OUTPUT_DIR/blank_issues.json"

# Identify duplicate or similar titles
echo -e "\n${BLUE}Analyzing for duplicate titles...${NC}"
duplicate_titles=$(echo "$all_issues" | jq '[.[] | .title] | group_by(.) | map(select(length > 1)) | map({title: .[0], count: length})')
duplicate_count=$(echo "$duplicate_titles" | jq '. | length')
echo -e "${YELLOW}Potentially duplicate titles:${NC} $duplicate_count"
echo "$duplicate_titles" > "$OUTPUT_DIR/duplicate_titles.json"

# Analyze labels
echo -e "\n${BLUE}Analyzing labels...${NC}"
all_labels=$(echo "$all_issues" | jq '[.[] | .labels[]] | map(.name) | unique')
label_count=$(echo "$all_labels" | jq '. | length')
echo -e "${BLUE}Unique labels:${NC} $label_count"
echo "$all_labels" > "$OUTPUT_DIR/all_labels.json"

# Count issues without labels
no_labels=$(echo "$all_issues" | jq '[.[] | select(.labels | length == 0)]')
no_labels_count=$(echo "$no_labels" | jq '. | length')
echo -e "${YELLOW}Issues without labels:${NC} $no_labels_count"
echo "$no_labels" > "$OUTPUT_DIR/no_labels.json"

# Analyze milestones
echo -e "\n${BLUE}Analyzing milestones...${NC}"
all_milestones=$(echo "$all_issues" | jq '[.[] | .milestone] | map(select(. != null)) | unique_by(.title)')
milestone_count=$(echo "$all_milestones" | jq '. | length')
echo -e "${BLUE}Unique milestones:${NC} $milestone_count"
echo "$all_milestones" > "$OUTPUT_DIR/all_milestones.json"

# Count issues without milestones
no_milestone=$(echo "$all_issues" | jq '[.[] | select(.milestone == null)]')
no_milestone_count=$(echo "$no_milestone" | jq '. | length')
echo -e "${YELLOW}Issues without milestones:${NC} $no_milestone_count"
echo "$no_milestone" > "$OUTPUT_DIR/no_milestone.json"

# Identify issues with very short titles
short_titles=$(echo "$all_issues" | jq '[.[] | select((.title | length) < 15)]')
short_titles_count=$(echo "$short_titles" | jq '. | length')
echo -e "${YELLOW}Issues with very short titles:${NC} $short_titles_count"
echo "$short_titles" > "$OUTPUT_DIR/short_titles.json"

# Identify old untouched issues (not updated in last 30 days)
thirty_days_ago=$(date -d "30 days ago" +%s)
old_issues=$(echo "$all_issues" | jq --arg date "$thirty_days_ago" '[.[] | select(.state == "OPEN" and ((.updatedAt | fromdateiso8601 | strftime("%s")) | tonumber) < ($date | tonumber))]')
old_issues_count=$(echo "$old_issues" | jq '. | length')
echo -e "${YELLOW}Open issues not updated in 30+ days:${NC} $old_issues_count"
echo "$old_issues" > "$OUTPUT_DIR/old_issues.json"

# Count issues by author
echo -e "\n${BLUE}Analyzing authors...${NC}"
authors=$(echo "$all_issues" | jq 'group_by(.author.login) | map({author: .[0].author.login, count: length})')
echo "$authors" > "$OUTPUT_DIR/authors.json"

# Generate issue quality scores (basic heuristic)
echo -e "\n${BLUE}Generating issue quality scores...${NC}"
quality_issues=$(echo "$all_issues" | jq '[.[] | {
  number: .number,
  title: .title,
  quality_score: (
    (if .body then (if (.body | length) > 500 then 3 elif (.body | length) > 100 then 2 elif (.body | length) > 0 then 1 else 0 end) else 0 end) +
    (if (.labels | length) > 0 then 1 else 0 end) +
    (if .milestone != null then 1 else 0 end) +
    (if (.assignees | length) > 0 then 1 else 0 end) +
    (if (.title | length) > 50 then 2 elif (.title | length) > 20 then 1 else 0 end)
  )
}] | sort_by(.quality_score) | reverse')
echo "$quality_issues" > "$OUTPUT_DIR/quality_issues.json"

# Count issues by quality score
high_quality=$(echo "$quality_issues" | jq '[.[] | select(.quality_score >= 5)] | length')
medium_quality=$(echo "$quality_issues" | jq '[.[] | select(.quality_score >= 3 and .quality_score < 5)] | length')
low_quality=$(echo "$quality_issues" | jq '[.[] | select(.quality_score < 3)] | length')

echo -e "${GREEN}High quality issues (score >= 5):${NC} $high_quality"
echo -e "${BLUE}Medium quality issues (score 3-4):${NC} $medium_quality"
echo -e "${YELLOW}Low quality issues (score < 3):${NC} $low_quality"

# Create action lists based on quality assessment
echo -e "\n${BLUE}Creating action recommendation lists...${NC}"

# Issues to keep as-is (high quality)
echo "$quality_issues" | jq '[.[] | select(.quality_score >= 5) | .number]' > "$OUTPUT_DIR/keep_issues.json"

# Issues to improve (medium quality)
echo "$quality_issues" | jq '[.[] | select(.quality_score >= 3 and .quality_score < 5) | .number]' > "$OUTPUT_DIR/improve_issues.json"

# Issues to review for possible archiving or deletion (low quality)
echo "$quality_issues" | jq '[.[] | select(.quality_score < 3) | .number]' > "$OUTPUT_DIR/review_issues.json"

# Generate summary report
echo -e "\n${BLUE}Generating summary report...${NC}"
cat > "$OUTPUT_DIR/issue_audit_summary.md" << EOL
# GitHub Issue Audit Summary for $repo_owner/$repo_name

## Overview
- **Total Issues:** $total_issues
- **Open Issues:** $open_issues
- **Closed Issues:** $closed_issues

## Quality Assessment
- **High Quality Issues:** $high_quality
- **Medium Quality Issues:** $medium_quality
- **Low Quality Issues:** $low_quality

## Issue Problems
- **Blank/Short Issues:** $blank_count
- **Duplicate Titles:** $duplicate_count
- **Issues Without Labels:** $no_labels_count
- **Issues Without Milestones:** $no_milestone_count
- **Issues with Very Short Titles:** $short_titles_count
- **Stale Issues (30+ days):** $old_issues_count

## Recommendations
1. **Keep As-Is:** $high_quality issues (see keep_issues.json)
2. **Improve:** $medium_quality issues (see improve_issues.json)
3. **Review for Cleanup:** $low_quality issues (see review_issues.json)

## Next Steps
1. Run the cleanup script with appropriate action for each issue category
2. Standardize remaining issues using the issue templates
3. Apply consistent labeling and milestone assignments
4. Close or merge duplicate issues
5. Archive issues that are no longer relevant

## Files Generated
- all_issues.json - Complete issue data
- blank_issues.json - Issues with little or no content
- duplicate_titles.json - Issues with similar titles
- all_labels.json - All labels used in the repository
- no_labels.json - Issues without any labels
- all_milestones.json - All milestones used
- no_milestone.json - Issues without milestones
- short_titles.json - Issues with very short titles
- old_issues.json - Stale issues not updated recently
- authors.json - Issue counts by author
- quality_issues.json - All issues with quality scores
- keep_issues.json - High quality issues to keep
- improve_issues.json - Medium quality issues to improve
- review_issues.json - Low quality issues to review or delete
EOL

echo -e "${GREEN}Audit complete! Reports saved to $OUTPUT_DIR${NC}"
echo -e "Review the summary at: $OUTPUT_DIR/issue_audit_summary.md"