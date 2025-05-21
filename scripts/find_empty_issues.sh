#\!/bin/bash
# Script to identify GitHub issues with missing or minimal descriptions

echo "Finding GitHub issues with minimal descriptions..."

# Get all open issues
gh issue list --limit 100 --json number,title,body > issues.json

# Count issues with "# Task Description" and nothing else substantial
echo "Issues with minimal '# Task Description' content:"
grep -F '"body":"# Task Description' issues.json  < /dev/null |  cut -d'"' -f8 | sed 's/^/Issue #/'
count1=$(grep -F '"body":"# Task Description' issues.json | wc -l)

# Count issues with very short descriptions (less than 100 characters)
echo -e "\nIssues with very short descriptions (less than 100 characters):"
jq -r '.[] | select(.body | length < 100) | .number, .title' issues.json | paste - - | sed 's/^/Issue #/'
count2=$(jq -r '.[] | select(.body | length < 100) | .number' issues.json | wc -l)

echo -e "\nTotal issues needing better descriptions: $count1"
echo -e "\nSample improved format for issue descriptions:"
echo "----------------------------------------"
echo "# [Issue Type] Title

## Description
Detailed explanation of the task or issue.

## Acceptance Criteria
- [ ] Criterion 1
- [ ] Criterion 2

## Dependencies
- Depends on #XX (if applicable)

## Technical Details
Any technical information or context needed.
"
echo "----------------------------------------"
