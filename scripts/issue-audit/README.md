# GitHub Issue Management Tools

This directory contains tools for managing and organizing GitHub issues for the CryoProtect project. These scripts help clean up the repository, standardize issues, and maintain a consistent workflow.

## Overview of Scripts

1. **audit_github_issues.sh** - Analyzes repository issues and generates a quality report
2. **cleanup_github_issues.sh** - Processes issues based on audit results
3. **standardize_github_repository.sh** - Sets up standard labels, milestones, and templates
4. **consolidate_github_issues.sh** - Helps merge duplicate issues and maintain organization

## Requirements

- GitHub CLI (`gh`) - Must be authenticated with appropriate permissions
- `jq` - For JSON processing
- `python3` - For content similarity analysis
- `bc` - For numerical comparisons

## Usage Workflow

### 1. Audit the Repository

Start by running the audit to analyze the current state of your issues:

```bash
./audit_github_issues.sh
```

This will:
- Analyze all GitHub issues in the repository
- Assign quality scores based on content, labels, etc.
- Identify potentially blank or duplicate issues
- Generate a comprehensive report in the `output` directory

### 2. Standardize the Repository

Next, set up standard labels, milestones, and issue templates:

```bash
./standardize_github_repository.sh
```

Choose from options to:
- Create standard labels (component, priority, type, status)
- Create project milestones aligned with your phases
- Create issue templates for consistent formatting
- Clean up existing labels

### 3. Clean Up Issues

Process issues based on the audit results:

```bash
./cleanup_github_issues.sh
```

Options include:
- Archive low-quality issues
- Standardize issue content with templates
- Fix issues missing labels or milestones
- Close empty issues
- Batch update issues

### 4. Consolidate Issues

Manage duplicate or related issues:

```bash
./consolidate_github_issues.sh
```

Features include:
- Finding duplicate issues by title or content similarity
- Consolidating related issues (merging content)
- Moving issues between milestones
- Mass closing issues with standardized comments
- Generating detailed consolidation reports

## Example Usage Scenarios

### Cleaning Up Empty Issues

```bash
# First, run the audit to identify empty issues
./audit_github_issues.sh

# Then run cleanup with option 4 to close empty issues
./cleanup_github_issues.sh
# Select option 4 when prompted
```

### Standardizing Issue Labels

```bash
# Run the standardization script
./standardize_github_repository.sh
# Select option 1 to create standard labels
# Select option 8 to delete unnecessary labels
```

### Finding and Merging Duplicate Issues

```bash
# Run the consolidation script
./consolidate_github_issues.sh
# Select option 1 to find duplicate issues by title
# Then select option 3 to consolidate specific issues
```

## Output Files

All scripts generate output in the `output` directory:
- JSON data files with issue information
- Markdown reports with analysis and recommendations
- Logs of actions taken during cleanup

## Best Practices

1. **Always start with an audit** before making changes
2. **Back up important data** before mass-closing or consolidating issues
3. **Handle high-quality issues with care** - focus cleanup on low-quality issues
4. **Move issues between milestones** rather than closing and recreating
5. **Use standard templates** for all new issues going forward
6. **Run periodic audits** to maintain repository cleanliness

## GitHub Issue Structure Guidelines

For optimal organization, issues should follow this structure:

1. **Title**: Clear, specific, and descriptive
2. **Description**: Detailed explanation of the issue or feature
3. **Requirements**: Specific functionality or acceptance criteria
4. **Technical Context**: Information about the existing system/architecture
5. **Constraints**: Any limitations or requirements to consider
6. **Acceptance Criteria**: Clear testing criteria for completion

## GitHub Issue Metadata Guidelines

All issues should include:

1. **Labels**:
   - Component (e.g., database, api, ui)
   - Priority (high, medium, low)
   - Type (bug, feature, enhancement, etc.)
   - Status (ready, in-progress, blocked, etc.)

2. **Milestone**: Assign to the appropriate project phase

3. **Assignee**: The person responsible for implementation

## Cross-Environment Integration with TaskMaster

These tools integrate with the project's TaskMaster system:

1. When consolidating issues, use `/github-prd update-mapping` to update TaskMaster mappings
2. After standardizing issues, use `/sync-tasks commit` to ensure consistency
3. Reference TaskMaster tasks in issue comments for cross-environment traceability

Remember that TaskMaster provides task management for Cursor IDE while GitHub issues remain the source of truth for the project.