# RooRoo GitHub Issue Manager for CryoProtect Fedora Migration

This document explains how to use the RooRoo GitHub Issue Manager to track tasks for the CryoProtect Fedora migration project.

## Setup

1. First, make sure you have authenticated with GitHub:

```bash
gh auth login
```

2. Run the setup script to prepare RooRoo:

```bash
./setup_rooroo.sh
```

3. To verify that RooRoo is properly working on your Fedora system, run:

```bash
./test_rooroo_fedora.sh
```

4. To set up the initial issue structure for the Fedora migration project:

```bash
./setup_fedora_migration_issues.sh
```

## Issue Structure

The Fedora migration project uses the following issue structure:

1. **Epics** - High-level categories of work:
   - Epic 1: Initial Setup and Analysis
   - Epic 2: Fedora Core Setup
   - Epic 3: Database Migration
   - Epic 4: Environment Configuration
   - Epic 5: Containerization
   - Epic 6: Testing and Verification

2. **Tasks** - Specific work items that belong to an Epic.

3. **Subtasks** - Smaller units of work that belong to a Task.

## Common Commands

Here are some common commands for working with the RooRoo GitHub Issue Manager:

### Create a New Task

```bash
.github/rooroo/issue-manager.sh create-issue \
  --repo "yourusername/CryoProtect" \
  --title "Task: Implement SELinux policy for PostgreSQL" \
  --body "# Task Description\n\nImplement SELinux policy that allows PostgreSQL to work properly with CryoProtect.\n\n## Requirements\n\n- Create policy module\n- Set contexts for database files\n- Test with SELinux enforcing\n\n## Acceptance Criteria\n\n- PostgreSQL works with SELinux enforcing\n- No SELinux denials during operation" \
  --labels "fedora-migration,selinux"
```

### Link a Task to an Epic

```bash
.github/rooroo/issue-manager.sh add-subissue \
  --repo "yourusername/CryoProtect" \
  --parent 2 \
  --child 8
```

This links issue #8 as a subtask of issue #2 (assuming #2 is an Epic).

### Add a Comment to an Issue

```bash
.github/rooroo/issue-manager.sh add-comment \
  --repo "yourusername/CryoProtect" \
  --issue 8 \
  --body "I've started work on this task. The initial SELinux policy module has been created."
```

### Close an Issue

```bash
.github/rooroo/issue-manager.sh close-issue \
  --repo "yourusername/CryoProtect" \
  --issue 8 \
  --comment "This task is now complete. SELinux policy has been implemented and tested."
```

### List Issues with a Specific Label

```bash
.github/rooroo/issue-manager.sh list-issues \
  --repo "yourusername/CryoProtect" \
  --labels "selinux"
```

### Get Subissues of an Epic

```bash
.github/rooroo/issue-manager.sh get-subissues \
  --repo "yourusername/CryoProtect" \
  --parent 2
```

## Creating Additional Tasks

When creating new tasks, it's recommended to follow this template:

```
# Task: [Task Title]

## Description
[Detailed description of the task]

## Requirements
- [Requirement 1]
- [Requirement 2]
- [Requirement 3]

## Acceptance Criteria
- [Criterion 1]
- [Criterion 2]
- [Criterion 3]

## Related Issues
- [Related Issue or Epic]
```

Then link the task to the appropriate Epic using the `add-subissue` command.

## Labels

The following labels are used in the Fedora migration project:

- `fedora-migration` - All issues related to the Fedora migration
- `selinux` - Issues related to SELinux configuration
- `database` - Issues related to database migration and configuration
- `docker` - Issues related to Docker configuration
- `environment` - Issues related to environment setup
- `testing` - Issues related to testing and verification
- `documentation` - Issues related to documentation
- `p0-blocker` - Highest priority issues that block progress
- `p1-critical` - Critical issues for project success
- `p2-important` - Important issues that should be addressed
- `p3-enhancement` - Nice-to-have enhancements

## Integration with Other Tools

RooRoo GitHub Issue Manager can be integrated with other tools in the CryoProtect workflow:

1. **Strategic Planner** - Create a plan and implement it as GitHub issues
2. **Workflow Coordinator** - Track task status using GitHub issues
3. **Verification Agent** - Report verification results as comments on issues

## Troubleshooting

If you encounter issues with RooRoo:

1. Make sure GitHub CLI is properly authenticated:
```bash
gh auth status
```

2. Check if RooRoo scripts are executable:
```bash
chmod +x .github/rooroo/issue-manager.sh
chmod +x .github/rooroo/scripts/*.sh
```

3. If using SELinux in enforcing mode, you may need to set appropriate contexts:
```bash
sudo chcon -t bin_t .github/rooroo/issue-manager.sh
sudo chcon -t bin_t .github/rooroo/scripts/*.sh
```

4. For more detailed troubleshooting, run the test script:
```bash
./test_rooroo_fedora.sh
```