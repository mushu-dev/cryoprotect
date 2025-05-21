# Claude Code Custom Commands

This directory contains custom commands for Claude Code to use when working with the CryoProtect project.

## Available Commands

When you type a command prefixed with `/` in Claude Code, it will automatically execute the corresponding script in this directory:

- `/github-issues` - Manage GitHub issues
- `/github-prs` - Manage GitHub pull requests
- `/github-consolidate` - Organize issues, setup milestones, labels, link issues
- `/github-project` - View project status, milestones, components, priorities
- `/github-workflow-template` - Create a standardized workflow issue template
- `/github-prd` - Manage Product Requirements Documents for TaskMaster integration

## Command Usage

Commands are prefixed with `/` when used in Claude Code and are automatically executed:

```
/github-prd create "Feature Name" 123
/github-project status
/github-issues list
```

All commands have a `help` option that displays usage information:

```
/github-prd help
```

## TaskMaster Integration

The `/github-prd` command facilitates integration between GitHub issues and TaskMaster:

1. Create a GitHub issue with detailed requirements
2. Extract it as a PRD with `/github-prd extract 123`
3. Generate TaskMaster tasks with `/github-prd generate-tasks scripts/prd_file.txt`
4. Implement tasks in Cursor IDE
5. Return to Claude Code for review

This workflow bridges GitHub (project tracking) and TaskMaster (implementation tasks).