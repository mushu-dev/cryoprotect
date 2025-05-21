# CryoProtect PRD Directory

This directory contains Product Requirements Documents (PRDs) for CryoProtect features. These PRDs are used to generate structured tasks via TaskMaster for implementation in Cursor IDE.

## Environment-Specific Workflow

The CryoProtect project uses a dual task management system:
- **GitHub Issues**: Used as the source of truth and central workspace (only visible to Claude Code)
- **TaskMaster**: Used for ongoing implementation tasks (visible to both Claude Code and Cursor IDE)

Since Cursor IDE cannot directly access GitHub issues, we need a system to synchronize information between these two environments.

## PRD Workflow

1. **Claude Code**: Create GitHub issue to track a new feature or major component
2. **Claude Code**: Structure the issue according to PRD template with all required sections
3. **Claude Code**: Extract the PRD from the issue using: `/github-prd extract issue-number`
4. **Claude Code**: Generate tasks with TaskMaster: `/github-prd generate-tasks scripts/prd_file.txt`
5. **Claude Code**: Update GitHub-TaskMaster mapping with: `/github-prd update-mapping prd-file issue-number`
6. **Claude Code**: Commit task files to Git with: `/sync-tasks commit`
7. **Cursor IDE**: Pull latest tasks from Git with: `/sync-tasks pull`
8. **Cursor IDE**: Implement the tasks in Cursor IDE with TaskMaster guidance (no GitHub access)
9. **Cursor IDE**: Update task status as work progresses and commit with: `/sync-tasks commit`
10. **Claude Code**: Pull task updates and update GitHub issues with implementation progress

## Available PRDs

- `prd_example.txt`: Example PRD for connection pool optimization
- Add additional PRDs as they are created

## PRD Structure

Each PRD should include the following sections:

1. **Overview**: High-level description of the feature/component
2. **Requirements**: Detailed list of must-have functionality
3. **Technical Context**: Existing architecture and integration points
4. **Constraints**: Performance, security, or compatibility requirements
5. **Acceptance Criteria**: How to verify successful implementation

## Helper Commands

### GitHub PRD Management

The custom `/github-prd` helper script provides several commands to manage PRDs:

```bash
# Create a PRD from a GitHub issue
/github-prd create "Feature Name" 123

# Extract PRD sections from an existing GitHub issue
/github-prd extract 123

# Generate TaskMaster tasks from a PRD file
/github-prd generate-tasks scripts/prd_file.txt

# Create a blank PRD template
/github-prd template "New Feature"

# Display a PRD file with syntax highlighting
/github-prd view scripts/prd_file.txt

# Link a PRD file to a GitHub issue for tracking
/github-prd link scripts/prd_file.txt 123

# List all PRD files
/github-prd list

# Update GitHub issue to TaskMaster task mapping
/github-prd update-mapping scripts/prd_file.txt 123

# Show GitHub issue to TaskMaster task mapping
/github-prd mapping [issue_number]
```

### TaskMaster Synchronization

Since both Claude Code and Cursor IDE need to access the same TaskMaster tasks, we use Git to synchronize task files between environments:

```bash
# Check task synchronization status
/sync-tasks status

# Pull latest task changes from Git
/sync-tasks pull

# Commit and push task changes to Git
/sync-tasks push

# List all tasks
/sync-tasks list

# Commit task files to Git with automatic message
/sync-tasks commit
```

### Consistency Verification

To ensure that GitHub issues and TaskMaster tasks remain consistent, you can run:

```bash
/scripts/verify_github_taskmaster_consistency.sh
```

This script checks for:
- Missing PRD files referenced in mappings
- Missing TaskMaster tasks referenced in mappings
- PRD files without GitHub issue references
- PRD files with unmapped GitHub issues
- Orphaned TaskMaster tasks (not referenced in any mapping)
```

## TaskMaster Integration

These PRDs are designed to work with [TaskMaster AI](https://github.com/eyaltoledano/claude-task-master), an AI-powered task management system for structured development with Claude and Cursor IDE.

### TaskMaster Commands

After generating tasks, use these commands to manage them:

```bash
# View all tasks
task-master list

# See next task to work on
task-master next

# View task details
task-master show 1

# Expand a task into subtasks
task-master expand --id=2 --num=5

# Analyze task complexity
task-master analyze-complexity

# Update task status
task-master set-status --id=3 --status=done

# Generate task files
task-master generate

# View complexity report
task-master complexity-report
```

### Environment-Specific Usage

**In Claude Code CLI:**
- Use direct bash commands: `task-master <command>`
- Focus on task generation, planning, and tracking
- Translate GitHub issues to TaskMaster tasks
- Maintain mapping between GitHub issues and TaskMaster tasks
- Update GitHub issues based on TaskMaster progress

**In Cursor IDE:**
- Use TaskMaster via natural language or direct command
- Focus on implementation and task execution
- Ask for task details: "Tell me about task 3"
- Request implementation guidance: "Help me implement task 2.1"
- Update task status: "Mark task 1.2 as completed"

### Synchronization Between Environments

1. **Before switching environments**:
   - Commit task changes with `/sync-tasks commit`
   
2. **When beginning work in a new environment**:
   - Pull latest tasks with `/sync-tasks pull`
   
3. **Regularly check synchronization status**:
   - Monitor with `/sync-tasks status`

Both environments share the same task files in the `/tasks` directory, maintained through Git synchronization.

## Best Practices

1. **Clear Environment Separation**:
   - GitHub Issues: Only in Claude Code - for project tracking and source of truth
   - TaskMaster Tasks: In both Claude Code and Cursor - for implementation details
   - Claude must translate between the two systems

2. **PRD Quality Standards**:
   - Always include GitHub issue reference in PRDs
   - Use detailed requirements with clear acceptance criteria
   - Include technical context with file references
   - Make PRDs self-contained with all necessary information for Cursor

3. **Task Management Discipline**:
   - Commit task changes before switching environments
   - Pull latest tasks when starting work in an environment
   - Regularly verify GitHub-TaskMaster consistency
   - Keep mapping file up-to-date with `/github-prd update-mapping`
   - Mark task dependencies clearly for proper sequencing

4. **Clear Task Communication**:
   - In task descriptions, include specific file paths and line numbers
   - Break complex tasks into smaller, manageable subtasks
   - Add notes about implementation decisions to tasks
   - Reference any constraints or limitations in task details
   - Include verification steps in task descriptions

5. **Regular Verification**:
   - Run the verification script periodically to check consistency
   - Immediately fix any detected inconsistencies
   - Ensure all GitHub issues have corresponding TaskMaster tasks
   - Confirm PRDs are properly linked to issues