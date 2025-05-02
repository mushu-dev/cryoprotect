# TaskMaster AI for CryoProtect v2

This document provides instructions for setting up and using TaskMaster AI with the CryoProtect v2 project.

## Setup

### Option 1: Using the Setup Scripts

#### Windows
1. Run `setup_taskmaster.bat`

#### Linux/macOS
1. Run `./setup_taskmaster.sh`

### Option 2: Manual Setup

1. Install TaskMaster AI:
   ```
   npm install -g task-master-ai
   ```

2. Set environment variables:
   ```
   export ANTHROPIC_API_KEY="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"
   export PERPLEXITY_API_KEY="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"
   export MODEL="claude-3-7-sonnet-20250219"
   export PERPLEXITY_MODEL="sonar-pro"
   export MAX_TOKENS="64000"
   export TEMPERATURE="0.2"
   export DEFAULT_SUBTASKS="5"
   export DEFAULT_PRIORITY="medium"
   ```

3. Run TaskMaster AI:
   ```
   task-master-ai
   ```

### Option 3: Using npm Script

Run the following command:
```
npm run taskmaster
```

## Configuration

The TaskMaster AI configuration is stored in `.taskmaster.json`. You can modify this file to change the configuration.

## Usage

TaskMaster AI will help manage and track tasks for the CryoProtect v2 project. It integrates with the project state tracking in project_state.json.

### Common Commands

- Create a new task: `task-master-ai create "Task description"`
- List all tasks: `task-master-ai list`
- Mark a task as complete: `task-master-ai complete <task-id>`
- Update task status: `task-master-ai update <task-id> --status=Running`

For more commands and options, run:
```
task-master-ai --help
```