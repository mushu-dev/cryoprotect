@echo off
ECHO Installing TaskMaster AI...

:: Install TaskMaster AI
npm install -g task-master-ai

:: Set environment variables
SET ANTHROPIC_API_KEY=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
SET PERPLEXITY_API_KEY=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
SET MODEL=claude-3-7-sonnet-20250219
SET PERPLEXITY_MODEL=sonar-pro
SET MAX_TOKENS=64000
SET TEMPERATURE=0.2
SET DEFAULT_SUBTASKS=5
SET DEFAULT_PRIORITY=medium

:: Run TaskMaster AI
task-master-ai

PAUSE