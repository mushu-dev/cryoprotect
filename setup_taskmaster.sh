#!/bin/bash

# Install TaskMaster AI
npm install -g task-master-ai

# Set environment variables
export ANTHROPIC_API_KEY="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"
export PERPLEXITY_API_KEY="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"
export MODEL="claude-3-7-sonnet-20250219"
export PERPLEXITY_MODEL="sonar-pro"
export MAX_TOKENS="64000"
export TEMPERATURE="0.2"
export DEFAULT_SUBTASKS="5"
export DEFAULT_PRIORITY="medium"

# Run TaskMaster AI
task-master-ai