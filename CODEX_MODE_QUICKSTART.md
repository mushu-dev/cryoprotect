# Codex CLI Mode for Roocode - Quickstart Guide

This guide shows you how to use the new Codex CLI mode in Roocode to analyze your codebase, create improvement plans, and implement changes.

## Prerequisites

1. Make sure you have the necessary tools installed:
   ```bash
   setup_codex_cli.bat
   ```

## Using the Codex Mode

The `codex` mode is now fully integrated into Roocode. You can use it directly with various operations:

### Analysis Operations

```bash
roo "analyze: Perform a deep analysis of the database schema" --mode codex
```

```bash
roo "analyze: Review the codebase for security vulnerabilities" --mode codex
```

```bash
roo "analyze: Identify performance bottlenecks in the API endpoints" --mode codex
```

### Planning Operations

```bash
roo "plan: Create an improvement plan based on the analysis" --mode codex
```

```bash
roo "plan: Design a refactoring strategy for the authentication module" --mode codex
```

### Implementation Operations

```bash
roo "implement: Refactor the database connection pooling" --mode codex
```

```bash
roo "implement: Add error handling to API endpoints" --mode codex
```

### Testing Operations

```bash
roo "test: Create comprehensive tests for the auth module" --mode codex
```

```bash
roo "test: Add integration tests for the API endpoints" --mode codex
```

## Using the Workflow

For a complete end-to-end experience, use the Codex Analysis Workflow:

```bash
roo @workflow:codex-analysis-workflow
```

This will:
1. Analyze your codebase
2. Create an improvement plan
3. Coordinate with other agents
4. Implement critical improvements
5. Test the changes
6. Update documentation

## Tips for Best Results

1. **Be Specific**: Specify exactly what you want to analyze or implement
2. **Focus Areas**: Start with focused areas before analyzing the entire codebase
3. **Review Results**: Always review analysis reports and plans before implementing
4. **Iterative Approach**: Use multiple analysis passes for complex codebases

## Common Operations Examples

### Database Analysis

```bash
roo "analyze: Examine the database schema and identify normalization issues" --mode codex
```

### API Improvement

```bash
roo "analyze: Review API endpoint structure and suggest RESTful improvements" --mode codex
```

### Security Audit

```bash
roo "analyze: Perform a security audit of the authentication implementation" --mode codex
```

### Performance Optimization

```bash
roo "plan: Create a performance optimization plan for database queries" --mode codex
```

## Troubleshooting

If you encounter issues:

1. Verify your API key is set in the .env file
2. Check that the Codex CLI is installed (`npm install -g openai-codex`)
3. Ensure Node.js and npm are properly installed
4. Check terminal output for specific error messages

For more details, see the full `CODEX_CLI_GUIDE.md` documentation.
