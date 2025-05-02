# Practical Guide to Using Codex CLI with Roocode

This guide explains how to use the Codex CLI integration with Roocode for deep codebase analysis, improvement planning, and task delegation in your CryoProtect v2 project.

## Prerequisites

1. **Install Node.js and npm** (if not already installed)
   - Download from: https://nodejs.org/

2. **Install Codex CLI**
   ```bash
   npm install -g openai-codex
   ```

3. **Set up your API key**
   - Add your OpenAI API key to your .env file:
   ```
   OPENAI_API_KEY=sk-your-key-here
   ```

4. **Verify setup**
   - Run the setup script:
   ```bash
   setup_codex_cli.bat
   ```

## Approach 1: Using the Helper Script (Recommended)

The `codex_analysis_helper.py` script provides a structured approach to using Codex CLI for your project:

### Step 1: Analyze Your Codebase

```bash
python codex_analysis_helper.py analyze-codebase --output analysis_report.md
```

This will:
- Identify important files in your project
- Use Codex CLI to analyze code structure and patterns
- Produce a comprehensive analysis report in Markdown format

For focused analysis:
```bash
python codex_analysis_helper.py analyze-codebase --output database_analysis.md --focus database
```

### Step 2: Create an Improvement Plan

```bash
python codex_analysis_helper.py create-plan --input analysis_report.md --output improvement_plan.json
```

This will:
- Read the analysis report
- Create a structured improvement plan with specific tasks
- Assign appropriate agent modes to each task
- Save the plan in JSON format compatible with Roocode

### Step 3: Delegate Tasks to Roocode

```bash
python codex_analysis_helper.py delegate-tasks --input improvement_plan.json
```

This will:
- Convert the improvement plan to Roocode format
- Save it as `.roo/codex_plan.json`
- Generate instructions for running the plan in Roocode

### Step 4: Execute with Roocode

```bash
roo @/.roo/codex_plan.json
```

This will start the implementation of the plan with the appropriate specialized agents.

## Approach 2: Using the Codex Analysis Mode

The newly added `codex-analysis` mode in Roocode provides an integrated approach:

```bash
roo "Perform a deep analysis of the codebase structure and suggest improvements" --mode codex-analysis
```

Use this approach for:
- Interactive analysis sessions
- Quick insights without running the full workflow
- Delegating specific tasks to other modes

### Example Prompts for the Codex Analysis Mode

- `"Analyze the database schema and identify optimization opportunities"`
- `"Review the authentication implementation and suggest security improvements"`
- `"Find performance bottlenecks in the API endpoints"`
- `"Create a plan to refactor the codebase into a more modular structure"`

## Approach 3: Direct Use of Codex CLI

For quick, one-off tasks:

```bash
codex "Analyze this code and suggest improvements: [paste code here]" > suggestions.txt
```

## Best Practices

1. **Start Small**: Begin with analyzing specific components before tackling the entire codebase
2. **Be Specific**: Provide focused prompts rather than overly general ones
3. **Review Before Implementing**: Always review analysis reports and plans before execution
4. **Combine Approaches**: Use both Codex CLI and other Roocode modes as needed
5. **Iterate**: Run multiple analysis passes with different focuses

## Troubleshooting

- **API Key Issues**: Ensure your `OPENAI_API_KEY` is correctly set in the .env file or environment
- **Installation Problems**: Make sure Node.js and npm are properly installed and in your PATH
- **Output Format Issues**: If JSON is malformed, try regenerating with a more specific prompt

## Using with Existing Project Structure

This integration works seamlessly with your existing Roocode setup, allowing you to:

1. Use Codex for deep analysis beyond Roocode's standard capabilities
2. Generate detailed improvement plans based on code-aware analysis
3. Delegate specific implementation tasks to specialized Roocode modes
4. Track progress through your existing workflows

For questions or issues, refer to the documentation in the `codex_cli_wrapper.py` and `codex_analysis_helper.py` files.
