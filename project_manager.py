"""
project_manager.py

CLI Project Manager for the roo codebase, orchestrating project management tasks via Codex CLI.

Features:
- generate-code: Generate code for a specified module or feature.
- update-docs: Generate or update documentation.
- create-task: Generate a project task or plan.
- run-tests: (Optionally) Generate and run tests.

Each command uses codex_cli_wrapper.run_codex() to invoke Codex CLI, applies the result (writes files, updates docs, runs tests), and prints a summary.

Extensible for future project management actions.

Usage Examples:
    python project_manager.py generate-code --feature "user authentication" --output "auth_module.py"
    python project_manager.py update-docs --module "api" --output "README_API.md"
    python project_manager.py create-task --description "Refactor database layer for performance"
    python project_manager.py run-tests --module "api"

Requirements:
- codex_cli_wrapper.py must be present in the same directory or Python path.

How to extend:
- Add a new function for the action.
- Register it as a subparser command in main().
- Follow the pattern of argument parsing, Codex CLI invocation, output handling, and summary printing.

Author: Roo (AI)
"""

import argparse
import sys
import os
import traceback

# Load environment variables from .env so they are available to this script and any subprocesses (e.g., Codex CLI).
from dotenv import load_dotenv
load_dotenv()
print("DEBUG: OPENAI_API_KEY loaded:", "OPENAI_API_KEY" in os.environ)

try:
    from codex_cli_wrapper import run_codex
except ImportError:
    print("ERROR: codex_cli_wrapper.py not found or import failed.", file=sys.stderr)
    sys.exit(1)

def generate_code(feature, output):
    """
    Generate code for a specified feature using Codex CLI and write to output file.

    Args:
        feature (str): Description of the feature/module to generate.
        output (str): Output file path for generated code.

    Returns:
        dict: Result summary.
    """
    prompt = f"Generate Python code for the following feature/module in the roo codebase:\n\nFeature: {feature}\n\nProvide only the code, no explanations."
    try:
        code = run_codex(prompt)
        if not code or not isinstance(code, str):
            raise ValueError("Codex CLI returned empty or invalid output.")
        with open(output, "w", encoding="utf-8") as f:
            f.write(code)
        return {
            "status": "SUCCESS",
            "summary": f"Code for feature '{feature}' generated and written to '{output}'.",
            "output_file": output
        }
    except Exception as e:
        return {
            "status": "ERROR",
            "summary": f"Failed to generate code for feature '{feature}': {e}",
            "traceback": traceback.format_exc()
        }

def update_docs(module, output):
    """
    Generate or update documentation for a module using Codex CLI.

    Args:
        module (str): Module or feature to document.
        output (str): Output file path for documentation.

    Returns:
        dict: Result summary.
    """
    prompt = f"Generate or update Markdown documentation for the following module/feature in the roo codebase:\n\nModule: {module}\n\nProvide only the documentation content."
    try:
        docs = run_codex(prompt)
        if not docs or not isinstance(docs, str):
            raise ValueError("Codex CLI returned empty or invalid output.")
        with open(output, "w", encoding="utf-8") as f:
            f.write(docs)
        return {
            "status": "SUCCESS",
            "summary": f"Documentation for module '{module}' generated/updated and written to '{output}'.",
            "output_file": output
        }
    except Exception as e:
        return {
            "status": "ERROR",
            "summary": f"Failed to update documentation for module '{module}': {e}",
            "traceback": traceback.format_exc()
        }

def create_task(description):
    """
    Generate a project task or plan using Codex CLI.

    Args:
        description (str): Description of the task or plan.

    Returns:
        dict: Result summary.
    """
    prompt = f"Create a detailed project management task or plan for the following description in the roo codebase:\n\nDescription: {description}\n\nProvide a Markdown-formatted task or plan."
    try:
        task = run_codex(prompt)
        if not task or not isinstance(task, str):
            raise ValueError("Codex CLI returned empty or invalid output.")
        print(task)
        return {
            "status": "SUCCESS",
            "summary": "Project task/plan generated and printed above.",
            "task": task
        }
    except Exception as e:
        return {
            "status": "ERROR",
            "summary": f"Failed to create project task/plan: {e}",
            "traceback": traceback.format_exc()
        }

def run_tests(module=None):
    """
    Generate and (optionally) run tests for a module using Codex CLI.

    Args:
        module (str): Module or feature to test (optional).

    Returns:
        dict: Result summary.
    """
    prompt = "Generate Python tests"
    if module:
        prompt += f" for the following module/feature in the roo codebase:\n\nModule: {module}\n"
    else:
        prompt += " for the roo codebase.\n"
    prompt += "Provide only the test code, no explanations."
    try:
        test_code = run_codex(prompt)
        if not test_code or not isinstance(test_code, str):
            raise ValueError("Codex CLI returned empty or invalid output.")
        # Write tests to a temporary file
        test_file = f"test_{module or 'roo'}.py"
        with open(test_file, "w", encoding="utf-8") as f:
            f.write(test_code)
        # Attempt to run the tests
        import subprocess
        result = subprocess.run([sys.executable, test_file], capture_output=True, text=True)
        test_output = result.stdout + "\n" + result.stderr
        print(test_output)
        return {
            "status": "SUCCESS" if result.returncode == 0 else "ERROR",
            "summary": f"Tests {'passed' if result.returncode == 0 else 'failed'} for '{module or 'roo'}'. See output above.",
            "test_file": test_file,
            "test_output": test_output
        }
    except Exception as e:
        return {
            "status": "ERROR",
            "summary": f"Failed to generate or run tests: {e}",
            "traceback": traceback.format_exc()
        }

def main():
    parser = argparse.ArgumentParser(
        description="Codex CLI Project Manager for the roo codebase.",
        epilog="See usage examples in the module docstring."
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Project management commands")

    # generate-code
    parser_code = subparsers.add_parser("generate-code", help="Generate code for a specified module or feature.")
    parser_code.add_argument("--feature", required=True, help="Feature or module description to generate code for.")
    parser_code.add_argument("--output", required=True, help="Output file path for generated code.")

    # update-docs
    parser_docs = subparsers.add_parser("update-docs", help="Generate or update documentation for a module.")
    parser_docs.add_argument("--module", required=True, help="Module or feature to document.")
    parser_docs.add_argument("--output", required=True, help="Output file path for documentation.")

    # create-task
    parser_task = subparsers.add_parser("create-task", help="Generate a project task or plan.")
    parser_task.add_argument("--description", required=True, help="Description of the task or plan.")

    # run-tests
    parser_tests = subparsers.add_parser("run-tests", help="Generate and run tests for a module or the codebase.")
    parser_tests.add_argument("--module", required=False, help="Module or feature to test (optional).")

    args = parser.parse_args()

    if args.command == "generate-code":
        result = generate_code(args.feature, args.output)
    elif args.command == "update-docs":
        result = update_docs(args.module, args.output)
    elif args.command == "create-task":
        result = create_task(args.description)
    elif args.command == "run-tests":
        result = run_tests(args.module)
    else:
        parser.print_help()
        sys.exit(1)

    # Print summary/result
    print("\n=== RESULT ===")
    for k, v in result.items():
        if k != "traceback":
            print(f"{k}: {v}")
    if result.get("status") == "ERROR" and "traceback" in result:
        print("\nTraceback:\n", result["traceback"])

if __name__ == "__main__":
    main()

"""
How to extend project_manager.py:

1. Define a new function for your project management action, following the pattern of argument parsing, Codex CLI invocation, output handling, and summary printing.
2. Register the new function as a subparser command in main(), specifying required/optional arguments.
3. Add a conditional branch in main() to call your new function when the command is invoked.
4. Document the new command in the module docstring and usage examples.

This design ensures the CLI is easily extensible for future project management actions.
"""