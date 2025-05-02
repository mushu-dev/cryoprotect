"""
codex_cli_wrapper.py

A Python wrapper for the OpenAI Codex CLI tool.

This module provides a function to invoke the OpenAI Codex CLI (a Node.js-based code generation/completion tool)
from Python, capturing its output and handling errors robustly. The Codex CLI must be installed via npm and
requires the OPENAI_API_KEY environment variable to be set.

Platform: Windows, macOS, Linux (platform-agnostic)
Dependencies: Node.js, npm, OpenAI Codex CLI (installed globally or available in PATH)

Usage Example:
--------------
from codex_cli_wrapper import run_codex

try:
    result = run_codex("Write a Python function to add two numbers.", output_format="json")
    print("Codex output:", result)
except Exception as e:
    print("Error:", e)

Note:
- The OPENAI_API_KEY environment variable must be set before using this module.
- The Codex CLI must be installed and available in your system PATH.
"""

import os
import subprocess
import sys
from typing import Optional
import shutil
import json

class CodexCLIError(Exception):
    """Custom exception for Codex CLI wrapper errors."""
    pass

def run_codex(prompt: str, output_format: str = "text") -> str:
    """
    Invoke the OpenAI Codex CLI with the given prompt and output format.

    Args:
        prompt (str): The prompt to send to Codex.
        output_format (str): Output format ("text" or "json"). Default is "text".

    Returns:
        str: The output from the Codex CLI.

    Raises:
        CodexCLIError: If the CLI fails, is not installed, or dependencies are missing.
        EnvironmentError: If OPENAI_API_KEY is not set in the environment.
    """
    # Ensure the API key is set
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise EnvironmentError("OPENAI_API_KEY environment variable is not set. Please set it before using Codex CLI.")

    # Resolve the full path to the Codex CLI executable
    codex_path = shutil.which("codex")
    if codex_path is None:
        raise CodexCLIError(
            "Codex CLI is not installed or not found in your PATH. "
            "Ensure Node.js, npm, and the Codex CLI are installed and available."
        )
    
    # Load user config for auto-approve and default model, if provided
    autoapprove_flag = False
    default_model = None
    cfg_path = os.path.join(os.getcwd(), 'config.json')
    if os.path.isfile(cfg_path):
        try:
            with open(cfg_path, 'r', encoding='utf-8') as cfg_file:
                cfg = json.load(cfg_file)
            # Support both camelCase and snake_case keys
            autoapprove_flag = bool(cfg.get('autoapprove', cfg.get('autoApprove', False)))
            default_model = cfg.get('model', cfg.get('modelName', None))
        except Exception:
            # Ignore config errors and proceed with defaults
            pass
    # Build the CLI command using the Codex CLI executable name
    # Insert prompt first, then optional flags
    cli_cmd = ["codex", prompt]
    # Apply default model override
    if default_model:
        cli_cmd.extend(["--model", default_model])
    # Apply auto-approve flag
    if autoapprove_flag:
        cli_cmd.append("--auto-approve")
    # JSON output flag
    if output_format == "json":
        cli_cmd.append("--json")

    # Prepare environment (inherit, but ensure API key is present)
    env = os.environ.copy()
    env["OPENAI_API_KEY"] = api_key

    try:
        # Run the CLI, capture stdout and stderr
        result = subprocess.run(
            cli_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            text=True,
            shell=False  # Platform-agnostic: do not use shell
        )
    except FileNotFoundError as e:
        # codex CLI or Node.js/npm not installed or not in PATH
        raise CodexCLIError(
            "Codex CLI is not installed or not found in your PATH. "
            "Ensure Node.js, npm, and the Codex CLI are installed and available."
        ) from e
    except Exception as e:
        raise CodexCLIError(f"Failed to invoke Codex CLI: {e}") from e

    # Ensure stdout and stderr are strings
    stdout_content = result.stdout or ""
    stderr_content = result.stderr or ""
    if result.returncode != 0:
        # Include stderr and stdout in the error message
        raise CodexCLIError(
            f"Codex CLI failed with exit code {result.returncode}.\n"
            f"STDERR:\n{stderr_content.strip()}\n"
            f"STDOUT:\n{stdout_content.strip()}"
        )

    return stdout_content.strip()

if __name__ == "__main__":
    # Example usage
    prompt = "Write a Python function to add two numbers."
    try:
        output = run_codex(prompt, output_format="text")
        print("Codex CLI output:\n", output)
    except Exception as err:
        print("Error invoking Codex CLI:", err, file=sys.stderr)