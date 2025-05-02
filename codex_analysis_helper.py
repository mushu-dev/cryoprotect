"""
codex_analysis_helper.py

A helper script for performing codebase analysis using Codex CLI.
This script integrates with your Roocode workflow by providing specialized functions for analyzing
your codebase, creating improvement plans, and delegating tasks to Roocode agents.

Usage:
    python codex_analysis_helper.py analyze-codebase --output analysis_report.md
    python codex_analysis_helper.py create-plan --input analysis_report.md --output improvement_plan.json
    python codex_analysis_helper.py delegate-tasks --input improvement_plan.json

Requirements:
    - codex_cli_wrapper.py must be in the same directory or Python path
    - OPENAI_API_KEY environment variable must be set
"""

import argparse
import json
import os
import sys
import traceback
from typing import Dict, List, Any

try:
    from codex_cli_wrapper import run_codex
except ImportError:
    print("ERROR: codex_cli_wrapper.py not found or import failed.", file=sys.stderr)
    sys.exit(1)

def analyze_codebase(output_file: str, focus_area: str = None) -> Dict:
    """
    Perform deep analysis of the codebase using Codex CLI.
    
    Args:
        output_file: Path to save the analysis report
        focus_area: Optional area to focus on (database, api, security, etc.)
    
    Returns:
        Dict with analysis results and status
    """
    focus_text = f" with special focus on {focus_area}" if focus_area else ""
    
    # Get a list of important files to analyze
    file_list_prompt = f"Analyze this project directory and identify the most important files to examine for a codebase analysis{focus_text}. Return only a JSON array of file paths, no explanation."
    
    try:
        file_list_str = run_codex(file_list_prompt)
        file_list = json.loads(file_list_str)
        
        # Read content of key files
        file_contents = {}
        for file_path in file_list[:15]:  # Limit to 15 files to avoid token limits
            try:
                if os.path.exists(file_path) and os.path.isfile(file_path):
                    with open(file_path, 'r', encoding='utf-8') as f:
                        file_contents[file_path] = f.read()
            except Exception as e:
                print(f"Warning: Could not read file {file_path}: {e}")
        
        # Create analysis prompt with file contents
        analysis_prompt = f"""
        Perform a comprehensive analysis of this codebase{focus_text}. Here are the key files:
        
        {json.dumps(file_contents, indent=2)}
        
        Provide a detailed analysis report including:
        1. Overview of codebase structure and architecture
        2. Identified issues, technical debt, and potential risks
        3. Improvement opportunities and recommendations
        4. Prioritized action items for remediation
        
        Format your response as Markdown.
        """
        
        # Run analysis with Codex CLI
        analysis_result = run_codex(analysis_prompt)
        
        # Save to output file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(analysis_result)
        
        return {
            "status": "SUCCESS",
            "message": f"Codebase analysis completed and saved to {output_file}",
            "file_count": len(file_list),
            "output_file": output_file
        }
    
    except Exception as e:
        return {
            "status": "ERROR",
            "message": f"Failed to analyze codebase: {str(e)}",
            "traceback": traceback.format_exc()
        }

def create_improvement_plan(input_file: str, output_file: str) -> Dict:
    """
    Create a detailed improvement plan based on analysis report.
    
    Args:
        input_file: Path to analysis report
        output_file: Path to save the improvement plan
    
    Returns:
        Dict with plan creation results and status
    """
    try:
        # Read analysis report
        with open(input_file, 'r', encoding='utf-8') as f:
            analysis_content = f.read()
        
        # Create plan prompt
        plan_prompt = f"""
        Based on this codebase analysis report, create a detailed improvement plan:
        
        {analysis_content}
        
        Create a JSON plan with this structure:
        {{
            "parent": {{
                "title": "Improvement Plan Title",
                "description": "Overall plan description"
            }},
            "tasks": [
                {{
                    "id": 1,
                    "title": "Task Title",
                    "description": "Detailed task description",
                    "ownerMode": "appropriate_agent_mode",
                    "priority": "high/medium/low",
                    "effort": 1-5,
                    "dependencies": [list of task IDs this depends on]
                }}
            ]
        }}
        
        Return only the valid JSON object, no explanation.
        """
        
        # Run plan creation with Codex CLI
        plan_result = run_codex(plan_prompt)
        
        # Try to parse as JSON to validate
        improvement_plan = json.loads(plan_result)
        
        # Save to output file
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(improvement_plan, f, indent=2)
        
        return {
            "status": "SUCCESS",
            "message": f"Improvement plan created and saved to {output_file}",
            "task_count": len(improvement_plan.get("tasks", [])),
            "output_file": output_file
        }
    
    except json.JSONDecodeError:
        # If not valid JSON, try to fix it with Codex
        try:
            fix_prompt = f"Fix this JSON so it's valid:\n\n{plan_result}\n\nReturn only the fixed JSON, no explanation."
            fixed_json = run_codex(fix_prompt)
            improvement_plan = json.loads(fixed_json)
            
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(improvement_plan, f, indent=2)
            
            return {
                "status": "SUCCESS",
                "message": f"Improvement plan created (with JSON fixes) and saved to {output_file}",
                "task_count": len(improvement_plan.get("tasks", [])),
                "output_file": output_file
            }
        except Exception as e:
            return {
                "status": "ERROR",
                "message": f"Failed to create a valid JSON plan: {str(e)}",
                "raw_output": plan_result,
                "traceback": traceback.format_exc()
            }
    
    except Exception as e:
        return {
            "status": "ERROR",
            "message": f"Failed to create improvement plan: {str(e)}",
            "traceback": traceback.format_exc()
        }

def delegate_tasks(input_file: str) -> Dict:
    """
    Delegate tasks from the improvement plan to Roocode agents.
    
    Args:
        input_file: Path to improvement plan JSON
    
    Returns:
        Dict with delegation results and status
    """
    try:
        # Read improvement plan
        with open(input_file, 'r', encoding='utf-8') as f:
            plan = json.load(f)
        
        # Save as Roocode-compatible plan
        roo_plan_path = os.path.join(".roo", "codex_plan.json")
        
        # Ensure task structure is compatible with Roocode
        tasks = plan.get("tasks", [])
        for task in tasks:
            # Map ownerMode to a valid Roocode mode if needed
            mode_mapping = {
                "Architect": "architect",
                "PM": "pm",
                "DBA": "db-engineer",
                "Backend": "backend",
                "Frontend": "frontend",
                "QA": "qa",
                "DevOps": "devops",
                "Ops": "devops",
                "Security": "backend",  # Assuming no dedicated security mode
                "Debug": "debug",
                "Documentation": "docs"
            }
            
            if task.get("ownerMode") in mode_mapping:
                task["ownerMode"] = mode_mapping[task.get("ownerMode")]
            
            # Remove fields not supported by Roocode if present
            for field in ["priority", "effort", "dependencies"]:
                if field in task:
                    del task[field]
        
        # Save the Roocode plan
        with open(roo_plan_path, 'w', encoding='utf-8') as f:
            json.dump(plan, f, indent=2)
        
        # Create instructions for Roocode invocation
        instructions = f"""
        # Task Delegation Instructions
        
        The improvement plan has been saved to `.roo/codex_plan.json` with {len(tasks)} tasks.
        
        ## How to Execute in Roocode
        
        Run the following command to start the implementation:
        
        ```
        roo @/.roo/codex_plan.json
        ```
        
        Alternatively, you can execute specific tasks by their ID:
        
        ```
        roo @/.roo/codex_plan.json#3  # Executes task with ID 3
        ```
        
        ## Task Overview
        
        {', '.join([f"#{task['id']}: {task['title']} ({task['ownerMode']})" for task in tasks[:5]])}
        {...more tasks}
        """
        
        # Save instructions
        instructions_path = "codex_delegation_instructions.md"
        with open(instructions_path, 'w', encoding='utf-8') as f:
            f.write(instructions)
        
        return {
            "status": "SUCCESS",
            "message": f"Tasks delegated to Roocode. Plan saved to {roo_plan_path}",
            "instructions": instructions_path,
            "task_count": len(tasks),
            "command": "roo @/.roo/codex_plan.json"
        }
    
    except Exception as e:
        return {
            "status": "ERROR",
            "message": f"Failed to delegate tasks: {str(e)}",
            "traceback": traceback.format_exc()
        }

def main():
    parser = argparse.ArgumentParser(
        description="Codebase analysis and improvement planning using Codex CLI",
        epilog="Integrates with Roocode for task execution"
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Analysis commands")
    
    # analyze-codebase command
    analyze_parser = subparsers.add_parser("analyze-codebase", 
                                          help="Perform deep analysis of the codebase")
    analyze_parser.add_argument("--output", required=True, 
                               help="Output file for analysis report (Markdown)")
    analyze_parser.add_argument("--focus", required=False,
                               help="Optional area to focus on (database, api, security, etc.)")
    
    # create-plan command
    plan_parser = subparsers.add_parser("create-plan",
                                       help="Create improvement plan based on analysis")
    plan_parser.add_argument("--input", required=True,
                            help="Input analysis report file")
    plan_parser.add_argument("--output", required=True,
                            help="Output file for improvement plan (JSON)")
    
    # delegate-tasks command
    delegate_parser = subparsers.add_parser("delegate-tasks",
                                          help="Delegate tasks to Roocode agents")
    delegate_parser.add_argument("--input", required=True,
                               help="Input improvement plan JSON file")
    
    args = parser.parse_args()
    
    if args.command == "analyze-codebase":
        result = analyze_codebase(args.output, args.focus)
    elif args.command == "create-plan":
        result = create_improvement_plan(args.input, args.output)
    elif args.command == "delegate-tasks":
        result = delegate_tasks(args.input)
    else:
        parser.print_help()
        sys.exit(1)
    
    # Print result
    print(f"\n=== {args.command.upper()} RESULT ===")
    for key, value in result.items():
        if key != "traceback":
            print(f"{key}: {value}")
    
    if result.get("status") == "ERROR" and "traceback" in result:
        print("\nTraceback:\n", result["traceback"])
    
    # Exit with appropriate code
    sys.exit(0 if result.get("status") == "SUCCESS" else 1)

if __name__ == "__main__":
    main()
