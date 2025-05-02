"""
codex_mode_handler.py

A dedicated handler for the Codex CLI mode in Roocode.
This script is designed to be executed by Roocode when the codex mode is activated.

Usage (from Roocode):
    roo "Your prompt" --mode codex

Environment:
    - OPENAI_API_KEY must be set
    - codex_cli_wrapper.py must be in the Python path
"""

import os
import sys
import json
import subprocess
import traceback

# Add the parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from codex_cli_wrapper import run_codex
except ImportError:
    print("ERROR: codex_cli_wrapper.py not found or import failed.")
    sys.exit(1)

def handle_codex_request(prompt, workspace_dir=None, context=None):
    """
    Handle a request to the Codex mode.
    
    Args:
        prompt (str): The user's prompt
        workspace_dir (str): The workspace directory
        context (dict): Additional context
        
    Returns:
        dict: The response
    """
    try:
        # Determine the operation type from the prompt
        operation = determine_operation(prompt)
        
        if operation == "analyze":
            return handle_analysis(prompt, workspace_dir, context)
        elif operation == "plan":
            return handle_planning(prompt, workspace_dir, context)
        elif operation == "implement":
            return handle_implementation(prompt, workspace_dir, context)
        elif operation == "test":
            return handle_testing(prompt, workspace_dir, context)
        else:
            # For other operations, just pass directly to Codex CLI
            result = run_codex(prompt)
            return {
                "status": "success",
                "result": result,
                "next_steps": [
                    "Review the results above",
                    "For more specific operations, try prefixing your prompt with 'analyze:', 'plan:', 'implement:', or 'test:'"
                ]
            }
    except Exception as e:
        return {
            "status": "error",
            "message": str(e),
            "traceback": traceback.format_exc()
        }

def determine_operation(prompt):
    """
    Determine the operation type from the prompt.
    
    Args:
        prompt (str): The user's prompt
        
    Returns:
        str: The operation type
    """
    prompt_lower = prompt.lower()
    
    if prompt_lower.startswith("analyze:") or "analyze" in prompt_lower:
        return "analyze"
    elif prompt_lower.startswith("plan:") or "plan" in prompt_lower:
        return "plan"
    elif prompt_lower.startswith("implement:") or "implement" in prompt_lower:
        return "implement"
    elif prompt_lower.startswith("test:") or "test" in prompt_lower:
        return "test"
    else:
        return "general"

def handle_analysis(prompt, workspace_dir, context):
    """
    Handle an analysis request.
    
    Args:
        prompt (str): The user's prompt
        workspace_dir (str): The workspace directory
        context (dict): Additional context
        
    Returns:
        dict: The response
    """
    # Extract focus area from prompt if any
    focus_area = None
    if "focus on" in prompt.lower():
        focus_area = prompt.lower().split("focus on")[1].split()[0].strip("., ")
    
    # Create an enhanced prompt for Codex
    enhanced_prompt = f"""
    Perform a comprehensive analysis of the codebase based on this request:
    
    {prompt}
    
    Focus your analysis on:
    1. Architecture and structure
    2. Code quality and patterns
    3. Potential issues and improvements
    4. Implementation recommendations
    
    Format your response as a structured Markdown report with sections.
    """
    
    # Run Codex CLI
    result = run_codex(enhanced_prompt)
    
    # Save to file
    output_file = os.path.join(workspace_dir or ".", "codex_analysis_report.md")
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(result)
    
    return {
        "status": "success",
        "result": result,
        "output_file": output_file,
        "next_steps": [
            f"Analysis saved to {output_file}",
            "To create an improvement plan based on this analysis, use: plan: create plan from analysis",
            "To focus on a specific area, use: analyze: focus on [area]"
        ]
    }

def handle_planning(prompt, workspace_dir, context):
    """
    Handle a planning request.
    
    Args:
        prompt (str): The user's prompt
        workspace_dir (str): The workspace directory
        context (dict): Additional context
        
    Returns:
        dict: The response
    """
    # Check if planning from analysis
    from_analysis = "from analysis" in prompt.lower()
    
    if from_analysis:
        # Use the latest analysis report
        analysis_file = os.path.join(workspace_dir or ".", "codex_analysis_report.md")
        if os.path.exists(analysis_file):
            with open(analysis_file, "r", encoding="utf-8") as f:
                analysis = f.read()
            
            enhanced_prompt = f"""
            Based on the following analysis report, create a detailed improvement plan:
            
            {analysis}
            
            Create a JSON plan with the following structure:
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
                        "ownerMode": "appropriate_roocode_mode"
                    }}
                ]
            }}
            
            Assign each task to an appropriate Roocode mode from: pm, architect, db-engineer, backend, frontend, qa, devops, docs, debug
            
            Return only valid JSON, no explanation or extra text.
            """
        else:
            return {
                "status": "error",
                "message": "No analysis report found. Please run an analysis first."
            }
    else:
        # Create plan directly from prompt
        enhanced_prompt = f"""
        Create a detailed improvement plan based on this request:
        
        {prompt}
        
        Create a JSON plan with the following structure:
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
                    "ownerMode": "appropriate_roocode_mode"
                }}
            ]
        }}
        
        Assign each task to an appropriate Roocode mode from: pm, architect, db-engineer, backend, frontend, qa, devops, docs, debug
        
        Return only valid JSON, no explanation or extra text.
        """
    
    # Run Codex CLI
    result = run_codex(enhanced_prompt)
    
    # Try to parse as JSON
    try:
        plan = json.loads(result)
        
        # Save to file
        output_file = os.path.join(workspace_dir or ".", ".roo", "codex_plan.json")
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(plan, f, indent=2)
        
        return {
            "status": "success",
            "result": f"Created improvement plan with {len(plan.get('tasks', []))} tasks",
            "output_file": output_file,
            "next_steps": [
                f"Plan saved to {output_file}",
                "To execute this plan in Roocode, run: roo @/.roo/codex_plan.json",
                "To execute a specific task, run: roo @/.roo/codex_plan.json#[task_id]"
            ]
        }
    except json.JSONDecodeError:
        # Try to fix JSON with Codex
        fix_prompt = f"Fix this JSON so it's valid:\n\n{result}\n\nReturn only the fixed JSON, no explanation."
        try:
            fixed_result = run_codex(fix_prompt)
            plan = json.loads(fixed_result)
            
            # Save to file
            output_file = os.path.join(workspace_dir or ".", ".roo", "codex_plan.json")
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, "w", encoding="utf-8") as f:
                json.dump(plan, f, indent=2)
            
            return {
                "status": "success",
                "result": f"Created improvement plan with {len(plan.get('tasks', []))} tasks (JSON required fixing)",
                "output_file": output_file,
                "next_steps": [
                    f"Plan saved to {output_file}",
                    "To execute this plan in Roocode, run: roo @/.roo/codex_plan.json",
                    "To execute a specific task, run: roo @/.roo/codex_plan.json#[task_id]"
                ]
            }
        except:
            return {
                "status": "error",
                "message": "Failed to create a valid JSON plan",
                "raw_result": result
            }

def handle_implementation(prompt, workspace_dir, context):
    """
    Handle an implementation request.
    
    Args:
        prompt (str): The user's prompt
        workspace_dir (str): The workspace directory
        context (dict): Additional context
        
    Returns:
        dict: The response
    """
    # Enhanced prompt for implementation
    enhanced_prompt = f"""
    Implement code based on this request:
    
    {prompt}
    
    Provide complete, well-documented code that is ready to use.
    Include any necessary imports and error handling.
    """
    
    # Run Codex CLI
    result = run_codex(enhanced_prompt)
    
    # Try to determine the appropriate file name and extension
    file_name = "codex_implementation"
    extension = ".py"  # Default to Python
    
    # Look for language hints in the prompt
    languages = {
        "python": ".py",
        "javascript": ".js",
        "typescript": ".ts",
        "html": ".html",
        "css": ".css",
        "sql": ".sql",
        "bash": ".sh",
        "batch": ".bat"
    }
    
    for lang, ext in languages.items():
        if lang in prompt.lower():
            extension = ext
            break
    
    # Save to file
    output_file = os.path.join(workspace_dir or ".", f"{file_name}{extension}")
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(result)
    
    return {
        "status": "success",
        "result": result,
        "output_file": output_file,
        "next_steps": [
            f"Implementation saved to {output_file}",
            "To test this implementation, use: test: test the implementation",
            "To modify this code, use: implement: modify [description of changes]"
        ]
    }

def handle_testing(prompt, workspace_dir, context):
    """
    Handle a testing request.
    
    Args:
        prompt (str): The user's prompt
        workspace_dir (str): The workspace directory
        context (dict): Additional context
        
    Returns:
        dict: The response
    """
    # Enhanced prompt for testing
    enhanced_prompt = f"""
    Create comprehensive tests based on this request:
    
    {prompt}
    
    Provide complete test code with:
    1. Test setup and teardown
    2. Test cases for normal operation
    3. Test cases for edge cases
    4. Test cases for error conditions
    
    Use appropriate testing frameworks and best practices.
    """
    
    # Run Codex CLI
    result = run_codex(enhanced_prompt)
    
    # Save to file
    output_file = os.path.join(workspace_dir or ".", "codex_tests.py")
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(result)
    
    return {
        "status": "success",
        "result": result,
        "output_file": output_file,
        "next_steps": [
            f"Tests saved to {output_file}",
            "To run these tests, use: roo 'Run the Codex tests' --mode qa"
        ]
    }

def main():
    """Main entry point when script is executed directly"""
    if len(sys.argv) < 2:
        print("Usage: python codex_mode_handler.py 'prompt' [workspace_dir]")
        sys.exit(1)
    
    prompt = sys.argv[1]
    workspace_dir = sys.argv[2] if len(sys.argv) > 2 else os.getcwd()
    
    result = handle_codex_request(prompt, workspace_dir)
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
