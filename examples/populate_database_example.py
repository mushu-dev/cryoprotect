#!/usr/bin/env python3
"""
Example script demonstrating how to use the populate_database.py script.

This script shows different ways to run the database population process
with various command-line arguments.
"""

import os
import sys
import subprocess
import argparse
from datetime import datetime

def run_command(command, description):
    """
    Run a command and print its output.
    
    Args:
        command: Command to run
        description: Description of the command
    """
    print(f"\n{'=' * 80}")
    print(f"EXAMPLE: {description}")
    print(f"COMMAND: {' '.join(command)}")
    print(f"{'=' * 80}\n")
    
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )
        
        # Print output in real-time
        for line in process.stdout:
            print(line, end='')
        
        # Wait for process to complete
        process.wait()
        
        print(f"\n{'=' * 80}")
        print(f"EXIT CODE: {process.returncode}")
        print(f"{'=' * 80}\n")
        
        return process.returncode
    except Exception as e:
        print(f"Error running command: {str(e)}")
        return 1

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Run database population examples")
    parser.add_argument("--example", type=int, choices=range(1, 7), help="Run a specific example (1-6)")
    parser.add_argument("--all", action="store_true", help="Run all examples")
    args = parser.parse_args()
    
    # Get the path to the populate_database.py script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    populate_script = os.path.join(os.path.dirname(script_dir), "populate_database.py")
    
    # Create timestamp for unique output directories
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Define examples
    examples = [
        {
            "id": 1,
            "description": "Run the full database population process",
            "command": [
                sys.executable, populate_script,
                "--checkpoint-dir", f"./checkpoints/example1_{timestamp}",
                "--report-dir", f"./reports/example1_{timestamp}"
            ]
        },
        {
            "id": 2,
            "description": "Run only the reference compounds import step",
            "command": [
                sys.executable, populate_script,
                "--steps", "reference",
                "--checkpoint-dir", f"./checkpoints/example2_{timestamp}",
                "--report-dir", f"./reports/example2_{timestamp}"
            ]
        },
        {
            "id": 3,
            "description": "Run all steps except performance optimization",
            "command": [
                sys.executable, populate_script,
                "--skip", "performance",
                "--checkpoint-dir", f"./checkpoints/example3_{timestamp}",
                "--report-dir", f"./reports/example3_{timestamp}"
            ]
        },
        {
            "id": 4,
            "description": "Run with smaller batch size for better memory usage",
            "command": [
                sys.executable, populate_script,
                "--batch-size", "5",
                "--checkpoint-dir", f"./checkpoints/example4_{timestamp}",
                "--report-dir", f"./reports/example4_{timestamp}"
            ]
        },
        {
            "id": 5,
            "description": "Run with verbose logging",
            "command": [
                sys.executable, populate_script,
                "--verbose",
                "--steps", "reference",  # Only run reference step to keep output manageable
                "--checkpoint-dir", f"./checkpoints/example5_{timestamp}",
                "--report-dir", f"./reports/example5_{timestamp}"
            ]
        },
        {
            "id": 6,
            "description": "Resume from checkpoint (requires a previous run)",
            "command": [
                sys.executable, populate_script,
                "--resume",
                "--checkpoint-dir", f"./checkpoints/example6_{timestamp}",
                "--report-dir", f"./reports/example6_{timestamp}"
            ]
        }
    ]
    
    # Run examples
    if args.example:
        # Run a specific example
        example = next((e for e in examples if e["id"] == args.example), None)
        if example:
            run_command(example["command"], example["description"])
        else:
            print(f"Example {args.example} not found")
    elif args.all:
        # Run all examples
        for example in examples:
            run_command(example["command"], example["description"])
    else:
        # Print usage
        print("Please specify an example to run with --example or run all with --all")
        print("\nAvailable examples:")
        for example in examples:
            print(f"  {example['id']}: {example['description']}")

if __name__ == "__main__":
    main()