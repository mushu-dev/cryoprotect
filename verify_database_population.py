"""
Verification script for the database population module.

This script tests the imports and CLI entry point for the database population module.
"""

import sys
import importlib
import subprocess
from typing import Dict, List, Tuple

def test_imports() -> List[Tuple[str, bool, str]]:
    """
    Test importing the required modules.
    
    Returns:
        List of tuples containing (module_name, success, error_message)
    """
    modules_to_test = [
        "database",
        "database.population",
        "database.population.runner",
        "database.population.molecules",
        "database.population.mixtures"
    ]
    
    results = []
    
    for module_name in modules_to_test:
        try:
            importlib.import_module(module_name)
            results.append((module_name, True, ""))
        except Exception as e:
            results.append((module_name, False, str(e)))
    
    return results

def test_cli_entry_point() -> Tuple[bool, str]:
    """
    Test the CLI entry point with --help flag.
    
    Returns:
        Tuple containing (success, output_or_error)
    """
    try:
        result = subprocess.run(
            ["python", "-m", "database.population.runner", "--help"],
            capture_output=True,
            text=True,
            check=True
        )
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        return False, f"Error: {e.stderr}"
    except Exception as e:
        return False, f"Unexpected error: {str(e)}"

def main():
    """
    Main function to run the verification tests.
    """
    print("Verifying database population module integration...")
    print("\n1. Testing imports:")
    
    import_results = test_imports()
    all_imports_successful = True
    
    for module_name, success, error_message in import_results:
        status = "SUCCESS" if success else "ERROR"
        print(f"  - {module_name}: {status}")
        if not success:
            all_imports_successful = False
            print(f"    Error: {error_message}")
    
    print("\n2. Testing CLI entry point:")
    cli_success, cli_output = test_cli_entry_point()
    cli_status = "SUCCESS" if cli_success else "ERROR"
    print(f"  - CLI entry point: {cli_status}")
    
    if not cli_success:
        print(f"    Error: {cli_output}")
    else:
        print("\nCLI Help Output:")
        print("-" * 40)
        print(cli_output)
        print("-" * 40)
    
    # Overall status
    overall_status = "SUCCESS" if all_imports_successful and cli_success else "ERROR"
    
    print("\nSummary:")
    print(f"  - Import tests: {'SUCCESS' if all_imports_successful else 'ERROR'}")
    print(f"  - CLI entry point test: {cli_status}")
    print(f"  - Overall status: {overall_status}")
    
    # Return structured result
    result = {
        "status": overall_status,
        "summary": {
            "imports": import_results,
            "cli_entry_point": cli_status
        },
        "deliverable": cli_output if cli_success else cli_output
    }
    
    print("\nStructured Result:")
    print(f"status: '{overall_status}'")
    print("summary: List of tests performed and any issues encountered")
    print("deliverable: Output of the CLI --help command and any error messages if present")
    
    return result

if __name__ == "__main__":
    main()