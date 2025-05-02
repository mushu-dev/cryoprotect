#!/usr/bin/env python3
"""
Unified Test Coverage Runner for CryoProtect v2

- Runs all tests using pytest and measures code coverage for the entire codebase.
- Generates an HTML coverage report in reports/htmlcov/.
- Prints a coverage summary by module to the terminal.
- Exits with nonzero status if tests fail.

Usage:
    python run_coverage.py
"""

import os
import sys
import subprocess

def main():
    # Ensure the reports/htmlcov directory exists
    os.makedirs("reports/htmlcov", exist_ok=True)

    # Build the pytest command
    # Add all main source directories and files to --cov
    cov_targets = [
        "api",
        "app.py",
        "auth_config.py",
        "config.py",
        "connection_pool_wrapper.py",
        "database_cli.py",
        "implement_connection_pooling.py",
        "logging_config.py",
        "notify.py",
        "rate_limit_config.py",
        "populate_database.py",
        "populate_experiments.py",
        "populate_mixtures.py",
        "populate_molecular_properties.py",
        "populate_predictions.py",
        "populate_molecules.py",
        "resources",
        "models"
    ]
    cov_args = []
    for target in cov_targets:
        if os.path.exists(target) or os.path.isdir(target):
            cov_args.extend(["--cov", target])

    # Run pytest with coverage
    cmd = [
        sys.executable, "-m", "pytest",
        *cov_args,
        "--cov-report=term",
        "--cov-report=html:reports/htmlcov"
    ]
    print("Running unified test coverage with command:")
    print(" ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print("Some tests failed or coverage run did not complete successfully.")
        sys.exit(result.returncode)
    else:
        print("\nCoverage HTML report generated in reports/htmlcov/index.html")
        print("Open this file in your browser to view detailed coverage by module and file.")

if __name__ == "__main__":
    main()