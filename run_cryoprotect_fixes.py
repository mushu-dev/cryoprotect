#!/usr/bin/env python3
"""
CryoProtect v2 - Master Integration Script

This script orchestrates the execution of all fixes for the CryoProtect project in the correct order:
1. Database backup (create_database_backup.py)
2. Schema standardization (standardize_schema.py)
3. Security implementation (implement_security.py)
4. Relationship fixes (fix_relationships.py)
5. API integration fixes (fix_api_integration.py)
6. Testing (run_tests.py)

The script provides options to:
- Run all fixes in sequence
- Run individual fixes
- Verify the state after each fix
- Rollback changes if needed

Usage:
    python run_cryoprotect_fixes.py [options]

Options:
    --all                  Run all fixes in sequence (default)
    --backup-only          Run only database backup
    --schema-only          Run only schema standardization
    --security-only        Run only security implementation
    --relationships-only   Run only relationship fixes
    --api-only             Run only API integration fixes
    --test-only            Run only tests
    --verify               Verify the state after each fix
    --rollback             Rollback changes if a fix fails
    --dry-run              Show what would be done without making changes
    --report               Generate a comprehensive final report
    --help                 Show this help message and exit
"""

import os
import sys
import json
import time
import argparse
import logging
import importlib.util
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Union, Optional, Tuple

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(f"cryoprotect_fixes_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Fix scripts
FIX_SCRIPTS = {
    "backup": "create_database_backup.py",
    "schema": "standardize_schema.py",
    "security": "implement_security.py",
    "relationships": "fix_relationships.py",
    "api": "fix_api_integration.py",
    "tests": "tests/run_tests.py"
}

# Results of each fix
fix_results = {}

# Backup directory for rollbacks
BACKUP_DIR = Path("backups") / f"master_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="CryoProtect v2 - Master Integration Script")
    
    # Fix selection options
    fix_group = parser.add_argument_group("Fix Selection")
    fix_group.add_argument("--all", action="store_true", help="Run all fixes in sequence (default)")
    fix_group.add_argument("--backup-only", action="store_true", help="Run only database backup")
    fix_group.add_argument("--schema-only", action="store_true", help="Run only schema standardization")
    fix_group.add_argument("--security-only", action="store_true", help="Run only security implementation")
    fix_group.add_argument("--relationships-only", action="store_true", help="Run only relationship fixes")
    fix_group.add_argument("--api-only", action="store_true", help="Run only API integration fixes")
    fix_group.add_argument("--test-only", action="store_true", help="Run only tests")
    
    # Execution options
    exec_group = parser.add_argument_group("Execution Options")
    exec_group.add_argument("--verify", action="store_true", help="Verify the state after each fix")
    exec_group.add_argument("--rollback", action="store_true", help="Rollback changes if a fix fails")
    exec_group.add_argument("--dry-run", action="store_true", help="Show what would be done without making changes")
    exec_group.add_argument("--report", action="store_true", help="Generate a comprehensive final report")
    
    args = parser.parse_args()
    
    # If no specific fix is selected, run all fixes
    if not any([args.backup_only, args.schema_only, args.security_only, 
                args.relationships_only, args.api_only, args.test_only]):
        args.all = True
    
    return args

def create_backup_directory():
    """Create a backup directory for rollbacks."""
    try:
        os.makedirs(BACKUP_DIR, exist_ok=True)
        logger.info(f"Created backup directory: {BACKUP_DIR}")
        return True
    except Exception as e:
        logger.error(f"Error creating backup directory: {str(e)}")
        return False

def backup_file(file_path):
    """Backup a file to the backup directory."""
    try:
        file_path = Path(file_path)
        if file_path.exists():
            # Create subdirectories in the backup directory if needed
            backup_path = BACKUP_DIR / file_path
            os.makedirs(backup_path.parent, exist_ok=True)
            
            # Copy the file
            import shutil
            shutil.copy2(file_path, backup_path)
            logger.info(f"Backed up file: {file_path}")
            return True
        else:
            logger.warning(f"File does not exist, cannot backup: {file_path}")
            return False
    except Exception as e:
        logger.error(f"Error backing up file {file_path}: {str(e)}")
        return False

def run_python_script(script_path, args=None, dry_run=False):
    """Run a Python script with the given arguments."""
    script_path = Path(script_path)
    if not script_path.exists():
        logger.error(f"Script does not exist: {script_path}")
        return False, None
    
    cmd = [sys.executable, str(script_path)]
    if args:
        cmd.extend(args)
    
    if dry_run:
        logger.info(f"[DRY RUN] Would execute: {' '.join(cmd)}")
        return True, "Dry run - no output"
    
    try:
        logger.info(f"Executing: {' '.join(cmd)}")
        start_time = time.time()
        result = subprocess.run(cmd, check=False, capture_output=True, text=True)
        end_time = time.time()
        
        execution_time = round(end_time - start_time, 2)
        
        if result.returncode == 0:
            logger.info(f"Script executed successfully in {execution_time}s: {script_path}")
            return True, result.stdout
        else:
            logger.error(f"Script execution failed in {execution_time}s: {script_path}")
            logger.error(f"Error output: {result.stderr}")
            return False, result.stderr
    except Exception as e:
        logger.error(f"Error executing script {script_path}: {str(e)}")
        return False, str(e)

def run_database_backup(args):
    """Run the database backup script."""
    logger.info("=== Running Database Backup ===")
    
    script_path = FIX_SCRIPTS["backup"]
    script_args = []
    
    if args.dry_run:
        script_args.append("--dry-run")
    
    success, output = run_python_script(script_path, script_args, args.dry_run)
    
    fix_results["backup"] = {
        "success": success,
        "timestamp": datetime.now().isoformat(),
        "output": output[:500] + "..." if output and len(output) > 500 else output,
        "script": script_path
    }
    
    if success:
        logger.info("Database backup completed successfully")
    else:
        logger.error("Database backup failed")
        if args.rollback:
            logger.info("No rollback needed for database backup")
    
    return success

def run_schema_standardization(args):
    """Run the schema standardization script."""
    logger.info("=== Running Schema Standardization ===")
    
    script_path = FIX_SCRIPTS["schema"]
    script_args = []
    
    if args.dry_run:
        script_args.append("--dry-run")
    
    if args.verify:
        script_args.append("--verify")
    
    success, output = run_python_script(script_path, script_args, args.dry_run)
    
    fix_results["schema"] = {
        "success": success,
        "timestamp": datetime.now().isoformat(),
        "output": output[:500] + "..." if output and len(output) > 500 else output,
        "script": script_path
    }
    
    if success:
        logger.info("Schema standardization completed successfully")
    else:
        logger.error("Schema standardization failed")
        if args.rollback:
            logger.info("Running schema rollback...")
            rollback_args = ["--rollback"]
            rollback_success, rollback_output = run_python_script(script_path, rollback_args, args.dry_run)
            if rollback_success:
                logger.info("Schema rollback completed successfully")
            else:
                logger.error("Schema rollback failed")
    
    return success

def run_security_implementation(args):
    """Run the security implementation script."""
    logger.info("=== Running Security Implementation ===")
    
    script_path = FIX_SCRIPTS["security"]
    script_args = []
    
    if args.dry_run:
        script_args.append("--dry-run")
    
    if args.verify:
        script_args.append("--verify-only")
    
    success, output = run_python_script(script_path, script_args, args.dry_run)
    
    fix_results["security"] = {
        "success": success,
        "timestamp": datetime.now().isoformat(),
        "output": output[:500] + "..." if output and len(output) > 500 else output,
        "script": script_path
    }
    
    if success:
        logger.info("Security implementation completed successfully")
    else:
        logger.error("Security implementation failed")
        if args.rollback:
            logger.info("Running security rollback...")
            rollback_args = ["--rollback"]
            rollback_success, rollback_output = run_python_script(script_path, rollback_args, args.dry_run)
            if rollback_success:
                logger.info("Security rollback completed successfully")
            else:
                logger.error("Security rollback failed")
    
    return success

def run_relationship_fixes(args):
    """Run the relationship fixes script."""
    logger.info("=== Running Relationship Fixes ===")
    
    script_path = FIX_SCRIPTS["relationships"]
    script_args = []
    
    if args.dry_run:
        script_args.append("--dry-run")
    
    if args.verify:
        script_args.append("--verify-only")
    
    success, output = run_python_script(script_path, script_args, args.dry_run)
    
    fix_results["relationships"] = {
        "success": success,
        "timestamp": datetime.now().isoformat(),
        "output": output[:500] + "..." if output and len(output) > 500 else output,
        "script": script_path
    }
    
    if success:
        logger.info("Relationship fixes completed successfully")
    else:
        logger.error("Relationship fixes failed")
        if args.rollback:
            logger.info("Running relationship fixes rollback...")
            rollback_args = ["--rollback"]
            rollback_success, rollback_output = run_python_script(script_path, rollback_args, args.dry_run)
            if rollback_success:
                logger.info("Relationship fixes rollback completed successfully")
            else:
                logger.error("Relationship fixes rollback failed")
    
    return success

def run_api_integration_fixes(args):
    """Run the API integration fixes script."""
    logger.info("=== Running API Integration Fixes ===")
    
    script_path = FIX_SCRIPTS["api"]
    script_args = []
    
    if args.dry_run:
        script_args.append("--dry-run")
    
    success, output = run_python_script(script_path, script_args, args.dry_run)
    
    fix_results["api"] = {
        "success": success,
        "timestamp": datetime.now().isoformat(),
        "output": output[:500] + "..." if output and len(output) > 500 else output,
        "script": script_path
    }
    
    if success:
        logger.info("API integration fixes completed successfully")
    else:
        logger.error("API integration fixes failed")
        if args.rollback:
            logger.info("Running API integration fixes rollback...")
            rollback_args = ["--rollback"]
            rollback_success, rollback_output = run_python_script(script_path, rollback_args, args.dry_run)
            if rollback_success:
                logger.info("API integration fixes rollback completed successfully")
            else:
                logger.error("API integration fixes rollback failed")
    
    return success

def run_tests(args):
    """Run the test suite."""
    logger.info("=== Running Tests ===")
    
    script_path = FIX_SCRIPTS["tests"]
    script_args = ["--report"]
    
    if args.verify:
        script_args.append("--verbose")
    
    success, output = run_python_script(script_path, script_args, args.dry_run)
    
    fix_results["tests"] = {
        "success": success,
        "timestamp": datetime.now().isoformat(),
        "output": output[:500] + "..." if output and len(output) > 500 else output,
        "script": script_path
    }
    
    if success:
        logger.info("Tests completed successfully")
    else:
        logger.error("Tests failed")
        # No rollback needed for tests
    
    return success

def verify_state():
    """Verify the state of the system after all fixes."""
    logger.info("=== Verifying System State ===")
    
    # Check if all fixes were successful
    all_successful = all(result.get("success", False) for result in fix_results.values())
    
    if all_successful:
        logger.info("All fixes were applied successfully")
    else:
        logger.warning("Some fixes failed to apply")
        for fix, result in fix_results.items():
            if not result.get("success", False):
                logger.warning(f"Fix '{fix}' failed")
    
    # Additional verification could be added here
    
    return all_successful

def generate_report():
    """Generate a comprehensive final report."""
    logger.info("=== Generating Final Report ===")
    
    report_dir = Path("reports")
    os.makedirs(report_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_path = report_dir / f"cryoprotect_fixes_report_{timestamp}.json"
    text_report_path = report_dir / f"cryoprotect_fixes_report_{timestamp}.txt"
    
    # Create the report data
    report_data = {
        "timestamp": datetime.now().isoformat(),
        "fixes": fix_results,
        "summary": {
            "total_fixes": len(fix_results),
            "successful_fixes": sum(1 for result in fix_results.values() if result.get("success", False)),
            "failed_fixes": sum(1 for result in fix_results.values() if not result.get("success", False))
        }
    }
    
    # Save JSON report
    with open(report_path, "w") as f:
        json.dump(report_data, f, indent=2)
    
    # Create text report
    with open(text_report_path, "w") as f:
        f.write("CryoProtect v2 - Fixes Report\n")
        f.write("==========================\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Summary
        f.write("Summary\n")
        f.write("-------\n")
        f.write(f"Total fixes: {report_data['summary']['total_fixes']}\n")
        f.write(f"Successful fixes: {report_data['summary']['successful_fixes']}\n")
        f.write(f"Failed fixes: {report_data['summary']['failed_fixes']}\n\n")
        
        # Details for each fix
        f.write("Fix Details\n")
        f.write("-----------\n")
        
        for fix, result in fix_results.items():
            f.write(f"\n{fix.upper()}\n")
            f.write(f"{'=' * len(fix)}\n")
            f.write(f"Status: {'SUCCESS' if result.get('success', False) else 'FAILED'}\n")
            f.write(f"Timestamp: {result.get('timestamp', 'N/A')}\n")
            f.write(f"Script: {result.get('script', 'N/A')}\n")
            
            # Include a snippet of the output
            output = result.get("output", "")
            if output:
                f.write("\nOutput snippet:\n")
                f.write("--------------\n")
                lines = output.split("\n")
                for line in lines[:20]:  # Show first 20 lines
                    f.write(f"{line}\n")
                if len(lines) > 20:
                    f.write("...\n")
            
            f.write("\n")
        
        # Before/After Comparisons
        f.write("\nBefore/After Comparisons\n")
        f.write("=======================\n")
        f.write("Database Schema: Changes applied to standardize table names and relationships\n")
        f.write("Security: Row-level security policies and app-specific roles implemented\n")
        f.write("API Integration: Updated to handle Supabase v2.x client responses\n")
        
        # Recommendations
        f.write("\nRecommendations\n")
        f.write("===============\n")
        
        if report_data['summary']['failed_fixes'] > 0:
            f.write("- Address failed fixes before deploying to production\n")
        
        f.write("- Regularly run the test suite to ensure system integrity\n")
        f.write("- Consider implementing automated CI/CD pipeline for future updates\n")
        f.write("- Maintain regular database backups\n")
    
    logger.info(f"Report saved to: {report_path}")
    logger.info(f"Text report saved to: {text_report_path}")
    
    return report_path, text_report_path

def main():
    """Main function to orchestrate the fixes."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Master Integration Script")
    print("=" * 80 + "\n")
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Create backup directory for rollbacks
    if args.rollback and not args.dry_run:
        if not create_backup_directory():
            logger.error("Failed to create backup directory, aborting")
            return 1
    
    # Track overall success
    overall_success = True
    
    # Run the selected fixes
    try:
        # Database backup
        if args.all or args.backup_only:
            backup_success = run_database_backup(args)
            overall_success = overall_success and backup_success
            if not backup_success and not args.all:
                logger.error("Database backup failed, aborting")
                return 1
        
        # Schema standardization
        if (args.all or args.schema_only) and (overall_success or args.all):
            schema_success = run_schema_standardization(args)
            overall_success = overall_success and schema_success
            if not schema_success and not args.all:
                logger.error("Schema standardization failed, aborting")
                return 1
        
        # Security implementation
        if (args.all or args.security_only) and (overall_success or args.all):
            security_success = run_security_implementation(args)
            overall_success = overall_success and security_success
            if not security_success and not args.all:
                logger.error("Security implementation failed, aborting")
                return 1
        
        # Relationship fixes
        if (args.all or args.relationships_only) and (overall_success or args.all):
            relationships_success = run_relationship_fixes(args)
            overall_success = overall_success and relationships_success
            if not relationships_success and not args.all:
                logger.error("Relationship fixes failed, aborting")
                return 1
        
        # API integration fixes
        if (args.all or args.api_only) and (overall_success or args.all):
            api_success = run_api_integration_fixes(args)
            overall_success = overall_success and api_success
            if not api_success and not args.all:
                logger.error("API integration fixes failed, aborting")
                return 1
        
        # Tests
        if (args.all or args.test_only) and (overall_success or args.all):
            test_success = run_tests(args)
            overall_success = overall_success and test_success
        
        # Verify state
        if args.verify:
            verify_success = verify_state()
            overall_success = overall_success and verify_success
        
        # Generate report
        if args.report:
            report_path, text_report_path = generate_report()
        
        # Print summary
        print("\n" + "=" * 60)
        print("CryoProtect v2 - Fix Summary")
        print("=" * 60)
        
        for fix, result in fix_results.items():
            status = "SUCCESS" if result.get("success", False) else "FAILED"
            print(f"{fix.ljust(15)}: {status}")
        
        print("\nOverall Status:", "SUCCESS" if overall_success else "FAILED")
        
        if args.report:
            print(f"\nDetailed report saved to: {text_report_path}")
        
        return 0 if overall_success else 1
    
    except KeyboardInterrupt:
        logger.warning("Process interrupted by user")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())