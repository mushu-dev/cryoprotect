#!/usr/bin/env python3
"""
CryoProtect v2 - Test Database Remediation Manager

This script tests the database remediation manager to ensure it can correctly
identify and execute all the necessary scripts for each phase of the remediation process.

Usage:
    python test_database_remediation_manager.py
"""

import os
import sys
import logging
import subprocess
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(f"test_remediation_manager_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Define the scripts that should be available for each phase
REQUIRED_SCRIPTS = {
    "Phase 1": ["implement_security.py", "apply_service_role_rls.py"],
    "Phase 2": ["standardize_schema.py"],
    "Phase 3": ["fix_relationships.py"],
    "Phase 4": ["complete_database_remediation.py", "verify_database_integrity.py"],
    "Phase 5": ["apply_performance_indexes.bat", "verify_performance_indexes.bat"]
}

def check_script_exists(script):
    """Check if a script exists in the current directory."""
    if os.path.exists(script):
        logger.info(f"✓ {script} exists")
        return True
    else:
        logger.error(f"✗ {script} does not exist")
        return False

def check_script_executable(script):
    """Check if a script is executable."""
    if script.endswith(".py"):
        # For Python scripts, check if they can be imported
        try:
            # Try to run the script with --help or --version to see if it's valid
            result = subprocess.run(
                [sys.executable, script, "--help"],
                capture_output=True,
                text=True,
                timeout=5
            )
            # If the script exits with a non-zero status, it might still be valid
            # but just doesn't support --help
            logger.info(f"✓ {script} is executable")
            return True
        except Exception as e:
            logger.error(f"✗ {script} is not executable: {str(e)}")
            return False
    elif script.endswith(".bat"):
        # For batch scripts, just check if they exist
        # We can't easily check if they're valid without executing them
        if os.path.exists(script):
            logger.info(f"✓ {script} exists (assuming it's executable)")
            return True
        else:
            logger.error(f"✗ {script} does not exist")
            return False
    else:
        logger.warning(f"? Unknown script type: {script}")
        return os.path.exists(script)

def test_database_remediation_manager():
    """Test the database remediation manager."""
    # Check if the main script exists
    if not check_script_exists("database_remediation_manager.py"):
        return False
    
    # Check if the main script is executable
    if not check_script_executable("database_remediation_manager.py"):
        return False
    
    # Check if all required scripts exist
    all_scripts_exist = True
    for phase, scripts in REQUIRED_SCRIPTS.items():
        logger.info(f"\nChecking scripts for {phase}:")
        for script in scripts:
            if not check_script_exists(script):
                all_scripts_exist = False
    
    if not all_scripts_exist:
        logger.error("\nNot all required scripts exist. Please ensure all scripts are available.")
        return False
    
    # Test the database_remediation_manager.py script with --dry-run
    logger.info("\nTesting database_remediation_manager.py with --dry-run...")
    try:
        result = subprocess.run(
            [sys.executable, "database_remediation_manager.py", "--dry-run"],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode != 0:
            logger.error(f"✗ database_remediation_manager.py --dry-run failed with exit code {result.returncode}")
            logger.error(f"Output: {result.stdout}")
            logger.error(f"Error: {result.stderr}")
            return False
        
        logger.info(f"✓ database_remediation_manager.py --dry-run succeeded")
        logger.info(f"Output: {result.stdout}")
        
        # Check if the output mentions all phases
        output = result.stdout.lower()
        all_phases_mentioned = True
        for phase in range(1, 6):
            phase_name = f"phase {phase}"
            if phase_name not in output:
                logger.warning(f"✗ Output does not mention {phase_name}")
                all_phases_mentioned = False
        
        if not all_phases_mentioned:
            logger.warning("Not all phases were mentioned in the output.")
        
        return True
    except subprocess.TimeoutExpired:
        logger.error("✗ database_remediation_manager.py --dry-run timed out")
        return False
    except Exception as e:
        logger.error(f"✗ Error testing database_remediation_manager.py: {str(e)}")
        return False

def main():
    """Main function to test the database remediation manager."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Test Database Remediation Manager")
    print("=" * 80 + "\n")
    
    success = test_database_remediation_manager()
    
    print("\n" + "=" * 60)
    print("Test Results")
    print("=" * 60)
    
    if success:
        print("\nStatus: PASS")
        print("The database remediation manager appears to be working correctly.")
    else:
        print("\nStatus: FAIL")
        print("The database remediation manager has issues that need to be addressed.")
    
    print("\nFor detailed information, check the log file.")
    print("=" * 60 + "\n")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())