#!/usr/bin/env python3
"""
CryoProtect Analyzer API - Standardization Runner

This script runs the API audit and applies standardization to all API endpoints.
It provides a convenient way to standardize the API in one step.

Usage:
    python run_api_standardization.py [--dry-run] [--verbose] [--skip-audit] [--skip-apply]

Options:
    --dry-run       Show changes without applying them
    --verbose       Enable verbose output
    --skip-audit    Skip the audit step
    --skip-apply    Skip the apply step
"""

import os
import sys
import argparse
import logging
import subprocess
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def run_audit(output_file: str, verbose: bool = False) -> bool:
    """
    Run the API audit.
    
    Args:
        output_file: Output file for the audit report
        verbose: Whether to enable verbose output
        
    Returns:
        True if the audit was successful, False otherwise
    """
    logger.info("Running API audit...")
    
    cmd = [sys.executable, 'api_audit.py', '--output-file', output_file]
    
    if verbose:
        cmd.append('--verbose')
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"API audit completed successfully. Report saved to {output_file}")
        
        if verbose:
            logger.info(result.stdout)
        
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"API audit failed: {e}")
        logger.error(e.stdout)
        logger.error(e.stderr)
        return False

def apply_standardization(audit_file: str, dry_run: bool = False, verbose: bool = False) -> bool:
    """
    Apply standardization to API endpoints.
    
    Args:
        audit_file: Audit report file
        dry_run: Whether to show changes without applying them
        verbose: Whether to enable verbose output
        
    Returns:
        True if the standardization was successful, False otherwise
    """
    logger.info("Applying API standardization...")
    
    cmd = [sys.executable, 'apply_api_standardization.py', '--audit-file', audit_file]
    
    if dry_run:
        cmd.append('--dry-run')
    
    if verbose:
        cmd.append('--verbose')
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        if dry_run:
            logger.info("API standardization dry run completed successfully")
        else:
            logger.info("API standardization applied successfully")
        
        if verbose:
            logger.info(result.stdout)
        
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"API standardization failed: {e}")
        logger.error(e.stdout)
        logger.error(e.stderr)
        return False

def generate_openapi_docs(verbose: bool = False) -> bool:
    """
    Generate OpenAPI documentation.
    
    Args:
        verbose: Whether to enable verbose output
        
    Returns:
        True if the documentation generation was successful, False otherwise
    """
    logger.info("Generating OpenAPI documentation...")
    
    # Create docs directory if it doesn't exist
    os.makedirs('docs', exist_ok=True)
    
    # Import the OpenAPI generator
    try:
        from api.openapi import generate_openapi_files
        
        # Generate OpenAPI files
        generate_openapi_files()
        
        logger.info("OpenAPI documentation generated successfully")
        return True
    except Exception as e:
        logger.error(f"OpenAPI documentation generation failed: {e}")
        return False

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Run API endpoint standardization')
    parser.add_argument('--dry-run', action='store_true', help='Show changes without applying them')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('--skip-audit', action='store_true', help='Skip the audit step')
    parser.add_argument('--skip-apply', action='store_true', help='Skip the apply step')
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Create a timestamp for the audit report
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    audit_file = f'api_audit_report_{timestamp}.md'
    
    logger.info("Starting API endpoint standardization process")
    
    # Run the audit
    if not args.skip_audit:
        if not run_audit(audit_file, args.verbose):
            logger.error("API audit failed. Aborting standardization process.")
            return 1
    else:
        logger.info("Skipping audit step")
        # Use the most recent audit file
        audit_files = [f for f in os.listdir('.') if f.startswith('api_audit_report_') and f.endswith('.md')]
        if audit_files:
            audit_files.sort(reverse=True)
            audit_file = audit_files[0]
            logger.info(f"Using existing audit file: {audit_file}")
        else:
            logger.error("No existing audit file found. Please run the audit first.")
            return 1
    
    # Apply standardization
    if not args.skip_apply:
        if not apply_standardization(audit_file, args.dry_run, args.verbose):
            logger.error("API standardization failed.")
            return 1
    else:
        logger.info("Skipping apply step")
    
    # Generate OpenAPI documentation
    if not args.skip_apply and not args.dry_run:
        if not generate_openapi_docs(args.verbose):
            logger.error("OpenAPI documentation generation failed.")
            return 1
    
    logger.info("API endpoint standardization process completed successfully")
    
    # Print next steps
    if not args.dry_run and not args.skip_apply:
        logger.info("\nNext steps:")
        logger.info("1. Review the changes made to the API endpoints")
        logger.info("2. Run the tests to ensure everything works correctly")
        logger.info("3. Update any client code to handle the standardized responses")
        logger.info("4. Check the OpenAPI documentation at /api/v1/docs/")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())