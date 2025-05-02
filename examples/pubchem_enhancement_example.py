#!/usr/bin/env python3
"""
Example script demonstrating how to use the enhanced PubChem property script.

This example shows how to:
1. Run the property enhancement in dry run mode
2. Run with checkpointing for resumable operations
3. Generate a detailed report
4. Process with custom batch size
"""

import os
import sys
import logging
from datetime import datetime

# Add parent directory to path to import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from enhance_pubchem_properties import enhance_pubchem_properties

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def run_dry_mode_example():
    """Run the enhancer in dry run mode to see what would be updated."""
    logger.info("Running PubChem property enhancement in dry run mode")
    
    # Create a timestamp for the report file
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_file = f"reports/pubchem_dry_run_{timestamp}.json"
    
    # Run the enhancement in dry run mode
    results = enhance_pubchem_properties(
        output_report=report_file,
        batch_size=5,  # Process 5 molecules in parallel
        dry_run=True   # Don't actually update the database
    )
    
    logger.info(f"Dry run completed: {results['molecules_enhanced']} molecules would be enhanced")
    logger.info(f"Report saved to: {report_file}")
    
    return results

def run_with_checkpointing():
    """Run the enhancer with checkpointing for resumable operations."""
    logger.info("Running PubChem property enhancement with checkpointing")
    
    # Create a timestamp for the report file
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_file = f"reports/pubchem_enhancement_{timestamp}.json"
    checkpoint_file = "checkpoints/pubchem_enhancement_example.json"
    
    # Run the enhancement with checkpointing
    results = enhance_pubchem_properties(
        output_report=report_file,
        batch_size=10,  # Process 10 molecules in parallel
        dry_run=False,  # Actually update the database
        checkpoint_file=checkpoint_file  # Use custom checkpoint file
    )
    
    logger.info(f"Enhancement completed: {results['molecules_enhanced']} molecules enhanced")
    logger.info(f"Report saved to: {report_file}")
    logger.info(f"Checkpoint saved to: {checkpoint_file}")
    
    return results

def run_full_enhancement():
    """Run a full enhancement with optimized settings."""
    logger.info("Running full PubChem property enhancement")
    
    # Create a timestamp for the report file
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_file = f"reports/pubchem_full_enhancement_{timestamp}.json"
    
    # Run the enhancement with optimized settings
    results = enhance_pubchem_properties(
        output_report=report_file,
        batch_size=20,  # Process 20 molecules in parallel for better performance
        dry_run=False   # Actually update the database
    )
    
    logger.info(f"Full enhancement completed: {results['molecules_enhanced']} molecules enhanced")
    logger.info(f"Report saved to: {report_file}")
    
    return results

if __name__ == "__main__":
    # Create necessary directories
    os.makedirs("reports", exist_ok=True)
    os.makedirs("checkpoints", exist_ok=True)
    
    # Choose which example to run
    example_type = "dry_run"  # Options: "dry_run", "checkpoint", "full"
    
    if example_type == "dry_run":
        run_dry_mode_example()
    elif example_type == "checkpoint":
        run_with_checkpointing()
    elif example_type == "full":
        run_full_enhancement()
    else:
        logger.error(f"Unknown example type: {example_type}")
        sys.exit(1)
    
    logger.info("Example completed successfully")