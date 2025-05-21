#!/usr/bin/env python3
"""
CryoProtect v2 - Master Database Population Script

This script orchestrates the entire database population process using the connection factory.
It sequentially executes the following steps:
1. Import reference compounds
2. Import ChEMBL data
3. Reconcile ChEMBL and PubChem cross-references
4. Enhance PubChem properties
5. Add performance indexes

Features:
- Command-line arguments for selective execution
- Comprehensive logging and reporting
- Checkpoint management for resumable operations
- Progress tracking and status updates
- Error handling with detailed reporting
"""

import os
import sys
import argparse
import logging
import json
import time
import importlib
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Set

# Import database connection utilities
from database.connection import get_db_connection, get_db_connection_info
import sql_executor

# Configure logging
os.makedirs("logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[
        logging.FileHandler(f"logs/database_population_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("database_population")

# Define the population steps
POPULATION_STEPS = [
    {
        "id": "reference",
        "name": "Reference Compounds Import",
        "module": "import_reference_compounds",
        "function": "import_reference_compounds",
        "description": "Import reference cryoprotectant compounds with complete property data",
        "enabled": True,
        "required": True,
        "dependencies": []
    },
    {
        "id": "chembl",
        "name": "ChEMBL Data Import",
        "module": "import_full_chembl",
        "function": "main",
        "description": "Import cryoprotectant-related compounds from ChEMBL",
        "enabled": True,
        "required": True,
        "dependencies": ["reference"]
    },
    {
        "id": "reconcile",
        "name": "Cross-Reference Reconciliation",
        "module": "reconcile_chembl_properties",
        "function": "reconcile_cross_references",
        "description": "Reconcile cross-references between ChEMBL and PubChem",
        "enabled": True,
        "required": True,
        "dependencies": ["reference", "chembl"]
    },
    {
        "id": "pubchem",
        "name": "PubChem Property Enhancement",
        "module": "enhance_pubchem_properties",
        "function": "enhance_pubchem_properties",
        "description": "Enhance molecules with additional properties from PubChem",
        "enabled": True,
        "required": True,
        "dependencies": ["reference", "chembl", "reconcile"]
    },
    {
        "id": "performance",
        "name": "Database Performance Optimization",
        "module": "add_performance_indexes",
        "function": "apply_performance_indexes",
        "description": "Add performance indexes to optimize database queries",
        "enabled": True,
        "required": True,
        "dependencies": ["reference", "chembl", "reconcile", "pubchem"]
    }
]

class PopulationManager:
    """
    Manages the database population process.
    Handles step execution, checkpointing, and reporting.
    """
    
    def __init__(self, args):
        """
        Initialize the population manager.
        
        Args:
            args: Command-line arguments
        """
        self.args = args
        self.checkpoint_dir = args.checkpoint_dir
        self.report_dir = args.report_dir
        self.steps = POPULATION_STEPS
        self.results = {
            "start_time": datetime.now().isoformat(),
            "end_time": None,
            "overall_status": "pending",
            "steps": {}
        }
        
        # Create directories
        os.makedirs(self.checkpoint_dir, exist_ok=True)
        os.makedirs(self.report_dir, exist_ok=True)
        
        # Load checkpoint if resuming
        self.checkpoint_file = os.path.join(self.checkpoint_dir, "population_checkpoint.json")
        if args.resume and os.path.exists(self.checkpoint_file):
            self._load_checkpoint()
        
        # Apply step filters from command line
        self._apply_step_filters()
        
        # Validate dependencies
        self._validate_dependencies()
    
    def _load_checkpoint(self):
        """Load checkpoint data from file."""
        try:
            with open(self.checkpoint_file, 'r') as f:
                checkpoint_data = json.load(f)
                
            # Restore results from checkpoint
            if "results" in checkpoint_data:
                self.results = checkpoint_data["results"]
                logger.info(f"Loaded checkpoint with {len(self.results['steps'])} completed steps")
                
                # Update start time if resuming
                if not self.args.restart:
                    self.results["start_time"] = datetime.now().isoformat()
        except Exception as e:
            logger.warning(f"Failed to load checkpoint file: {str(e)}")
    
    def _save_checkpoint(self):
        """Save checkpoint data to file."""
        try:
            checkpoint_data = {
                "timestamp": datetime.now().isoformat(),
                "results": self.results
            }
            
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
                
            logger.debug(f"Saved checkpoint with {len(self.results['steps'])} steps")
        except Exception as e:
            logger.error(f"Failed to save checkpoint: {str(e)}")
    
    def _apply_step_filters(self):
        """Apply step filters from command line arguments."""
        # If specific steps are specified, disable all others
        if self.args.steps:
            step_ids = self.args.steps.split(',')
            for step in self.steps:
                step["enabled"] = step["id"] in step_ids
        
        # If steps to skip are specified, disable them
        if self.args.skip:
            skip_ids = self.args.skip.split(',')
            for step in self.steps:
                if step["id"] in skip_ids:
                    step["enabled"] = False
        
        # If restart flag is set, clear previous results
        if self.args.restart:
            self.results = {
                "start_time": datetime.now().isoformat(),
                "end_time": None,
                "overall_status": "pending",
                "steps": {}
            }
    
    def _validate_dependencies(self):
        """Validate that all dependencies are satisfied."""
        enabled_steps = {step["id"] for step in self.steps if step["enabled"]}
        
        for step in self.steps:
            if step["enabled"]:
                for dependency in step["dependencies"]:
                    if dependency not in enabled_steps:
                        if self.args.force:
                            logger.warning(f"Step '{step['id']}' depends on disabled step '{dependency}', but --force is set")
                        else:
                            logger.error(f"Step '{step['id']}' depends on disabled step '{dependency}'")
                            step["enabled"] = False
                            enabled_steps.remove(step["id"])
                            break
    
    def _get_step_status(self, step_id):
        """Get the status of a step from results."""
        if step_id in self.results["steps"]:
            return self.results["steps"][step_id]["status"]
        return "pending"
    
    def _should_run_step(self, step):
        """Determine if a step should be run."""
        # Skip disabled steps
        if not step["enabled"]:
            return False
        
        # Skip completed steps unless restart is specified
        if self.args.restart:
            return True
        
        status = self._get_step_status(step["id"])
        if status == "success" and not self.args.force:
            return False
        
        # Check dependencies
        for dependency in step["dependencies"]:
            dep_status = self._get_step_status(dependency)
            if dep_status != "success":
                if self.args.force:
                    logger.warning(f"Dependency '{dependency}' has status '{dep_status}', but --force is set")
                else:
                    logger.error(f"Cannot run step '{step['id']}' because dependency '{dependency}' has status '{dep_status}'")
                    return False
        
        return True
    
    def _import_module_function(self, module_name, function_name):
        """
        Import a function from a module.
        
        Args:
            module_name: Name of the module
            function_name: Name of the function
            
        Returns:
            The imported function
        """
        try:
            module = importlib.import_module(module_name)
            return getattr(module, function_name)
        except ImportError:
            logger.error(f"Failed to import module '{module_name}'")
            raise
        except AttributeError:
            logger.error(f"Function '{function_name}' not found in module '{module_name}'")
            raise
    
    def _execute_step(self, step):
        """
        Execute a single population step.
        
        Args:
            step: Step configuration dictionary
            
        Returns:
            Tuple of (success, result)
        """
        step_id = step["id"]
        step_name = step["name"]
        
        logger.info(f"Starting step: {step_name}")
        
        # Initialize step result
        step_result = {
            "id": step_id,
            "name": step_name,
            "start_time": datetime.now().isoformat(),
            "end_time": None,
            "status": "running",
            "error": None,
            "details": {}
        }
        
        # Save to results and checkpoint
        self.results["steps"][step_id] = step_result
        self._save_checkpoint()
        
        try:
            # Import the module function
            func = self._import_module_function(step["module"], step["function"])
            
            # Prepare arguments based on step
            kwargs = {}
            
            # Common arguments for all steps
            if step_id == "reference":
                # Reference compounds import
                report_path = os.path.join(self.report_dir, f"reference_import_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
                checkpoint_path = os.path.join(self.checkpoint_dir, "reference_import_checkpoint.json")
                kwargs = {
                    "output_report": report_path,
                    "checkpoint_path": checkpoint_path,
                    "resume": self.args.resume,
                    "batch_size": self.args.batch_size
                }
            
            elif step_id == "chembl":
                # ChEMBL import - this uses argparse internally, so we need to modify sys.argv
                # Save original argv
                original_argv = sys.argv.copy()
                
                # Build new argv for the ChEMBL import
                chembl_argv = [
                    "import_full_chembl.py",
                    f"--batch-size={self.args.batch_size}",
                    f"--checkpoint-dir={self.checkpoint_dir}",
                    "--mode=resume" if self.args.resume else "--mode=full"
                ]
                
                if self.args.verbose:
                    chembl_argv.append("--verbose")
                
                # Set the new argv
                sys.argv = chembl_argv
                
                # The main function handles its own arguments
                result = func()
                
                # Restore original argv
                sys.argv = original_argv
                
                # Update step result
                step_result["status"] = "success" if result == 0 else "error"
                step_result["end_time"] = datetime.now().isoformat()
                step_result["details"]["exit_code"] = result
                
                # Save checkpoint
                self._save_checkpoint()
                
                return step_result["status"] == "success", step_result
            
            elif step_id == "reconcile":
                # Reconciliation
                report_path = os.path.join(self.report_dir, f"reconciliation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
                kwargs = {
                    "output_report": report_path,
                    "dry_run": False
                }
            
            elif step_id == "pubchem":
                # PubChem property enhancement
                report_path = os.path.join(self.report_dir, f"pubchem_enhancement_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
                checkpoint_path = os.path.join(self.checkpoint_dir, "pubchem_enhancement.json")
                kwargs = {
                    "output_report": report_path,
                    "batch_size": self.args.batch_size,
                    "dry_run": False,
                    "checkpoint_file": checkpoint_path
                }
            
            elif step_id == "performance":
                # Performance optimization
                kwargs = {
                    "schema": "public",
                    "test_only": False
                }
            
            # Execute the function with appropriate arguments
            if step_id != "chembl":  # ChEMBL is handled separately above
                result = func(**kwargs)
                
                # Update step result
                step_result["status"] = "success"
                step_result["details"]["result"] = result
            
            logger.info(f"Step '{step_name}' completed successfully")
            
        except Exception as e:
            logger.error(f"Error executing step '{step_name}': {str(e)}")
            logger.error(traceback.format_exc())
            
            # Update step result
            step_result["status"] = "error"
            step_result["error"] = str(e)
            step_result["details"]["traceback"] = traceback.format_exc()
            
            # Save checkpoint
            step_result["end_time"] = datetime.now().isoformat()
            self._save_checkpoint()
            
            return False, step_result
        
        # Update step result
        step_result["end_time"] = datetime.now().isoformat()
        self._save_checkpoint()
        
        return True, step_result
    
    def run(self):
        """
        Run the database population process.
        
        Returns:
            Overall success status (True/False)
        """
        logger.info("Starting database population process")
        
        # Verify database connection using connection factory
        try:
            logger.info("Verifying database connection using connection factory")
            db = get_db_connection()
            success, message = db.test_connection()
            if success:
                logger.info(f"Successfully connected to database: {message}")
                
                # Get connection info
                conn_info = get_db_connection_info()
                logger.info(f"Connection info: {conn_info}")
            else:
                logger.error(f"Database connection test failed: {message}")
                return False
        except Exception as e:
            logger.error(f"Database connection verification failed: {str(e)}")
            return False
        
        # Execute each step
        all_success = True
        for step in self.steps:
            if self._should_run_step(step):
                success, _ = self._execute_step(step)
                if not success:
                    all_success = False
                    if step["required"] and not self.args.continue_on_error:
                        logger.error(f"Required step '{step['name']}' failed, stopping process")
                        break
            else:
                logger.info(f"Skipping step: {step['name']}")
        
        # Update overall status
        self.results["end_time"] = datetime.now().isoformat()
        self.results["overall_status"] = "success" if all_success else "error"
        
        # Save final checkpoint
        self._save_checkpoint()
        
        # Generate final report
        self._generate_final_report()
        
        # Log completion
        if all_success:
            logger.info("Database population process completed successfully")
        else:
            logger.warning("Database population process completed with errors")
        
        return all_success
    
    def _generate_final_report(self):
        """Generate a final report of the population process."""
        # Calculate duration
        try:
            start_time = datetime.fromisoformat(self.results["start_time"])
            end_time = datetime.fromisoformat(self.results["end_time"])
            duration_seconds = (end_time - start_time).total_seconds()
            
            # Format duration
            hours, remainder = divmod(duration_seconds, 3600)
            minutes, seconds = divmod(remainder, 60)
            duration_formatted = f"{int(hours)}h {int(minutes)}m {int(seconds)}s"
            
            # Add to results
            self.results["duration_seconds"] = duration_seconds
            self.results["duration_formatted"] = duration_formatted
        except Exception as e:
            logger.error(f"Error calculating duration: {str(e)}")
        
        # Count steps by status
        status_counts = {}
        for step_id, step in self.results["steps"].items():
            status = step.get("status", "unknown")
            if status not in status_counts:
                status_counts[status] = 0
            status_counts[status] += 1
        
        self.results["status_counts"] = status_counts
        
        # Save final report
        report_path = os.path.join(self.report_dir, f"population_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
        try:
            with open(report_path, 'w') as f:
                json.dump(self.results, f, indent=2)
            logger.info(f"Final report saved to {report_path}")
            
            # Also generate a markdown report
            markdown_path = os.path.join(self.report_dir, f"population_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md")
            self._generate_markdown_report(markdown_path)
        except Exception as e:
            logger.error(f"Error saving final report: {str(e)}")
    
    def _generate_markdown_report(self, path):
        """
        Generate a markdown report of the population process.
        
        Args:
            path: Path to save the markdown report
        """
        try:
            with open(path, 'w') as f:
                # Title
                f.write("# CryoProtect v2 Database Population Report\n\n")
                
                # Summary
                f.write("## Summary\n\n")
                f.write(f"- **Status**: {self.results['overall_status'].upper()}\n")
                f.write(f"- **Start Time**: {self.results['start_time']}\n")
                f.write(f"- **End Time**: {self.results['end_time']}\n")
                if "duration_formatted" in self.results:
                    f.write(f"- **Duration**: {self.results['duration_formatted']}\n")
                f.write("\n")
                
                # Status counts
                if "status_counts" in self.results:
                    f.write("## Status Counts\n\n")
                    for status, count in self.results["status_counts"].items():
                        f.write(f"- **{status.capitalize()}**: {count}\n")
                    f.write("\n")
                
                # Step details
                f.write("## Step Details\n\n")
                for step in self.steps:
                    step_id = step["id"]
                    if step_id in self.results["steps"]:
                        step_result = self.results["steps"][step_id]
                        status = step_result.get("status", "unknown")
                        status_emoji = "✅" if status == "success" else "❌" if status == "error" else "⏳"
                        
                        f.write(f"### {status_emoji} {step['name']} ({step_id})\n\n")
                        f.write(f"- **Status**: {status.upper()}\n")
                        f.write(f"- **Start Time**: {step_result.get('start_time', 'N/A')}\n")
                        f.write(f"- **End Time**: {step_result.get('end_time', 'N/A')}\n")
                        
                        if "error" in step_result and step_result["error"]:
                            f.write(f"- **Error**: {step_result['error']}\n")
                        
                        f.write("\n")
                    else:
                        f.write(f"### ⏹️ {step['name']} ({step_id})\n\n")
                        f.write("- **Status**: NOT RUN\n\n")
                
                # Footer
                f.write("---\n")
                f.write(f"Report generated at {datetime.now().isoformat()}\n")
            
            logger.info(f"Markdown report saved to {path}")
        except Exception as e:
            logger.error(f"Error generating markdown report: {str(e)}")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="CryoProtect v2 Database Population")
    
    # Step selection
    parser.add_argument("--steps", help="Comma-separated list of steps to run (reference,chembl,reconcile,pubchem,performance)")
    parser.add_argument("--skip", help="Comma-separated list of steps to skip")
    
    # Execution control
    parser.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    parser.add_argument("--restart", action="store_true", help="Restart from beginning, ignoring checkpoints")
    parser.add_argument("--force", action="store_true", help="Force execution even if dependencies are not satisfied")
    parser.add_argument("--continue-on-error", action="store_true", help="Continue execution even if a required step fails")
    
    # Performance options
    parser.add_argument("--batch-size", type=int, default=10, help="Batch size for processing")
    
    # Output options
    parser.add_argument("--checkpoint-dir", default="./checkpoints", help="Directory for checkpoint files")
    parser.add_argument("--report-dir", default="./reports", help="Directory for report files")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    
    return parser.parse_args()

def main():
    """Main function."""
    args = parse_arguments()
    
    # Set log level based on verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize the population manager
    manager = PopulationManager(args)
    
    # Run the population process
    success = manager.run()
    
    # Return exit code
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
