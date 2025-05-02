#!/usr/bin/env python3
"""
Wrapper script to run add_performance_indexes.py with correct database credentials.
"""

import os
import sys
import subprocess
import json
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    """Set environment variables and run add_performance_indexes.py"""
    # Map Supabase environment variables to the ones expected by get_db_config
    env_mapping = {
        'DB_HOST': os.environ.get('SUPABASE_DB_HOST', 'db.tsdlmynydfuypiugmkev.supabase.co'),
        'DB_PORT': os.environ.get('SUPABASE_DB_PORT', '5432'),
        'DB_NAME': os.environ.get('SUPABASE_DB_NAME', 'postgres'),
        'DB_USER': os.environ.get('SUPABASE_DB_USER', 'postgres'),
        'DB_PASSWORD': os.environ.get('SUPABASE_DB_PASSWORD', 'postgres'),
    }
    
    # Set environment variables
    for key, value in env_mapping.items():
        os.environ[key] = value
        logger.info(f"Set environment variable {key}={value}")
    
    # Run add_performance_indexes.py
    logger.info("Running add_performance_indexes.py...")
    result = subprocess.run([sys.executable, 'add_performance_indexes.py'], 
                           capture_output=True, text=True)
    
    # Log the output
    if result.stdout:
        logger.info(f"Output: {result.stdout}")
    if result.stderr:
        logger.error(f"Error: {result.stderr}")
    
    # Create a report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/performance_indexes_{timestamp}.json"
    
    # Ensure reports directory exists
    os.makedirs("reports", exist_ok=True)
    
    report = {
        "timestamp": timestamp,
        "success": result.returncode == 0,
        "exit_code": result.returncode,
        "output": result.stdout,
        "error": result.stderr
    }
    
    # Write report to file
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Report saved to {report_path}")
    
    # Update project_state.json
    try:
        with open('project_state.json', 'r') as f:
            project_state = json.load(f)
        
        # Update task status
        if result.returncode == 0:
            project_state['tasks']['task-dbpop-05']['status'] = 'Done'
            project_state['tasks']['task-dbpop-05']['outputs'].append(report_path)
            project_state['tasks']['task-dbpop-05']['log'].append(
                f"{datetime.now().isoformat()}Z - Performance optimization completed successfully. Added indexes to database tables."
            )
            # Update overall phase status
            for phase in project_state['highLevelPlan']:
                if phase.get('taskId') == 'task-dbpop-05':
                    phase['status'] = 'Done'
        else:
            project_state['tasks']['task-dbpop-05']['status'] = 'Error'
            project_state['tasks']['task-dbpop-05']['log'].append(
                f"{datetime.now().isoformat()}Z - Performance optimization failed: {result.stderr}"
            )
        
        # Write updated project state
        with open('project_state.json', 'w') as f:
            json.dump(project_state, f, indent=2)
        
        logger.info("Updated project_state.json")
    except Exception as e:
        logger.error(f"Failed to update project_state.json: {str(e)}")
    
    return result.returncode

if __name__ == "__main__":
    sys.exit(main())