#!/usr/bin/env python3
"""
Update project state for task-dbpop-05.

This script updates the project_state.json file to mark the performance optimization task as completed.
"""

import json
import os
import sys
from datetime import datetime

def update_project_state():
    """Update project_state.json to mark task-dbpop-05 as completed."""
    try:
        # Read current project state
        with open('project_state.json', 'r') as f:
            project_state = json.load(f)
        
        # Get the latest report file
        report_files = [f for f in os.listdir('reports') if f.startswith('performance_indexes_')]
        if not report_files:
            print("No performance indexes report found.")
            return 1
        
        latest_report = sorted(report_files)[-1]
        report_path = os.path.join('reports', latest_report)
        
        # Read the report
        with open(report_path, 'r') as f:
            report = json.load(f)
        
        # Update task status
        project_state['tasks']['task-dbpop-05']['status'] = 'Done'
        if report_path not in project_state['tasks']['task-dbpop-05']['outputs']:
            project_state['tasks']['task-dbpop-05']['outputs'].append(report_path)
        
        # Add log entry
        timestamp = datetime.now().isoformat() + 'Z'
        log_entry = f"{timestamp} - Performance optimization completed successfully. Added {report['indexes_added']} indexes to {len(report['tables_modified'])} tables."
        project_state['tasks']['task-dbpop-05']['log'].append(log_entry)
        
        # Update phase status
        for phase in project_state['highLevelPlan']:
            if phase.get('taskId') == 'task-dbpop-05':
                phase['status'] = 'Done'
        
        # Write updated project state
        with open('project_state.json', 'w') as f:
            json.dump(project_state, f, indent=2)
        
        print(f"Updated project_state.json - Task task-dbpop-05 marked as Done")
        return 0
    
    except Exception as e:
        print(f"Error updating project state: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(update_project_state())