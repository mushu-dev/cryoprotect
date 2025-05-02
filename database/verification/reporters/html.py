"""
HTML report generator for verification results.

This module provides functions for generating HTML reports.
"""

import logging
from datetime import datetime
from typing import Dict

logger = logging.getLogger(__name__)

def generate_html_report(results: Dict) -> str:
    """
    Generate HTML report from verification results.

    Args:
        results: Verification results dictionary

    Returns:
        HTML report as a string
    """
    timestamp = results.get('timestamp', datetime.now().isoformat())
    level = results.get('level', 'unknown')
    success = results.get('success', False)

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Database Verification Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            color: #333;
        }}
        h1, h2, h3 {{
            color: #444;
        }}
        .header {{
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 1px solid #eee;
        }}
        .timestamp {{
            color: #888;
            font-size: 0.9em;
        }}
        .success {{
            color: green;
        }}
        .failure {{
            color: red;
        }}
        .module {{
            margin-bottom: 20px;
            padding: 15px;
            background-color: #f9f9f9;
            border-radius: 5px;
        }}
        .issue {{
            margin: 10px 0;
            padding: 10px;
            background-color: #fff;
            border-radius: 3px;
            border-left: 4px solid #ddd;
        }}
        .issue.error {{
            border-left-color: #e74c3c;
        }}
        .issue.warning {{
            border-left-color: #f39c12;
        }}
        .summary {{
            margin-top: 20px;
            padding: 15px;
            background-color: #f0f0f0;
            border-radius: 5px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Database Verification Report</h1>
        <div class="timestamp">Generated on {timestamp}</div>
    </div>

    <h2>Overview</h2>
    <p>Verification level: <strong>{level}</strong></p>
    <p>Overall status: <span class="{('success' if success else 'failure')}">
        {('PASSED' if success else 'FAILED')}
    </span></p>

    <h2>Module Results</h2>
"""

    # Add module results
    for module, module_results in results.get('results', {}).items():
        module_success = module_results.get('success', False)
        module_status = 'PASSED' if module_success else 'FAILED'
        status_class = 'success' if module_success else 'failure'

        html += f"""
    <div class="module">
        <h3>{module} <span class="{status_class}">{module_status}</span></h3>
"""

        # Add issues
        issues = module_results.get('issues', [])
        if issues:
            html += f"<p>Found {len(issues)} issues:</p>"

            for issue in issues:
                severity = issue.get('severity', 'info')
                message = issue.get('message', 'Unknown issue')

                html += f"""
        <div class="issue {severity}">
            <strong>{severity.upper()}:</strong> {message}
        </div>
"""
        else:
            html += "<p>No issues found.</p>"

        html += "</div>"

    # Add summary
    total_issues = sum(
        len(module_results.get('issues', []))
        for module_results in results.get('results', {}).values()
    )

    error_count = sum(
        sum(1 for issue in module_results.get('issues', []) if issue.get('severity') == 'error')
        for module_results in results.get('results', {}).values()
    )

    warning_count = sum(
        sum(1 for issue in module_results.get('issues', []) if issue.get('severity') == 'warning')
        for module_results in results.get('results', {}).values()
    )

    html += f"""
    <div class="summary">
        <h2>Summary</h2>
        <p>Total issues: {total_issues}</p>
        <p>Errors: {error_count}</p>
        <p>Warnings: {warning_count}</p>
    </div>
</body>
</html>
"""

    return html