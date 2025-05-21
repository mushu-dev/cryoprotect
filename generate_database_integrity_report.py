#!/usr/bin/env python3
"""
generate_database_integrity_report.py

Generate a visual HTML report from database integrity verification results.
This script takes the JSON report from verify_database_integrity_enhanced.py
and creates a visual, interactive HTML report for easier analysis.

Usage:
    python generate_database_integrity_report.py [--input=report.json] [--output=report.html]
"""

import json
import argparse
import os
from datetime import datetime
import base64
import sys
from typing import Dict, List, Any, Optional

# HTML template for the report
HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CryoProtect Database Integrity Report</title>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            margin: 0;
            padding: 20px;
            background-color: #f5f7fa;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 0 20px rgba(0,0,0,0.1);
        }
        h1, h2, h3 {
            color: #2c3e50;
        }
        h1 {
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }
        .summary-card {
            background-color: #ecf0f1;
            border-radius: 8px;
            padding: 15px;
            margin-bottom: 20px;
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
        }
        .stat-box {
            flex: 1;
            min-width: 200px;
            background-color: white;
            border-radius: 6px;
            padding: 12px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            text-align: center;
        }
        .stat-value {
            font-size: 24px;
            font-weight: bold;
            color: #2c3e50;
        }
        .stat-label {
            font-size: 14px;
            color: #7f8c8d;
        }
        .status {
            font-size: 18px;
            font-weight: bold;
            padding: 8px 16px;
            border-radius: 4px;
            display: inline-block;
            margin-bottom: 20px;
        }
        .status-pass {
            background-color: #e6ffed;
            color: #22863a;
            border: 1px solid #b7eb8f;
        }
        .status-fail {
            background-color: #ffeef0;
            color: #cb2431;
            border: 1px solid #f5c6cb;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 20px;
        }
        th, td {
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        th {
            background-color: #f8f9fa;
            font-weight: 600;
        }
        tr:hover {
            background-color: #f8f9fa;
        }
        .issue-card {
            background-color: white;
            border-radius: 8px;
            padding: 15px;
            margin-bottom: 15px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.05);
            border-left: 4px solid #3498db;
        }
        .issue-error {
            border-left-color: #e74c3c;
        }
        .issue-warning {
            border-left-color: #f39c12;
        }
        .issue-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 10px;
        }
        .issue-title {
            margin: 0;
            font-size: 16px;
            font-weight: 600;
        }
        .issue-severity {
            font-size: 12px;
            padding: 4px 8px;
            border-radius: 4px;
            text-transform: uppercase;
        }
        .severity-error {
            background-color: #ffeef0;
            color: #cb2431;
        }
        .severity-warning {
            background-color: #fff5e6;
            color: #b36b00;
        }
        .issue-description {
            font-size: 14px;
            color: #555;
            margin-bottom: 0;
        }
        
        .chart-container {
            width: 100%;
            max-width: 800px;
            margin: 0 auto 30px auto;
            height: 300px;
        }
        
        @media (max-width: 768px) {
            .stat-box {
                min-width: 100%;
            }
        }
        
        .tabs {
            overflow: hidden;
            margin-bottom: 20px;
        }
        .tab-button {
            background-color: #f1f1f1;
            border: none;
            outline: none;
            cursor: pointer;
            padding: 10px 16px;
            transition: 0.3s;
            font-size: 14px;
            border-radius: 4px 4px 0 0;
        }
        .tab-button:hover {
            background-color: #ddd;
        }
        .tab-button.active {
            background-color: #3498db;
            color: white;
        }
        .tab-content {
            display: none;
            padding: 10px;
            border-top: none;
        }
        
        .search-container {
            margin-bottom: 20px;
        }
        #issueSearch {
            width: 100%;
            padding: 10px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }
        
        /* Chart legend */
        .legend {
            display: flex;
            justify-content: center;
            margin-top: 10px;
            flex-wrap: wrap;
        }
        .legend-item {
            display: flex;
            align-items: center;
            margin-right: 20px;
            margin-bottom: 5px;
        }
        .legend-color {
            width: 15px;
            height: 15px;
            margin-right: 5px;
            border-radius: 3px;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>CryoProtect Database Integrity Report</h1>
        <p>Generated on: {timestamp}</p>
        
        <div class="status {status_class}">Status: {status}</div>
        
        <div class="summary-card">
            <div class="stat-box">
                <div class="stat-value">{tables_checked}</div>
                <div class="stat-label">Tables Checked</div>
            </div>
            <div class="stat-box">
                <div class="stat-value">{rows_checked}</div>
                <div class="stat-label">Rows Checked</div>
            </div>
            <div class="stat-box">
                <div class="stat-value">{foreign_keys_checked}</div>
                <div class="stat-label">FK Constraints Verified</div>
            </div>
            <div class="stat-box">
                <div class="stat-value">{issue_count}</div>
                <div class="stat-label">Issues Found</div>
            </div>
            <div class="stat-box">
                <div class="stat-value">{execution_time}s</div>
                <div class="stat-label">Execution Time</div>
            </div>
        </div>
        
        <!-- Charts Section -->
        <h2>Data Overview</h2>
        <div class="chart-container">
            <canvas id="tableRowsChart"></canvas>
        </div>
        
        <div class="chart-container">
            <canvas id="issuesChart"></canvas>
            <div class="legend">
                <div class="legend-item">
                    <div class="legend-color" style="background-color: #e74c3c;"></div>
                    <span>Errors</span>
                </div>
                <div class="legend-item">
                    <div class="legend-color" style="background-color: #f39c12;"></div>
                    <span>Warnings</span>
                </div>
            </div>
        </div>
        
        <!-- Data Tables Section -->
        <h2>Database Tables</h2>
        
        <div class="tabs">
            <button class="tab-button active" onclick="openTab(event, 'tablesData')">Tables Data</button>
            <button class="tab-button" onclick="openTab(event, 'issuesData')">Issues</button>
        </div>
        
        <div id="tablesData" class="tab-content" style="display: block;">
            <table>
                <thead>
                    <tr>
                        <th>Table Name</th>
                        <th>Row Count</th>
                    </tr>
                </thead>
                <tbody>
                    {tables_rows}
                </tbody>
            </table>
        </div>
        
        <div id="issuesData" class="tab-content">
            <div class="search-container">
                <input type="text" id="issueSearch" placeholder="Search issues..." onkeyup="filterIssues()">
            </div>
            
            <div id="issues-container">
                {issues_html}
            </div>
        </div>
    </div>
    
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <script>
        // Table row count chart
        const tableRowsCtx = document.getElementById('tableRowsChart');
        new Chart(tableRowsCtx, {
            type: 'bar',
            data: {
                labels: {table_names_json},
                datasets: [{
                    label: 'Row Count',
                    data: {table_counts_json},
                    backgroundColor: 'rgba(52, 152, 219, 0.7)',
                    borderColor: 'rgba(52, 152, 219, 1)',
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    title: {
                        display: true,
                        text: 'Table Row Counts',
                        font: {
                            size: 16
                        }
                    },
                    legend: {
                        display: false
                    }
                },
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Number of Rows'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Table Name'
                        }
                    }
                }
            }
        });
        
        // Issues chart
        const issuesCtx = document.getElementById('issuesChart');
        new Chart(issuesCtx, {
            type: 'pie',
            data: {
                labels: ['Errors', 'Warnings'],
                datasets: [{
                    data: [{error_count}, {warning_count}],
                    backgroundColor: [
                        '#e74c3c',
                        '#f39c12'
                    ],
                    borderColor: [
                        '#c0392b',
                        '#d35400'
                    ],
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    title: {
                        display: true,
                        text: 'Issues by Severity',
                        font: {
                            size: 16
                        }
                    },
                    legend: {
                        display: false
                    }
                }
            }
        });
        
        // Tabs functionality
        function openTab(evt, tabName) {
            var i, tabcontent, tablinks;
            
            tabcontent = document.getElementsByClassName("tab-content");
            for (i = 0; i < tabcontent.length; i++) {
                tabcontent[i].style.display = "none";
            }
            
            tablinks = document.getElementsByClassName("tab-button");
            for (i = 0; i < tablinks.length; i++) {
                tablinks[i].className = tablinks[i].className.replace(" active", "");
            }
            
            document.getElementById(tabName).style.display = "block";
            evt.currentTarget.className += " active";
        }
        
        // Issue filtering
        function filterIssues() {
            const search = document.getElementById('issueSearch').value.toLowerCase();
            const issues = document.querySelectorAll('.issue-card');
            
            issues.forEach(issue => {
                const title = issue.querySelector('.issue-title').textContent.toLowerCase();
                const description = issue.querySelector('.issue-description').textContent.toLowerCase();
                
                if (title.includes(search) || description.includes(search)) {
                    issue.style.display = '';
                } else {
                    issue.style.display = 'none';
                }
            });
        }
    </script>
</body>
</html>
"""

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate visual database integrity report")
    parser.add_argument("--input", default="database_integrity_report.json", 
                        help="Input JSON report file")
    parser.add_argument("--output", default="database_integrity_report.html", 
                        help="Output HTML report file")
    return parser.parse_args()

def format_number(num: int) -> str:
    """Format large numbers with commas."""
    return f"{num:,}"

def generate_report(input_file: str, output_file: str) -> bool:
    """Generate HTML report from JSON data."""
    try:
        # Load the JSON report
        with open(input_file, 'r') as f:
            report = json.load(f)
        
        # Format timestamp
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if "timestamp" in report:
            try:
                timestamp = datetime.fromisoformat(report["timestamp"]).strftime("%Y-%m-%d %H:%M:%S")
            except ValueError:
                pass
        
        # Status class
        status = report.get("status", "UNKNOWN")
        status_class = "status-pass" if status == "PASS" else "status-fail"
        
        # Summary data
        summary = report.get("summary", {})
        tables_checked = format_number(summary.get("tables_checked", 0))
        rows_checked = format_number(summary.get("rows_checked", 0))
        foreign_keys_checked = format_number(summary.get("foreign_keys_checked", 0))
        issue_count = format_number(summary.get("issues_found", 0))
        execution_time = report.get("execution_time_seconds", 0)
        
        # Count issues by severity
        issues = report.get("issues", [])
        error_count = sum(1 for issue in issues if issue.get("severity") == "error")
        warning_count = sum(1 for issue in issues if issue.get("severity") == "warning")
        
        # Tables data
        data_counts = report.get("data_counts", {})
        tables_rows = ""
        
        # Sort tables by row count (descending)
        sorted_tables = sorted(data_counts.items(), key=lambda x: x[1], reverse=True)
        
        for table, count in sorted_tables:
            tables_rows += f"""
                <tr>
                    <td>{table}</td>
                    <td>{format_number(count)}</td>
                </tr>
            """
        
        # Issues HTML
        issues_html = ""
        for issue in issues:
            severity = issue.get("severity", "warning")
            severity_class = f"severity-{severity}"
            issue_class = f"issue-{severity}"
            
            issues_html += f"""
                <div class="issue-card {issue_class}">
                    <div class="issue-header">
                        <h3 class="issue-title">{issue.get("title", "Unknown Issue")}</h3>
                        <span class="issue-severity {severity_class}">{severity}</span>
                    </div>
                    <p class="issue-description">{issue.get("description", "No description provided.")}</p>
                </div>
            """
        
        # Prepare table data for charts
        table_names = [table for table, _ in sorted_tables[:20]]  # Limit to top 20 for readability
        table_counts = [count for _, count in sorted_tables[:20]]
        
        # Prepare for JSON embedding in JavaScript
        table_names_json = json.dumps(table_names)
        table_counts_json = json.dumps(table_counts)
        
        # Replace placeholders in template
        html_content = HTML_TEMPLATE.format(
            timestamp=timestamp,
            status=status,
            status_class=status_class,
            tables_checked=tables_checked,
            rows_checked=rows_checked,
            foreign_keys_checked=foreign_keys_checked,
            issue_count=issue_count,
            execution_time=execution_time,
            tables_rows=tables_rows,
            issues_html=issues_html,
            table_names_json=table_names_json,
            table_counts_json=table_counts_json,
            error_count=error_count,
            warning_count=warning_count
        )
        
        # Write the HTML report
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        print(f"Report generated successfully: {output_file}")
        return True
        
    except Exception as e:
        print(f"Error generating report: {e}", file=sys.stderr)
        return False

def main():
    """Main entry point."""
    args = parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found")
        return 1
    
    success = generate_report(args.input, args.output)
    return 0 if success else 1

if __name__ == "__main__":
    exit(main())