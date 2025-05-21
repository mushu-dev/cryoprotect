#!/bin/bash
#
# Production Testing Script for CryoProtect
# This script runs the production tests inside the test container

set -e  # Exit on any error

echo "Running CryoProtect production tests..."

# Make sure the containers are running
if ! podman ps | grep -q cryoprotect-test; then
    echo "Error: Test container is not running."
    echo "Please run ./production_test_container.sh first."
    exit 1
fi

if ! podman ps | grep -q cryoprotect-app; then
    echo "Error: App container is not running."
    echo "Please run ./production_test_container.sh first."
    exit 1
fi

if ! podman ps | grep -q cryoprotect-rdkit; then
    echo "Error: RDKit container is not running."
    echo "Please run ./production_test_container.sh first."
    exit 1
fi

# Check if the containers can communicate
echo "Verifying container communication..."
if ! podman exec cryoprotect-test curl -s http://cryoprotect-app:5000/health &>/dev/null; then
    echo "Error: Cannot connect to app container. Check network configuration."
    exit 1
fi

if ! podman exec cryoprotect-test curl -s http://cryoprotect-rdkit:5000/health &>/dev/null; then
    echo "Error: Cannot connect to RDKit container. Check network configuration."
    exit 1
fi

# Ensure test data directory exists
echo "Ensuring test data directory exists..."
podman exec cryoprotect-test mkdir -p /app/tests/test_data

# Copy test data files if needed
echo "Checking for test data files..."
if ! podman exec cryoprotect-test ls /app/tests/test_data/core_cryoprotectants.json; then
    echo "Creating sample test data..."
    podman exec cryoprotect-test sh -c "cd /app && python -c \"
import json
import os

# Create test data directory
os.makedirs('/app/tests/test_data', exist_ok=True)

# Create core cryoprotectants data
core_data = {
    'molecules': [
        {
            'name': 'Dimethyl sulfoxide',
            'smiles': 'CS(=O)C',
            'formula': 'C2H6OS',
            'properties': {
                'molecular_weight': 78.13,
                'logp': -1.35,
                'tpsa': 17.07,
                'h_bond_donors': 0,
                'h_bond_acceptors': 1
            }
        },
        {
            'name': 'Glycerol',
            'smiles': 'C(C(CO)O)O',
            'formula': 'C3H8O3',
            'properties': {
                'molecular_weight': 92.09,
                'logp': -1.76,
                'tpsa': 60.69,
                'h_bond_donors': 3,
                'h_bond_acceptors': 3
            }
        }
    ]
}

# Create edge cases data
edge_data = {
    'molecules': [
        {
            'name': 'Sucrose',
            'smiles': 'C(C1C(C(C(C(O1)O)O)O)O)OC2C(C(C(C(O2)CO)O)O)O',
            'formula': 'C12H22O11',
            'properties': {
                'molecular_weight': 342.3,
                'logp': -3.76,
                'tpsa': 189.53,
                'h_bond_donors': 8,
                'h_bond_acceptors': 11
            }
        }
    ]
}

# Create mixtures data
mixtures_data = {
    'mixtures': [
        {
            'name': 'DMSO-Glycerol Mix',
            'description': 'Common cryoprotectant mixture',
            'components': [
                {
                    'molecule_name': 'Dimethyl sulfoxide',
                    'concentration': 70,
                    'concentration_unit': '%v/v'
                },
                {
                    'molecule_name': 'Glycerol',
                    'concentration': 30,
                    'concentration_unit': '%v/v'
                }
            ]
        }
    ]
}

# Write to files
with open('/app/tests/test_data/core_cryoprotectants.json', 'w') as f:
    json.dump(core_data, f, indent=2)

with open('/app/tests/test_data/edge_cases.json', 'w') as f:
    json.dump(edge_data, f, indent=2)

with open('/app/tests/test_data/mixtures.json', 'w') as f:
    json.dump(mixtures_data, f, indent=2)

print('Created sample test data files.')
\""
fi

# Load test data
echo "Loading test data into the database..."
podman exec cryoprotect-test python /app/tests/test_data/load_test_data.py --dataset all || echo "Warning: Test data loading failed, but continuing anyway."

# Run the production workflow tests
echo "Running production workflow tests..."
podman exec cryoprotect-test python /app/production_workflow_test.py --app http://cryoprotect-app:5000 --rdkit http://cryoprotect-rdkit:5000 --output /app/production_workflow_results.json

# Copy the test results
echo "Copying test results..."
podman cp cryoprotect-test:/app/production_workflow_results.json ./production_workflow_results.json
podman cp cryoprotect-test:/app/production_workflow_test.log ./production_workflow_test.log

# Copy the test results
podman cp cryoprotect-test:/app/real_data_test_results.json ./
podman cp cryoprotect-test:/app/real_data_test.log ./

echo "Tests completed. Results saved to real_data_test_results.json and real_data_test.log"

# Generate HTML report
echo "Generating HTML report..."
cat > ./generate_report.py << EOL
#!/usr/bin/env python3
import json
import datetime

# Load test results
with open('real_data_test_results.json', 'r') as f:
    results = json.load(f)

# Generate HTML report
html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>CryoProtect Production Test Results</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1, h2, h3 {{ color: #333; }}
        .summary {{ margin: 20px 0; padding: 10px; background-color: #f5f5f5; border-radius: 5px; }}
        .passed {{ color: green; }}
        .failed {{ color: red; }}
        .skipped {{ color: orange; }}
        .test-case {{ margin: 10px 0; padding: 10px; border: 1px solid #ddd; border-radius: 5px; }}
        .test-case.passed {{ border-left: 5px solid green; }}
        .test-case.failed {{ border-left: 5px solid red; }}
        .test-case.skipped {{ border-left: 5px solid orange; }}
        .performance {{ margin: 20px 0; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        th {{ background-color: #4CAF50; color: white; }}
    </style>
</head>
<body>
    <h1>CryoProtect Production Test Results</h1>
    <p>Report generated on {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    
    <div class="summary">
        <h2>Test Summary</h2>
        <p>Status: <span class="{results['status'].lower()}">{results['status']}</span></p>
        <p>Total Tests: {results['total_tests']}</p>
        <p>Passed: <span class="passed">{results['passed_tests']}</span></p>
        <p>Failed: <span class="failed">{results['failed_tests']}</span></p>
        <p>Skipped: <span class="skipped">{results['skipped_tests']}</span></p>
        <p>Start Time: {results.get('start_time', 'N/A')}</p>
        <p>End Time: {results.get('end_time', 'N/A')}</p>
    </div>
    
    <h2>Test Results</h2>
"""

# Group test cases by workflow
workflows = {}
for tc in results['test_cases']:
    workflow_id = tc['id'].split('.')[0] if '.' in tc['id'] else tc['id']
    if workflow_id not in workflows:
        workflows[workflow_id] = []
    workflows[workflow_id].append(tc)

# Sort workflows by ID
sorted_workflows = sorted(workflows.items())

# Add each workflow and its test cases
for workflow_id, test_cases in sorted_workflows:
    # Find the workflow summary test case
    workflow_summary = None
    for tc in test_cases:
        if tc['id'] == workflow_id:
            workflow_summary = tc
            break
    
    if workflow_summary:
        html += f"""
        <h3>{workflow_summary['name']}</h3>
        <div class="test-case {workflow_summary['status'].lower()}">
            <p>Status: <span class="{workflow_summary['status'].lower()}">{workflow_summary['status']}</span></p>
            <p>Message: {workflow_summary['message']}</p>
        </div>
        """
    else:
        html += f"<h3>Workflow {workflow_id}</h3>"
    
    # Add the individual test cases
    for tc in test_cases:
        if tc['id'] == workflow_id:
            continue  # Skip the summary test case
        
        html += f"""
        <div class="test-case {tc['status'].lower()}">
            <h4>{tc['name']} ({tc['id']})</h4>
            <p>Status: <span class="{tc['status'].lower()}">{tc['status']}</span></p>
            <p>Message: {tc['message']}</p>
        </div>
        """

# Add performance metrics
if 'performance_metrics' in results:
    html += """
    <h2>Performance Metrics</h2>
    <div class="performance">
        <table>
            <tr>
                <th>Metric</th>
                <th>Value</th>
            </tr>
    """
    
    for name, value in results['performance_metrics'].items():
        if isinstance(value, list):
            # Show average for lists
            avg = sum(value) / len(value) if value else 0
            html += f"""
            <tr>
                <td>{name}</td>
                <td>{avg:.4f}s (avg of {len(value)} measurements)</td>
            </tr>
            """
        elif isinstance(value, float):
            html += f"""
            <tr>
                <td>{name}</td>
                <td>{value:.4f}s</td>
            </tr>
            """
        else:
            html += f"""
            <tr>
                <td>{name}</td>
                <td>{value}</td>
            </tr>
            """
    
    html += """
        </table>
    </div>
    """

# Add performance summary if available
if 'performance_summary' in results:
    html += """
    <h2>Performance Summary</h2>
    <div class="performance">
        <table>
            <tr>
                <th>Metric</th>
                <th>Value</th>
            </tr>
    """
    
    for name, value in results['performance_summary'].items():
        if isinstance(value, dict):
            # Show statistics
            html += f"""
            <tr>
                <td>{name}</td>
                <td>
                    <ul>
            """
            for stat_name, stat_value in value.items():
                html += f"<li>{stat_name}: {stat_value:.4f}s</li>" if isinstance(stat_value, float) else f"<li>{stat_name}: {stat_value}</li>"
            
            html += """
                    </ul>
                </td>
            </tr>
            """
        elif isinstance(value, float):
            html += f"""
            <tr>
                <td>{name}</td>
                <td>{value:.4f}s</td>
            </tr>
            """
        else:
            html += f"""
            <tr>
                <td>{name}</td>
                <td>{value}</td>
            </tr>
            """
    
    html += """
        </table>
    </div>
    """

html += """
</body>
</html>
"""

# Write HTML report
with open('production_test_report.html', 'w') as f:
    f.write(html)

print("HTML report generated: production_test_report.html")
EOL

# Make the script executable
chmod +x ./generate_report.py

# Generate the HTML report
python ./generate_report.py

echo "HTML report generated: production_test_report.html"
echo "Production testing completed."