#!/usr/bin/env python3
"""
CryoProtect v2 - API Monitoring Dashboard

This script provides a real-time monitoring dashboard for the CryoProtect API:
1. Monitors API endpoint health
2. Tracks database connection status
3. Displays performance metrics
4. Logs errors and warnings

Usage:
    python api_monitoring_dashboard.py [--port <port>]
"""

import os
import sys
import json
import time
import uuid
import logging
import argparse
import threading
import webbrowser
from datetime import datetime
from flask import Flask, render_template_string, jsonify
from dotenv import load_dotenv

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("api_monitoring.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Supabase connection
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

# Flask app
app = Flask(__name__)

# Global state
monitoring_data = {
    "api_status": "Unknown",
    "database_status": "Unknown",
    "endpoints": {},
    "performance": {
        "response_times": [],
        "error_rates": [],
        "timestamps": []
    },
    "errors": [],
    "last_updated": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
}

# HTML template for the dashboard
DASHBOARD_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <title>CryoProtect API Monitoring Dashboard</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
        }
        .header {
            background-color: #2c3e50;
            color: white;
            padding: 20px;
            border-radius: 5px;
            margin-bottom: 20px;
        }
        .card {
            background-color: white;
            border-radius: 5px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            padding: 20px;
            margin-bottom: 20px;
        }
        .status {
            display: inline-block;
            padding: 5px 10px;
            border-radius: 3px;
            font-weight: bold;
        }
        .status-ok {
            background-color: #2ecc71;
            color: white;
        }
        .status-warning {
            background-color: #f39c12;
            color: white;
        }
        .status-error {
            background-color: #e74c3c;
            color: white;
        }
        .status-unknown {
            background-color: #95a5a6;
            color: white;
        }
        table {
            width: 100%;
            border-collapse: collapse;
        }
        table, th, td {
            border: 1px solid #ddd;
        }
        th, td {
            padding: 10px;
            text-align: left;
        }
        th {
            background-color: #f2f2f2;
        }
        .error-log {
            max-height: 200px;
            overflow-y: auto;
            background-color: #f9f9f9;
            padding: 10px;
            border: 1px solid #ddd;
            border-radius: 3px;
        }
        .error-item {
            margin-bottom: 5px;
            padding: 5px;
            border-bottom: 1px solid #eee;
        }
        .error-item:last-child {
            border-bottom: none;
        }
        .refresh-button {
            background-color: #3498db;
            color: white;
            border: none;
            padding: 10px 15px;
            border-radius: 3px;
            cursor: pointer;
            font-size: 14px;
        }
        .refresh-button:hover {
            background-color: #2980b9;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>CryoProtect API Monitoring Dashboard</h1>
            <p>Last updated: <span id="last-updated">{{ data.last_updated }}</span></p>
            <button class="refresh-button" onclick="refreshData()">Refresh Data</button>
        </div>
        
        <div class="card">
            <h2>System Status</h2>
            <div>
                <p>API Status: 
                    <span class="status {% if data.api_status == 'OK' %}status-ok{% elif data.api_status == 'Warning' %}status-warning{% elif data.api_status == 'Error' %}status-error{% else %}status-unknown{% endif %}" id="api-status">
                        {{ data.api_status }}
                    </span>
                </p>
                <p>Database Status: 
                    <span class="status {% if data.database_status == 'OK' %}status-ok{% elif data.database_status == 'Warning' %}status-warning{% elif data.database_status == 'Error' %}status-error{% else %}status-unknown{% endif %}" id="db-status">
                        {{ data.database_status }}
                    </span>
                </p>
            </div>
        </div>
        
        <div class="card">
            <h2>API Endpoints</h2>
            <table id="endpoints-table">
                <thead>
                    <tr>
                        <th>Endpoint</th>
                        <th>Method</th>
                        <th>Status</th>
                        <th>Response Time (ms)</th>
                        <th>Last Checked</th>
                    </tr>
                </thead>
                <tbody>
                    {% for endpoint, data in data.endpoints.items() %}
                    <tr>
                        <td>{{ endpoint }}</td>
                        <td>{{ data.method }}</td>
                        <td>
                            <span class="status {% if data.status == 'OK' %}status-ok{% elif data.status == 'Warning' %}status-warning{% elif data.status == 'Error' %}status-error{% else %}status-unknown{% endif %}">
                                {{ data.status }}
                            </span>
                        </td>
                        <td>{{ data.response_time }}</td>
                        <td>{{ data.last_checked }}</td>
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
        </div>
        
        <div class="card">
            <h2>Performance Metrics</h2>
            <div id="performance-metrics">
                <p>Average Response Time: <span id="avg-response-time">{{ data.performance.avg_response_time if data.performance.avg_response_time else 'N/A' }}</span> ms</p>
                <p>Error Rate: <span id="error-rate">{{ data.performance.error_rate if data.performance.error_rate else 'N/A' }}</span>%</p>
            </div>
        </div>
        
        <div class="card">
            <h2>Error Log</h2>
            <div class="error-log" id="error-log">
                {% if data.errors %}
                    {% for error in data.errors %}
                    <div class="error-item">
                        <strong>{{ error.timestamp }}</strong>: {{ error.message }}
                    </div>
                    {% endfor %}
                {% else %}
                    <p>No errors recorded.</p>
                {% endif %}
            </div>
        </div>
    </div>
    
    <script>
        function refreshData() {
            fetch('/api/status')
                .then(response => response.json())
                .then(data => {
                    document.getElementById('last-updated').textContent = data.last_updated;
                    
                    // Update API status
                    const apiStatus = document.getElementById('api-status');
                    apiStatus.textContent = data.api_status;
                    apiStatus.className = 'status';
                    if (data.api_status === 'OK') apiStatus.classList.add('status-ok');
                    else if (data.api_status === 'Warning') apiStatus.classList.add('status-warning');
                    else if (data.api_status === 'Error') apiStatus.classList.add('status-error');
                    else apiStatus.classList.add('status-unknown');
                    
                    // Update DB status
                    const dbStatus = document.getElementById('db-status');
                    dbStatus.textContent = data.database_status;
                    dbStatus.className = 'status';
                    if (data.database_status === 'OK') dbStatus.classList.add('status-ok');
                    else if (data.database_status === 'Warning') dbStatus.classList.add('status-warning');
                    else if (data.database_status === 'Error') dbStatus.classList.add('status-error');
                    else dbStatus.classList.add('status-unknown');
                    
                    // Update endpoints table
                    const endpointsTable = document.getElementById('endpoints-table').getElementsByTagName('tbody')[0];
                    endpointsTable.innerHTML = '';
                    
                    for (const [endpoint, endpointData] of Object.entries(data.endpoints)) {
                        const row = endpointsTable.insertRow();
                        
                        const endpointCell = row.insertCell(0);
                        endpointCell.textContent = endpoint;
                        
                        const methodCell = row.insertCell(1);
                        methodCell.textContent = endpointData.method;
                        
                        const statusCell = row.insertCell(2);
                        const statusSpan = document.createElement('span');
                        statusSpan.textContent = endpointData.status;
                        statusSpan.className = 'status';
                        if (endpointData.status === 'OK') statusSpan.classList.add('status-ok');
                        else if (endpointData.status === 'Warning') statusSpan.classList.add('status-warning');
                        else if (endpointData.status === 'Error') statusSpan.classList.add('status-error');
                        else statusSpan.classList.add('status-unknown');
                        statusCell.appendChild(statusSpan);
                        
                        const responseTimeCell = row.insertCell(3);
                        responseTimeCell.textContent = endpointData.response_time;
                        
                        const lastCheckedCell = row.insertCell(4);
                        lastCheckedCell.textContent = endpointData.last_checked;
                    }
                    
                    // Update performance metrics
                    document.getElementById('avg-response-time').textContent = data.performance.avg_response_time || 'N/A';
                    document.getElementById('error-rate').textContent = data.performance.error_rate || 'N/A';
                    
                    // Update error log
                    const errorLog = document.getElementById('error-log');
                    errorLog.innerHTML = '';
                    
                    if (data.errors && data.errors.length > 0) {
                        for (const error of data.errors) {
                            const errorItem = document.createElement('div');
                            errorItem.className = 'error-item';
                            
                            const timestamp = document.createElement('strong');
                            timestamp.textContent = error.timestamp;
                            
                            errorItem.appendChild(timestamp);
                            errorItem.appendChild(document.createTextNode(': ' + error.message));
                            
                            errorLog.appendChild(errorItem);
                        }
                    } else {
                        const noErrors = document.createElement('p');
                        noErrors.textContent = 'No errors recorded.';
                        errorLog.appendChild(noErrors);
                    }
                })
                .catch(error => {
                    console.error('Error fetching data:', error);
                    alert('Failed to refresh data. See console for details.');
                });
        }
        
        // Auto-refresh every 30 seconds
        setInterval(refreshData, 30000);
    </script>
</body>
</html>
"""

def connect_to_supabase():
    """Connect to Supabase and check connection status."""
    try:
        from supabase import create_client, Client
        
        supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
        
        # Test connection by making a simple query
        response = supabase.from_("property_types").select("*").limit(1).execute()
        
        if hasattr(response, 'data'):
            return True, "Connected to Supabase successfully"
        else:
            return False, "Connected to Supabase but received invalid response"
    except Exception as e:
        return False, f"Error connecting to Supabase: {str(e)}"

def check_api_endpoint(endpoint, method="GET", data=None):
    """Check if an API endpoint is working."""
    try:
        # Import Flask app for local testing
        from app import create_app
        
        app = create_app(testing=True)
        client = app.test_client()
        
        start_time = time.time()
        
        with app.app_context():
            if method == "GET":
                response = client.get(endpoint)
            elif method == "POST":
                response = client.post(endpoint, json=data)
            elif method == "PUT":
                response = client.put(endpoint, json=data)
            elif method == "DELETE":
                response = client.delete(endpoint)
            else:
                return "Error", f"Unsupported method: {method}", 0
            
            response_time = int((time.time() - start_time) * 1000)  # in milliseconds
            
            if response.status_code < 300:
                status = "OK"
            elif response.status_code < 500:
                status = "Warning"
            else:
                status = "Error"
            
            return status, f"Status code: {response.status_code}", response_time
    except Exception as e:
        return "Error", str(e), 0

def monitor_api_endpoints():
    """Monitor all API endpoints and update the monitoring data."""
    global monitoring_data
    
    # Define endpoints to monitor
    endpoints = [
        {"endpoint": "/health", "method": "GET", "data": None},
        {"endpoint": "/api/v1/molecules", "method": "GET", "data": None},
        {"endpoint": "/api/v1/mixtures", "method": "GET", "data": None},
        {"endpoint": "/api/v1/rdkit/properties", "method": "POST", "data": {"molecule_data": "CCO", "input_format": "smiles"}},
        {"endpoint": "/api/v1/compare-properties", "method": "POST", "data": {"ids": [str(uuid.uuid4()), str(uuid.uuid4())]}},
    ]
    
    total_response_time = 0
    error_count = 0
    endpoint_count = len(endpoints)
    
    for endpoint_info in endpoints:
        endpoint = endpoint_info["endpoint"]
        method = endpoint_info["method"]
        data = endpoint_info["data"]
        
        status, message, response_time = check_api_endpoint(endpoint, method, data)
        
        if status == "Error":
            error_count += 1
            monitoring_data["errors"].append({
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "message": f"Error on {method} {endpoint}: {message}"
            })
            # Keep only the last 100 errors
            if len(monitoring_data["errors"]) > 100:
                monitoring_data["errors"] = monitoring_data["errors"][-100:]
        
        total_response_time += response_time
        
        monitoring_data["endpoints"][endpoint] = {
            "method": method,
            "status": status,
            "message": message,
            "response_time": response_time,
            "last_checked": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
    
    # Update performance metrics
    avg_response_time = total_response_time / endpoint_count if endpoint_count > 0 else 0
    error_rate = (error_count / endpoint_count) * 100 if endpoint_count > 0 else 0
    
    monitoring_data["performance"]["response_times"].append(avg_response_time)
    monitoring_data["performance"]["error_rates"].append(error_rate)
    monitoring_data["performance"]["timestamps"].append(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    
    # Keep only the last 100 data points
    if len(monitoring_data["performance"]["response_times"]) > 100:
        monitoring_data["performance"]["response_times"] = monitoring_data["performance"]["response_times"][-100:]
        monitoring_data["performance"]["error_rates"] = monitoring_data["performance"]["error_rates"][-100:]
        monitoring_data["performance"]["timestamps"] = monitoring_data["performance"]["timestamps"][-100:]
    
    # Calculate averages
    monitoring_data["performance"]["avg_response_time"] = sum(monitoring_data["performance"]["response_times"]) / len(monitoring_data["performance"]["response_times"])
    monitoring_data["performance"]["error_rate"] = sum(monitoring_data["performance"]["error_rates"]) / len(monitoring_data["performance"]["error_rates"])
    
    # Update overall API status
    if error_count == 0:
        monitoring_data["api_status"] = "OK"
    elif error_count < endpoint_count / 2:
        monitoring_data["api_status"] = "Warning"
    else:
        monitoring_data["api_status"] = "Error"
    
    # Update last updated timestamp
    monitoring_data["last_updated"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def check_database_status():
    """Check the database status and update the monitoring data."""
    global monitoring_data
    
    success, message = connect_to_supabase()
    
    if success:
        monitoring_data["database_status"] = "OK"
    else:
        monitoring_data["database_status"] = "Error"
        monitoring_data["errors"].append({
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "message": f"Database error: {message}"
        })
        # Keep only the last 100 errors
        if len(monitoring_data["errors"]) > 100:
            monitoring_data["errors"] = monitoring_data["errors"][-100:]

def update_monitoring_data():
    """Update all monitoring data."""
    try:
        check_database_status()
        monitor_api_endpoints()
        logger.info("Monitoring data updated successfully")
    except Exception as e:
        logger.error(f"Error updating monitoring data: {str(e)}")
        monitoring_data["errors"].append({
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "message": f"Monitoring error: {str(e)}"
        })

def monitoring_thread():
    """Thread function to periodically update monitoring data."""
    while True:
        update_monitoring_data()
        time.sleep(30)  # Update every 30 seconds

@app.route('/')
def dashboard():
    """Render the dashboard."""
    return render_template_string(DASHBOARD_TEMPLATE, data=monitoring_data)

@app.route('/api/status')
def api_status():
    """Return the current monitoring data as JSON."""
    return jsonify(monitoring_data)

def main():
    """Main function to run the monitoring dashboard."""
    parser = argparse.ArgumentParser(description='API Monitoring Dashboard')
    parser.add_argument('--port', type=int, default=5001, help='Port to run the dashboard on')
    args = parser.parse_args()
    
    print("\n" + "=" * 80)
    print("CryoProtect v2 - API Monitoring Dashboard")
    print("=" * 80)
    
    try:
        # Start the monitoring thread
        monitor_thread = threading.Thread(target=monitoring_thread)
        monitor_thread.daemon = True
        monitor_thread.start()
        
        # Initial update
        update_monitoring_data()
        
        # Open browser
        webbrowser.open(f"http://localhost:{args.port}")
        
        # Run the Flask app
        app.run(host='0.0.0.0', port=args.port, debug=False)
        
        return 0
    except KeyboardInterrupt:
        print("\nMonitoring dashboard stopped by user")
        return 0
    except Exception as e:
        logger.error(f"Error in monitoring dashboard: {str(e)}")
        print(f"\nError in monitoring dashboard: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())