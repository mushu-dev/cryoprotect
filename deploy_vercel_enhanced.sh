#!/bin/bash
# Enhanced Vercel deployment script with API connectivity

set -e  # Exit on any error

BACKEND_URL="https://cryoprotect-8030e4025428.herokuapp.com"
echo "🚀 Starting enhanced Vercel deployment with backend: $BACKEND_URL"

# Clean up any temporary files
echo "🧹 Cleaning up temporary files..."
rm -rf .vercel/output 2>/dev/null || true
rm -rf deploy_temp 2>/dev/null || true

# Create a temporary project for deployment
echo "🔧 Creating deployment structure..."
mkdir -p deploy_temp/api deploy_temp/static/js deploy_temp/static/css

# Create essential API files
echo "📝 Creating serverless functions..."

# Create index.py
cat > deploy_temp/api/index.py << EOF
import json
import os
from datetime import datetime
import requests

# API connectivity configuration
BACKEND_URL = "${BACKEND_URL}"

def handler(request):
    """
    Main API handler for Vercel serverless function
    """
    try:
        # Extract path and headers from request
        path = request.get('path', '')
        headers = request.get('headers', {})
        origin = headers.get('origin', 'unknown')
        
        # Health check endpoint
        if '/api/v1/health' in path:
            return {
                'statusCode': 200,
                'headers': {
                    'Content-Type': 'application/json',
                    'Access-Control-Allow-Origin': '*'
                },
                'body': json.dumps({
                    'status': 'ok',
                    'version': '1.0.0',
                    'timestamp': datetime.now().isoformat(),
                    'environment': 'vercel-production'
                })
            }
            
        # Backend connectivity test
        if '/api/v1/connect-test' in path:
            try:
                # Forward request to backend connectivity endpoint
                response = requests.get(
                    f"{BACKEND_URL}/api/connect",
                    headers={
                        'Origin': 'https://frontend-cryoprotect.vercel.app'
                    }
                )
                
                # Get data from backend
                data = response.json()
                
                # Add verification data
                data['vercel_verified'] = True
                data['vercel_origin'] = origin
                data['backend_url'] = BACKEND_URL
                
                return {
                    'statusCode': 200,
                    'headers': {
                        'Content-Type': 'application/json',
                        'Access-Control-Allow-Origin': '*'
                    },
                    'body': json.dumps(data)
                }
            except Exception as e:
                return {
                    'statusCode': 500,
                    'headers': {
                        'Content-Type': 'application/json',
                        'Access-Control-Allow-Origin': '*'
                    },
                    'body': json.dumps({
                        'status': 'error',
                        'message': f'Backend connectivity error: {str(e)}',
                        'backend_url': BACKEND_URL
                    })
                }
        
        # Default API response
        return {
            'statusCode': 404,
            'headers': {
                'Content-Type': 'application/json',
                'Access-Control-Allow-Origin': '*'
            },
            'body': json.dumps({
                'error': 'Unknown API route',
                'path': path,
                'timestamp': datetime.now().isoformat()
            })
        }
    except Exception as e:
        return {
            'statusCode': 500,
            'headers': {
                'Content-Type': 'application/json',
                'Access-Control-Allow-Origin': '*'
            },
            'body': json.dumps({
                'error': str(e),
                'timestamp': datetime.now().isoformat()
            })
        }
EOF

# Create connectivity test page
cat > deploy_temp/api/connectivity.py << EOF
from datetime import datetime

def handler(request):
    """
    Connectivity test page handler
    """
    backend_url = "${BACKEND_URL}"
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>CryoProtect - API Connectivity Test</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {{ 
            font-family: -apple-system, system-ui, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; 
            line-height: 1.6; 
            color: #333; 
            max-width: 800px; 
            margin: 0 auto; 
            padding: 20px; 
            background-color: #f8f9fa;
        }}
        h1, h2 {{ 
            color: #2c3e50; 
            margin-bottom: 0.5em;
        }}
        .container {{ 
            border: 1px solid #ddd; 
            border-radius: 8px; 
            padding: 30px; 
            margin-top: 40px; 
            background-color: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }}
        .success {{ 
            color: #28a745; 
            font-weight: bold;
        }}
        .error {{
            color: #dc3545;
            font-weight: bold;
        }}
        .pending {{
            color: #ffc107;
            font-weight: bold;
        }}
        pre {{
            background: #f8f9fa;
            border-radius: 4px;
            padding: 10px;
            overflow: auto;
            font-size: 14px;
            border: 1px solid #ddd;
        }}
        .api-link {{
            display: inline-block;
            margin-top: 1em;
            padding: 8px 16px;
            background-color: #f8f9fa;
            border-radius: 4px;
            text-decoration: none;
            color: #007bff;
            border: 1px solid #ddd;
        }}
        .api-link:hover {{
            background-color: #e9ecef;
            color: #0056b3;
        }}
        .header {{
            margin-bottom: 1.5em;
            padding-bottom: 1em;
            border-bottom: 1px solid #eee;
        }}
        .footer {{
            margin-top: 2em;
            padding-top: 1em;
            border-top: 1px solid #eee;
            font-size: 0.9em;
            color: #6c757d;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 16px 0;
        }}
        th, td {{
            text-align: left;
            padding: 8px;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #f8f9fa;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>CryoProtect API Connectivity Test</h1>
            <p>Verify connectivity between Vercel frontend and Heroku backend</p>
        </div>
        
        <div>
            <h2>Configuration</h2>
            <table>
                <tr>
                    <th>Backend URL</th>
                    <td id="backendUrl">{backend_url}</td>
                </tr>
                <tr>
                    <th>Frontend URL</th>
                    <td id="frontendUrl">Loading...</td>
                </tr>
                <tr>
                    <th>Current Time</th>
                    <td>{datetime.now().strftime("%Y-%m-%d %H:%M:%S UTC")}</td>
                </tr>
            </table>
        </div>
        
        <div>
            <h2>Connectivity Test</h2>
            <p>Status: <span id="status" class="pending">⏳ Testing connectivity...</span></p>
            <pre id="results">Running tests...</pre>
            
            <p><strong>Test Endpoints:</strong></p>
            <ul>
                <li><a href="/api/v1/health" target="_blank" class="api-link">Vercel Health Check</a></li>
                <li><a href="/api/v1/connect-test" target="_blank" class="api-link">Backend Connectivity Test</a></li>
                <li><a href="{backend_url}/health" target="_blank" class="api-link">Heroku Health Check</a></li>
                <li><a href="{backend_url}/api/connect" target="_blank" class="api-link">Heroku Connectivity API</a></li>
            </ul>
        </div>
        
        <div class="footer">
            <p>© 2025 CryoProtect Team</p>
        </div>
    </div>
    
    <script>
        document.getElementById('frontendUrl').textContent = window.location.origin;
        
        async function runConnectivityTests() {
            const statusEl = document.getElementById('status');
            const resultsEl = document.getElementById('results');
            let results = {};
            
            try {
                // Test Vercel Health
                try {
                    const vercelHealthResponse = await fetch('/api/v1/health');
                    const vercelHealth = await vercelHealthResponse.json();
                    results.vercelHealth = {
                        status: vercelHealthResponse.status,
                        data: vercelHealth
                    };
                } catch (error) {
                    results.vercelHealth = {
                        status: 'error',
                        message: error.message
                    };
                }
                
                // Test Backend Connectivity
                try {
                    const connectivityResponse = await fetch('/api/v1/connect-test');
                    const connectivity = await connectivityResponse.json();
                    results.backendConnectivity = {
                        status: connectivityResponse.status,
                        data: connectivity
                    };
                } catch (error) {
                    results.backendConnectivity = {
                        status: 'error',
                        message: error.message
                    };
                }
                
                // Direct test to Backend (might be blocked by CORS)
                try {
                    const directResponse = await fetch('{backend_url}/health', {
                        mode: 'cors'
                    });
                    const directData = await directResponse.json();
                    results.directBackend = {
                        status: directResponse.status,
                        data: directData
                    };
                } catch (error) {
                    results.directBackend = {
                        status: 'error',
                        message: error.message
                    };
                }
                
                // Determine overall status
                const isSuccess = results.backendConnectivity?.status === 200;
                
                // Update status
                if (isSuccess) {
                    statusEl.textContent = '✅ Success! Backend is reachable';
                    statusEl.className = 'success';
                } else {
                    statusEl.textContent = '❌ Failed! Cannot connect to backend';
                    statusEl.className = 'error';
                }
                
                // Display results
                resultsEl.textContent = JSON.stringify(results, null, 2);
                
            } catch (error) {
                statusEl.textContent = '❌ Error running tests: ' + error.message;
                statusEl.className = 'error';
                resultsEl.textContent = error.stack;
            }
        }
        
        // Run the tests when page loads
        window.onload = runConnectivityTests;
    </script>
</body>
</html>"""
    
    return {
        'statusCode': 200,
        'headers': {
            'Content-Type': 'text/html'
        },
        'body': html
    }
EOF

# Create app.py (landing page)
cat > deploy_temp/api/app.py << 'EOF'
from datetime import datetime

def handler(request):
    """
    Landing page handler
    """
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>CryoProtect</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {{ 
            font-family: -apple-system, system-ui, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; 
            line-height: 1.6; 
            color: #333; 
            max-width: 800px; 
            margin: 0 auto; 
            padding: 20px; 
            background-color: #f8f9fa;
        }}
        h1 {{ 
            color: #2c3e50; 
            margin-bottom: 0.5em;
        }}
        .container {{ 
            border: 1px solid #ddd; 
            border-radius: 8px; 
            padding: 30px; 
            margin-top: 40px; 
            background-color: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }}
        .success {{ 
            color: #28a745; 
            font-weight: bold;
        }}
        .api-link {{
            display: inline-block;
            margin-top: 1em;
            padding: 8px 16px;
            background-color: #f8f9fa;
            border-radius: 4px;
            text-decoration: none;
            color: #007bff;
            border: 1px solid #ddd;
        }}
        .api-link:hover {{
            background-color: #e9ecef;
            color: #0056b3;
        }}
        .header {{
            margin-bottom: 1.5em;
            padding-bottom: 1em;
            border-bottom: 1px solid #eee;
        }}
        .footer {{
            margin-top: 2em;
            padding-top: 1em;
            border-top: 1px solid #eee;
            font-size: 0.9em;
            color: #6c757d;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>CryoProtect</h1>
            <p>Molecular cryoprotectant database and analysis platform</p>
        </div>
        
        <div>
            <p><strong>Status:</strong> <span class="success">✅ Deployment successful!</span></p>
            <p><strong>Environment:</strong> Vercel Production</p>
            <p><strong>Server time:</strong> {current_time}</p>
            <p><strong>Version:</strong> 1.0.0</p>
        </div>
        
        <div>
            <a href="/api/v1/health" class="api-link">Test API Health</a>
            <a href="/connectivity" class="api-link">Test Backend Connectivity</a>
        </div>
        
        <div class="footer">
            <p>© 2025 CryoProtect Team</p>
        </div>
    </div>
</body>
</html>"""
    
    return {
        'statusCode': 200,
        'headers': {
            'Content-Type': 'text/html'
        },
        'body': html
    }
EOF

# Create minimal requirements.txt
cat > deploy_temp/requirements.txt << 'EOF'
requests==2.32.3
EOF

# Create more explicit vercel.json 
cat > deploy_temp/vercel.json << EOF
{
  "version": 2,
  "installCommand": "pip install -r requirements.txt",
  "buildCommand": "echo 'Skipping build'",
  "outputDirectory": ".",
  "framework": null,
  "env": {
    "BACKEND_URL": "${BACKEND_URL}"
  },
  "functions": {
    "api/index.py": {
      "runtime": "python@3.9",
      "memory": 256,
      "maxDuration": 10
    },
    "api/app.py": {
      "runtime": "python@3.9",
      "memory": 128,
      "maxDuration": 5
    },
    "api/connectivity.py": {
      "runtime": "python@3.9",
      "memory": 128,
      "maxDuration": 10
    }
  },
  "routes": [
    { 
      "src": "/api/v1/(.*)", 
      "dest": "/api/index.py"
    },
    { 
      "src": "/static/(.*)", 
      "dest": "/static/$1",
      "headers": {
        "cache-control": "public, max-age=86400"
      }
    },
    {
      "src": "/connectivity",
      "dest": "/api/connectivity.py",
      "headers": {
        "cache-control": "no-cache, no-store, must-revalidate"
      }
    },
    { 
      "src": "/(.*)", 
      "dest": "/api/app.py"
    }
  ],
  "regions": ["iad1"]
}
EOF

# Create static placeholder
cat > deploy_temp/static/css/placeholder.css << 'EOF'
/* Placeholder CSS file */
body {
  font-family: system-ui, sans-serif;
}
EOF

# Create a simple package.json
cat > deploy_temp/package.json << 'EOF'
{
  "name": "cryoprotect-enhanced",
  "version": "1.0.0",
  "private": true
}
EOF

# Create essential __init__.py files to ensure Python packaging works
touch deploy_temp/api/__init__.py

# Switch to deployment directory
cd deploy_temp

# Deploy to Vercel
echo "🚀 Deploying enhanced build to Vercel..."
vercel --prod --yes

echo "✅ Deployment process completed!"
echo "Note: Check the Vercel dashboard for final deployment status"

# Return to original directory and clean up
cd ..
rm -rf deploy_temp

echo "🎉 Done!"