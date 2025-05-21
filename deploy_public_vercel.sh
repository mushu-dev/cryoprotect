#!/bin/bash
# Deployment script for Vercel with public access enabled

set -e  # Exit on any error

BACKEND_URL="https://cryoprotect-8030e4025428.herokuapp.com"
echo "🚀 Starting Vercel deployment with backend: $BACKEND_URL"

# Clean up any temporary files
echo "🧹 Cleaning up temporary files..."
rm -rf .vercel/output 2>/dev/null || true
rm -rf deploy_temp 2>/dev/null || true

# Create a temporary project for deployment
echo "🔧 Creating deployment structure..."
mkdir -p deploy_temp/api deploy_temp/static/js deploy_temp/static/css

# Create index.py for API endpoints
cat > deploy_temp/api/index.py << EOF
import json
from datetime import datetime

def handler(request):
    """
    Main API handler for Vercel serverless function
    """
    try:
        # Extract path from request
        path = request.get('path', '')
        headers = request.get('headers', {})
        origin = headers.get('origin', 'Unknown')
        
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
            
        # Backend connectivity info endpoint
        if '/api/v1/backend-info' in path:
            return {
                'statusCode': 200,
                'headers': {
                    'Content-Type': 'application/json',
                    'Access-Control-Allow-Origin': '*'
                },
                'body': json.dumps({
                    'status': 'success',
                    'backend_url': '${BACKEND_URL}',
                    'api_connect_url': '${BACKEND_URL}/api/connect',
                    'health_url': '${BACKEND_URL}/health',
                    'frontend_origin': origin,
                    'timestamp': datetime.now().isoformat()
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

# Create app.py landing page
cat > deploy_temp/api/app.py << EOF
from datetime import datetime

def handler(request):
    """
    Landing page handler
    """
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    backend_url = "${BACKEND_URL}"
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>CryoProtect API</title>
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
            margin-right: 10px;
            margin-bottom: 10px;
        }}
        .api-link:hover {{
            background-color: #e9ecef;
            color: #0056b3;
        }}
        .test-button {{
            background-color: #007bff;
            color: white;
            border: none;
            padding: 8px 16px;
            border-radius: 4px;
            cursor: pointer;
            font-weight: 500;
            margin-top: 10px;
        }}
        .test-button:hover {{
            background-color: #0069d9;
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
            text-align: center;
        }}
        pre {{
            background: #f8f9fa;
            border-radius: 4px;
            padding: 10px;
            overflow: auto;
        }}
        #results {{
            height: 200px;
            overflow: auto;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>CryoProtect API</h1>
            <p>Molecular cryoprotectant database and analysis platform</p>
        </div>
        
        <div>
            <p><strong>Status:</strong> <span class="success">✅ Deployment successful!</span></p>
            <p><strong>Environment:</strong> Vercel Production</p>
            <p><strong>Server time:</strong> {current_time}</p>
            <p><strong>Backend URL:</strong> <a href="{backend_url}" target="_blank">{backend_url}</a></p>
        </div>
        
        <div>
            <h2>API Endpoints</h2>
            <a href="/api/v1/health" class="api-link">Vercel Health Check</a>
            <a href="/api/v1/backend-info" class="api-link">Backend Info</a>
            <a href="{backend_url}/health" class="api-link">Heroku Health Check</a>
            <a href="{backend_url}/api/connect" class="api-link">Heroku Connectivity API</a>
            <a href="{backend_url}/api/molecules?limit=5" class="api-link">List Molecules (5)</a>
        </div>
        
        <div>
            <h2>Connection Test</h2>
            <p>Click the button to test connectivity to the backend:</p>
            <button id="testButton" class="test-button">Test Backend Connection</button>
            <pre id="results">No test run yet...</pre>
        </div>
        
        <div class="footer">
            <p>© 2025 CryoProtect Team</p>
        </div>
    </div>
    
    <script>
        document.getElementById('testButton').addEventListener('click', async function() {
            const resultsEl = document.getElementById('results');
            resultsEl.textContent = 'Running test...';
            
            try {{
                const response = await fetch('{backend_url}/api/connect', {{
                    headers: {{
                        'Origin': window.location.origin
                    }}
                }});
                
                const data = await response.json();
                
                resultsEl.textContent = JSON.stringify({{
                    status: response.status,
                    statusText: response.statusText,
                    data: data,
                    timestamp: new Date().toISOString()
                }}, null, 2);
            }} catch (error) {{
                resultsEl.textContent = 'Error: ' + error.message;
            }}
        }});
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

# Create minimal requirements.txt
cat > deploy_temp/requirements.txt << EOF
# No external dependencies for this deployment
EOF

# Create custom vercel.json with public flag
cat > deploy_temp/vercel.json << EOF
{
  "version": 2,
  "public": true,
  "github": {
    "enabled": false
  },
  "installCommand": "echo 'Skipping install'",
  "buildCommand": "echo 'Skipping build'",
  "outputDirectory": ".",
  "framework": null,
  "functions": {
    "api/index.py": {
      "memory": 256,
      "maxDuration": 10
    },
    "api/app.py": {
      "memory": 128,
      "maxDuration": 5
    }
  },
  "routes": [
    { 
      "src": "/api/v1/(.*)", 
      "dest": "/api/index.py",
      "headers": {
        "Access-Control-Allow-Origin": "*"
      }
    },
    { 
      "src": "/static/(.*)", 
      "dest": "/static/$1",
      "headers": {
        "cache-control": "public, max-age=86400"
      }
    },
    { 
      "src": "/(.*)", 
      "dest": "/api/app.py"
    }
  ],
  "regions": ["iad1"],
  "env": {
    "VERCEL_PROTECTION_BYPASS": "1"
  }
}
EOF

# Create a placeholder package.json
cat > deploy_temp/package.json << EOF
{
  "name": "cryoprotect-public",
  "version": "1.0.0",
  "private": true
}
EOF

# Create necessary folders and files
touch deploy_temp/api/__init__.py
mkdir -p deploy_temp/static/css
cat > deploy_temp/static/css/placeholder.css << EOF
/* Placeholder CSS file */
body {
  font-family: system-ui, sans-serif;
}
EOF

# Deploy to Vercel
cd deploy_temp

echo "🚀 Deploying to Vercel with public access..."
vercel --public --yes

echo "✅ Deployment completed."
echo "Note: Check the URL above for your new public deployment"
echo "If authentication is still enabled, try using the deployment command with team:"
echo "vercel --scope team_name --public --yes"

# Return to original directory and clean up
cd ..