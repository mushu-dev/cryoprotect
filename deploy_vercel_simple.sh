#!/bin/bash
# Ultra-minimal Vercel deployment script that requires no Python in build step

set -e  # Exit on any error

BACKEND_URL="https://cryoprotect-8030e4025428.herokuapp.com"
echo "🚀 Starting minimal Vercel deployment with backend: $BACKEND_URL..."

# Clean up any temporary files
echo "🧹 Cleaning up temporary files..."
rm -rf .vercel/output 2>/dev/null || true
rm -rf deploy_temp 2>/dev/null || true

# Create a temporary minimal project for deployment
echo "🔧 Creating an ultra-minimal deployment structure..."
mkdir -p deploy_temp/api deploy_temp/static
mkdir -p deploy_temp/static/css deploy_temp/static/js

# Create essential API files from scratch directly in shell
echo "📝 Creating minimal serverless functions..."

# Create index.py
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

# Create connectivity.py
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
                <li><a href="/api/v1/backend-info" target="_blank" class="api-link">Backend Info</a></li>
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
                
                // Get Backend Info
                try {
                    const backendInfoResponse = await fetch('/api/v1/backend-info');
                    const backendInfo = await backendInfoResponse.json();
                    results.backendInfo = {
                        status: backendInfoResponse.status,
                        data: backendInfo
                    };
                } catch (error) {
                    results.backendInfo = {
                        status: 'error',
                        message: error.message
                    };
                }
                
                // Direct test to Backend (might be blocked by CORS)
                try {
                    const directResponse = await fetch('{backend_url}/api/connect', {
                        headers: {
                            'Origin': window.location.origin
                        }
                    });
                    const directData = await directResponse.json();
                    results.directBackend = {
                        status: directResponse.status,
                        data: directData
                    };
                    
                    // If we got here, connectivity is working!
                    statusEl.textContent = '✅ Success! Backend is reachable';
                    statusEl.className = 'success';
                    
                } catch (error) {
                    results.directBackend = {
                        status: 'error',
                        message: error.message
                    };
                    
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
            <p><strong>Backend:</strong> <a href="{backend_url}" target="_blank">{backend_url}</a></p>
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
# No external dependencies for the minimal deployment
EOF

# Create connection-test.js
cat > deploy_temp/static/js/connection-test.js << EOF
// API Connection test utility
document.addEventListener('DOMContentLoaded', () => {
  const backendUrl = '${BACKEND_URL}';
  
  async function testConnection() {
    try {
      const response = await fetch(backendUrl + '/api/connect', {
        headers: {
          'Origin': window.location.origin
        }
      });
      
      const data = await response.json();
      return {
        success: true,
        data: data
      };
    } catch (error) {
      return {
        success: false,
        error: error.message
      };
    }
  }
  
  window.testApiConnection = testConnection;
});
EOF

# Create more explicit vercel.json 
cat > deploy_temp/vercel.json << 'EOF'
{
  "version": 2,
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
    },
    "api/connectivity.py": {
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
  "name": "cryoprotect-minimal",
  "version": "1.0.0",
  "private": true
}
EOF

# Create essential __init__.py files to ensure Python packaging works
touch deploy_temp/api/__init__.py

# Switch to deployment directory
cd deploy_temp

# Deploy to Vercel
echo "🚀 Deploying minimal build to Vercel..."
vercel --prod --yes

echo "✅ Deployment process completed!"
echo "Note: Check the Vercel dashboard for final deployment status"

# Return to original directory and clean up
cd ..
rm -rf deploy_temp

echo "🎉 Done!"