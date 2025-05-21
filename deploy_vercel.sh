#!/bin/bash
# Ultra-minimal Vercel deployment script that requires no Python in build step

set -e  # Exit on any error

echo "🚀 Starting minimal Vercel deployment..."

# Clean up any temporary files
echo "🧹 Cleaning up temporary files..."
rm -rf .vercel/output 2>/dev/null || true
rm -rf deploy_temp 2>/dev/null || true

# Create a temporary minimal project for deployment
echo "🔧 Creating an ultra-minimal deployment structure..."
mkdir -p deploy_temp/api deploy_temp/static
mkdir -p deploy_temp/static/css

# Create essential API files from scratch directly in shell
echo "📝 Creating minimal serverless functions..."

# Create index.py
cat > deploy_temp/api/index.py << 'EOF'
import json
from datetime import datetime

def handler(request):
    """
    Main API handler for Vercel serverless function
    """
    try:
        # Extract path from request
        path = request.get('path', '')
        
        # Health check endpoint
        if '/api/v1/health' in path:
            return {
                'statusCode': 200,
                'headers': {
                    'Content-Type': 'application/json'
                },
                'body': json.dumps({
                    'status': 'ok',
                    'version': '1.0.0',
                    'timestamp': datetime.now().isoformat(),
                    'environment': 'vercel-production'
                })
            }
        
        # Default API response
        return {
            'statusCode': 404,
            'headers': {
                'Content-Type': 'application/json'
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
                'Content-Type': 'application/json'
            },
            'body': json.dumps({
                'error': str(e),
                'timestamp': datetime.now().isoformat()
            })
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
            <a href="/api/v1/health" class="api-link">Test API Health Endpoint</a>
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
      "runtime": "python3",
      "memory": 256,
      "maxDuration": 10
    },
    "api/app.py": {
      "runtime": "python3",
      "memory": 128,
      "maxDuration": 5
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
vercel --prod

echo "✅ Deployment process completed!"
echo "Note: Check the Vercel dashboard for final deployment status"

# Return to original directory and clean up
cd ..
rm -rf deploy_temp

echo "🎉 Done!"