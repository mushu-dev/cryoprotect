import sys
import os
import json
from datetime import datetime

# Add the parent directory to sys.path to import app.py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Simple response for main app page
def handler(request):
    """Handler for root path - serves HTML landing page"""
    # Check for protection bypass token
    headers = request.get('headers', {})
    query_params = request.get('queryStringParameters', {})
    
    bypass_header = headers.get('x-protection-bypass')
    bypass_query = query_params.get('bypass')
    expected_bypass = os.environ.get('PROTECTION_BYPASS', 'TAt23KbtFE8dkZobJU3hpgTP4L5ja07V')
    
    # If bypass token doesn't match, return access denied
    if not (bypass_header == expected_bypass or bypass_query == expected_bypass):
        return {
            'statusCode': 403,
            'headers': {
                'Content-Type': 'text/html'
            },
            'body': '<h1>Access Denied</h1><p>This deployment requires a valid protection bypass token.</p>'
        }
    
    html = """
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>CryoProtect - Vercel Deployment</title>
        <style>
            body { font-family: Arial, sans-serif; line-height: 1.6; margin: 40px; }
            h1 { color: #333; }
            .container { max-width: 800px; margin: 0 auto; }
            .info { background: #f9f9f9; padding: 20px; border-radius: 5px; }
            .links a { margin-right: 15px; color: #0070f3; text-decoration: none; }
            .links a:hover { text-decoration: underline; }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>CryoProtect Analyzer API</h1>
            <div class="info">
                <p>This is the Vercel serverless deployment of CryoProtect Analyzer API.</p>
                <p>The API endpoints are available at <code>/api/v1/...</code></p>
            </div>
            <div class="links">
                <h3>API Resources:</h3>
                <p>
                    <a href="/api/v1/health">Health Check</a>
                    <a href="/api/v1/docs">API Documentation</a>
                </p>
            </div>
            <div class="info">
                <p>Deployment Time: """ + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + """</p>
            </div>
        </div>
    </body>
    </html>
    """
    
    return {
        'statusCode': 200,
        'headers': {
            'Content-Type': 'text/html'
        },
        'body': html
    }