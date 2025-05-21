from datetime import datetime

def handler(request):
    """
    Landing page handler
    """
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    backend_url = "https://cryoprotect-8030e4025428.herokuapp.com"
    
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
