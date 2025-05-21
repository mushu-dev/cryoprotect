#!/bin/bash
# Script to try multiple Vercel authentication bypass methods

set -e  # Exit on any error

echo "🔍 Attempting multiple Vercel authentication bypass methods..."

# Create a clean deployment directory
echo "🧹 Creating clean deployment directory..."
rm -rf bypass_temp
mkdir -p bypass_temp

# Create a minimal index.html
cat > bypass_temp/index.html << EOF
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>CryoProtect - Authentication Bypass Test</title>
  <style>
    body {
      font-family: system-ui, sans-serif;
      max-width: 800px;
      margin: 0 auto;
      padding: 20px;
      line-height: 1.6;
    }
    h1 { color: #2c3e50; }
    .container {
      border: 1px solid #ddd;
      border-radius: 8px;
      padding: 30px;
      margin-top: 40px;
      background-color: white;
      box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
    .success { color: #28a745; font-weight: bold; }
    pre {
      background: #f8f9fa;
      padding: 10px;
      border-radius: 4px;
      overflow: auto;
    }
    .btn {
      display: inline-block;
      padding: 8px 16px;
      background: #0070f3;
      color: white;
      border-radius: 4px;
      text-decoration: none;
      margin-top: 10px;
    }
  </style>
</head>
<body>
  <div class="container">
    <h1>CryoProtect API Authentication Bypass Test</h1>
    <p>If you can see this page, the authentication bypass was successful! 🎉</p>

    <h2>Backend Connection</h2>
    <p>Backend URL: <code>https://cryoprotect-8030e4025428.herokuapp.com</code></p>
    <button onclick="testConnection()" class="btn">Test Connection</button>
    <pre id="result">Click button to test connection...</pre>

    <script>
      function testConnection() {
        const resultEl = document.getElementById('result');
        resultEl.innerHTML = 'Testing connection...';
        
        fetch('https://cryoprotect-8030e4025428.herokuapp.com/health')
          .then(response => response.json())
          .then(data => {
            resultEl.innerHTML = JSON.stringify(data, null, 2);
          })
          .catch(error => {
            resultEl.innerHTML = 'Error: ' + error.message;
          });
      }
    </script>
  </div>
</body>
</html>
EOF

# Create a minimal vercel.json with all bypass methods
cat > bypass_temp/vercel.json << EOF
{
  "version": 2,
  "public": true,
  "github": {
    "enabled": false
  },
  "build": {
    "env": {
      "VERCEL_PROTECTION_BYPASS": "1"
    }
  },
  "env": {
    "VERCEL_PROTECTION_BYPASS": "1"
  }
}
EOF

# Create a package.json with specific settings
cat > bypass_temp/package.json << EOF
{
  "name": "cryoprotect-bypass",
  "version": "1.0.0",
  "private": true,
  "scripts": {
    "vercel-build": "echo 'Skip build step'"
  }
}
EOF

# Move to the temp directory
cd bypass_temp

# Try deploying with special flags to bypass protection
echo "🚀 Deploying with authentication bypass methods..."

# Method 1: Try with protection bypass env
echo "Method 1: Protection bypass environment variable"
vercel --env VERCEL_PROTECTION_BYPASS=1 --prod --confirm

# Method 2: Output deployment URL
echo
echo "✅ Deployment complete! Test the URL above to see if authentication is bypassed."
echo
echo "If authentication is still required, try these manual steps:"
echo "1. Go to https://vercel.com/dashboard"
echo "2. Click on your project"
echo "3. Go to Settings → Security"
echo "4. Look for 'Password Protection' and disable it"
echo "5. Redeploy with: vercel --prod"
echo
echo "Alternative: Consider Netlify deployment with this simple command:"
echo "npx netlify-cli deploy --dir=. --prod"