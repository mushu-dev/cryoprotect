#!/bin/bash
# Simpler deployment to Netlify without npm install

set -e  # Exit on error

# Create a sample directory with a basic page
mkdir -p test-deploy
cat > test-deploy/index.html << EOF
<!DOCTYPE html>
<html>
<head>
  <title>CryoProtect on Netlify</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      margin: 0;
      padding: 20px;
      display: flex;
      justify-content: center;
      align-items: center;
      min-height: 100vh;
      background-color: #f5f8fa;
    }
    .container {
      background-color: white;
      border-radius: 8px;
      box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
      padding: 30px;
      max-width: 600px;
      text-align: center;
    }
    h1 {
      color: #2c3e50;
    }
    .api-status {
      margin-top: 20px;
      padding: 15px;
      background-color: #f1f1f1;
      border-radius: 4px;
    }
    button {
      background-color: #3498db;
      color: white;
      border: none;
      padding: 10px 20px;
      border-radius: 4px;
      cursor: pointer;
      margin-top: 20px;
    }
    button:hover {
      background-color: #2980b9;
    }
  </style>
</head>
<body>
  <div class="container">
    <h1>CryoProtect Successfully Deployed to Netlify! 🚀</h1>
    <p>This is a test deployment page for verification purposes.</p>
    
    <div class="api-status">
      <h3>API Connection Status</h3>
      <p id="status-message">Click the button below to test connection to the Heroku backend.</p>
      <button id="test-connection">Test Connection</button>
    </div>
    
    <script>
      document.getElementById('test-connection').addEventListener('click', async () => {
        const statusMessage = document.getElementById('status-message');
        statusMessage.textContent = 'Testing connection...';
        
        try {
          const response = await fetch('https://cryoprotect-8030e4025428.herokuapp.com/v1/health');
          
          if (response.ok) {
            statusMessage.textContent = '✅ Successfully connected to the Heroku backend!';
            statusMessage.style.color = 'green';
          } else {
            statusMessage.textContent = '❌ Failed to connect to the backend. Status: ' + response.status;
            statusMessage.style.color = 'red';
          }
        } catch (error) {
          statusMessage.textContent = '❌ Error connecting to backend: ' + error.message;
          statusMessage.style.color = 'red';
        }
      });
    </script>
  </div>
</body>
</html>
EOF

# Create a new Netlify site
SITE_NAME="cryoprotect-$(date +%s)"
echo "Creating new Netlify site: $SITE_NAME"

# Create a new site and extract site ID using a direct approach
NETLIFY_RESPONSE=$(netlify sites:create --name "$SITE_NAME" --json || echo "{}")
SITE_ID=$(echo "$NETLIFY_RESPONSE" | grep -o '"id":"[^"]*"' | cut -d'"' -f4)

if [ -z "$SITE_ID" ]; then
  echo "Failed to extract site ID. Using manual deployment instead..."
  netlify deploy --dir=test-deploy
else
  echo "Site created with ID: $SITE_ID"
  
  # Link the site
  echo "Linking to site..."
  netlify link --id "$SITE_ID" --force
  
  # Deploy to Netlify
  echo "Deploying to Netlify..."
  netlify deploy --prod --dir=test-deploy
  
  # Print site URL
  echo "Deployment complete!"
  echo "Your site is available at: https://$SITE_NAME.netlify.app"
fi