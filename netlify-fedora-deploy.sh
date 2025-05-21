#!/bin/bash
set -e

echo "🚀 Running direct Netlify deployment for Fedora systems..."

# Check if necessary commands are available
if ! command -v zip &> /dev/null || ! command -v curl &> /dev/null; then
  echo "⚠️ This script requires zip and curl. Installing..."
  sudo dnf install -y zip curl
fi

# Create a clean deployment directory
DEPLOY_DIR="netlify-direct-deploy"
mkdir -p $DEPLOY_DIR
rm -rf $DEPLOY_DIR/*

# Copy pre-built static files
echo "📦 Copying existing static files..."
cp -r frontend/out/* $DEPLOY_DIR/

# Create _redirects file
echo "📝 Creating _redirects file..."
cat > $DEPLOY_DIR/_redirects << EOL
# API routes redirect to Heroku backend
/api/*  https://cryoprotect-8030e4025428.herokuapp.com/api/:splat  200

# Handle dynamic routes for molecules
/molecules/*  /molecules/[id].html  200

# Handle dynamic routes for mixtures
/mixtures/*  /mixtures/[id].html  200

# SPA fallback for all other routes
/*  /index.html  200
EOL

echo "📁 Packaging files for direct deployment..."
DEPLOY_ZIP="deploy.zip"
cd $DEPLOY_DIR
zip -r -q ../$DEPLOY_ZIP *
cd ..

echo "⚠️ IMPORTANT: To deploy directly, we need a Netlify personal access token."
echo "Please generate one at https://app.netlify.com/user/applications/personal if you don't have one."
read -p "Enter your Netlify personal access token: " TOKEN

# Get site ID
SITE_ID=$(curl -s -H "Authorization: Bearer $TOKEN" https://api.netlify.com/api/v1/sites | grep -o '"id":"[^"]*' | head -1 | cut -d'"' -f4)

if [[ -z "$SITE_ID" ]]; then
  echo "⚠️ Could not find site ID. Please enter manually."
  read -p "Enter your Netlify site ID: " SITE_ID
fi

echo "🔄 Deploying to Netlify site ID: $SITE_ID"
echo "⏳ This may take a few minutes..."

# Direct API deployment
DEPLOY_RESPONSE=$(curl -s -H "Authorization: Bearer $TOKEN" \
  -F "file=@$DEPLOY_ZIP" \
  -F "function_paths=[]" \
  -F "draft=false" \
  https://api.netlify.com/api/v1/sites/$SITE_ID/deploys)

# Extract deploy URL
DEPLOY_URL=$(echo $DEPLOY_RESPONSE | grep -o '"deploy_ssl_url":"[^"]*' | cut -d'"' -f4)

if [[ -n "$DEPLOY_URL" ]]; then
  echo "✅ Deployment successful!"
  echo "🌐 Deployed to: $DEPLOY_URL"
else
  echo "❌ Deployment failed. Response:"
  echo $DEPLOY_RESPONSE
fi

# Cleanup
rm -f $DEPLOY_ZIP