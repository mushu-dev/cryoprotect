#!/bin/bash
set -e

echo "🚀 Preparing Netlify deployment package..."

# Create a clean deployment directory
DEPLOY_DIR="netlify-deploy-package"
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

# Create a README for deployment
echo "📝 Creating deployment instructions..."
cat > $DEPLOY_DIR/DEPLOY_INSTRUCTIONS.txt << EOL
# CryoProtect Netlify Deployment Instructions

1. Go to https://app.netlify.com/drop
2. Drag and drop this entire folder (netlify-deploy-package)
3. Wait for deployment to complete
4. After deployment, go to Site settings > Build & deploy > Environment
5. Add the following environment variables:
   - NEXT_PUBLIC_NETLIFY=true
   - NEXT_PUBLIC_API_URL=https://cryoprotect-8030e4025428.herokuapp.com/api
   - NEXTAUTH_URL=https://[your-site-name].netlify.app
   - NEXTAUTH_SECRET=[your-auth-secret]
   - NEXT_PUBLIC_SUPABASE_URL=[your-supabase-url]
   - NEXT_PUBLIC_SUPABASE_ANON_KEY=[your-supabase-anon-key]
6. Go to Site settings > Build & deploy > Continuous deployment
7. Connect your GitHub repository to enable continuous deployment for future updates
EOL

echo "✅ Deployment package prepared successfully in $DEPLOY_DIR."
echo "📤 Instructions for manual drag-and-drop deployment are in $DEPLOY_DIR/DEPLOY_INSTRUCTIONS.txt"
echo "🌐 Go to https://app.netlify.com/drop to deploy this package."