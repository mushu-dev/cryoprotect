#!/bin/bash
set -e

echo "🚀 Running ultra-minimal Netlify deploy script..."

# Create a clean deployment directory
DEPLOY_DIR="netlify-minimal-deploy"
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

echo "✅ Deployment files prepared successfully!"
echo "📂 Your deployment files are in the '$DEPLOY_DIR' directory"

# Try using npx netlify-cli to avoid installation issues
echo "🚀 Deploying to Netlify using npx..."
cd $DEPLOY_DIR
npx netlify-cli deploy --dir=. --prod

echo "✨ Deployment complete!"