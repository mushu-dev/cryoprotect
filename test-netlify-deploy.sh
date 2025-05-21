#\!/bin/bash
# Test script for Netlify deployment
set -e

# Print banner
echo "=========================================="
echo "📦 Testing Netlify Deployment Configuration"
echo "=========================================="

# Check if we're in the right directory
if [ \! -d "frontend" ]; then
  echo "❌ Error: Must run from the project root directory"
  exit 1
fi

# Go to frontend directory
cd frontend

# Validate build configuration
echo "🔍 Validating build configuration..."
if [ -f "./validate-netlify-build.js" ]; then
  node ./validate-netlify-build.js
else
  echo "⚠️ Warning: validation script not found"
fi

# Clean build directory
echo "🧹 Cleaning build directory..."
rm -rf .next out .api-routes-backup

# Create test environment variables
echo "🌍 Creating test environment variables..."
cat > .env.local << EOL
NEXT_PUBLIC_API_URL=https://cryoprotect-8030e4025428.herokuapp.com/v1
NEXT_PUBLIC_NETLIFY=true
NEXT_PUBLIC_USE_MOCK_DATA=false
NEXT_PUBLIC_ENABLE_API_LOGGING=true
NEXT_PUBLIC_ENVIRONMENT=test
EOL

# Exclude API routes
echo "🚫 Excluding API routes..."
./exclude-api-routes.sh

# Test the build
echo "🏗️ Testing build process..."
yarn build

# Create redirects file
echo "🔄 Creating Netlify redirects file..."
mkdir -p out
cat > out/_redirects << EOL
# Handle API requests to external backend
/api/*  https://cryoprotect-8030e4025428.herokuapp.com/api/:splat  200
# Handle dynamic routes for molecules
/molecules/*  /molecules/placeholder/index.html  200
# Handle dynamic routes for mixtures
/mixtures/*  /mixtures/placeholder/index.html  200
# SPA fallback
/*    /index.html   200
EOL

# Restore API routes
echo "♻️ Restoring API routes..."
./restore-api-routes.sh

# Check for dynamic routes
echo "🔍 Checking for generated dynamic routes..."
if [ -d "./out/molecules/962" ]; then
  echo "✅ Pre-rendered molecule route found"
else
  echo "❌ Error: Pre-rendered molecule route not found"
fi

# Check if redirects are in the _redirects file
echo "🔀 Checking redirects..."
if [ -f "./out/_redirects" ]; then
  echo "✅ _redirects file exists"
  cat ./out/_redirects
else
  echo "⚠️ Warning: _redirects file not found"
fi

# Check the analytics implementation
echo "📊 Checking analytics implementation..."
if grep -q "trackPageView" ./out/main*.js; then
  echo "✅ Analytics tracking found in main.js"
else
  echo "⚠️ Warning: Analytics tracking not found in compiled JS"
fi

# Final message
echo ""
echo "=========================================="
echo "✅ Test build completed\!"
echo "👉 Next step: Deploy to Netlify using netlify deploy --prod"
echo "=========================================="

# Return to original directory
cd ..
