#!/bin/bash
# Minimal build script for Netlify
set -ex

# Clean install with yarn which handles conflicts better
echo "Setting up environment..."
npm install -g yarn

# Show versions
echo "Node version: $(node -v)"
echo "Yarn version: $(yarn -v)"

# Install dependencies with yarn
echo "Installing dependencies with yarn..."
yarn install

# Build Next.js as static export
echo "Building Next.js app as static export..."

# Create env file for build
echo "Creating .env.local file..."
cat > .env.local << EOL
NEXT_PUBLIC_API_URL=${NEXT_PUBLIC_API_URL:-https://cryoprotect-8030e4025428.herokuapp.com/v1}
NEXT_PUBLIC_NETLIFY=true
NEXT_PUBLIC_USE_MOCK_DATA=false
NEXT_PUBLIC_ENABLE_API_LOGGING=true
NEXT_PUBLIC_ENVIRONMENT=production
EOL

# Exclude API routes before build
./exclude-api-routes.sh

# Build the app
yarn build

# Create _redirects file for Netlify
echo "Creating Netlify _redirects file..."
cat > out/_redirects << EOL
# Handle API requests to external backend
/api/*  https://cryoprotect-8030e4025428.herokuapp.com/api/:splat  200
# Handle dynamic routes for molecules
/molecules/*  /molecules/[id]/index.html  200
# Handle dynamic routes for mixtures
/mixtures/*  /mixtures/[id]/index.html  200
# SPA fallback
/*    /index.html   200
EOL

# Restore API routes after build
./restore-api-routes.sh

echo "Build and configuration completed successfully!"