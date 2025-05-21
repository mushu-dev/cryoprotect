#!/bin/bash
set -e

echo "🚀 Running simplified Netlify deploy script..."

# Create a clean deployment directory
mkdir -p netlify-deploy
rm -rf netlify-deploy/*

# Copy the existing static output files
echo "📦 Copying existing static files..."
cp -r frontend/out/* netlify-deploy/

# Create Netlify configuration file
echo "📝 Creating netlify.toml file..."
cat > netlify-deploy/netlify.toml << EOL
[build]
  publish = "/"

# Explicitly set Node.js version
[build.environment]
  NODE_VERSION = "18"

# Security headers for all pages
[[headers]]
  for = "/*"
  [headers.values]
    Content-Security-Policy = "default-src 'self'; script-src 'self' 'unsafe-eval' 'unsafe-inline' https://cdn.jsdelivr.net https://plausible.io https://*.netlify.app https://*.netlify.com; style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; img-src 'self' data: https:; font-src 'self' https://fonts.gstatic.com; connect-src 'self' https://api.cryoprotect.app https://*.supabase.co https://plausible.io https://*.netlify.app https://*.netlify.com;"
    X-Frame-Options = "DENY"
    X-Content-Type-Options = "nosniff"
    Referrer-Policy = "strict-origin-when-cross-origin"
    Permissions-Policy = "camera=(), microphone=(), geolocation=()"
EOL

# Create _redirects file
echo "📝 Creating _redirects file..."
cat > netlify-deploy/_redirects << EOL
# API routes redirect to Heroku backend
/api/*  https://cryoprotect-8030e4025428.herokuapp.com/api/:splat  200

# API connectivity endpoint - CORS enabled
/api/v1/health/connectivity  https://cryoprotect-8030e4025428.herokuapp.com/api/v1/health/connectivity  200!
/api/*  https://cryoprotect-8030e4025428.herokuapp.com/api/:splat  200!

# Handle dynamic routes for molecules
/molecules/*  /molecules/[id].html  200

# Handle dynamic routes for mixtures
/mixtures/*  /mixtures/[id].html  200

# SPA fallback for all other routes
/*  /index.html  200
/*  /404.html  404
EOL

echo "✅ Deployment files prepared successfully!"
echo "📂 Your deployment files are in the 'netlify-deploy' directory"
echo "🌐 You can now deploy these files to Netlify"

# Deploy using netlify-cli if installed
if command -v netlify &> /dev/null; then
  echo "🚀 Deploying to Netlify..."
  cd netlify-deploy
  netlify deploy --dir=. --prod
else
  echo "ℹ️ To deploy to Netlify, run:"
  echo "cd netlify-deploy && netlify deploy --dir=. --prod"
fi