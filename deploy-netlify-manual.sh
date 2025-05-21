#!/bin/bash
set -e

echo "🚀 Running manual Netlify deploy script..."

# Create a clean deployment directory
mkdir -p netlify-manual-deploy
rm -rf netlify-manual-deploy/*

# Copy the existing static output files
echo "📦 Copying existing static files..."
cp -r frontend/out/* netlify-manual-deploy/

# Create _redirects file
echo "📝 Creating _redirects file..."
cat > netlify-manual-deploy/_redirects << EOL
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

# Create netlify.json config for direct deployment
echo "📝 Creating netlify.json for manual deployment..."
cat > netlify-manual-deploy/netlify.json << EOL
{
  "headers": [
    {
      "for": "/*",
      "values": {
        "Content-Security-Policy": "default-src 'self'; script-src 'self' 'unsafe-eval' 'unsafe-inline' https://cdn.jsdelivr.net https://plausible.io https://*.netlify.app https://*.netlify.com; style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; img-src 'self' data: https:; font-src 'self' https://fonts.gstatic.com; connect-src 'self' https://api.cryoprotect.app https://*.supabase.co https://plausible.io https://*.netlify.app https://*.netlify.com;",
        "X-Frame-Options": "DENY",
        "X-Content-Type-Options": "nosniff",
        "Referrer-Policy": "strict-origin-when-cross-origin",
        "Permissions-Policy": "camera=(), microphone=(), geolocation=()"
      }
    }
  ]
}
EOL

echo "✅ Manual deployment files prepared successfully!"
echo "📂 Your deployment files are in the 'netlify-manual-deploy' directory"

# Install netlify-cli if not present
if ! command -v netlify &> /dev/null; then
    echo "🔧 Installing netlify-cli..."
    npm install -g netlify-cli
fi

echo "🔑 Authenticating with Netlify (you may need to log in)..."
netlify status || netlify login

echo "🚀 Deploying to Netlify manually without build step..."
echo "This will bypass any build settings in the Netlify UI"

# Deploy site directly
cd netlify-manual-deploy
netlify deploy --dir=. --prod --message "Manual deployment with static files"

echo "✨ Deployment complete!"