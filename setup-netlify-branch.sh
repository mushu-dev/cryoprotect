#!/bin/bash
set -e

echo "🚀 Setting up Netlify deployment branch..."

# Create a new branch specifically for Netlify deployment
git checkout -b netlify-deployment

# Create simplified netlify.toml
cat > netlify.toml << EOL
[build]
  publish = "out"
  command = "echo 'Using pre-built files, no build needed'"

[[redirects]]
  from = "/api/*"
  to = "https://cryoprotect-8030e4025428.herokuapp.com/api/:splat"
  status = 200

[[redirects]]
  from = "/molecules/*"
  to = "/molecules/[id].html"
  status = 200

[[redirects]]
  from = "/mixtures/*"
  to = "/mixtures/[id].html"
  status = 200

[[redirects]]
  from = "/*"
  to = "/index.html"
  status = 200

[[headers]]
  for = "/*"
  [headers.values]
    Content-Security-Policy = "default-src 'self'; script-src 'self' 'unsafe-eval' 'unsafe-inline' https://cdn.jsdelivr.net https://plausible.io https://*.netlify.app https://*.netlify.com; style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; img-src 'self' data: https:; font-src 'self' https://fonts.gstatic.com; connect-src 'self' https://api.cryoprotect.app https://*.supabase.co https://plausible.io https://*.netlify.app https://*.netlify.com;"
    X-Frame-Options = "DENY"
    X-Content-Type-Options = "nosniff"
    Referrer-Policy = "strict-origin-when-cross-origin"
    Permissions-Policy = "camera=(), microphone=(), geolocation=()"
EOL

# Create a README specifically for Netlify deployment
cat > NETLIFY_DEPLOYMENT_README.md << EOL
# CryoProtect Netlify Deployment

This branch contains pre-built static files for direct deployment to Netlify.

## Features
- Uses pre-built static files in the /out directory
- Configured with proper redirects for dynamic routes
- Includes security headers
- API routes forwarded to the Heroku backend

## How to deploy

Connect this branch to Netlify and it will automatically deploy the static site.

1. In the Netlify UI, go to "Site settings" > "Build & deploy" > "Continuous deployment"
2. Connect to your GitHub repository and select this branch (netlify-deployment)
3. The deployment should happen automatically without running a build command

## Environment Variables

Make sure to set these environment variables in the Netlify UI:

- NEXT_PUBLIC_NETLIFY=true
- NEXT_PUBLIC_API_URL=https://cryoprotect-8030e4025428.herokuapp.com/api
- NEXTAUTH_URL=https://[your-site-name].netlify.app
- NEXTAUTH_SECRET=[your-auth-secret]
- NEXT_PUBLIC_SUPABASE_URL=[your-supabase-url]
- NEXT_PUBLIC_SUPABASE_ANON_KEY=[your-supabase-anon-key]
EOL

# Commit the changes
git add netlify.toml NETLIFY_DEPLOYMENT_README.md
git commit -m "Setup Netlify deployment configuration

This commit adds a simplified Netlify configuration that:
- Uses pre-built static files
- Configures redirects for dynamic routes
- Sets up security headers
- Forwards API routes to Heroku backend

🤖 Generated with [Claude Code](https://claude.ai/code)

Co-Authored-By: Claude <noreply@anthropic.com>"

echo "✅ Netlify deployment branch created successfully!"
echo "🌐 Next steps:"
echo "1. Push this branch to GitHub with: git push origin netlify-deployment"
echo "2. In the Netlify UI, connect to this branch for deployment"
echo "3. Set up the environment variables listed in NETLIFY_DEPLOYMENT_README.md"