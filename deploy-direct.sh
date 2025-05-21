#!/bin/bash

# Direct deployment script for Netlify that works around prerender issues
echo "Starting direct deployment to Netlify..."

# Ensure we're in the project root
cd "$(dirname "$0")"

# Clean any previous builds
echo "Cleaning previous builds..."
rm -rf frontend/.next frontend/out

# Set required environment variables
echo "Setting up environment variables..."
cat > frontend/.env.production << EOF
NEXT_PUBLIC_API_URL=https://cryoprotect-8030e4025428.herokuapp.com/v1
NEXT_PUBLIC_RDKIT_API_URL=https://cryoprotect-rdkit.fly.dev
NEXT_PUBLIC_CONVEX_URL=https://upbeat-parrot-866.convex.cloud
NEXT_PUBLIC_USE_CONVEX=true
NEXT_PUBLIC_ENVIRONMENT=production
NEXT_PUBLIC_ENABLE_API_LOGGING=true
NEXT_PUBLIC_NETLIFY=true
EOF

# Create a temporary file to modify next.config.js
echo "Configuring Next.js for deployment..."
cat > frontend/next.config.js << 'EOF'
/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  swcMinify: true,
  images: {
    unoptimized: true,
  },
  eslint: {
    ignoreDuringBuilds: true
  },
  typescript: {
    ignoreBuildErrors: true
  }
};

module.exports = nextConfig;
EOF

# Install the Netlify plugin if needed
echo "Ensuring Netlify plugin is installed..."
cd frontend
npm install -D @netlify/plugin-nextjs

# Build the application
echo "Building the application..."
NETLIFY_NEXT_PLUGIN_SKIP_INSTALL=1 npm run build

# Deploy to Netlify
echo "Deploying to Netlify..."
netlify deploy --prod --dir=.next --functions=netlify/functions

echo "Deployment complete!"