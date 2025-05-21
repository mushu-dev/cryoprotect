#!/bin/bash
set -e

echo "🚀 Running fixed Netlify deployment script using Pages Router..."

# Set up environment variables
NEXT_PUBLIC_API_URL=${NEXT_PUBLIC_API_URL:-https://cryoprotect-8030e4025428.herokuapp.com/v1}
NEXT_PUBLIC_RDKIT_API_URL=${NEXT_PUBLIC_RDKIT_API_URL:-https://cryoprotect-rdkit.fly.dev}
NEXT_PUBLIC_CONVEX_URL=${NEXT_PUBLIC_CONVEX_URL:-https://upbeat-parrot-866.convex.cloud}
NEXT_PUBLIC_USE_CONVEX=${NEXT_PUBLIC_USE_CONVEX:-true}
NEXT_PUBLIC_ENVIRONMENT=${NEXT_PUBLIC_ENVIRONMENT:-production}
NEXT_PUBLIC_ENABLE_API_LOGGING=${NEXT_PUBLIC_ENABLE_API_LOGGING:-true}
NEXT_PUBLIC_NETLIFY=${NEXT_PUBLIC_NETLIFY:-true}

# Create environment file
echo "📝 Creating .env.local file..."
cat > .env.local << EOL
NEXT_PUBLIC_API_URL=${NEXT_PUBLIC_API_URL}
NEXT_PUBLIC_RDKIT_API_URL=${NEXT_PUBLIC_RDKIT_API_URL}
NEXT_PUBLIC_CONVEX_URL=${NEXT_PUBLIC_CONVEX_URL}
NEXT_PUBLIC_USE_CONVEX=${NEXT_PUBLIC_USE_CONVEX}
NEXT_PUBLIC_ENVIRONMENT=${NEXT_PUBLIC_ENVIRONMENT}
NEXT_PUBLIC_ENABLE_API_LOGGING=${NEXT_PUBLIC_ENABLE_API_LOGGING}
NEXT_PUBLIC_NETLIFY=${NEXT_PUBLIC_NETLIFY}
EOL

# Install dependencies
echo "📦 Installing dependencies..."
npm install

# Install TypeScript babel dependencies
echo "📦 Installing TypeScript babel dependencies..."
npm install --save-dev @babel/preset-typescript @babel/plugin-transform-typescript

# Configure to use Pages Router only
echo "🔧 Configuring for Pages Router only..."
./use-pages-router.sh

# Disable TypeScript components that cause issues
echo "🔧 Disabling TypeScript components for compatibility..."
./disable-typescript-components.sh

# Update navigation component for build
echo "🔧 Updating navigation component for build..."
./update-navigation-for-build.sh

# Create simplified next.config.js for static export with Pages Router
echo "📝 Creating optimized next.config.js for Pages Router static export..."
cat > next.config.js << EOL
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
  },
  
  // React 18 specific settings
  compiler: {
    styledComponents: true
  },
  
  // Static export configuration for Pages Router
  output: 'export',
  trailingSlash: true,
  distDir: 'out',
  
  // Explicitly define routes for static export (Pages Router only)
  exportPathMap: async function() {
    return {
      '/': { page: '/' },
      '/molecules': { page: '/molecules' },
      '/mixtures': { page: '/mixtures' },
      '/experiments': { page: '/experiments' },
      '/protocols': { page: '/protocols' },
      '/properties': { page: '/properties' },
      '/protocols/create': { page: '/protocols/create' },
      '/protocols/[id]': { page: '/protocols/[id]' },
      '/dashboard': { page: '/dashboard' },
      '/components': { page: '/components' },
    };
  }
};

module.exports = nextConfig;
EOL

# Build the project
echo "🏗️ Building project for Netlify..."
npm run build

# Create Netlify redirects file
echo "📝 Creating _redirects file for Netlify..."
cat > out/_redirects << EOL
# API routes should redirect to the backend
/api/*  ${NEXT_PUBLIC_API_URL}/:splat  200

# Handle dynamic routes 
/molecules/*  /molecules/index.html  200
/mixtures/*  /mixtures/index.html  200
/experiments/*  /experiments/index.html  200
/protocols/*  /protocols/index.html  200

# SPA fallback for all other routes
/*  /index.html  200
EOL

# Restore TypeScript components
echo "🔄 Restoring TypeScript components..."
./restore-typescript-components.sh

# Restore original navigation component
echo "🔄 Restoring original navigation component..."
./restore-navigation.sh

# Restore App Router for continued development
echo "🔄 Restoring App Router for continued development..."
./restore-app-router.sh

echo "✅ Build and configuration completed successfully!"
echo "🚀 Ready to deploy to Netlify!"
echo ""
echo "To deploy to Netlify, run:"
echo "netlify deploy --prod"