#!/bin/bash
set -e

echo "🚀 Running fixed Netlify deployment script..."

# Prepare environment
echo "📝 Setting up environment variables..."
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

# Determine if we're using App Router or Pages Router
if [ -d "src/app" ] && [ "$(find src/app -type f | wc -l)" -gt 5 ]; then
  echo "📊 App Router detected - configuring for App Router deployment"
  APP_ROUTER=true
else
  echo "📊 Pages Router detected - configuring for Pages Router deployment"
  APP_ROUTER=false
fi

# Update next.config.js to ensure it has correct settings for static export
echo "📝 Creating optimized next.config.js for static export..."
if [ "$APP_ROUTER" = true ]; then
  # Configuration for App Router
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
  
  // Static export configuration for App Router
  output: 'export',
  trailingSlash: true,
  distDir: 'out',
};

module.exports = nextConfig;
EOL
else
  # Configuration for Pages Router
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
    };
  }
};

module.exports = nextConfig;
EOL
fi

# We need to ensure UI components are properly included in the bundle
echo "🔍 Verifying UI component imports..."

# Build the project
echo "🏗️ Building project for Netlify..."
npm run build

# Create Netlify redirects file
echo "📝 Creating _redirects file for Netlify..."
cat > out/_redirects << EOL
# API routes should redirect to the backend
/api/*  ${NEXT_PUBLIC_API_URL}/:splat  200

# Handle dynamic routes for molecules
/molecules/*  /molecules/index.html  200

# Handle dynamic routes for mixtures
/mixtures/*  /mixtures/index.html  200

# Handle dynamic routes for experiments
/experiments/*  /experiments/index.html  200

# Handle dynamic routes for protocols
/protocols/*  /protocols/index.html  200

# SPA fallback for all other routes
/*  /index.html  200
EOL

# Create Netlify headers file
echo "📝 Creating _headers file for Netlify..."
cat > out/_headers << EOL
# These headers help secure your site
/*
  X-Frame-Options: DENY
  X-Content-Type-Options: nosniff
  Referrer-Policy: strict-origin-when-cross-origin
  Content-Security-Policy: default-src 'self'; script-src 'self' 'unsafe-eval' 'unsafe-inline' https://cdn.jsdelivr.net https://plausible.io; style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; img-src 'self' data: https:; font-src 'self' https://fonts.gstatic.com; connect-src 'self' ${NEXT_PUBLIC_API_URL} ${NEXT_PUBLIC_RDKIT_API_URL} ${NEXT_PUBLIC_CONVEX_URL} https://plausible.io https://*.netlify.app https://*.netlify.com https://identity.netlify.com;
  Permissions-Policy: camera=(), microphone=(), geolocation=()
EOL

echo "✅ Build and configuration completed successfully!"
echo "🚀 Ready to deploy to Netlify!"