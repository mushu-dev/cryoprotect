#!/bin/bash
set -e

echo "🚀 Running fixed Netlify deployment script (Pages Router only)..."

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

# Hide the App Router directory
echo "🔧 Hiding App Router directory..."
if [ -d "src/app" ]; then
  mkdir -p .app-router-backup
  cp -r src/app/* .app-router-backup/
  rm -rf src/app
fi

# Create simplified next.config.js for static export (Pages Router)
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
      '/dashboard': { page: '/dashboard' },
    };
  }
};

module.exports = nextConfig;
EOL

# Fix package.json to use correct build commands
echo "📝 Updating package.json build script..."
sed -i 's/"build": "next build"/"build": "next build"/g' package.json

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

# Create a simple index.html for non-JavaScript fallback
echo "📝 Creating non-JavaScript fallback page..."
cat > out/noscript.html << EOL
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>CryoProtect - JavaScript Required</title>
  <style>
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
      line-height: 1.6;
      color: #333;
      max-width: 650px;
      margin: 0 auto;
      padding: 20px;
    }
    h1 { color: #0070f3; }
    .card {
      border: 1px solid #eaeaea;
      border-radius: 5px;
      padding: 20px;
      margin: 20px 0;
      box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    .btn {
      display: inline-block;
      background-color: #0070f3;
      color: white;
      padding: 8px 16px;
      border-radius: 4px;
      text-decoration: none;
    }
  </style>
</head>
<body>
  <h1>CryoProtect</h1>
  <div class="card">
    <h2>JavaScript Required</h2>
    <p>CryoProtect requires JavaScript to function properly. Please enable JavaScript in your browser settings to access the full application.</p>
    <p>Once JavaScript is enabled, you can:</p>
    <ul>
      <li>Explore cryoprotectant molecules</li>
      <li>Analyze molecular properties</li>
      <li>Manage protocols and experiments</li>
      <li>Access the full suite of features</li>
    </ul>
    <a href="/" class="btn">Reload Page</a>
  </div>
</body>
</html>
EOL

# Restore App Router if it was backed up
echo "🔄 Restoring App Router directory..."
if [ -d ".app-router-backup" ] && [ "$(ls -A .app-router-backup)" ]; then
  mkdir -p src/app
  cp -r .app-router-backup/* src/app/
  rm -rf .app-router-backup
fi

echo "✅ Build and configuration completed successfully!"
echo "🚀 Ready to deploy to Netlify!"
echo ""
echo "To deploy to Netlify, run:"
echo "netlify deploy --prod"