#!/bin/bash
set -e

echo "Running Netlify build script..."

# Create environment variables file
echo "Creating .env.local file for Netlify..."
cat > .env.local << EOL
NEXT_PUBLIC_API_URL=${NEXT_PUBLIC_API_URL}
NEXT_PUBLIC_USE_MOCK_DATA=${NEXT_PUBLIC_USE_MOCK_DATA:-false}
NEXT_PUBLIC_ENABLE_API_LOGGING=${NEXT_PUBLIC_ENABLE_API_LOGGING:-true}
NEXT_PUBLIC_ENVIRONMENT=${NEXT_PUBLIC_ENVIRONMENT:-production}
NEXT_PUBLIC_NETLIFY=${NEXT_PUBLIC_NETLIFY:-true}
NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS=${NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS}
NEXTAUTH_URL=${NEXTAUTH_URL}
NEXTAUTH_SECRET=${NEXTAUTH_SECRET}
PROTECTION_BYPASS=${PROTECTION_BYPASS}
EOL

# Create simplified next.config.js for static export
echo "Creating simplified next.config.js for static export..."
cat > next.config.js << EOL
/** @type {import('next').NextConfig} */
const nextConfig = {
  output: 'export',
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
EOL

# Move API routes out of the way for static build
echo "Temporarily excluding API routes for static build..."
if [ -d "src/app/api" ]; then
  mkdir -p .api-routes-backup
  mv src/app/api .api-routes-backup/
  echo "API routes moved to backup"
fi

# Build the Next.js app
echo "Building Next.js app..."
npx next build

# Create Netlify redirects file
echo "Creating _redirects file for Netlify..."
cat > out/_redirects << EOL
# API routes should redirect to the backend
/api/*  https://cryoprotect-8030e4025428.herokuapp.com/api/:splat  200

# Handle dynamic routes 
/molecules/*  /molecules/[id].html  200
/mixtures/*  /mixtures/[id].html  200

# SPA fallback
/*  /index.html  200
EOL

# Restore API routes
echo "Restoring API routes..."
if [ -d ".api-routes-backup/api" ]; then
  mkdir -p src/app
  mv .api-routes-backup/api src/app/
  rm -rf .api-routes-backup
fi

echo "Build completed successfully!"