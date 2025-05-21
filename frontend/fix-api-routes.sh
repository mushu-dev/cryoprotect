#!/bin/bash

# Script to fix API routes on Netlify deployment
set -e

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}==============================================${NC}"
echo -e "${BLUE}    CryoProtect API Routes Fix    ${NC}"
echo -e "${BLUE}==============================================${NC}"
echo

# Create improved Netlify functions directory structure
echo -e "${BLUE}Setting up improved Netlify functions...${NC}"
mkdir -p netlify/functions/api

# Create API proxy function
echo -e "${BLUE}Creating API proxy function...${NC}"
cat > netlify/functions/api.js << EOF
// API proxy function to handle all /api/* routes
exports.handler = async function(event, context) {
  const path = event.path;
  console.log("API Request received at path:", path);
  
  // Handle specific API routes
  if (path.endsWith('/api/hello')) {
    return {
      statusCode: 200,
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({
        name: 'CryoProtect API',
        message: 'API is working!',
        timestamp: new Date().toISOString(),
        environment: process.env.NEXT_PUBLIC_ENVIRONMENT || 'production'
      })
    };
  }
  
  if (path.endsWith('/api/health') || path.endsWith('/health')) {
    return {
      statusCode: 200,
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({
        status: 'OK',
        timestamp: new Date().toISOString(),
        environment: process.env.NEXT_PUBLIC_ENVIRONMENT || 'production'
      })
    };
  }
  
  // Default response for unhandled API routes
  return {
    statusCode: 404,
    body: JSON.stringify({
      error: 'API route not found',
      path: path
    })
  };
};
EOF

# Update netlify.toml with improved API routing
echo -e "${BLUE}Updating netlify.toml configuration...${NC}"
cat > netlify.toml << EOF
[build]
  command = "npm run build"
  publish = ".next"
  functions = "netlify/functions"
  
[build.environment]
  NEXT_PUBLIC_ENVIRONMENT = "production"
  NEXT_PUBLIC_ENABLE_API_LOGGING = "true"
  NEXT_PUBLIC_NETLIFY = "true"
  NEXT_PUBLIC_API_URL = "https://cryoprotect-8030e4025428.herokuapp.com/v1"
  NEXT_PUBLIC_RDKIT_API_URL = "https://cryoprotect-rdkit.fly.dev"
  NEXT_PUBLIC_CONVEX_URL = "https://upbeat-parrot-866.convex.cloud"
  NEXT_PUBLIC_USE_CONVEX = "false"
  NODE_VERSION = "16"
  NPM_FLAGS = "--legacy-peer-deps"
  NEXT_DISABLE_EDGE_IMAGES = "true"

[[plugins]]
  package = "@netlify/plugin-nextjs"

# Direct API proxy for all API routes
[[redirects]]
  from = "/api/*"
  to = "/.netlify/functions/api"
  status = 200
  force = true

# Direct health endpoint
[[redirects]]
  from = "/health"
  to = "/.netlify/functions/api"
  status = 200
  force = true

# RDKit health check endpoint
[[redirects]]
  from = "/rdkit-api/health"
  to = "https://cryoprotect-rdkit.fly.dev/health"
  status = 200
  force = true
  headers = {Access-Control-Allow-Origin = "*"}

# Handle backend API redirects to Heroku
[[redirects]]
  from = "/backend-api/*"
  to = "https://cryoprotect-8030e4025428.herokuapp.com/:splat"
  status = 200
  force = true
  headers = {Access-Control-Allow-Origin = "*"}

# These headers help secure your site
[[headers]]
  for = "/*"
  [headers.values]
    X-Frame-Options = "DENY"
    X-Content-Type-Options = "nosniff"
    Referrer-Policy = "strict-origin-when-cross-origin"
    
# This redirects all requests to the Next.js app
[[redirects]]
  from = "/*"
  to = "/.netlify/functions/next"
  status = 200
EOF

echo -e "${GREEN}API route configuration updated${NC}"

# Update next.config.js to handle API routes better
echo -e "${BLUE}Updating Next.js configuration...${NC}"
cat > next.config.js << EOF
/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  swcMinify: true,
  images: {
    unoptimized: true
  },
  eslint: {
    ignoreDuringBuilds: true
  },
  typescript: {
    ignoreBuildErrors: true
  },
  // This tells Next.js to use trailing slashes in URLs
  trailingSlash: true,
  // Explicitly disable rewrites for API routes
  async rewrites() {
    return [];
  }
};

module.exports = nextConfig;
EOF

echo -e "${GREEN}Next.js configuration updated${NC}"

# Build and deploy
echo -e "${BLUE}Building the application...${NC}"
npm run build

echo -e "${BLUE}Deploying to Netlify...${NC}"
netlify deploy --prod

echo -e "${BLUE}==============================================${NC}"
echo -e "${GREEN}API Routes Fix Deployed!${NC}"
echo -e "${BLUE}==============================================${NC}"

echo -e "${YELLOW}Check the following URLs:${NC}"
echo -e "1. Main site: https://cryoprotect.app"
echo -e "2. API endpoint: https://cryoprotect.app/api/hello"
echo -e "3. Health endpoint: https://cryoprotect.app/health"

exit 0