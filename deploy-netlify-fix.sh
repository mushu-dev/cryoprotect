#!/bin/bash

# Fix Netlify Deployment Script for Convex Integration
# This script fixes the Netlify deployment to support Next.js with Convex

set -e

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}==============================================${NC}"
echo -e "${BLUE}    CryoProtect Netlify Deployment Fix    ${NC}"
echo -e "${BLUE}==============================================${NC}"
echo

# Check if Netlify CLI is installed
if ! command -v netlify &> /dev/null; then
    echo -e "${YELLOW}Netlify CLI not found. Installing...${NC}"
    npm install -g netlify-cli
fi

# Create a temporary directory
DEPLOY_TEMP=$(mktemp -d)
echo -e "${BLUE}Created temporary directory: ${DEPLOY_TEMP}${NC}"

# Function to clean up on exit
cleanup() {
    echo -e "${BLUE}Cleaning up temporary files...${NC}"
    rm -rf "$DEPLOY_TEMP"
    echo -e "${BLUE}Cleanup complete.${NC}"
}

# Register cleanup function to run on exit
trap cleanup EXIT

# Step 1: Update netlify.toml with correct settings
echo -e "${BLUE}Updating netlify.toml configuration...${NC}"

# Use absolute paths to ensure files are found
PROJECT_ROOT="/home/mushu/Projects/cryoprotect"
FRONTEND_DIR="${PROJECT_ROOT}/frontend"

cat > "${FRONTEND_DIR}/netlify.toml" << EOF
[build]
  command = "npm run build"
  publish = ".next"
  
[build.environment]
  NEXT_PUBLIC_ENVIRONMENT = "production"
  NEXT_PUBLIC_ENABLE_API_LOGGING = "true"
  NEXT_PUBLIC_NETLIFY = "true"
  NEXT_PUBLIC_API_URL = "https://cryoprotect-8030e4025428.herokuapp.com/v1"
  NEXT_PUBLIC_RDKIT_API_URL = "https://cryoprotect-rdkit.fly.dev"
  NEXT_PUBLIC_CONVEX_URL = "https://upbeat-parrot-866.convex.cloud"
  NEXT_PUBLIC_USE_CONVEX = "true"
  # Other environment variables will be set through the Netlify UI or CLI

# API health check endpoint
[[redirects]]
  from = "/api/health"
  to = "https://cryoprotect-8030e4025428.herokuapp.com/health"
  status = 200
  force = true
  headers = {Access-Control-Allow-Origin = "*"}

# API CORS test endpoint
[[redirects]]
  from = "/api/test-cors"
  to = "https://cryoprotect-8030e4025428.herokuapp.com/test-cors"
  status = 200
  force = true
  headers = {Access-Control-Allow-Origin = "*"}

# RDKit health check endpoint
[[redirects]]
  from = "/rdkit-api/health"
  to = "https://cryoprotect-rdkit.fly.dev/health"
  status = 200
  force = true
  headers = {Access-Control-Allow-Origin = "*"}

# RDKit CORS test endpoint
[[redirects]]
  from = "/rdkit-api/test-cors"
  to = "https://cryoprotect-rdkit.fly.dev/test-cors"
  status = 200
  force = true
  headers = {Access-Control-Allow-Origin = "*"}

# Handle API redirects to Heroku backend
[[redirects]]
  from = "/api/*"
  to = "https://cryoprotect-8030e4025428.herokuapp.com/:splat"
  status = 200
  force = true
  headers = {Access-Control-Allow-Origin = "*"}

# These headers help secure your site
[[headers]]
  for = "/*"
  [headers.values]
    Content-Security-Policy = "default-src 'self'; script-src 'self' 'unsafe-eval' 'unsafe-inline' https://cdn.jsdelivr.net https://plausible.io https://*.netlify.app https://*.netlify.com; style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; img-src 'self' data: https:; font-src 'self' https://fonts.gstatic.com; connect-src 'self' https://*.cryoprotect.app https://cryoprotect-8030e4025428.herokuapp.com https://cryoprotect-rdkit.fly.dev https://*.convex.cloud https://upbeat-parrot-866.convex.cloud https://*.supabase.co https://plausible.io https://*.netlify.app https://*.netlify.com;"
    X-Frame-Options = "DENY"
    X-Content-Type-Options = "nosniff"
    Referrer-Policy = "strict-origin-when-cross-origin"
    Permissions-Policy = "camera=(), microphone=(), geolocation=()"

# This redirects all requests to the Next.js app
[[redirects]]
  from = "/*"
  to = "/.netlify/functions/next"
  status = 200
EOF

echo -e "${GREEN}Updated netlify.toml with correct configuration${NC}"

# Step 2: Update next.config.js to work with Netlify
echo -e "${BLUE}Updating next.config.js...${NC}"

cat > "${FRONTEND_DIR}/next.config.js" << EOF
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

echo -e "${GREEN}Updated next.config.js - removed 'output: export' to enable SSR${NC}"

# Step 3: Update package.json build scripts
echo -e "${BLUE}Updating package.json build scripts...${NC}"

# Save current directory
CURRENT_DIR=$(pwd)

# Change to frontend directory
cd "${FRONTEND_DIR}"

# Use a temporary file to update package.json
jq '.scripts.build = "next build"' package.json > "$DEPLOY_TEMP/package.json.tmp"
mv "$DEPLOY_TEMP/package.json.tmp" package.json

echo -e "${GREEN}Updated package.json build script${NC}"

# Step 4: Install the Netlify Next.js plugin
echo -e "${BLUE}Installing Netlify Next.js plugin...${NC}"
npm install -D @netlify/plugin-nextjs

# Add the plugin to package.json
jq '.plugins = [{"package": "@netlify/plugin-nextjs"}]' package.json > "$DEPLOY_TEMP/package.json.tmp"
mv "$DEPLOY_TEMP/package.json.tmp" package.json

echo -e "${GREEN}Installed Netlify Next.js plugin${NC}"

# Step 5: Deploy to Netlify
echo -e "${BLUE}Deploying to Netlify...${NC}"

# Check if user is logged in to Netlify
if ! netlify status 2>&1 | grep -q "Logged in"; then
    echo -e "${YELLOW}Not logged in to Netlify. Please log in:${NC}"
    netlify login
fi

# Check if site is linked
if ! netlify status 2>&1 | grep -q "cryoprotect"; then
    echo -e "${YELLOW}Site not linked. Linking to cryoprotect site...${NC}"
    netlify unlink 2>/dev/null
    netlify link --name cryoprotect
fi

# Deploy to Netlify
echo -e "${BLUE}Starting Netlify deployment...${NC}"
netlify deploy --build --prod

# Go back to original directory
cd "$CURRENT_DIR"

echo -e "${BLUE}==============================================${NC}"
echo -e "${GREEN}Netlify deployment completed successfully!${NC}"
echo -e "${BLUE}==============================================${NC}"

echo -e "${YELLOW}Note: It may take a few minutes for the changes to propagate.${NC}"
echo -e "${YELLOW}Visit your Netlify site to verify that the deployment was successful.${NC}"

exit 0