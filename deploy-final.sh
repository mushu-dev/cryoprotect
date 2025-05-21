#!/bin/bash
# Deploy to Vercel with temporary package.json modifications

# Environment variables
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V" 
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"
FRONTEND_URL="https://cryoprotect.vercel.app"
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

echo "====================================================="
echo "CryoProtect Vercel Deployment (Final Solution)"
echo "====================================================="

# Make sure we're in the right directory
if [ ! -d "frontend" ]; then
  echo "Error: frontend directory not found. Please run from project root."
  exit 1
fi

# Backup package.json
echo "Backing up original package.json..."
cp frontend/package.json frontend/package.json.original

# Create a clean package.json without overrides
echo "Creating modified package.json without overrides..."
cat frontend/package.json | jq 'del(.overrides) | del(.resolutions)' > frontend/package.json.tmp
mv frontend/package.json.tmp frontend/package.json

# Remove existing vercel.json files
echo "Temporarily moving vercel.json files..."
if [ -f "vercel.json" ]; then
  mv vercel.json vercel.json.bak
fi

if [ -f "frontend/vercel.json" ]; then
  mv frontend/vercel.json frontend/vercel.json.bak
fi

# Set up .npmrc
echo "legacy-peer-deps=true" > frontend/.npmrc

# Clean previous builds
echo "Cleaning previous builds..."
rm -rf frontend/.next
rm -rf frontend/node_modules/.cache

# Create a package-lock.json file if it doesn't exist
if [ ! -f "frontend/package-lock.json" ]; then
  echo "Creating minimal package-lock.json..."
  echo '{"name":"cryoprotect-frontend","lockfileVersion":2,"requires":true,"packages":{}}' > frontend/package-lock.json
fi

# Deploy to Vercel
echo "Deploying to Vercel..."
(cd frontend && 
 FORCE_COLOR=1 vercel deploy --prod --yes \
  -e NEXT_PUBLIC_API_URL=$API_URL \
  -e NEXT_PUBLIC_USE_MOCK_DATA=false \
  -e NEXT_PUBLIC_ENABLE_API_LOGGING=true \
  -e NEXT_PUBLIC_ENVIRONMENT=production \
  -e NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  -e NEXTAUTH_URL=$FRONTEND_URL \
  -e NEXTAUTH_SECRET="$NEXTAUTH_SECRET" \
  -e PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  -e VERCEL_ANALYTICS_ID=true \
  -e VERCEL_SPEED_INSIGHTS=true)

RESULT=$?

# Restore original files
echo "Restoring original files..."
if [ -f "frontend/package.json.original" ]; then
  mv frontend/package.json.original frontend/package.json
fi

if [ -f "vercel.json.bak" ]; then
  mv vercel.json.bak vercel.json
fi

if [ -f "frontend/vercel.json.bak" ]; then
  mv frontend/vercel.json.bak frontend/vercel.json
fi

if [ $RESULT -eq 0 ]; then
  echo "✅ Deployment successful!"
  echo "Your site should now be live with Vercel Analytics and Speed Insights."
else
  echo "❌ Deployment failed."
  echo "Please check the error logs above for more details."
fi