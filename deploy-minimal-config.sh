#!/bin/bash
# Deploy with minimal config, avoiding configuration conflicts

# Environment variables
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V" 
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"
FRONTEND_URL="https://cryoprotect.vercel.app"
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

echo "====================================================="
echo "CryoProtect Vercel Deployment (Minimal Config)"
echo "====================================================="

# Make sure we're in the right directory
if [ ! -d "frontend" ]; then
  echo "Error: frontend directory not found. Please run from project root."
  exit 1
fi

# Rename both vercel.json files temporarily to avoid conflicts
echo "Temporarily renaming vercel.json files to avoid conflicts..."
if [ -f "vercel.json" ]; then
  mv vercel.json vercel.json.bak
fi

if [ -f "frontend/vercel.json" ]; then
  mv frontend/vercel.json frontend/vercel.json.bak
fi

# Make sure .npmrc is set up
echo "Configuring .npmrc..."
echo "legacy-peer-deps=true" > frontend/.npmrc

# Clean previous builds
echo "Cleaning previous builds..."
rm -rf frontend/.next
rm -rf frontend/node_modules/.cache

# Deploy with minimal configuration via command line
echo "Deploying to Vercel..."
vercel deploy --prod --yes \
  --cwd frontend \
  --build-env NEXT_PUBLIC_API_URL=$API_URL \
  --build-env NEXT_PUBLIC_USE_MOCK_DATA=false \
  --build-env NEXT_PUBLIC_ENABLE_API_LOGGING=true \
  --build-env NEXT_PUBLIC_ENVIRONMENT=production \
  --build-env NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  --build-env NEXTAUTH_URL=$FRONTEND_URL \
  --build-env NEXTAUTH_SECRET="$NEXTAUTH_SECRET" \
  --build-env PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  --build-env VERCEL_ANALYTICS_ID=true \
  --build-env VERCEL_SPEED_INSIGHTS=true \
  -e NEXT_PUBLIC_API_URL=$API_URL \
  -e NEXT_PUBLIC_USE_MOCK_DATA=false \
  -e NEXT_PUBLIC_ENABLE_API_LOGGING=true \
  -e NEXT_PUBLIC_ENVIRONMENT=production \
  -e NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  -e NEXTAUTH_URL=$FRONTEND_URL \
  -e NEXTAUTH_SECRET="$NEXTAUTH_SECRET" \
  -e PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  -e VERCEL_ANALYTICS_ID=true \
  -e VERCEL_SPEED_INSIGHTS=true

RESULT=$?

# Restore the original vercel.json files
echo "Restoring original vercel.json files..."
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