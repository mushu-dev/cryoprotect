#!/bin/bash
# Deploy to Vercel with React dependency conflict fix

# Environment variables
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V" 
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"
FRONTEND_URL="https://cryoprotect.vercel.app"
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

echo "====================================================="
echo "CryoProtect Vercel Deployment (Fixed for React Conflict)"
echo "====================================================="

# Make sure we're in the right directory
if [ ! -d "frontend" ]; then
  echo "Error: frontend directory not found. Please run from project root."
  exit 1
fi

# Make sure the .npmrc file is correct
echo "Setting up .npmrc for legacy peer dependencies..."
echo "fund=false" > frontend/.npmrc
echo "audit=false" >> frontend/.npmrc
echo "legacy-peer-deps=true" >> frontend/.npmrc

# Create a vercel.json that explicitly uses legacy-peer-deps
echo "Updating vercel.json..."
cat > vercel.json << 'EOF'
{
  "version": 2,
  "installCommand": "echo 'Skipping root install' && cd frontend && npm install --legacy-peer-deps",
  "buildCommand": "cd frontend && npm run build",
  "outputDirectory": "frontend/.next",
  "framework": "nextjs",
  "git": {
    "deploymentEnabled": {
      "main": true
    }
  },
  "github": {
    "silent": true
  },
  "cleanUrls": true,
  "regions": ["iad1"]
}
EOF

# Clean up any caches/previous builds
echo "Cleaning up previous builds..."
rm -rf frontend/.next
rm -rf frontend/node_modules/.cache

# Deploy to Vercel with proper environment variables
echo "Deploying to Vercel..."
vercel deploy --prod --yes \
  --archive=tgz \
  -e NEXT_PUBLIC_API_URL=$API_URL \
  -e NEXT_PUBLIC_USE_MOCK_DATA=false \
  -e NEXT_PUBLIC_ENABLE_API_LOGGING=true \
  -e NEXT_PUBLIC_ENVIRONMENT=production \
  -e NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  -e NEXTAUTH_URL=$FRONTEND_URL \
  -e NEXTAUTH_SECRET="$NEXTAUTH_SECRET" \
  -e PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  -e VERCEL_ANALYTICS_ID=true \
  -e VERCEL_SPEED_INSIGHTS=true \
  -e NPM_FLAGS="--legacy-peer-deps"

if [ $? -eq 0 ]; then
  echo "✅ Deployment successful!"
else
  echo "❌ Deployment failed."
  echo "Try with --debug flag for more information."
fi