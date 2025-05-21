#!/bin/bash
# Script to prepare and deploy CryoProtect to Vercel

echo "====================================================="
echo "CryoProtect Vercel Deployment Helper"
echo "====================================================="

# Environment variables
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V" 
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"
FRONTEND_URL="https://cryoprotect.vercel.app"
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

# Step 1: Clean up previous builds
echo "Cleaning up previous builds..."
rm -rf frontend/.next
rm -rf frontend/node_modules/.cache

# Step 2: Make sure .npmrc exists in both root and frontend
echo "Setting up .npmrc files..."
cat > .npmrc << 'EOF'
legacy-peer-deps=true
strict-peer-dependencies=false
fund=false
audit=false
loglevel=error
EOF

cat > frontend/.npmrc << 'EOF'
legacy-peer-deps=true
strict-peer-dependencies=false
fund=false
audit=false
loglevel=error
EOF

# Step 3: Update vercel.json
echo "Updating vercel.json..."
cat > vercel.json << 'EOF'
{
  "version": 2,
  "buildCommand": "npm --prefix frontend run build",
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

# Step 4: Update .vercelignore to minimize files
echo "Updating .vercelignore..."
cat > .vercelignore << 'EOF'
# Exclude everything by default
*

# Include only the frontend directory
!frontend/**

# Exclude node_modules from frontend
frontend/node_modules/**

# Include necessary config files
!package.json
!.npmrc
!frontend/.npmrc
!vercel.json
!.vercelignore
EOF

# Step 5: Create a minimal root package.json
echo "Creating minimal root package.json..."
cat > package.json << 'EOF'
{
  "name": "cryoprotect",
  "version": "1.0.0",
  "private": true,
  "engines": {
    "node": "18.x",
    "npm": "9.x"
  },
  "scripts": {
    "build": "cd frontend && npm run build"
  }
}
EOF

# Step 6: Deploy to Vercel
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
  echo "Try running vercel deploy --debug for more information."
fi