#!/bin/bash
# Deploy directly to Vercel by modifying the frontend package.json

echo "====================================================="
echo "CryoProtect Vercel Deployment - Direct Analytics"
echo "====================================================="

# Environment variables
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V" 
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"
FRONTEND_URL="https://cryoprotect.vercel.app"
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

# Make sure we're in the project root
if [ ! -d "frontend" ]; then
  echo "Error: frontend directory not found. Please run from project root."
  exit 1
fi

# 1. Backup the original package.json
echo "Backing up original package.json..."
cp frontend/package.json frontend/package.json.bak

# 2. Remove the overrides section from package.json to prevent the conflict
echo "Removing overrides from package.json..."
jq 'del(.overrides) | del(.resolutions)' frontend/package.json > frontend/package.json.tmp
mv frontend/package.json.tmp frontend/package.json

# 3. Make sure we have necessary packages in the root package.json
echo "Updating root package.json..."
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
  },
  "dependencies": {
    "@vercel/analytics": "^1.5.0",
    "@vercel/speed-insights": "^1.2.0"
  }
}
EOF

# 4. Update vercel.json to use the correct build approach
echo "Updating vercel.json..."
cat > vercel.json << 'EOF'
{
  "version": 2,
  "installCommand": "echo 'Skipping root install' && cd frontend && npm install --no-fund --no-audit",
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

# 5. Make sure we have .npmrc (but without legacy-peer-deps now)
echo "Creating .npmrc files..."
echo "fund=false" > .npmrc
echo "audit=false" >> .npmrc

echo "fund=false" > frontend/.npmrc
echo "audit=false" >> frontend/.npmrc

# 6. Verify Analytics components are in layout.tsx
if [ -f "frontend/src/app/layout.tsx" ]; then
  echo "Checking for Analytics components in layout..."
  
  if ! grep -q "@vercel/analytics/react" frontend/src/app/layout.tsx || ! grep -q "@vercel/speed-insights/next" frontend/src/app/layout.tsx; then
    echo "Adding Analytics imports to layout.tsx..."
    cp frontend/src/app/layout.tsx frontend/src/app/layout.tsx.bak
    
    if ! grep -q "@vercel/analytics/react" frontend/src/app/layout.tsx; then
      sed -i "1i import { Analytics } from '@vercel/analytics/react'" frontend/src/app/layout.tsx
    fi
    
    if ! grep -q "@vercel/speed-insights/next" frontend/src/app/layout.tsx; then
      sed -i "1i import { SpeedInsights } from '@vercel/speed-insights/next'" frontend/src/app/layout.tsx
    fi
  fi
  
  if ! grep -q "<Analytics" frontend/src/app/layout.tsx || ! grep -q "<SpeedInsights" frontend/src/app/layout.tsx; then
    echo "Adding Analytics components to layout body..."
    
    if ! grep -q "<Analytics" frontend/src/app/layout.tsx; then
      sed -i "s/<\/body>/  <Analytics \/>\n  <\/body>/" frontend/src/app/layout.tsx
    fi
    
    if ! grep -q "<SpeedInsights" frontend/src/app/layout.tsx; then
      sed -i "s/<\/body>/  <SpeedInsights \/>\n  <\/body>/" frontend/src/app/layout.tsx
    fi
  fi
fi

# 7. Update .vercelignore
echo "Updating .vercelignore..."
cat > .vercelignore << 'EOF'
# Directories to exclude
**/node_modules
**/.git
**/.github
**/__pycache__
**/venv
**/alembic
**/api
**/tests
**/migrations
**/cache
**/backups
**/checkpoints
**/chembl
**/reports

# Allow only what's needed
!frontend
!package.json
!.npmrc
!frontend/.npmrc
!vercel.json
!.vercelignore
EOF

# 8. Deploy to Vercel
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
  -e VERCEL_SPEED_INSIGHTS=true

RESULT=$?

# 9. Restore the original package.json
echo "Restoring original package.json..."
if [ -f "frontend/package.json.bak" ]; then
  mv frontend/package.json.bak frontend/package.json
fi

if [ $RESULT -eq 0 ]; then
  echo "✅ Deployment successful!"
  echo "Your site should now be live with Vercel Analytics and Speed Insights."
else
  echo "❌ Deployment failed."
  echo "Please check the error logs above for more details."
fi