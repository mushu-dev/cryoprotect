#!/bin/bash
# Final version of deployment script that fixes React version conflicts

# Environment variables
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V" 
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"
FRONTEND_URL="https://cryoprotect.vercel.app"
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

echo "====================================================="
echo "CryoProtect Final Vercel Deployment Fix"
echo "====================================================="

# Make sure we're in the right directory
if [ ! -d "frontend" ]; then
  echo "Error: frontend directory not found. Please run from project root."
  exit 1
fi

# Backup original package.json
echo "Backing up original package.json..."
cp frontend/package.json frontend/package.json.original

# Replace with version without overrides
echo "Replacing package.json with version that fixes React conflicts..."
cp frontend/package.json.nooverrides frontend/package.json

# Update .npmrc
echo "Setting up .npmrc..."
echo "legacy-peer-deps=true" > frontend/.npmrc

# Update vercel.json
echo "Updating vercel.json..."
cat > vercel.json << 'EOF'
{
  "version": 2,
  "buildCommand": "cd frontend && npm install --legacy-peer-deps && npm run build",
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

# Clean up any previous builds
echo "Cleaning previous builds..."
rm -rf frontend/.next
rm -rf frontend/node_modules/.cache

# Deploy to Vercel
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

# Restore original package.json
echo "Restoring original package.json..."
mv frontend/package.json.original frontend/package.json

if [ $RESULT -eq 0 ]; then
  echo "✅ Deployment successful!"
  echo "Your site should now be live with Vercel Analytics and Speed Insights."
else
  echo "❌ Deployment failed."
  echo "Please check the error logs above for more details."
fi