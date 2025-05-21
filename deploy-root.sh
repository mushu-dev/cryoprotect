#!/bin/bash
# Deploy directly to Vercel specifying the root directory via command line

# Environment variables
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V" 
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"
FRONTEND_URL="https://cryoprotect.vercel.app"
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

echo "====================================================="
echo "CryoProtect Vercel Deployment (Root Directory via CLI)"
echo "====================================================="

# Make sure we're in the right directory
if [ ! -d "frontend" ]; then
  echo "Error: frontend directory not found. Please run from project root."
  exit 1
fi

# Update vercel.json
echo "Updating vercel.json..."
cat > vercel.json << 'EOF'
{
  "version": 2,
  "buildCommand": "npm install --legacy-peer-deps && npm run build",
  "outputDirectory": ".next",
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

# Set up .npmrc
echo "Setting up .npmrc..."
echo "legacy-peer-deps=true" > frontend/.npmrc

# Clean up any previous builds
echo "Cleaning previous builds..."
rm -rf frontend/.next
rm -rf frontend/node_modules/.cache

# Deploy to Vercel with root directory specified
echo "Deploying to Vercel (specifying root directory)..."
vercel deploy --prod --yes \
  --cwd frontend \
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

if [ $? -eq 0 ]; then
  echo "✅ Deployment successful!"
  echo "Your site should now be live with Vercel Analytics and Speed Insights."
else
  echo "❌ Deployment failed."
  echo "Please check the error logs above for more details."
fi