#!/bin/bash
# Script for deploying CryoProtect to Vercel with proper project detection

# Environment variables
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V" 
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"
FRONTEND_URL="https://cryoprotect.vercel.app"
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

echo "====================================================="
echo "CryoProtect Vercel Deployment (Fixed)"
echo "====================================================="

# Create a temporary deployment directory
TEMP_DIR=$(mktemp -d)
echo "Created temporary directory: $TEMP_DIR"

# Copy only what's needed for deployment
echo "Copying frontend files..."
mkdir -p $TEMP_DIR/frontend/src
cp -r /home/mushu/Projects/cryoprotect/frontend/src $TEMP_DIR/frontend/
cp -r /home/mushu/Projects/cryoprotect/frontend/public $TEMP_DIR/frontend/ 2>/dev/null || mkdir $TEMP_DIR/frontend/public
cp /home/mushu/Projects/cryoprotect/frontend/package.json $TEMP_DIR/frontend/
cp /home/mushu/Projects/cryoprotect/frontend/next.config.js $TEMP_DIR/frontend/ 2>/dev/null
cp /home/mushu/Projects/cryoprotect/frontend/tsconfig.json $TEMP_DIR/frontend/ 2>/dev/null
cp /home/mushu/Projects/cryoprotect/frontend/tailwind.config.js $TEMP_DIR/frontend/ 2>/dev/null
cp /home/mushu/Projects/cryoprotect/frontend/postcss.config.js $TEMP_DIR/frontend/ 2>/dev/null

# Create .npmrc file in frontend directory
echo "Creating .npmrc in frontend directory..."
cat > $TEMP_DIR/frontend/.npmrc << 'EOF'
legacy-peer-deps=true
strict-peer-dependencies=false
EOF

# Create package.json in root - this is important for Vercel to recognize Next.js
echo "Creating root package.json..."
cat > $TEMP_DIR/package.json << 'EOF'
{
  "name": "cryoprotect",
  "version": "1.0.0",
  "private": true,
  "dependencies": {
    "next": "14.0.4"
  },
  "devDependencies": {},
  "scripts": {
    "dev": "cd frontend && npm run dev",
    "build": "cd frontend && npm run build",
    "start": "cd frontend && npm start"
  },
  "workspaces": ["frontend"]
}
EOF

# Create vercel.json with specialized config
echo "Creating vercel.json..."
cat > $TEMP_DIR/vercel.json << 'EOF'
{
  "version": 2,
  "buildCommand": "npm --prefix frontend install --legacy-peer-deps && npm --prefix frontend run build",
  "installCommand": "echo 'Skipping install - handled in buildCommand'",
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

# Deploy from the temporary directory
echo "Deploying to Vercel..."
cd $TEMP_DIR
vercel deploy --prod --yes \
  --name cryoprotect \
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
  -e NODE_OPTIONS="--max_old_space_size=3072"

RESULT=$?

# Clean up
echo "Cleaning up temporary directory..."
cd /home/mushu/Projects/cryoprotect
rm -rf $TEMP_DIR

if [ $RESULT -eq 0 ]; then
  echo "✅ Deployment successful!"
else
  echo "❌ Deployment failed."
  echo "Try running with vercel --debug for more information."
fi