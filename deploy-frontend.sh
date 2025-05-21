#!/bin/bash
# Master deployment script that fixes dependency conflicts and deploys to Vercel

# Environment variables
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V"
FRONTEND_URL="https://cryoprotect.vercel.app"
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

echo "==========================================="
echo "CryoProtect Frontend Deployment to Vercel"
echo "==========================================="

# Ensure .vercelignore is set up to exclude non-frontend files
if [ ! -f ".vercelignore" ]; then
  echo "Creating .vercelignore file to minimize deployment size..."
  cat > .vercelignore << 'EOF'
# Exclude everything by default
*

# Include only the frontend directory
!frontend/**
frontend/node_modules
frontend/.next

# Include necessary config files
!vercel.json
!.vercelignore

# Allow only the minimum needed to make Vercel deploy work
!package.json
!README.md
EOF
fi

# Create package.json in the root to help Vercel identify project type
if [ ! -f "package.json" ]; then
  echo "Creating minimal package.json in project root..."
  cat > package.json << 'EOF'
{
  "name": "cryoprotect",
  "private": true,
  "version": "0.1.0",
  "workspaces": ["frontend"],
  "scripts": {
    "build": "cd frontend && npm run build",
    "dev": "cd frontend && npm run dev",
    "start": "cd frontend && npm run start",
    "deploy": "vercel --prod"
  }
}
EOF
fi

# Check if frontend directory exists
if [ ! -d "frontend" ]; then
  echo "Error: frontend directory not found. Please run this script from the project root."
  exit 1
fi

# Check analytics packages are installed
cd frontend

if ! grep -q "@vercel/analytics" package.json || ! grep -q "@vercel/speed-insights" package.json; then
  echo "Installing Vercel Analytics and Speed Insights..."
  npm install @vercel/analytics@1.1.1 @vercel/speed-insights@1.0.2 --save
fi

# Check if the layout includes analytics components
if [ -f "src/app/layout.tsx" ]; then
  if ! grep -q "Analytics" src/app/layout.tsx || ! grep -q "SpeedInsights" src/app/layout.tsx; then
    echo "Adding Analytics components to layout.tsx..."
    
    # Backup the file
    cp src/app/layout.tsx src/app/layout.tsx.bak
    
    # Check imports
    if ! grep -q "@vercel/analytics/react" src/app/layout.tsx; then
      sed -i '1i import { Analytics } from "@vercel/analytics/react"' src/app/layout.tsx
    fi
    
    if ! grep -q "@vercel/speed-insights/next" src/app/layout.tsx; then
      sed -i '1i import { SpeedInsights } from "@vercel/speed-insights/next"' src/app/layout.tsx
    fi
    
    # Check if components are in the body
    if ! grep -q "<Analytics" src/app/layout.tsx; then
      sed -i 's/<\/body>/  <Analytics \/>\n  <\/body>/' src/app/layout.tsx
    fi
    
    if ! grep -q "<SpeedInsights" src/app/layout.tsx; then
      sed -i 's/<\/body>/  <SpeedInsights \/>\n  <\/body>/' src/app/layout.tsx
    fi
  fi
fi

# Navigate back to project root
cd ..

echo "Preparing for deployment..."
echo "Using protection bypass token for API access"
echo "API URL: $API_URL"
echo "Frontend URL: $FRONTEND_URL"

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

# Check deployment result
if [ $? -eq 0 ]; then
  echo "✅ Deployment successful!"
  echo "Your site is now live with Vercel Analytics and Speed Insights enabled."
  
  # Notify Heroku if API key exists
  if [ -n "$HEROKU_API_KEY" ] && [ -n "$HEROKU_APP_NAME" ]; then
    echo "Notifying Heroku backend about frontend URL..."
    curl -X PATCH \
      -H "Content-Type: application/json" \
      -H "Accept: application/vnd.heroku+json; version=3" \
      -H "Authorization: Bearer $HEROKU_API_KEY" \
      -d "{\"config\": {\"VERCEL_FRONTEND_URL\": \"$FRONTEND_URL\"}}" \
      "https://api.heroku.com/apps/${HEROKU_APP_NAME}/config-vars"
    echo "Heroku configuration updated with frontend URL."
  fi
else
  echo "❌ Deployment failed."
  echo "Try running with --debug for more information."
fi