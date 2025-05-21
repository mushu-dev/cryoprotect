#!/bin/bash
# Script for integrating Analytics into the existing CryoProtect project
# This script directly modifies the actual frontend files to add Analytics

# Set up environment variables
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V"
FRONTEND_URL="https://cryoprotect.vercel.app"
API_URL="https://your-api-endpoint.com/api"
NEXTAUTH_SECRET="your-nextauth-secret-key"

# Check if we're in the right directory
if [ ! -d "frontend" ]; then
    echo "Error: Must be run from the project root directory"
    exit 1
fi

cd frontend || exit 1

# Install Vercel Analytics and Speed Insights
echo "Installing Vercel Analytics and Speed Insights..."
npm install @vercel/analytics@1.1.1 @vercel/speed-insights@1.0.2 --save

# Check for existing layout file (App Router)
if [ -f "src/app/layout.tsx" ]; then
    echo "Found App Router layout file, adding Analytics components..."
    # Backup the layout file
    cp src/app/layout.tsx src/app/layout.tsx.bak
    
    # Check if Analytics is already imported
    if ! grep -q "@vercel/analytics/react" src/app/layout.tsx; then
        # Add Analytics import at the top
        sed -i '1i import { Analytics } from "@vercel/analytics/react";' src/app/layout.tsx
    fi
    
    # Check if SpeedInsights is already imported
    if ! grep -q "@vercel/speed-insights/next" src/app/layout.tsx; then
        # Add SpeedInsights import at the top
        sed -i '1i import { SpeedInsights } from "@vercel/speed-insights/next";' src/app/layout.tsx
    fi
    
    # Add the components to the layout body (if not already there)
    if ! grep -q "<Analytics />" src/app/layout.tsx; then
        sed -i 's/<\/body>/  <Analytics \/>\n  <\/body>/' src/app/layout.tsx
    fi
    
    if ! grep -q "<SpeedInsights />" src/app/layout.tsx; then
        sed -i 's/<\/body>/  <SpeedInsights \/>\n  <\/body>/' src/app/layout.tsx
    fi
# Check for _app.js/tsx (Pages Router)
elif [ -f "src/pages/_app.tsx" ] || [ -f "src/pages/_app.js" ]; then
    APP_FILE=""
    if [ -f "src/pages/_app.tsx" ]; then
        APP_FILE="src/pages/_app.tsx"
    else
        APP_FILE="src/pages/_app.js"
    fi
    
    echo "Found Pages Router app file, adding Analytics components..."
    # Backup the app file
    cp "$APP_FILE" "${APP_FILE}.bak"
    
    # Check if Analytics is already imported
    if ! grep -q "@vercel/analytics/react" "$APP_FILE"; then
        # Add Analytics import at the top
        sed -i '1i import { Analytics } from "@vercel/analytics/react";' "$APP_FILE"
    fi
    
    # Check if SpeedInsights is already imported
    if ! grep -q "@vercel/speed-insights/next" "$APP_FILE"; then
        # Add SpeedInsights import at the top
        sed -i '1i import { SpeedInsights } from "@vercel/speed-insights/next";' "$APP_FILE"
    fi
    
    # Add the components to the app return (if not already there)
    if ! grep -q "<Analytics />" "$APP_FILE"; then
        sed -i 's/<\/.*>/  <Analytics \/>\n&/' "$APP_FILE"
    fi
    
    if ! grep -q "<SpeedInsights />" "$APP_FILE"; then
        sed -i 's/<\/.*>/  <SpeedInsights \/>\n&/' "$APP_FILE"
    fi
else
    echo "Error: Could not find layout.tsx or _app.js/tsx. Please add Analytics components manually."
    exit 1
fi

# Update next.config.js to simplify the build
echo "Updating next.config.js..."
if [ -f "next.config.js" ]; then
    cp next.config.js next.config.js.bak
    
    # Check if the file already contains analytics settings
    if ! grep -q "analyticsId" next.config.js; then
        # Check if there's a module.exports line to modify
        if grep -q "module.exports" next.config.js; then
            # Modify existing config
            sed -i 's/module.exports = {/module.exports = {\n  analyticsId: true,\n  speedInsights: {\n    enabled: true,\n  },/' next.config.js
        else
            # Create new config
            cat > next.config.js << 'EOF'
/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  analyticsId: true,
  speedInsights: {
    enabled: true,
  },
}

module.exports = nextConfig
EOF
        fi
    fi
fi

# Update package.json to fix React version if needed
echo "Checking package.json for potential React version conflicts..."
if [ -f "package.json" ]; then
    # Check if React version is compatible
    if grep -q "\"react\": \"19" package.json; then
        echo "Warning: Found React 19.x in package.json. This may cause conflicts."
        echo "Would you like to fix the React version to 18.2.0? (y/n)"
        read -r response
        if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
            sed -i 's/"react": "19[^"]*"/"react": "18.2.0"/' package.json
            sed -i 's/"react-dom": "19[^"]*"/"react-dom": "18.2.0"/' package.json
            npm install
        fi
    fi
fi

# Create/update deployment script
echo "Updating deploy-to-vercel.sh to include Analytics environment variables..."
cat > deploy-to-vercel.sh << 'EOF'
#!/bin/bash
# Deploy the frontend to Vercel with Analytics enabled

# Load environment variables if .env exists
if [ -f ../.env ]; then
  source ../.env
fi

# Set protection bypass token if not already set
if [ -z "$PROTECTION_BYPASS" ]; then
  PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V"
fi

# Set frontend URL if not already set
if [ -z "$FRONTEND_URL" ]; then
  FRONTEND_URL="https://cryoprotect.vercel.app"
fi

# Set API URL if not already set
if [ -z "$API_URL" ]; then
  API_URL="https://your-api-endpoint.com/api"
fi

echo "Deploying frontend to Vercel..."
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
EOF

chmod +x deploy-to-vercel.sh

echo "Done! Analytics components have been added to your CryoProtect project."
echo "To deploy to Vercel with analytics enabled, run:"
echo "cd frontend && ./deploy-to-vercel.sh"

echo "Done!"